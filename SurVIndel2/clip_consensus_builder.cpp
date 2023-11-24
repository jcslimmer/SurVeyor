#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <sstream>
#include <unordered_map>

#include "utils.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
#include "../libs/IntervalTree.h"
#include "sam_utils.h"
#include "stat_tests.h"
#include "extend_1sr_consensus.h"
#include "../src/sw_utils.h"

config_t config;
stats_t stats;
std::string workdir, complete_bam_fname, reference_fname;
std::mutex out_mtx, log_mtx;

sv2_chr_seqs_map_t chr_seqs;
sv2_contig_map_t contig_map;

std::mutex bam_pool_mtx;
std::queue<open_samFile_t*> bam_pool;
open_samFile_t* get_bam_reader(std::string bam_fname) {
    open_samFile_t* o;
    bam_pool_mtx.lock();
    if (!bam_pool.empty()) {
        o = bam_pool.front();
        bam_pool.pop();
    } else {
        o = open_samFile(bam_fname.c_str());
        hts_set_fai_filename(o->file, fai_path(reference_fname.c_str()));
    }
    bam_pool_mtx.unlock();
    return o;
}
void release_bam_reader(open_samFile_t* reader) {
    bam_pool_mtx.lock();
    bam_pool.push(reader);
    bam_pool_mtx.unlock();
}

std::unordered_map<std::string, std::vector<sv2_deletion_t*> > deletions_by_chr;
std::unordered_map<std::string, std::vector<sv2_duplication_t*> > duplications_by_chr;
std::unordered_map<std::string, std::vector<consensus_t*> > rc_sr_consensuses_by_chr, lc_sr_consensuses_by_chr, rc_hsr_consensuses_by_chr, lc_hsr_consensuses_by_chr;
std::unordered_map<std::string, std::vector<consensus_t*> > unpaired_consensuses_by_chr;
std::mutex mtx, indel_out_mtx, up_consensus_mtx;


struct pair_w_score_t {
    int rc_idx, lc_idx;
    int score;

    pair_w_score_t(int rc_idx, int lc_idx, int score) : rc_idx(rc_idx), lc_idx(lc_idx), score(score) {}
};

void extend_consensuses(int id, std::vector<consensus_t*>* consensuses, std::string contig_name,
	std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq, int start_idx, int end_idx) {

	std::vector<consensus_t*> consensuses_to_consider(consensuses->begin()+start_idx, consensuses->begin()+end_idx);

	open_samFile_t* bam_file = open_samFile(complete_bam_fname);
	std::vector<ext_read_t*> candidate_reads_for_extension = get_extension_reads_from_consensuses(consensuses_to_consider, contig_name, chr_seqs.get_len(contig_name), *mateseqs_w_mapq, stats, bam_file);

	std::vector<Interval<ext_read_t*>> it_ivals;
	for (ext_read_t* ext_read : candidate_reads_for_extension) {
		Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
		it_ivals.push_back(it_ival);
	}
	IntervalTree<ext_read_t*> candidate_reads_for_extension_itree(it_ivals);
	for (consensus_t* consensus : consensuses_to_consider) {
		if (consensus->left_clipped) {
			extend_consensus_to_right(consensus, candidate_reads_for_extension_itree, consensus->right_ext_target_start(stats.max_is, stats.read_len), 
				consensus->right_ext_target_end(stats.max_is, stats.read_len), contig_name, chr_seqs.get_len(contig_name), config.high_confidence_mapq, stats, *mateseqs_w_mapq);
		} else {
			extend_consensus_to_left(consensus, candidate_reads_for_extension_itree, consensus->left_ext_target_start(stats.max_is, stats.read_len), 
				consensus->left_ext_target_end(stats.max_is, stats.read_len), contig_name, chr_seqs.get_len(contig_name), config.high_confidence_mapq, stats, *mateseqs_w_mapq);
		}
	}
	for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

	close_samFile(bam_file);
}

void find_indels_from_rc_lc_pairs(std::string contig_name, std::vector<consensus_t*>& rc_consensuses, std::vector<consensus_t*>& lc_consensuses,
		std::vector<sv2_deletion_t*>& contig_deletions, std::vector<sv2_duplication_t*>& contig_duplications, StripedSmithWaterman::Aligner& aligner,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	// build interval tree of left-clipped consensuses (for quick search)
	std::vector<Interval<int>> consensus_iv;
	for (int i = 0; i < lc_consensuses.size(); i++) {
		consensus_t* lc_anchor = lc_consensuses[i];
		consensus_iv.push_back(Interval<int>(lc_anchor->breakpoint, lc_anchor->breakpoint+1, i));
	}
	IntervalTree<int> consensus_ivtree = IntervalTree<int>(consensus_iv);

	std::vector<pair_w_score_t> rc_lc_scored_pairs;
	for (int i = 0; i < rc_consensuses.size(); i++) {
		consensus_t* rc_anchor = rc_consensuses[i];

		// let us query the interval tree for left-clipped clusters in compatible positions
		// the compatible positions are calculated based on the mates of the right-clipped cluster
		hts_pos_t query_start, query_end;
		const int MAX_SEARCH_DIST = 10000;
		if (rc_anchor->remap_boundary == INT32_MAX) {
			query_start = rc_anchor->breakpoint - MAX_SEARCH_DIST;
			query_end = rc_anchor->breakpoint + MAX_SEARCH_DIST;
		} else {
			query_start = rc_anchor->remap_boundary - stats.max_is;
			query_end = rc_anchor->remap_boundary;
		}
		std::vector<Interval<int>> compatible_lc_idxs = consensus_ivtree.findOverlapping(query_start, query_end);

		for (auto& iv : compatible_lc_idxs) {
			consensus_t* lc_anchor = lc_consensuses[iv.value];
			if (lc_anchor->max_mapq < config.high_confidence_mapq && rc_anchor->max_mapq < config.high_confidence_mapq) continue;

			int min_overlap = (lc_anchor->is_hsr && rc_anchor->is_hsr) ? 50 : std::min(rc_anchor->clip_len, lc_anchor->clip_len)+config.min_clip_len;
			double max_mm_rate = (lc_anchor->is_hsr && rc_anchor->is_hsr) ? 0 : config.max_seq_error;

			suffix_prefix_aln_t spa = aln_suffix_prefix(rc_anchor->sequence, lc_anchor->sequence, 1, -4, max_mm_rate, min_overlap);
			if (spa.overlap > 0 && !is_homopolymer(lc_anchor->sequence.c_str(), spa.overlap)) {
				rc_lc_scored_pairs.push_back(pair_w_score_t(i, iv.value, spa.score));
			} else { // trim low quality (i.e., supported by less than 2 reads) bases
				// TODO: investigate if we can use base qualities for this
				std::string rc_consensus_trim = rc_anchor->sequence.substr(0, rc_anchor->sequence.length()-rc_anchor->lowq_clip_portion);
				std::string lc_consensus_trim = lc_anchor->sequence.substr(lc_anchor->lowq_clip_portion);

				suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus_trim, lc_consensus_trim, 1, -4, max_mm_rate, min_overlap);
				if (spa.overlap > 0 && !is_homopolymer(lc_consensus_trim.c_str(), spa.overlap)) {
					rc_lc_scored_pairs.push_back(pair_w_score_t(i, iv.value, spa.score));
				}
			}
		}
	}

	std::sort(rc_lc_scored_pairs.begin(), rc_lc_scored_pairs.end(),
	            [](const pair_w_score_t& ps1, const pair_w_score_t& ps2) {return ps1.score > ps2.score;});

	std::vector<bool> used_consensus_rc(rc_consensuses.size(), false), used_consensus_lc(lc_consensuses.size(), false);
	for (pair_w_score_t& ps : rc_lc_scored_pairs) {
		if (used_consensus_rc[ps.rc_idx] || used_consensus_lc[ps.lc_idx]) continue;

		consensus_t* rc_consensus = rc_consensuses[ps.rc_idx];
		consensus_t* lc_consensus = lc_consensuses[ps.lc_idx];

		int min_overlap = (lc_consensus->is_hsr && rc_consensus->is_hsr) ? 50 : std::min(rc_consensus->clip_len, lc_consensus->clip_len)+config.min_clip_len;
		double max_mm_rate = (lc_consensus->is_hsr && rc_consensus->is_hsr) ? 0 : config.max_seq_error;
		std::vector<sv_t*> svs = detect_svs(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), rc_consensus, lc_consensus, 
			aligner, min_overlap, config.min_clip_len, max_mm_rate);
		if (svs.empty()) continue;

		used_consensus_rc[ps.rc_idx] = used_consensus_lc[ps.lc_idx] = true;

		for (sv_t* sv : svs) {
			if (sv->svtype() == "DEL") {
				sv2_deletion_t* del = new sv2_deletion_t(sv);
				contig_deletions.push_back(del);
			} else {
				sv2_duplication_t* dup = new sv2_duplication_t(sv);
				contig_duplications.push_back(dup);
			}
		}
	}
	
	remove_marked_consensuses(rc_consensuses, used_consensus_rc);
	remove_marked_consensuses(lc_consensuses, used_consensus_lc);
}

void find_indels_from_unpaired_consensuses(int id, std::string contig_name, std::vector<consensus_t*>* consensuses,
		int start, int end, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq) {

	std::vector<sv2_deletion_t*> local_dels;
	std::vector<sv2_duplication_t*> local_dups;

	char* contig_seq = chr_seqs.get_seq(contig_name);
	hts_pos_t contig_len = chr_seqs.get_len(contig_name);

	StripedSmithWaterman::Aligner aligner(1,4,6,1,true);
	open_samFile_t* bam_file = open_samFile(complete_bam_fname);

	std::vector<consensus_t*> consensuses_to_consider(consensuses->begin()+start, consensuses->begin()+end);
	std::vector<ext_read_t*> candidate_reads_for_extension = get_extension_reads_from_consensuses(consensuses_to_consider, contig_name, contig_len, *mateseqs_w_mapq, stats, bam_file);

	std::vector<Interval<ext_read_t*>> it_ivals;
	for (ext_read_t* ext_read : candidate_reads_for_extension) {
		Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
		it_ivals.push_back(it_ival);
	}
	IntervalTree<ext_read_t*> candidate_reads_for_extension_itree(it_ivals);

	for (int i = start; i < consensuses->size() && i < end; i++) {
		consensus_t* consensus = consensuses->at(i);

		extend_consensus_to_left(consensus, candidate_reads_for_extension_itree, consensus->left_ext_target_start(stats.max_is, stats.read_len), 
			consensus->left_ext_target_end(stats.max_is, stats.read_len), contig_name, contig_len, config.high_confidence_mapq, stats, *mateseqs_w_mapq);
		extend_consensus_to_right(consensus, candidate_reads_for_extension_itree, consensus->right_ext_target_start(stats.max_is, stats.read_len), 
			consensus->right_ext_target_end(stats.max_is, stats.read_len), contig_name, contig_len, config.high_confidence_mapq, stats, *mateseqs_w_mapq);

		std::vector<indel_t*> indels;
		if (!consensus->left_clipped) {
			std::vector<sv_t*> svs = detect_svs(contig_name, contig_seq, contig_len, consensus, NULL, aligner, 0, config.min_clip_len, 0.0);
			for (sv_t* sv : svs) {
				indels.push_back(sv_to_indel(sv));
			}
		} else if (consensus->left_clipped) {
			std::vector<sv_t*> svs = detect_svs(contig_name, contig_seq, contig_len, NULL, consensus, aligner, 0, config.min_clip_len, 0.0);
			for (sv_t* sv : svs) {
				indels.push_back(sv_to_indel(sv));
			}
		}

		if (indels.empty()) {
			delete consensus;
			continue;
		}

		for (indel_t* indel : indels) {
			if (indel->sv->svtype() == "DEL") {
				local_dels.push_back((sv2_deletion_t*) indel);
			} else {
				local_dups.push_back((sv2_duplication_t*) indel);
			}
		}
	}
	close_samFile(bam_file);
	for (ext_read_t* read : candidate_reads_for_extension) delete read;

	indel_out_mtx.lock();
	deletions_by_chr[contig_name].insert(deletions_by_chr[contig_name].end(), local_dels.begin(), local_dels.end());
	duplications_by_chr[contig_name].insert(duplications_by_chr[contig_name].end(), local_dups.begin(), local_dups.end());
	indel_out_mtx.unlock();
}

// remove HSR clusters that overlap with a clipped position, i.e., the clipped position of a clipped cluster is contained in the HSR cluster
// direction of clip must be the same 
// clipped_consensuses must be sorted by breakpoint
void remove_hsr_overlapping_clipped(std::vector<consensus_t*>& hsr_consensuses, std::vector<consensus_t*>& clipped_consensuses) {
	std::sort(hsr_consensuses.begin(), hsr_consensuses.end(),
			[](const consensus_t* c1, const consensus_t* c2) { return c1->start < c2->start; });
	std::vector<consensus_t*> kept_consensus;
	int i = 0;
	for (consensus_t* c : hsr_consensuses) {
		while (i < clipped_consensuses.size() && clipped_consensuses[i]->breakpoint < c->start) i++;

		// clipped_consensuses[i] must be same clip direction as cluster and breakpoint must be the smallest s.t. >= cluster.start
		if (i >= clipped_consensuses.size() || clipped_consensuses[i]->breakpoint > c->end) {
			kept_consensus.push_back(c);
		} else {
			delete c;
		}
	}
	kept_consensus.swap(hsr_consensuses);
}

void read_consensuses(int id, int contig_id, std::string contig_name) {
	std::string chr, dir, seq;
    hts_pos_t start, end, breakpoint;
    int fwd_clipped, rev_clipped;
    int max_mapq, lowq_clip_portion;
    hts_pos_t remap_boundary;

	mtx.lock();
	std::vector<consensus_t*>& rc_sr_consensuses = rc_sr_consensuses_by_chr[contig_name];
	std::vector<consensus_t*>& lc_sr_consensuses = lc_sr_consensuses_by_chr[contig_name];
	std::vector<consensus_t*>& rc_hsr_consensuses = rc_hsr_consensuses_by_chr[contig_name];
	std::vector<consensus_t*>& lc_hsr_consensuses = lc_hsr_consensuses_by_chr[contig_name];
	mtx.unlock();

	std::string sr_consensuses_fname = workdir + "/workspace/sr_consensuses/" + std::to_string(contig_id) + ".txt";
	if (file_exists(sr_consensuses_fname)) {
		std::ifstream sr_consensuses_fin(sr_consensuses_fname);
		while (sr_consensuses_fin >> chr >> start >> end >> breakpoint >> dir >> seq >> 
			fwd_clipped >> rev_clipped >> max_mapq >> remap_boundary >> lowq_clip_portion) {
			if (dir == "L") {
				lc_sr_consensuses.push_back(new consensus_t(true, start, breakpoint, end, seq, fwd_clipped, rev_clipped, 
					breakpoint-start, max_mapq, remap_boundary, lowq_clip_portion));
			} else {
				rc_sr_consensuses.push_back(new consensus_t(false, start, breakpoint, end, seq, fwd_clipped, rev_clipped, 
					end-breakpoint, max_mapq, remap_boundary, lowq_clip_portion));
			}
		}
	}

	std::string hsr_consensuses_fname = workdir + "/workspace/hsr_consensuses/" + std::to_string(contig_id) + ".txt";
	if (file_exists(hsr_consensuses_fname)) {
		std::ifstream hsr_consensuses_fin(hsr_consensuses_fname);
		while (hsr_consensuses_fin >> chr >> start >> end >> breakpoint >> dir >> seq >> 
			fwd_clipped >> rev_clipped >> max_mapq >> remap_boundary >> lowq_clip_portion) {
			if (dir == "L") {
				consensus_t* consensus = new consensus_t(true, start, breakpoint, end, seq, fwd_clipped, rev_clipped, 
					breakpoint-start, max_mapq, remap_boundary, lowq_clip_portion);
				consensus->is_hsr = true;
				consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
				lc_hsr_consensuses.push_back(consensus);
			} else {
				consensus_t* consensus = new consensus_t(false, start, breakpoint, end, seq, fwd_clipped, rev_clipped, 
					end-breakpoint, max_mapq, remap_boundary, lowq_clip_portion);
				consensus->is_hsr = true;
				consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
				rc_hsr_consensuses.push_back(consensus);
			}
		}
	}

	// sort by breakpoint position
	auto consensus_cmp = [](consensus_t* c1, consensus_t* c2) {
		return c1->breakpoint < c2->breakpoint;
	};
    std::sort(rc_sr_consensuses.begin(), rc_sr_consensuses.end(), consensus_cmp);
    std::sort(lc_sr_consensuses.begin(), lc_sr_consensuses.end(), consensus_cmp);

	remove_hsr_overlapping_clipped(rc_hsr_consensuses, rc_sr_consensuses);
    remove_hsr_overlapping_clipped(lc_hsr_consensuses, lc_sr_consensuses);
}

void build_clip_consensuses(int id, int contig_id, std::string contig_name, std::vector<sv2_deletion_t*>& deletions,
                            std::vector<sv2_duplication_t*>& duplications, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	std::vector<sv2_deletion_t*> contig_deletions;
	std::vector<sv2_duplication_t*> contig_duplications;

	StripedSmithWaterman::Aligner aligner(1,4,6,1,true);
	StripedSmithWaterman::Filter filter;

	mtx.lock();
	std::vector<consensus_t*>& rc_sr_consensuses = rc_sr_consensuses_by_chr[contig_name];
	std::vector<consensus_t*>& lc_sr_consensuses = lc_sr_consensuses_by_chr[contig_name];
	std::vector<consensus_t*>& rc_hsr_consensuses = rc_hsr_consensuses_by_chr[contig_name];
	std::vector<consensus_t*>& lc_hsr_consensuses = lc_hsr_consensuses_by_chr[contig_name];
	mtx.unlock();

	find_indels_from_rc_lc_pairs(contig_name, rc_sr_consensuses, lc_sr_consensuses, contig_deletions, contig_duplications, aligner, mateseqs_w_mapq);
	find_indels_from_rc_lc_pairs(contig_name, rc_hsr_consensuses, lc_sr_consensuses, contig_deletions, contig_duplications, aligner, mateseqs_w_mapq);
	find_indels_from_rc_lc_pairs(contig_name, rc_sr_consensuses, lc_hsr_consensuses, contig_deletions, contig_duplications, aligner, mateseqs_w_mapq);
	find_indels_from_rc_lc_pairs(contig_name, rc_hsr_consensuses, lc_hsr_consensuses, contig_deletions, contig_duplications, aligner, mateseqs_w_mapq);

	deletions.insert(deletions.end(), contig_deletions.begin(), contig_deletions.end());
	duplications.insert(duplications.end(), contig_duplications.begin(), contig_duplications.end());

    /* == Deal with unpaired consensuses == */
	up_consensus_mtx.lock();
	std::vector<consensus_t*>& unpaired_consensuses = unpaired_consensuses_by_chr[contig_name];
	unpaired_consensuses.insert(unpaired_consensuses.end(), rc_sr_consensuses.begin(), rc_sr_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), lc_sr_consensuses.begin(), lc_sr_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), rc_hsr_consensuses.begin(), rc_hsr_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), lc_hsr_consensuses.begin(), lc_hsr_consensuses.end());
	up_consensus_mtx.unlock();
}


void build_consensuses(int id, int contig_id, std::string contig_name, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq) {
    mtx.lock();
    std::vector<sv2_deletion_t*>& deletions = deletions_by_chr[contig_name];
    std::vector<sv2_duplication_t*>& duplications = duplications_by_chr[contig_name];
    mtx.unlock();

	build_clip_consensuses(id, contig_id, contig_name, deletions, duplications, *mateseqs_w_mapq);
}

void size_and_depth_filtering(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    out_mtx.lock();
    std::vector<sv2_deletion_t*>& deletions = deletions_by_chr[contig_name];
    out_mtx.unlock();
    std::vector<double> temp1;
    std::vector<uint32_t> temp2;
    calculate_confidence_interval_size(contig_name, temp1, temp2, deletions, bam_file, stats, config.min_sv_size);
    depth_filter_del(contig_name, deletions, bam_file, config.min_size_for_depth_filtering, stats);

    out_mtx.lock();
    std::vector<sv2_duplication_t*>& duplications = duplications_by_chr[contig_name];
    out_mtx.unlock();
//    std::vector<duplication_t*> duplications_w_cleanup, duplications_wo_cleanup;
//    for (duplication_t* dup : duplications) {
//        if (dup->len() >= config.min_size_for_depth_filtering) duplications_w_cleanup.push_back(dup);
//        else duplications_wo_cleanup.push_back(dup);
//    }
//    depth_filter_dup_w_cleanup(contig_name, duplications_w_cleanup, bam_file, stats, config, max_allowed_frac_normalized, workdir);
//    depth_filter_dup(contig_name, duplications_wo_cleanup, bam_file, config.min_size_for_depth_filtering, config);
	depth_filter_dup(contig_name, duplications, bam_file, config.min_size_for_depth_filtering, stats);
    release_bam_reader(bam_file);
}


int main(int argc, char* argv[]) {

    complete_bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    reference_fname = argv[3];
    std::string sample_name = argv[4];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt");

    chr_seqs.read_fasta_into_map(reference_fname);

	ctpl::thread_pool read_consensuses_thread_pool(config.threads);
    std::vector<std::future<void> > futures;
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = read_consensuses_thread_pool.push(read_consensuses, contig_id, contig_name);
        futures.push_back(std::move(future));
	}
	read_consensuses_thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
	futures.clear();
	
	std::cout << "Extending consensuses." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();

	int block_size = 1000;
	std::vector<std::unordered_map<std::string, std::pair<std::string, int> > > mateseqs_w_mapq(contig_map.size());
	ctpl::thread_pool extend_consensuses_thread_pool(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
		std::ifstream fin(fname);
		std::string qname, read_seq, qual; int mapq;
		while (fin >> qname >> read_seq >> qual >> mapq) {
			mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
		}

		std::future<void> future;
		std::string contig_name = contig_map.get_name(contig_id);
		std::vector<consensus_t*>& rc_sr_consensuses = rc_sr_consensuses_by_chr[contig_name];
		for (int i = 0; i < rc_sr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) rc_sr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &rc_sr_consensuses, contig_name, &mateseqs_w_mapq[contig_id], start, end);
			futures.push_back(std::move(future));
		}
		std::vector<consensus_t*>& lc_sr_consensuses = lc_sr_consensuses_by_chr[contig_name];
		for (int i = 0; i < lc_sr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) lc_sr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &lc_sr_consensuses, contig_name, &mateseqs_w_mapq[contig_id], start, end);
			futures.push_back(std::move(future));
		}
		std::vector<consensus_t*>& rc_hsr_consensuses = rc_hsr_consensuses_by_chr[contig_name];
		for (int i = 0; i < rc_hsr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) rc_hsr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &rc_hsr_consensuses, contig_name, &mateseqs_w_mapq[contig_id], start, end);
			futures.push_back(std::move(future));
		}
		std::vector<consensus_t*>& lc_hsr_consensuses = lc_hsr_consensuses_by_chr[contig_name];
		for (int i = 0; i < lc_hsr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) lc_hsr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &lc_hsr_consensuses, contig_name, &mateseqs_w_mapq[contig_id], start, end);
			futures.push_back(std::move(future));
		}
	}
	extend_consensuses_thread_pool.stop(true);
	for (size_t i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Consensuses extended in " << elapsed_time << " seconds" << std::endl;

	std::cout << "Finding indels." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();

    ctpl::thread_pool thread_pool1(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool1.push(build_consensuses, contig_id, contig_name, &mateseqs_w_mapq[contig_id]);
        futures.push_back(std::move(future));
    }
    thread_pool1.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Indels found in " << elapsed_time << " seconds" << std::endl;

	std::cout << "Finding indels from unpaired consensuses." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();

    ctpl::thread_pool thread_pool2(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
    	std::vector<consensus_t*>& consensuses = unpaired_consensuses_by_chr[contig_name];
    	if (consensuses.empty()) continue;
    	for (int i = 0; i <= consensuses.size()/block_size; i++) {
    		int start = i * block_size;
    		int end = std::min(start+block_size, (int) consensuses.size());
			std::future<void> future = thread_pool2.push(find_indels_from_unpaired_consensuses,
					contig_name, &consensuses, i*block_size, end, &mateseqs_w_mapq[contig_id]);
			futures.push_back(std::move(future));
		}
    }
    thread_pool2.stop(true);
	for (size_t i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Indels found from unpaired consensuses in " << elapsed_time << " seconds" << std::endl;

	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		auto& deletions = deletions_by_chr[contig_name];
		for (int i = 0; i < deletions.size(); i++) {
			if (-deletions[i]->sv->svlen() < config.min_sv_size) {
				delete deletions[i];
				deletions[i] = NULL;
			}
		}
		deletions_by_chr[contig_name].erase(std::remove(deletions_by_chr[contig_name].begin(), deletions_by_chr[contig_name].end(), (sv2_deletion_t*) NULL), deletions_by_chr[contig_name].end());
		auto& duplications = duplications_by_chr[contig_name];
		for (int i = 0; i < duplications.size(); i++) {
			if (duplications[i]->sv->svlen() < config.min_sv_size) {
				delete duplications[i];
				duplications[i] = NULL;
			}
		}
		duplications_by_chr[contig_name].erase(std::remove(duplications_by_chr[contig_name].begin(), duplications_by_chr[contig_name].end(), (sv2_duplication_t*) NULL), duplications_by_chr[contig_name].end());
	}

    // create VCF out files
	bcf_hdr_t* out_vcf_header = sv2_generate_vcf_header(chr_seqs, sample_name, config, full_cmd_str);
    std::string out_vcf_fname = workdir + "/sr.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	std::cout << "Computing statistics for indels." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();

    ctpl::thread_pool thread_pool3(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool3.push(size_and_depth_filtering, contig_name);
        futures.push_back(std::move(future));
    }
    thread_pool3.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Statistics computed in " << elapsed_time << " seconds" << std::endl;

    std::unordered_map<std::string, std::vector<bcf1_t*>> bcf_entries;
    bcf1_t* bcf_entry = bcf_init();
    int del_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!deletions_by_chr.count(contig_name)) continue;
		std::vector<sv2_deletion_t*>& dels = deletions_by_chr[contig_name];
		std::sort(dels.begin(), dels.end(), [](sv2_deletion_t* del1, sv2_deletion_t* del2) {
			return std::tie(del1->sv->start, del1->sv->end) < std::tie(del2->sv->start, del2->sv->end);
		});
        for (sv2_deletion_t* del : dels) {
        	del->sv->id = "DEL_SR_" + std::to_string(del_id++);
            std::vector<std::string> filters;

            if (-del->sv->svlen() < config.min_sv_size) {
            	filters.push_back("SMALL");
            }
            if (-del->sv->svlen() < stats.max_is) {
            	if (-del->sv->svlen()/2 > del->max_conf_size) filters.push_back("SIZE_FILTER");
            }
            if (del->remap_boundary_lower > del->sv->start) {
                filters.push_back("REMAP_BOUNDARY_FILTER");
            } else if (del->remap_boundary_upper < del->sv->end) {
                filters.push_back("REMAP_BOUNDARY_FILTER");
            }
            if (-del->sv->svlen() >= config.min_size_for_depth_filtering &&
                    (del->med_left_flanking_cov*0.74<=del->med_indel_left_cov || del->med_right_flanking_cov*0.74<=del->med_indel_right_cov)) {
            	filters.push_back("DEPTH_FILTER");
            }
            if (del->med_left_flanking_cov > stats.max_depth || del->med_right_flanking_cov > stats.max_depth ||
            	del->med_left_flanking_cov < stats.min_depth || del->med_right_flanking_cov < stats.min_depth ||
				del->med_left_cluster_cov > stats.max_depth || del->med_right_cluster_cov > stats.max_depth) {
                filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (del->med_indel_left_cov > stats.max_depth || del->med_indel_right_cov > stats.max_depth) {
                filters.push_back("ANOMALOUS_DEL_DEPTH");
            }
            if (del->sv->full_junction_aln != NULL && del->sv->left_anchor_aln->best_score+del->sv->right_anchor_aln->best_score-del->sv->full_junction_aln->best_score < config.min_score_diff) {
            	filters.push_back("WEAK_SPLIT_ALIGNMENT");
            }
            if ((del->sv->lc_consensus == NULL || del->sv->lc_consensus->max_mapq < config.high_confidence_mapq) &&
            	(del->sv->rc_consensus == NULL || del->sv->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
            }
            if (del->sv->source == "1SR_RC" || del->sv->source == "1HSR_RC") {
            	if (del->sv->rc_consensus->right_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            } else if (del->sv->source == "1SR_LC" || del->sv->source == "1HSR_LC") {
            	if (del->sv->lc_consensus->left_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            }

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            sv2_del2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, del, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
            delete del;
        }
    }

    int dup_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!duplications_by_chr.count(contig_name)) continue;
    	std::vector<sv2_duplication_t*>& dups = duplications_by_chr[contig_name];
    	std::sort(dups.begin(), dups.end(), [](sv2_duplication_t* dup1, sv2_duplication_t* dup2) {
    		return std::tie(dup1->sv->start, dup1->sv->end) < std::tie(dup2->sv->start, dup2->sv->end);
    	});
        for (sv2_duplication_t* dup : dups) {
        	dup->sv->id = "DUP_SR_" + std::to_string(dup_id++);
            std::vector<std::string> filters;

            if (dup->sv->svlen() >= config.min_size_for_depth_filtering &&
                    (dup->med_left_flanking_cov*1.26>=dup->med_indel_left_cov || dup->med_indel_right_cov<=dup->med_right_flanking_cov*1.26)) {
                // note: using >= so that a 0 0 0 0 depth will not be accepted
                filters.push_back("DEPTH_FILTER");
            }

            if (dup->med_left_flanking_cov > stats.max_depth || dup->med_right_flanking_cov > stats.max_depth ||
                dup->med_left_flanking_cov < stats.min_depth || dup->med_right_flanking_cov < stats.min_depth) {
                filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (dup->sv->full_junction_aln != NULL && dup->sv->left_anchor_aln->best_score+dup->sv->right_anchor_aln->best_score-dup->sv->full_junction_aln->best_score < config.min_score_diff) {
				filters.push_back("WEAK_SPLIT_ALIGNMENT");
			}
            if ((dup->sv->lc_consensus == NULL || dup->sv->lc_consensus->max_mapq < config.high_confidence_mapq) &&
				(dup->sv->rc_consensus == NULL || dup->sv->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
			}
            if (dup->sv->svlen() >= config.min_size_for_depth_filtering && dup->disc_pairs < 3) {
                filters.push_back("NOT_ENOUGH_OW_PAIRS");
            }
            if (dup->sv->source == "1SR_RC" || dup->sv->source == "1HSR_RC") {
				if (dup->sv->rc_consensus->right_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            } else if (dup->sv->source == "1SR_LC" || dup->sv->source == "1HSR_LC") {
            	if (dup->sv->lc_consensus->left_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
			}

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            sv2_dup2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, dup, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
            delete dup;
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	auto& bcf_entries_contig = bcf_entries[contig_name];
    	std::sort(bcf_entries_contig.begin(), bcf_entries_contig.end(), [](const bcf1_t* b1, const bcf1_t* b2) {return b1->pos < b2->pos;});
		for (bcf1_t* bcf_entry : bcf_entries_contig) {
 			if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
			}
		}
	
		for (consensus_t* consensus : unpaired_consensuses_by_chr[contig_name]) delete consensus;
    }

    chr_seqs.clear();

    bcf_hdr_destroy(out_vcf_header);
    bcf_close(out_vcf_file);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
