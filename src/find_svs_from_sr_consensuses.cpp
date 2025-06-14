#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
#include "../libs/IntervalTree.h"
#include "htslib/hts.h"
#include "sam_utils.h"
#include "extend_1sr_consensus.h"
#include "sw_utils.h"
#include "types.h"
#include "vcf_utils.h"
#include "clustering_utils.h"

config_t config;
stats_t stats;
std::string workdir, complete_bam_fname, reference_fname;

chr_seqs_map_t chr_seqs;
contig_map_t contig_map;

std::unordered_map<std::string, std::vector<sv_t*> > svs_by_chr;
std::unordered_map<std::string, std::vector<std::shared_ptr<consensus_t>> > rc_sr_consensuses_by_chr, lc_sr_consensuses_by_chr, rc_hsr_consensuses_by_chr, lc_hsr_consensuses_by_chr;
std::unordered_map<std::string, std::vector<std::shared_ptr<consensus_t>> > unpaired_consensuses_by_chr;

std::unordered_map<std::string, std::vector<std::shared_ptr<cluster_t>> > ss_clusters_by_chr;
std::unordered_map<std::string, std::vector<breakend_t*> > dp_bnds_by_chr;

std::vector<std::unordered_map<std::string, std::pair<std::string, int> > > mateseqs_w_mapq;
std::vector<int> active_threads_per_chr;
std::vector<std::mutex> mutex_per_chr;

std::mutex mtx;


struct pair_w_score_t {
    int c1_idx, c2_idx;
	bool c1_lc, c2_lc;
	suffix_prefix_aln_t spa;
	std::shared_ptr<cluster_t> dp_cluster;

    pair_w_score_t(int c1_idx, bool c1_lc, int c2_idx, bool c2_lc, suffix_prefix_aln_t spa, std::shared_ptr<cluster_t> dp_cluster) : 
		c1_idx(c1_idx), c1_lc(c1_lc), c2_idx(c2_idx), c2_lc(c2_lc), spa(spa), dp_cluster(dp_cluster) {}
};

void extend_consensuses(int id, std::vector<std::shared_ptr<consensus_t>>* consensuses, std::string contig_name, int start_idx, int end_idx, bool extend_in_clip_direction) {

	int contig_id = contig_map.get_id(contig_name);
	mutex_per_chr[contig_id].lock();
	if (active_threads_per_chr[contig_id] == 0) {
		std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
		std::ifstream fin(fname);
		std::string qname, read_seq, qual; int mapq;
		while (fin >> qname >> read_seq >> qual >> mapq) {
			mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
		}
	}
	active_threads_per_chr[contig_id]++;
	mutex_per_chr[contig_id].unlock();

	std::vector<std::shared_ptr<consensus_t>> consensuses_to_consider(consensuses->begin()+start_idx, consensuses->begin()+end_idx);

	
	open_samFile_t* bam_file = open_samFile(complete_bam_fname);
	if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + reference_fname;
	}
	std::vector<ext_read_t*> candidate_reads_for_extension = get_extension_reads_from_consensuses(consensuses_to_consider, contig_name, chr_seqs.get_len(contig_name), stats, bam_file);
	if (!candidate_reads_for_extension.empty()) {
		std::vector<Interval<ext_read_t*>> it_ivals;
		for (ext_read_t* ext_read : candidate_reads_for_extension) {
			Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
			it_ivals.push_back(it_ival);
		}
		IntervalTree<ext_read_t*> candidate_reads_for_extension_itree(it_ivals);
		for (std::shared_ptr<consensus_t> consensus : consensuses_to_consider) {
			if (consensus->left_clipped != extend_in_clip_direction) {
				extend_consensus_to_right(consensus, candidate_reads_for_extension_itree, consensus->right_ext_target_start(stats.max_is, stats.read_len), 
					consensus->right_ext_target_end(stats.max_is, stats.read_len), chr_seqs.get_len(contig_name), config.high_confidence_mapq, stats, mateseqs_w_mapq[contig_id]);
			} else {
				extend_consensus_to_left(consensus, candidate_reads_for_extension_itree, consensus->left_ext_target_start(stats.max_is, stats.read_len), 
					consensus->left_ext_target_end(stats.max_is, stats.read_len), chr_seqs.get_len(contig_name), config.high_confidence_mapq, stats, mateseqs_w_mapq[contig_id]);
			}
		}
		for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;
	}

	close_samFile(bam_file);

	mutex_per_chr[contig_id].lock();
	active_threads_per_chr[contig_id]--;
	if (active_threads_per_chr[contig_id] == 0) {
		mateseqs_w_mapq[contig_id].clear();
	}
	mutex_per_chr[contig_id].unlock();
}

void remove_marked_consensuses(std::vector<std::shared_ptr<consensus_t>>& consensuses, std::vector<bool>& used) {
	for (int i = 0; i < consensuses.size(); i++) if (used[i]) consensuses[i] = nullptr;
	consensuses.erase(std::remove(consensuses.begin(), consensuses.end(), nullptr), consensuses.end());
}

void find_indels_from_rc_lc_pairs(std::string contig_name, std::vector<std::shared_ptr<consensus_t>>& rc_consensuses, std::vector<std::shared_ptr<consensus_t>>& lc_consensuses,
		StripedSmithWaterman::Aligner& aligner) {

	std::vector<sv_t*> local_svs;

	auto min_overlap_f = [](consensus_t* c1, consensus_t* c2) {
		return c1->is_hsr && c2->is_hsr ? 50 : config.min_clip_len;
	};
	auto max_seq_error_f = [](consensus_t* c1, consensus_t* c2) {
		return c1->is_hsr && c2->is_hsr ? 0 : config.max_seq_error;
	};

	int initial_rc_consensus_size = rc_consensuses.size();
	int initial_lc_consensus_size = lc_consensuses.size();

	// build interval tree of left-clipped consensuses (for quick search)
	std::vector<Interval<int>> lc_consensus_iv;
	for (int i = 0; i < lc_consensuses.size(); i++) {
		std::shared_ptr<consensus_t> lc_anchor = lc_consensuses[i];
		lc_consensus_iv.push_back(Interval<int>(lc_anchor->breakpoint, lc_anchor->breakpoint+1, i));
	}
	IntervalTree<int> lc_consensus_ivtree = IntervalTree<int>(lc_consensus_iv);

	std::vector<Interval<int>> rc_consensus_iv;
	for (int i = 0; i < rc_consensuses.size(); i++) {
		std::shared_ptr<consensus_t> rc_anchor = rc_consensuses[i];
		rc_consensus_iv.push_back(Interval<int>(rc_anchor->breakpoint, rc_anchor->breakpoint+1, i));
	}
	IntervalTree<int> rc_consensus_ivtree = IntervalTree<int>(rc_consensus_iv);

	std::vector<pair_w_score_t> consensuses_scored_pairs;
	mtx.lock();
	std::vector<std::shared_ptr<cluster_t>>& ss_clusters = ss_clusters_by_chr[contig_name];
	mtx.unlock();
	// find viable pairs of consensuses that form BND, and push them to consensuses_scored_pairs
	for (std::shared_ptr<cluster_t> c : ss_clusters) {

		std::vector<Interval<int>> compatible_la_idxs, compatible_ra_idxs;
		if (c->la_rev) {
			compatible_la_idxs = lc_consensus_ivtree.findOverlapping(c->la_start-stats.max_is, c->la_start+stats.read_len/2);
			compatible_ra_idxs = lc_consensus_ivtree.findOverlapping(c->ra_start-stats.max_is, c->ra_end+stats.read_len/2);
		} else {
			compatible_la_idxs = rc_consensus_ivtree.findOverlapping(c->la_end-stats.read_len/2, c->la_end+stats.max_is);
			compatible_ra_idxs = rc_consensus_ivtree.findOverlapping(c->ra_end-stats.read_len/2, c->ra_end+stats.max_is);	
		}

		// for SR consensues that could belong to both the left and the right anchor, assign them to the most suitable anchor
		for (Interval<int>& la_iv : compatible_la_idxs) {
			for (Interval<int>& ra_iv : compatible_ra_idxs) {
				if (la_iv.value == ra_iv.value) {
					if (c->la_rev) {
						std::shared_ptr<consensus_t> sr_cluster = lc_consensuses[la_iv.value];
						if (sr_cluster->breakpoint <= c->la_start+10) ra_iv.value = -1;
						else la_iv.value = -1;
					} else {
						std::shared_ptr<consensus_t> sr_cluster = rc_consensuses[la_iv.value];
						if (sr_cluster->breakpoint >= c->ra_end-10) la_iv.value = -1;
						else ra_iv.value = -1;
					}
					break;
				}
			}
		}
		compatible_la_idxs.erase(std::remove_if(compatible_la_idxs.begin(), compatible_la_idxs.end(),
			[](Interval<int>& iv) { return iv.value == -1; }), compatible_la_idxs.end());
		compatible_ra_idxs.erase(std::remove_if(compatible_ra_idxs.begin(), compatible_ra_idxs.end(),
			[](Interval<int>& iv) { return iv.value == -1; }), compatible_ra_idxs.end());

		if (compatible_la_idxs.empty() != compatible_ra_idxs.empty()) {
			std::shared_ptr<consensus_t> consensus;
			if (compatible_la_idxs.empty()) {
				hts_pos_t breakpoint = c->la_rev ? c->la_start : c->la_end;
				consensus = std::make_shared<consensus_t>(c->la_rev, c->la_start, breakpoint, c->la_end, c->la_furthermost_seq, !c->la_rev, c->la_rev, 0, c->la_max_mapq, 0, 0);
			} else {
				hts_pos_t breakpoint = c->ra_rev ? c->ra_start : c->ra_end;
				consensus = std::make_shared<consensus_t>(c->ra_rev, c->ra_start, breakpoint, c->ra_end, c->ra_furthermost_seq, !c->ra_rev, c->ra_rev, 0, c->ra_max_mapq, 0, 0);
			}
			consensus->is_hsr = true;
			int i;
			if (c->la_rev) {
				i = lc_consensuses.size();
				lc_consensuses.push_back(consensus);
			} else {
				i = rc_consensuses.size();
				rc_consensuses.push_back(consensus);
			}
			if (compatible_la_idxs.empty()) {
				compatible_la_idxs.push_back(Interval<int>(0, 0, i));
			} else {
				compatible_ra_idxs.push_back(Interval<int>(0, 0, i));
			}
		}
		
		for (Interval<int>& la_iv : compatible_la_idxs) {
			std::shared_ptr<consensus_t> la_consensus = (c->la_rev ? lc_consensuses[la_iv.value] : rc_consensuses[la_iv.value]);
			for (Interval<int>& ra_iv : compatible_ra_idxs) {
				if (la_iv.value == ra_iv.value) continue;
				
				std::shared_ptr<consensus_t> ra_consensus = (c->la_rev ? lc_consensuses[ra_iv.value] : rc_consensuses[ra_iv.value]);
				std::shared_ptr<consensus_t> leftmost_consensus = la_consensus->breakpoint < ra_consensus->breakpoint ? la_consensus : ra_consensus;
				std::shared_ptr<consensus_t> rightmost_consensus = la_consensus->breakpoint < ra_consensus->breakpoint ? ra_consensus : la_consensus;
				
				std::string lm_seq = leftmost_consensus->sequence, rm_seq = rightmost_consensus->sequence;
				if (c->la_rev) rc(lm_seq);
				else rc(rm_seq);

				suffix_prefix_aln_t spa = aln_suffix_prefix(lm_seq, rm_seq, 1, -4, config.max_seq_error, config.min_clip_len);
				if (spa.overlap > 0) {
					consensuses_scored_pairs.push_back(pair_w_score_t(la_iv.value, la_consensus->left_clipped, ra_iv.value, ra_consensus->left_clipped, spa, c));
				}
			}
		}
	}

	for (int i = 0; i < rc_consensuses.size(); i++) { // find viable pairs of consensuses that form deletions and tandem dups, and push them to consensuses_scored_pairs
		std::shared_ptr<consensus_t> rc_consensus = rc_consensuses[i];

		// let us query the interval tree for left-clipped clusters in compatible positions
		// the compatible positions are calculated based on the mates of the right-clipped cluster
		hts_pos_t query_start, query_end;
		const int MAX_SEARCH_DIST = 10000;
		if (rc_consensus->remap_boundary == INT32_MAX) {
			query_start = rc_consensus->breakpoint - MAX_SEARCH_DIST;
			query_end = rc_consensus->breakpoint + MAX_SEARCH_DIST;
		} else {
			query_start = rc_consensus->remap_boundary - stats.max_is;
			query_end = rc_consensus->remap_boundary;
		}
		std::vector<Interval<int>> compatible_lc_idxs = lc_consensus_ivtree.findOverlapping(query_start, query_end);

		for (auto& iv : compatible_lc_idxs) {
			std::shared_ptr<consensus_t> lc_consensus = lc_consensuses[iv.value];
			if (lc_consensus->max_mapq < config.high_confidence_mapq && rc_consensus->max_mapq < config.high_confidence_mapq) continue;

			int min_overlap = min_overlap_f(rc_consensus.get(), lc_consensus.get());
			double max_mm_rate = max_seq_error_f(rc_consensus.get(), lc_consensus.get());

			// workaround because this can happen for inversions, when the right and the left clipped consensuses inside the the inverted part,
			// i.e., the left clipped consensus at the beginning of the inverted part and the right clipped consensus at the end of the inverted part
			// these two are often extended to the same sequence, and they form a strong (perfect) match, instead of being used to call precise inversions
			if (rc_consensus->sequence == lc_consensus->sequence) continue; 

			suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus->sequence, lc_consensus->sequence, 1, -4, max_mm_rate, min_overlap);
			std::shared_ptr<consensus_t> c1_consensus = rc_consensus->left_clipped ? rc_consensus : lc_consensus;
			std::shared_ptr<consensus_t> c2_consensus = rc_consensus->left_clipped ? lc_consensus : rc_consensus;
			if (spa.overlap > 0 && !is_homopolymer(lc_consensus->sequence.c_str(), spa.overlap) && spa.overlap <= 2*stats.read_len) { 
				consensuses_scored_pairs.push_back(pair_w_score_t(i, false, iv.value, true, spa, NULL));
			} else { // trim low quality (i.e., supported by less than 2 reads) bases
				// TODO: investigate if we can use base qualities for this
				std::string rc_consensus_trim = rc_consensus->sequence.substr(0, rc_consensus->sequence.length()-rc_consensus->lowq_clip_portion);
				std::string lc_consensus_trim = lc_consensus->sequence.substr(lc_consensus->lowq_clip_portion);

				suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus_trim, lc_consensus_trim, 1, -4, max_mm_rate, min_overlap);
				if (spa.overlap > 0 && !is_homopolymer(lc_consensus_trim.c_str(), spa.overlap) && spa.overlap <= 2*stats.read_len) {
					consensuses_scored_pairs.push_back(pair_w_score_t(i, false, iv.value, true, spa, NULL));
				}
			}
		}
	}

	// same-strand get priority, because the LC and RC SR consensues within small inversions may end up artificially
	// forming a pair
	std::sort(consensuses_scored_pairs.begin(), consensuses_scored_pairs.end(),
			[](const pair_w_score_t& ps1, const pair_w_score_t& ps2) {
				bool ps1_is_ss = ps1.c1_lc == ps1.c2_lc;
				bool ps2_is_ss = ps2.c1_lc == ps2.c2_lc;
				if (ps1_is_ss && !ps2_is_ss) return true;
				else if (!ps1_is_ss && ps2_is_ss) return false;
				else return ps1.spa.score > ps2.spa.score;
			});

	mtx.lock();
	std::vector<breakend_t*> bnds_lf, bnds_rf;
	mtx.unlock();
	std::vector<bool> used_consensus_rc(rc_consensuses.size(), false), used_consensus_lc(lc_consensuses.size(), false);

	for (pair_w_score_t& ps : consensuses_scored_pairs) {
		std::shared_ptr<consensus_t> c1_consensus = ps.c1_lc ? lc_consensuses[ps.c1_idx] : rc_consensuses[ps.c1_idx]; // in rc/lc pairs, this is rc
		std::shared_ptr<consensus_t> c2_consensus = ps.c2_lc ? lc_consensuses[ps.c2_idx] : rc_consensuses[ps.c2_idx]; // in rc/lc pairs, this is lc
		bool used_c1 = ps.c1_lc ? used_consensus_lc[ps.c1_idx] : used_consensus_rc[ps.c1_idx];
		bool used_c2 = ps.c2_lc ? used_consensus_lc[ps.c2_idx] : used_consensus_rc[ps.c2_idx];
		if (used_c1 || used_c2) continue;

		std::vector<sv_t*> svs;
		if (c1_consensus->left_clipped == c2_consensus->left_clipped) {
			std::shared_ptr<consensus_t> leftmost_consensus = c1_consensus->breakpoint < c2_consensus->breakpoint ? c1_consensus : c2_consensus;
			std::shared_ptr<consensus_t> rightmost_consensus = c1_consensus->breakpoint < c2_consensus->breakpoint ? c2_consensus : c1_consensus;
			breakend_t* bnd = detect_bnd(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), leftmost_consensus, rightmost_consensus, ps.spa, aligner, config.min_clip_len);
			if (bnd == NULL || bnd->end-bnd->start < config.min_sv_size) {
				delete bnd;
				continue;
			}
			bnd->source = "2SR";
			svs.push_back(bnd);
			if (c1_consensus->left_clipped) {
				bnds_lf.push_back(bnd);
			} else {
				bnds_rf.push_back(bnd);
			}
		} else {
			int min_overlap = min_overlap_f(c1_consensus.get(), c2_consensus.get());
			double max_mm_rate = max_seq_error_f(c1_consensus.get(), c2_consensus.get());
			svs = detect_svs(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), c1_consensus, c2_consensus,
				aligner, min_overlap, config.min_clip_len, max_mm_rate);
		}

		if (svs.empty()) continue;

		if (ps.c1_lc) used_consensus_lc[ps.c1_idx] = true;
		else used_consensus_rc[ps.c1_idx] = true;
		if (ps.c2_lc) used_consensus_lc[ps.c2_idx] = true;
		else used_consensus_rc[ps.c2_idx] = true;

		local_svs.insert(local_svs.end(), svs.begin(), svs.end());
	}

	// detect inversions from ss clusters
	std::vector<std::shared_ptr<cluster_t>> lf_ss_clusters, rf_ss_clusters;
	for (std::shared_ptr<cluster_t> c : ss_clusters) {
		if (c->la_rev) {
			lf_ss_clusters.push_back(c);
		} else {
			rf_ss_clusters.push_back(c);
		}
	}

	int bnd_idx = 0;
	std::sort(bnds_lf.begin(), bnds_lf.end(), [](const breakend_t* i1, const breakend_t* i2) { return i1->start < i2->start; });
	for (std::shared_ptr<cluster_t> c : lf_ss_clusters) {	
		while (bnd_idx < bnds_lf.size() && bnds_lf[bnd_idx]->start < c->la_end-stats.max_is) bnd_idx++;

		bool has_bnd = false;
		for (int i = bnd_idx; i < bnds_lf.size() && bnds_lf[i]->start < c->la_end; i++) {
			if (bnds_lf[i]->end >= c->ra_end-stats.max_is && bnds_lf[i]->end <= c->ra_end) {
				has_bnd = true;
				break;
			}
		}

		if (!has_bnd) {
			std::string lm_seq = c->la_furthermost_seq, rm_seq = c->ra_furthermost_seq;
			rc(lm_seq);
			suffix_prefix_aln_t spa = aln_suffix_prefix(lm_seq, rm_seq, 1, -4, config.max_seq_error, config.min_clip_len);
			breakend_t* bnd = NULL;
			if (spa.overlap > 0) {
				std::shared_ptr<consensus_t> leftmost_consensus = std::make_shared<consensus_t>(true, c->la_start, c->la_start, c->la_end, c->la_furthermost_seq, 0, 1, 0, c->la_max_mapq, 0, 0);
				std::shared_ptr<consensus_t> rightmost_consensus = std::make_shared<consensus_t>(true, c->ra_start, c->ra_start, c->ra_end, c->ra_furthermost_seq, 0, 1, 0, c->ra_max_mapq, 0, 0);
				bnd = detect_bnd(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), leftmost_consensus, rightmost_consensus, spa, aligner, config.min_clip_len);
			} else {
				auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->la_start, c->la_end, c->la_end-c->la_start, 0);
				auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->ra_start, c->ra_end, c->ra_end-c->ra_start, 0);
				bnd = new breakend_t(contig_name, c->la_start, c->ra_start, "", NULL, NULL, left_anchor_aln, right_anchor_aln, '-');
				bnd->imprecise = true;
			}

			if (bnd->end-bnd->start < config.min_sv_size) {
				delete bnd;
				continue;
			}
			bnd->source = "DP";
			bnds_lf.push_back(bnd);
		}
	}

	std::sort(bnds_rf.begin(), bnds_rf.end(), [](const breakend_t* i1, const breakend_t* i2) { return i1->start < i2->start; });
	bnd_idx = 0;
	for (std::shared_ptr<cluster_t> c : rf_ss_clusters) {
		while (bnd_idx < bnds_rf.size() && bnds_rf[bnd_idx]->start < c->la_start) bnd_idx++;

		bool has_bnd = false;
		for (int i = bnd_idx; i < bnds_rf.size() && bnds_rf[i]->start <= c->la_start+stats.max_is; i++) {
			if (bnds_rf[i]->end >= c->ra_start && bnds_rf[i]->end <= c->ra_start+stats.max_is) {
				has_bnd = true;
				break;
			}
		}

		if (!has_bnd) {
			std::string lm_seq = c->la_furthermost_seq, rm_seq = c->ra_furthermost_seq;
			rc(rm_seq);
			suffix_prefix_aln_t spa = aln_suffix_prefix(lm_seq, rm_seq, 1, -4, config.max_seq_error, config.min_clip_len);
			breakend_t* bnd = NULL;
			if (spa.overlap) {
				std::shared_ptr<consensus_t> leftmost_consensus = std::make_shared<consensus_t>(false, c->la_start, c->la_end, c->la_end, c->la_furthermost_seq, 1, 0, 0, c->la_max_mapq, 0, 0);
				std::shared_ptr<consensus_t> rightmost_consensus = std::make_shared<consensus_t>(false, c->ra_start, c->ra_end, c->ra_end, c->ra_furthermost_seq, 1, 0, 0, c->ra_max_mapq, 0, 0);
				bnd = detect_bnd(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), leftmost_consensus, rightmost_consensus, spa, aligner, config.min_clip_len);
			} else {
				auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->la_start, c->la_end, c->la_end-c->la_start, 0);
				auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->ra_start, c->ra_end, c->ra_end-c->ra_start, 0);
				if (left_anchor_aln->end > right_anchor_aln->end) {
					std::swap(left_anchor_aln, right_anchor_aln);
				}
				bnd = new breakend_t(contig_name, left_anchor_aln->end, right_anchor_aln->end, "", NULL, NULL, left_anchor_aln, right_anchor_aln, '+');
				bnd->imprecise = true;
			}

			if (bnd == NULL || bnd->end-bnd->start < config.min_sv_size) {
				delete bnd;
				continue;
			}
			bnd->source = "DP";
			bnds_rf.push_back(bnd);
		}
	}

	std::unordered_set<breakend_t*> used_breakends;
	std::vector<std::tuple<hts_pos_t, breakend_t*, breakend_t*>> bnd_pairs;
	for (breakend_t* bnd_rf : bnds_rf) {
		for (breakend_t* bnd_lf : bnds_lf) {
			if (bnd_lf->start-bnd_rf->start > -stats.read_len && bnd_lf->end-bnd_rf->end > -stats.read_len &&
				(bnd_lf->start-bnd_rf->start <= stats.max_is || bnd_lf->end-bnd_rf->end <= stats.max_is)) {
				hts_pos_t dist = std::abs(bnd_rf->start-bnd_lf->start) + std::abs(bnd_rf->end-bnd_lf->end);
				bnd_pairs.push_back(std::make_tuple(dist, bnd_rf, bnd_lf));
			}
		}
	}

	std::sort(bnd_pairs.begin(), bnd_pairs.end());
	for (auto& bnd_pair : bnd_pairs) {
		breakend_t* bnd_rf = std::get<1>(bnd_pair);
		breakend_t* bnd_lf = std::get<2>(bnd_pair);
		if (used_breakends.count(bnd_rf) || used_breakends.count(bnd_lf)) continue;
		used_breakends.insert(bnd_rf);
		used_breakends.insert(bnd_lf);
		std::string seq = "";

		std::shared_ptr<consensus_t> rc_consensus = NULL, lc_consensus = NULL;
		if (bnd_rf->rc_consensus != NULL) {
			rc_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, seq, bnd_rf->rc_consensus->fwd_reads+bnd_rf->lc_consensus->fwd_reads, 
				bnd_rf->rc_consensus->rev_reads+bnd_rf->lc_consensus->rev_reads, 0, std::max(bnd_rf->rc_consensus->max_mapq, bnd_rf->lc_consensus->max_mapq), 0, 0);
			rc_consensus->max_mapq = std::max(bnd_rf->rc_consensus->max_mapq, bnd_rf->lc_consensus->max_mapq);
			rc_consensus->left_ext_reads = bnd_rf->rc_consensus->left_ext_reads;
			rc_consensus->right_ext_reads = bnd_rf->lc_consensus->left_ext_reads;
			rc_consensus->hq_left_ext_reads = bnd_rf->rc_consensus->hq_left_ext_reads;
			rc_consensus->hq_right_ext_reads = bnd_rf->lc_consensus->hq_left_ext_reads;
		}
		if (bnd_lf->lc_consensus != NULL) {
			lc_consensus = std::make_shared<consensus_t>(true, 0, 0, 0, seq, bnd_lf->rc_consensus->fwd_reads+bnd_lf->lc_consensus->fwd_reads, 
				bnd_lf->rc_consensus->rev_reads+bnd_lf->lc_consensus->rev_reads, 0, std::max(bnd_lf->rc_consensus->max_mapq, bnd_lf->lc_consensus->max_mapq), 0, 0);
			lc_consensus->max_mapq = std::max(bnd_lf->rc_consensus->max_mapq, bnd_lf->lc_consensus->max_mapq);
			lc_consensus->left_ext_reads = bnd_lf->rc_consensus->right_ext_reads;
			lc_consensus->right_ext_reads = bnd_lf->lc_consensus->right_ext_reads;
			lc_consensus->hq_left_ext_reads = bnd_lf->rc_consensus->hq_right_ext_reads;
			lc_consensus->hq_right_ext_reads = bnd_lf->lc_consensus->hq_right_ext_reads;
		}

		bool imprecise = bnd_rf->imprecise || bnd_lf->imprecise;
		inversion_t* inv = NULL;
		hts_pos_t inv_start = bnd_lf->start, inv_end = bnd_rf->end;
		hts_pos_t sv_start = inv_start, sv_end = inv_end;
		if (!imprecise && (bnd_rf->start+10 <= bnd_lf->start || bnd_rf->end+10 <= bnd_lf->end)) {
			sv_start = bnd_rf->start, sv_end = bnd_lf->end;
		}
		inv = new inversion_t(contig_name, sv_start, sv_end, "", rc_consensus, lc_consensus, bnd_rf->left_anchor_aln, bnd_lf->left_anchor_aln, bnd_rf->right_anchor_aln, bnd_lf->right_anchor_aln);
		inv->inv_start = inv_start;
		inv->inv_end = inv_end;
		inv->source = bnd_rf->source + "-" + bnd_lf->source;
		inv->imprecise = imprecise;
		if (inv->inv_end-inv->inv_start < config.min_sv_size) {
			delete inv;
		} else {
			local_svs.push_back(inv);
		}
	}

	mtx.lock();
	std::vector<breakend_t*>& dp_bnds = dp_bnds_by_chr[contig_name];
	mtx.unlock();
	local_svs.erase(std::remove_if(local_svs.begin(), local_svs.end(), [](sv_t* sv) { return sv->svtype() == "BND"; }), local_svs.end());
	for (breakend_t* bnd : bnds_lf) {
		if (!used_breakends.count(bnd)) {
			if (bnd->source == "DP") {
				dp_bnds.push_back(bnd);
			} else {
				local_svs.push_back(bnd);
			}
		}
	}
	for (breakend_t* bnd : bnds_rf) {
		if (!used_breakends.count(bnd)) {
			if (bnd->source == "DP") {
				dp_bnds.push_back(bnd);
			} else {
				local_svs.push_back(bnd);
			}
		}
	}

	rc_consensuses.resize(initial_rc_consensus_size);
	lc_consensuses.resize(initial_lc_consensus_size);
	used_consensus_rc.resize(initial_rc_consensus_size);
	used_consensus_lc.resize(initial_lc_consensus_size);
	remove_marked_consensuses(rc_consensuses, used_consensus_rc);
	remove_marked_consensuses(lc_consensuses, used_consensus_lc);

	mtx.lock();
	std::vector<sv_t*>& svs = svs_by_chr[contig_name];
	svs.insert(svs.end(), local_svs.begin(), local_svs.end());
	mtx.unlock();

}

// remove HSR clusters that overlap with a clipped position, i.e., the clipped position of a clipped cluster is contained in the HSR cluster
// direction of clip must be the same 
// clipped_consensuses must be sorted by breakpoint
void remove_hsr_overlapping_clipped(std::vector<std::shared_ptr<consensus_t>>& hsr_consensuses, std::vector<std::shared_ptr<consensus_t>>& clipped_consensuses) {
	std::sort(hsr_consensuses.begin(), hsr_consensuses.end(),
			[](std::shared_ptr<consensus_t> c1, std::shared_ptr<consensus_t> c2) { return c1->start < c2->start; });
	std::vector<std::shared_ptr<consensus_t>> kept_consensus;
	int i = 0;
	for (std::shared_ptr<consensus_t> c : hsr_consensuses) {
		while (i < clipped_consensuses.size() && clipped_consensuses[i]->breakpoint < c->start) i++;

		// clipped_consensuses[i] must be same clip direction as cluster and breakpoint must be the smallest s.t. >= cluster.start
		if (i >= clipped_consensuses.size() || clipped_consensuses[i]->breakpoint > c->end) {
			kept_consensus.push_back(c);
		}
	}
	kept_consensus.swap(hsr_consensuses);
}

void read_consensuses(int id, int contig_id, std::string contig_name) {
	std::string dir, seq;
    hts_pos_t start, end, breakpoint;
    int fwd_clipped, rev_clipped;
    int max_mapq, lowq_clip_portion;
    hts_pos_t remap_boundary;

	mtx.lock();
	std::vector<std::shared_ptr<consensus_t>>& rc_sr_consensuses = rc_sr_consensuses_by_chr[contig_name];
	std::vector<std::shared_ptr<consensus_t>>& lc_sr_consensuses = lc_sr_consensuses_by_chr[contig_name];
	std::vector<std::shared_ptr<consensus_t>>& rc_hsr_consensuses = rc_hsr_consensuses_by_chr[contig_name];
	std::vector<std::shared_ptr<consensus_t>>& lc_hsr_consensuses = lc_hsr_consensuses_by_chr[contig_name];
	mtx.unlock();

	std::string sr_consensuses_fname = workdir + "/workspace/sr_consensuses/" + std::to_string(contig_id) + ".txt";
	if (file_exists(sr_consensuses_fname)) {
		std::ifstream sr_consensuses_fin(sr_consensuses_fname);
		std::string line;
		while (std::getline(sr_consensuses_fin, line)) {
			std::shared_ptr<consensus_t> consensus = std::make_shared<consensus_t>(line, false);
			if (consensus->left_clipped) {
				lc_sr_consensuses.push_back(consensus);
			} else {
				rc_sr_consensuses.push_back(consensus);
			}
		}
	}

	std::string hsr_consensuses_fname = workdir + "/workspace/hsr_consensuses/" + std::to_string(contig_id) + ".txt";
	if (file_exists(hsr_consensuses_fname)) {
		std::ifstream hsr_consensuses_fin(hsr_consensuses_fname);
		std::string line;
		while (std::getline(hsr_consensuses_fin, line)) {
			std::shared_ptr<consensus_t> consensus = std::make_shared<consensus_t>(line, true);
			if (consensus->left_clipped) {
				lc_hsr_consensuses.push_back(consensus);
			} else {
				rc_hsr_consensuses.push_back(consensus);
			}
		}
	}

	// sort by breakpoint position
	auto consensus_cmp = [](std::shared_ptr<consensus_t> c1, std::shared_ptr<consensus_t> c2) {
		return c1->breakpoint < c2->breakpoint;
	};
    std::sort(rc_sr_consensuses.begin(), rc_sr_consensuses.end(), consensus_cmp);
    std::sort(lc_sr_consensuses.begin(), lc_sr_consensuses.end(), consensus_cmp);

	remove_hsr_overlapping_clipped(rc_hsr_consensuses, rc_sr_consensuses);
    remove_hsr_overlapping_clipped(lc_hsr_consensuses, lc_sr_consensuses);
}

void find_indels_from_paired_consensuses(int id, int contig_id, std::string contig_name, 
	std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq) {

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

	mtx.lock();
	std::vector<std::shared_ptr<consensus_t>>& rc_sr_consensuses = rc_sr_consensuses_by_chr[contig_name];
	std::vector<std::shared_ptr<consensus_t>>& lc_sr_consensuses = lc_sr_consensuses_by_chr[contig_name];
	std::vector<std::shared_ptr<consensus_t>>& rc_hsr_consensuses = rc_hsr_consensuses_by_chr[contig_name];
	std::vector<std::shared_ptr<consensus_t>>& lc_hsr_consensuses = lc_hsr_consensuses_by_chr[contig_name];
	mtx.unlock();

	std::vector<std::shared_ptr<consensus_t>> rc_consensuses, lc_consensuses;
	rc_consensuses.insert(rc_consensuses.end(), rc_sr_consensuses.begin(), rc_sr_consensuses.end());
	rc_consensuses.insert(rc_consensuses.end(), rc_hsr_consensuses.begin(), rc_hsr_consensuses.end());
	lc_consensuses.insert(lc_consensuses.end(), lc_sr_consensuses.begin(), lc_sr_consensuses.end());
	lc_consensuses.insert(lc_consensuses.end(), lc_hsr_consensuses.begin(), lc_hsr_consensuses.end());

	find_indels_from_rc_lc_pairs(contig_name, rc_consensuses, lc_consensuses, aligner);

    /* == Deal with unpaired consensuses == */
	mtx.lock();
	std::vector<std::shared_ptr<consensus_t>>& unpaired_consensuses = unpaired_consensuses_by_chr[contig_name];
	unpaired_consensuses.insert(unpaired_consensuses.end(), rc_consensuses.begin(), rc_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), lc_consensuses.begin(), lc_consensuses.end());
	mtx.unlock();
}

void find_indels_from_unpaired_consensuses(int id, std::string contig_name, std::vector<std::shared_ptr<consensus_t>>* consensuses, int start, int end) {

	std::vector<sv_t*> local_svs;

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	char* contig_seq = chr_seqs.get_seq(contig_name);
	hts_pos_t contig_len = chr_seqs.get_len(contig_name);

	for (int i = start; i < consensuses->size() && i < end; i++) {
		std::shared_ptr<consensus_t> consensus = consensuses->at(i);

		std::vector<sv_t*> svs;
		if (!consensus->left_clipped) {
			svs = detect_svs(contig_name, contig_seq, contig_len, consensus, NULL, aligner, 0, config.min_clip_len, 0.0);
		} else if (consensus->left_clipped) {
			svs = detect_svs(contig_name, contig_seq, contig_len, NULL, consensus, aligner, 0, config.min_clip_len, 0.0);
		}

		if (svs.empty()) {
			continue;
		}

		local_svs.insert(local_svs.end(), svs.begin(), svs.end());
	}

	mtx.lock();
	svs_by_chr[contig_name].insert(svs_by_chr[contig_name].end(), local_svs.begin(), local_svs.end());
	mtx.unlock();
}

void cluster_ss_dps(int id, int contig_id, std::string contig_name) {

	std::string dp_fname = workdir + "/workspace/same-strand/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(dp_fname)) return;

	std::unordered_map<std::string, std::string> qname_to_seq;
	std::unordered_map<std::string, int64_t> qname_to_mate_nm; 
	std::ifstream mateseqs_fin(workdir + "/workspace/same-strand/" + std::to_string(contig_id) + ".txt");
	std::string qname, seq;
	int64_t nm;
	while (mateseqs_fin >> qname >> seq >> nm) {
		qname_to_seq[qname] = seq;
		qname_to_mate_nm[qname] = nm;
	}

	std::vector<std::shared_ptr<cluster_t>> ss_clusters;
	open_samFile_t* dp_bam_file = open_samFile(dp_fname, true);
	hts_itr_t* iter = sam_itr_querys(dp_bam_file->idx, dp_bam_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(dp_bam_file->file, iter, read) >= 0) {
		std::string qname = bam_get_qname(read);
		std::shared_ptr<cluster_t> cluster = std::make_shared<cluster_t>(read, config.high_confidence_mapq);
		ss_clusters.push_back(cluster);
	}
	close_samFile(dp_bam_file);

	int min_cluster_size = std::max(3, int(stats.get_median_depth(contig_name)+5)/10);
	int max_cluster_size = (stats.get_max_depth(contig_name) * stats.max_is)/stats.read_len;
	std::vector<std::shared_ptr<bam1_t>> reads;
	if (!ss_clusters.empty()) cluster_clusters(ss_clusters, reads, stats.max_is, max_cluster_size);
	ss_clusters.erase(std::remove_if(ss_clusters.begin(), ss_clusters.end(), [min_cluster_size](std::shared_ptr<cluster_t> c) { return c == NULL || c->count < min_cluster_size; }), ss_clusters.end());

	for (std::shared_ptr<cluster_t> c : ss_clusters) {
		c->ra_furthermost_seq = qname_to_seq[c->ra_furthermost_seq];
	}

	mtx.lock();
	ss_clusters_by_chr[contig_name] = ss_clusters;
	mtx.unlock();
}

int main(int argc, char* argv[]) {

    complete_bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    reference_fname = argv[3];
    std::string sample_name = argv[4];

    std::string full_cmd_fname = workdir + "/call_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    chr_seqs.read_lens_into_map(reference_fname);

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

	mateseqs_w_mapq.resize(contig_map.size());
	active_threads_per_chr = std::vector<int>(contig_map.size());
	mutex_per_chr = std::vector<std::mutex>(contig_map.size());

	std::cout << "Clustering same-strand clusters." << std::endl;

	ctpl::thread_pool cluster_ss_dps_thread_pool(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = cluster_ss_dps_thread_pool.push(cluster_ss_dps, contig_id, contig_name);
		futures.push_back(std::move(future));
	}
	cluster_ss_dps_thread_pool.stop(true);
	for (size_t i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	std::cout << "Extending consensuses." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();

	int block_size = 100;
	ctpl::thread_pool extend_consensuses_thread_pool(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::future<void> future;
		std::string contig_name = contig_map.get_name(contig_id);
		std::vector<std::shared_ptr<consensus_t>>& rc_sr_consensuses = rc_sr_consensuses_by_chr[contig_name];
		for (int i = 0; i < rc_sr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) rc_sr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &rc_sr_consensuses, contig_name, start, end, false);
			futures.push_back(std::move(future));
		}
		std::vector<std::shared_ptr<consensus_t>>& lc_sr_consensuses = lc_sr_consensuses_by_chr[contig_name];
		for (int i = 0; i < lc_sr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) lc_sr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &lc_sr_consensuses, contig_name, start, end, false);
			futures.push_back(std::move(future));
		}
		std::vector<std::shared_ptr<consensus_t>>& rc_hsr_consensuses = rc_hsr_consensuses_by_chr[contig_name];
		for (int i = 0; i < rc_hsr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) rc_hsr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &rc_hsr_consensuses, contig_name, start, end, false);
			futures.push_back(std::move(future));
		}
		std::vector<std::shared_ptr<consensus_t>>& lc_hsr_consensuses = lc_hsr_consensuses_by_chr[contig_name];
		for (int i = 0; i < lc_hsr_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) lc_hsr_consensuses.size());
			future = extend_consensuses_thread_pool.push(extend_consensuses, &lc_hsr_consensuses, contig_name, start, end, false);
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

	chr_seqs.read_fasta_into_map(reference_fname);

    ctpl::thread_pool finding_indels_thread_pool(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = finding_indels_thread_pool.push(find_indels_from_paired_consensuses, contig_id, contig_name, &mateseqs_w_mapq[contig_id]);
        futures.push_back(std::move(future));
    }
    finding_indels_thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Indels found in " << elapsed_time << " seconds" << std::endl;

	std::cout << "Extending unpaired consensuses." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();

	chr_seqs.read_lens_into_map(reference_fname);

	ctpl::thread_pool extend_unpaired_consensuses_thread_pool(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::vector<std::shared_ptr<consensus_t>>& unpaired_consensuses = unpaired_consensuses_by_chr[contig_name];
		for (int i = 0; i <= unpaired_consensuses.size(); i += block_size) {
			int start = i, end = std::min(i+block_size, (int) unpaired_consensuses.size());
			std::future<void> future = extend_unpaired_consensuses_thread_pool.push(extend_consensuses, &unpaired_consensuses, contig_name, start, end, true);
			futures.push_back(std::move(future));
		}
	}
	extend_unpaired_consensuses_thread_pool.stop(true);
	for (size_t i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Unpaired consensuses extended in " << elapsed_time << " seconds" << std::endl;


	std::cout << "Finding indels from unpaired consensuses." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();

	chr_seqs.read_fasta_into_map(reference_fname);

    ctpl::thread_pool finding_indels_from_up_consenensus_thread_pool(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
    	std::vector<std::shared_ptr<consensus_t>>& consensuses = unpaired_consensuses_by_chr[contig_name];
    	if (consensuses.empty()) continue;
    	for (int i = 0; i <= consensuses.size()/block_size; i++) {
    		int start = i * block_size;
    		int end = std::min(start+block_size, (int) consensuses.size());
			std::future<void> future = finding_indels_from_up_consenensus_thread_pool.push(find_indels_from_unpaired_consensuses,
					contig_name, &consensuses, i*block_size, end);
			futures.push_back(std::move(future));
		}
    }
    finding_indels_from_up_consenensus_thread_pool.stop(true);
	for (size_t i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Indels found from unpaired consensuses in " << elapsed_time << " seconds" << std::endl;

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	StripedSmithWaterman::Filter filter;
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);

		auto& svs = svs_by_chr[contig_name];
		for (int i = 0; i < svs.size(); i++) {
			sv_t* sv = svs[i];
			if (!sv->imprecise && sv->end-sv->start >= config.min_sv_size && abs(sv->svlen()) < config.min_sv_size && sv->ins_seq.length() >= config.min_sv_size
			&& (sv->svtype() == "DEL" || sv->svtype() == "INS") && double(sv->ins_seq.length())/(sv->end-sv->start) >= 0.75) {
				char* contig_seq = chr_seqs.get_seq(contig_name);
				std::string inv_seq = svs[i]->ins_seq;
				rc(inv_seq);
				StripedSmithWaterman::Alignment aln;
				aligner.Align(inv_seq.c_str(), contig_seq+svs[i]->start, svs[i]->end-svs[i]->start, filter, &aln, 0);
				if (aln.sw_score >= 0.9*inv_seq.length()) {
					std::string ins_seq = sv->ins_seq;
					if (sv->svlen() == 0 && aln.sw_score == inv_seq.length()) ins_seq = "";
					inversion_t* inv = new inversion_t(contig_name, sv->start, sv->end, ins_seq, sv->rc_consensus, sv->lc_consensus, sv->left_anchor_aln, sv->left_anchor_aln, sv->right_anchor_aln, sv->right_anchor_aln);
					inv->source = sv->source;
					delete sv;
					svs[i] = inv;
				}
			} else if (svs[i]->svtype() == "DEL" && svs[i]->end-svs[i]->start < config.min_sv_size) {
				delete svs[i];
				svs[i] = NULL;
			} else if ((svs[i]->svtype() == "DUP" && svs[i]->svlen() < config.min_sv_size) || 
					   (svs[i]->svtype() == "INS" && svs[i]->ins_seq.length() < config.min_sv_size)) {
				delete svs[i];
				svs[i] = NULL;
			}
		}
		svs.erase(std::remove(svs.begin(), svs.end(), (sv_t*) NULL), svs.end());
	}

    // create VCF out files
	bcf_hdr_t* out_vcf_header = generate_vcf_header(chr_seqs, sample_name, config, full_cmd_str);
    std::string out_vcf_fname = workdir + "/intermediate_results/sr.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	std::string out_dp_inv_vcf_fname = workdir + "/intermediate_results/dp_invs.vcf.gz";
	htsFile* out_dp_inv_vcf_file = bcf_open(out_dp_inv_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_dp_inv_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_dp_inv_vcf_fname + ".");
	}

    bcf1_t* bcf_entry = bcf_init();
	std::unordered_map<std::string, int> svtype_id;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
		std::vector<sv_t*>& svs = svs_by_chr[contig_name];
		std::sort(svs.begin(), svs.end(), [](sv_t* sv1, sv_t* sv2) {
			return std::tie(sv1->start, sv1->end) < std::tie(sv2->start, sv2->end);
		});

        for (sv_t* sv : svs) {
			sv->id = sv->svtype() +  "_SR_" + std::to_string(svtype_id[sv->svtype()]++);

			// do some light filtering here - it helps merge_identical_calls not merge good calls with calls that will get filtered
			if (sv->svtype() == "DEL") {
				if (sv->remap_boundary_lower() > sv->start) {
					sv->sample_info.filters.push_back("REMAP_BOUNDARY_FILTER");
				} else if (sv->remap_boundary_upper() < sv->end) {
					sv->sample_info.filters.push_back("REMAP_BOUNDARY_FILTER");
				}
			} else if (sv->svtype() == "DUP") {
				if (sv->start > sv->remap_boundary_upper()) {
					sv->sample_info.filters.push_back("REMAP_BOUNDARY_FILTER");
				} else if (sv->end < sv->remap_boundary_lower()) {
					sv->sample_info.filters.push_back("REMAP_BOUNDARY_FILTER");
				}
			}

			sv->mh_len = 0; // current value is just a temporary approximation, reset it. genotype will calculate it correctly

			if (sv->source == "1SR_RC" || sv->source == "1HSR_RC") {
            	if (sv->rc_consensus->right_ext_reads < 3) sv->sample_info.filters.push_back("FAILED_TO_EXTEND");
            } else if (sv->source == "1SR_LC" || sv->source == "1HSR_LC") {
            	if (sv->lc_consensus->left_ext_reads < 3) sv->sample_info.filters.push_back("FAILED_TO_EXTEND");
            }

			sv2bcf(out_vcf_header, bcf_entry, sv, chr_seqs.get_seq(contig_name));
            if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
			}
            delete sv;
        }
    }

	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		auto dp_invs = dp_bnds_by_chr[contig_name];
		std::sort(dp_invs.begin(), dp_invs.end(), [](breakend_t* inv1, breakend_t* inv2) {
			return std::tie(inv1->start, inv1->end) < std::tie(inv2->start, inv2->end);
		});

		for (breakend_t* bnd : dp_invs) {
			bnd->id = "BND_DP_" + std::to_string(svtype_id["BND"]++);
			sv2bcf(out_vcf_header, bcf_entry, bnd, chr_seqs.get_seq(contig_name));
			if (bcf_write(out_dp_inv_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
			}
			delete bnd;
		}
	}

	bcf_close(out_dp_inv_vcf_file);

    chr_seqs.clear();

    bcf_hdr_destroy(out_vcf_header);
    bcf_close(out_vcf_file);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
