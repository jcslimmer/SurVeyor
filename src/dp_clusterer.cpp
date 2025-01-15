#include <fstream>
#include <string>
#include <mutex>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "../libs/cptl_stl.h"
#include "../libs/IntervalTree.h"
#include "sam_utils.h"
#include "sw_utils.h"
#include "types.h"
#include "utils.h"
#include "vcf_utils.h"
#include "clustering_utils.h"

std::string workdir;

config_t config;
stats_t stats;
chr_seqs_map_t chr_seqs;

std::unordered_map<std::string, std::vector<deletion_t*>> deletions_by_chr;
std::unordered_map<std::string, std::vector<sv_t*>> sr_entries_by_chr, sr_nonpass_entries_by_chr;
std::mutex maps_mtx;

struct cc_pair {
	cluster_t* c1, * c2;
	hts_pos_t dist;

	cc_pair(cluster_t* c1, cluster_t* c2) : c1(c1), c2(c2), dist(cluster_t::distance(c1, c2)) {}

	bool compatible() {
		return dist <= stats.max_is;
	}
};
bool operator < (const cc_pair& cc1, const cc_pair& cc2) { return cc1.dist > cc2.dist; }

void cluster_clusters(std::vector<cluster_t*>& clusters, int min_cluster_size, int max_cluster_size) {
	std::vector<bam1_t*> reads;
	cluster_clusters(clusters, reads, stats.max_is, max_cluster_size, false);

	clusters.erase(std::remove_if(clusters.begin(), clusters.end(), [min_cluster_size](cluster_t* c) {
		return c == NULL || c->count < min_cluster_size;
	}), clusters.end());

	std::sort(clusters.begin(), clusters.end());
}

void cluster_lp_dps(int contig_id, std::string contig_name, std::vector<deletion_t*>& deletions) {

	std::string dp_fname = workdir + "/workspace/long-pairs/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(dp_fname)) return;
	
	std::unordered_map<std::string, int64_t> qname_to_mate_nm; 
	std::ifstream mateseqs_fin(workdir + "/workspace/long-pairs/" + std::to_string(contig_id) + ".txt");
	std::string qname, seq;
	int64_t nm;
	while (mateseqs_fin >> qname >> seq >> nm) {
		qname_to_mate_nm[qname] = nm;
	}

	std::vector<cluster_t*> lp_clusters;
	open_samFile_t* dp_bam_file = open_samFile(dp_fname, true);
	hts_itr_t* iter = sam_itr_querys(dp_bam_file->idx, dp_bam_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(dp_bam_file->file, iter, read) >= 0) {
		std::string qname = bam_get_qname(read);
		cluster_t* cluster = new cluster_t(read, qname_to_mate_nm[qname], config.high_confidence_mapq);
		lp_clusters.push_back(cluster);
	}
	close_samFile(dp_bam_file);
	qname_to_mate_nm.clear();

	int min_cluster_size = std::max(3, int(stats.get_median_depth(contig_name)+5)/10);
	int max_cluster_size = (stats.get_max_depth(contig_name) * stats.max_is)/stats.read_len;
	if (!lp_clusters.empty()) cluster_clusters(lp_clusters, min_cluster_size, max_cluster_size);

	// set leftmost_rseq
	std::unordered_map<std::string, cluster_t*> mateseqs_to_retrieve;
	for (cluster_t* c : lp_clusters) {
		mateseqs_to_retrieve[c->ra_furthermost_seq] = c;
	}

	mateseqs_fin.clear();
	mateseqs_fin.seekg(0);
	while (mateseqs_fin >> qname >> seq >> nm) {
		if (!mateseqs_to_retrieve.count(qname)) continue;
		mateseqs_to_retrieve[qname]->ra_furthermost_seq = seq;
	}
	mateseqs_to_retrieve.clear();

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	for (cluster_t* c : lp_clusters) {
		if (c->used) continue;

		consensus_t* rc_consensus = new consensus_t(false, 0, c->la_end, 0, c->la_furthermost_seq, 0, 0, 0, 0, 0, 0);
		consensus_t* lc_consensus = new consensus_t(false, 0, c->ra_start, 0, c->ra_furthermost_seq, 0, 0, 0, 0, 0, 0);
		std::vector<sv_t*> svs = detect_svs(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), rc_consensus, lc_consensus, aligner, stats.read_len/3, config.min_clip_len, 0.0);

		deletion_t* del = NULL;
		if (!svs.empty() && svs[0]->svtype() == "DEL") { // first check if I can obtain precise coordinates for the deletion
			del = (deletion_t*) svs[0];
			del->original_range = std::to_string(c->la_end) + "-" + std::to_string(c->ra_start);
			del->remapped = true;
		} else if (c->la_end < c->ra_start) {
			sv_t::anchor_aln_t* left_anchor_aln = new sv_t::anchor_aln_t(c->la_start, c->la_end, c->la_end-c->la_start, 0, 0, "");
			sv_t::anchor_aln_t* right_anchor_aln = new sv_t::anchor_aln_t(c->ra_start, c->ra_end, c->ra_end-c->ra_start, 0, 0, "");
			del = new deletion_t(contig_name, c->la_end, c->ra_start, "", NULL, NULL, left_anchor_aln, right_anchor_aln, NULL);
			del->precompute_base_frequencies(chr_seqs.get_seq(contig_name));
			del->imprecise = true;
		}

		if (del && -del->svlen() >= config.min_sv_size) {
			del->sample_info.alt_bp1.supp_pairs = del->sample_info.alt_bp2.supp_pairs = c->count;
			del->disc_pairs_lf_maxmapq = c->la_max_mapq;
			del->disc_pairs_rf_maxmapq = c->ra_max_mapq;
			del->sample_info.alt_bp1.supp_pairs_high_mapq = c->la_confident_count;
			del->sample_info.alt_bp2.supp_pairs_high_mapq = c->ra_confident_count;
			del->disc_pairs_lf_avg_nm = double(c->la_cum_nm)/c->count;
			del->disc_pairs_rf_avg_nm = double(c->ra_cum_nm)/c->count;
			del->source = "DP";
			deletions.push_back(del);
		}

		delete c;
	}
}

void cluster_dps(int id, int contig_id, std::string contig_name, std::string bam_fname) {

	std::vector<deletion_t*> deletions;
	std::vector<inversion_t*> inversions;	
	cluster_lp_dps(contig_id, contig_name, deletions);

	maps_mtx.lock();
	deletions_by_chr[contig_name] = deletions;
	maps_mtx.unlock();
}

void merge_sr_dp(int id, int contig_id, std::string contig_name) {
	maps_mtx.lock();
	std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
	std::vector<sv_t*>& sv_sr_entries = sr_entries_by_chr[contig_name];
	maps_mtx.unlock();

	if (deletions.empty() || sv_sr_entries.empty()) return;

	std::vector<Interval<sv_t*> > del_intervals, dup_intervals;
	for (sv_t* sr_entry : sv_sr_entries) {
		if (sr_entry->svtype() == "DEL") {
			del_intervals.push_back(Interval<sv_t*>(sr_entry->start, sr_entry->end, sr_entry));
		}
	}
	IntervalTree<sv_t*> del_tree = IntervalTree<sv_t*>(del_intervals);

	for (int i = 0; i < deletions.size(); i++) {
		deletion_t* del = deletions[i];
		std::vector<Interval<sv_t*> > ivs = del_tree.findContained(del->left_anchor_aln->start, del->right_anchor_aln->end);

		std::vector<sv_t*> compatible_dels;
		for (Interval<sv_t*> iv : ivs) {
			sv_t* corr_sr_del = iv.value;
			if (corr_sr_del->start-del->left_anchor_aln->start <= stats.max_is && del->right_anchor_aln->end-corr_sr_del->end <= stats.max_is &&
				corr_sr_del->start >= del->start-stats.read_len/2 && corr_sr_del->end <= del->end+stats.read_len/2) {
				compatible_dels.push_back(corr_sr_del);
			}
		}

		if (compatible_dels.size() == 1) {
			sv_t* corr_sr_del = compatible_dels[0];
			corr_sr_del->sample_info.alt_bp1.supp_pairs = del->sample_info.alt_bp1.supp_pairs;
			corr_sr_del->sample_info.alt_bp2.supp_pairs = del->sample_info.alt_bp2.supp_pairs;
			corr_sr_del->disc_pairs_lf_maxmapq = del->disc_pairs_lf_maxmapq;
			corr_sr_del->disc_pairs_rf_maxmapq = del->disc_pairs_rf_maxmapq;
			corr_sr_del->sample_info.alt_bp1.supp_pairs_high_mapq = del->sample_info.alt_bp1.supp_pairs_high_mapq;
			corr_sr_del->sample_info.alt_bp2.supp_pairs_high_mapq = del->sample_info.alt_bp2.supp_pairs_high_mapq;
			corr_sr_del->disc_pairs_lf_avg_nm = del->disc_pairs_lf_avg_nm;
			corr_sr_del->disc_pairs_rf_avg_nm = del->disc_pairs_rf_avg_nm;
			deletions[i] = NULL;
		}
	}
	deletions.erase(std::remove(deletions.begin(), deletions.end(), (deletion_t*) NULL), deletions.end());
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];

    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    std::string reference_fname = argv[3];
    std::string sample_name = argv[4];

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");

	stats.parse(workdir + "/stats.txt", config.per_contig_stats);

	chr_seqs.read_fasta_into_map(reference_fname);

	ctpl::thread_pool thread_pool1(config.threads);
	std::vector<std::future<void> > futures;
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool1.push(cluster_dps, contig_id, contig_name, bam_fname);
		futures.push_back(std::move(future));
	}
	thread_pool1.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}
	futures.clear();

	std::string sr_vcf_fname = workdir + "/intermediate_results/sr.dedup.vcf.gz";
	htsFile* sr_vcf_file = bcf_open(sr_vcf_fname.c_str(), "r");
	bcf_hdr_t* hdr = bcf_hdr_read(sr_vcf_file);

	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(sr_vcf_file, hdr, bcf_entry) == 0) {
		sv_t* sv = bcf_to_sv(hdr, bcf_entry);
		if (sv->svtype() == "DEL" && sv->is_pass()) {
			sr_entries_by_chr[sv->chr].push_back(sv);
		} else {
			sr_nonpass_entries_by_chr[sv->chr].push_back(sv);
		}
	}
	bcf_close(sr_vcf_file);

	// merge SR and DISC
	ctpl::thread_pool thread_pool2(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool2.push(merge_sr_dp, contig_id, contig_name);
		futures.push_back(std::move(future));
	}
	thread_pool2.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}
	futures.clear();

	std::string dp_vcf_fname = workdir + "/intermediate_results/dp.vcf.gz";
	htsFile* dp_vcf_file = bcf_open(dp_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(dp_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + dp_vcf_fname + ".");
	}

	// read dp inversions dp_invs.vcf.gz
	std::string dp_inv_vcf_fname = workdir + "/intermediate_results/dp_invs.vcf.gz";
	htsFile* dp_inv_vcf_file = bcf_open(dp_inv_vcf_fname.c_str(), "r");
	bcf_hdr_t* hdr_inv = bcf_hdr_read(dp_inv_vcf_file);

	std::unordered_map<std::string, std::vector<sv_t*>> dp_svs_by_chr;
	while (bcf_read(dp_inv_vcf_file, hdr_inv, bcf_entry) == 0) {
		sv_t* sv = bcf_to_sv(hdr_inv, bcf_entry);
		dp_svs_by_chr[sv->chr].push_back(sv);
	}

	int del_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
		auto& deletions = deletions_by_chr[contig_name];
		for (deletion_t* del : deletions) del->id = "DEL_DP_" + std::to_string(del_id++);

		dp_svs_by_chr[contig_name].insert(dp_svs_by_chr[contig_name].end(), deletions.begin(), deletions.end());
		std::sort(dp_svs_by_chr[contig_name].begin(), dp_svs_by_chr[contig_name].end(), [](sv_t* sv1, sv_t* sv2) {
			return sv1->start < sv2->start;
		});
		for (sv_t* dp_sv : dp_svs_by_chr[contig_name]) {
			sv2bcf(hdr, bcf_entry, dp_sv, chr_seqs.get_seq(contig_name));
			if (bcf_write(dp_vcf_file, hdr, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + dp_vcf_fname + ".");
			}
		}
    }

	bcf_close(dp_vcf_file);

	tbx_index_build(dp_vcf_fname.c_str(), 0, &tbx_conf_vcf);

	std::string merged_vcf_fname = workdir + "/intermediate_results/sr_dp.vcf.gz";
	htsFile* merged_vcf_file = bcf_open(merged_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(merged_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + merged_vcf_fname + ".");
	}

	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		std::vector<sv_t*> all_entries;
		all_entries.insert(all_entries.end(), dp_svs_by_chr[contig_name].begin(), dp_svs_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), sr_entries_by_chr[contig_name].begin(), sr_entries_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), sr_nonpass_entries_by_chr[contig_name].begin(), sr_nonpass_entries_by_chr[contig_name].end());
		std::sort(all_entries.begin(), all_entries.end(), [](sv_t* sv1, sv_t* sv2) {
			return sv1->start < sv2->start;
		});
		for (sv_t* sv : all_entries) {
			sv2bcf(hdr, bcf_entry, sv, chr_seqs.get_seq(contig_name));
			if (bcf_write(merged_vcf_file, hdr, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + merged_vcf_fname + ".");
			}
		}
	}

	bcf_close(merged_vcf_file);

	tbx_index_build(merged_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
