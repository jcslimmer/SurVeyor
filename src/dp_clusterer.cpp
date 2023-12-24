#include <fstream>
#include <string>
#include <mutex>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "../libs/cptl_stl.h"
#include "../libs/IntervalTree.h"
#include "../libs/kdtree.h"
#include "sam_utils.h"
#include "sw_utils.h"
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
	std::sort(clusters.begin(), clusters.end(), [](cluster_t* c1, cluster_t* c2) {
		return std::make_tuple(c1->la_start, c1->la_end, c1->ra_start, c1->ra_end) < std::make_tuple(c2->la_start, c2->la_end, c2->ra_start, c2->ra_end);
	});

	std::vector<cluster_t*> dedup_clusters;
	dedup_clusters.push_back(clusters[0]);
	for (int i = 1; i < clusters.size(); i++) {
		if (*clusters[i-1] != *clusters[i]) {
			dedup_clusters.push_back(clusters[i]);
		}
	}
	clusters.swap(dedup_clusters);

	std::priority_queue<cc_pair> pq;
	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i]->used) continue;

		std::vector<cc_pair> ccps;
		for (int j = i+1; j < clusters.size(); j++) {
			if (clusters[j]->la_start-clusters[i]->la_start > stats.max_is || ccps.size() > max_cluster_size) break;
			cc_pair ccp = cc_pair(clusters[i], clusters[j]);
			if (ccp.compatible()) ccps.push_back(ccp);
		}
		if (ccps.size() <= max_cluster_size) {
			for (cc_pair& ccp : ccps) pq.push(ccp);
		} else { // very large cluster - probably abnormal region - mark all involved pairs to be discarded
			clusters[i]->used = true;
			for (cc_pair& ccp : ccps) {
				ccp.c2->used = true;
			}
		}
	}

	kdtree* kd_tree_endpoints = kd_create(2);
	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i]->used) continue;
		double p[2] = {double(clusters[i]->la_end), double(clusters[i]->ra_end)};
		kd_insert(kd_tree_endpoints, p, clusters[i]);
	}

	while (!pq.empty()) {
		cc_pair ccp = pq.top();
		pq.pop();

		if (ccp.c1->used || ccp.c2->used) continue;
		ccp.c1->used = ccp.c2->used = true;

		cluster_t* merged = cluster_t::merge(ccp.c1, ccp.c2);
		double ps[2] = {double(merged->la_start), double(merged->ra_start)};
		kdres* res = kd_nearest_range(kd_tree_endpoints, ps, 2*stats.max_is);
		while (!kd_res_end(res)) {
			cluster_t* c = (cluster_t*) kd_res_item_data(res);
			if (!c->used) {
				cc_pair ccp = cc_pair(merged, c);
				if (ccp.compatible()) pq.push(ccp);
			}
			kd_res_next(res);
		}
		kd_res_free(res);

		double pe[2] = {double(merged->la_end), double(merged->ra_end)};
		kd_insert(kd_tree_endpoints, pe, merged);
		clusters.push_back(merged);
	}
	kd_free(kd_tree_endpoints);

	clusters.erase(std::remove_if(clusters.begin(), clusters.end(), [min_cluster_size](cluster_t* c) {
		return c->count < min_cluster_size;
	}), clusters.end());

	std::sort(clusters.begin(), clusters.end());
}

void cluster_dps(int id, int contig_id, std::string contig_name, std::string bam_fname,
		std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq) {
	
	std::string dp_fname = workdir + "/workspace/long-pairs/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(dp_fname)) return;

	open_samFile_t* dp_bam_file = open_samFile(dp_fname, true);

	std::vector<cluster_t*> lp_clusters, ow_clusters;

	hts_itr_t* iter = sam_itr_querys(dp_bam_file->idx, dp_bam_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(dp_bam_file->file, iter, read) >= 0) {
		cluster_t* cluster = new cluster_t(read, config.high_confidence_mapq);
		if (read->core.isize > 0) {
			lp_clusters.push_back(cluster);
		} else {
			ow_clusters.push_back(cluster);
		}
	}
	close_samFile(dp_bam_file);

	int min_cluster_size = std::max(3, int(stats.get_median_depth(contig_name)+5)/10);
	int max_cluster_size = (stats.get_max_depth(contig_name) * stats.max_is)/stats.read_len;

	if (!lp_clusters.empty()) cluster_clusters(lp_clusters, min_cluster_size, max_cluster_size);
	if (!ow_clusters.empty()) cluster_clusters(ow_clusters, min_cluster_size, max_cluster_size);

	// set leftmost_rseq
	std::unordered_map<std::string, cluster_t*> mateseqs_to_retrieve;
	for (cluster_t* c : lp_clusters) {
		mateseqs_to_retrieve[c->leftmost_rseq] = c;
	}

	std::ifstream mateseqs_fin(workdir + "/workspace/sc_mateseqs/" + std::to_string(contig_id) + ".txt");
	std::string qname, seq;
	while (mateseqs_fin >> qname >> seq) {
		if (!mateseqs_to_retrieve.count(qname)) continue;
		mateseqs_to_retrieve[qname]->leftmost_rseq = seq;
	}
	mateseqs_to_retrieve.clear();

	// find deletions from clusters
	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	StripedSmithWaterman::Filter filter;
	
	std::vector<deletion_t*> deletions;
	for (cluster_t* c : lp_clusters) {
		if (c->used) continue;

		consensus_t* rc_consensus = new consensus_t(false, 0, c->la_end, 0, c->rightmost_lseq, 0, 0, 0, 0, 0, 0);
		consensus_t* lc_consensus = new consensus_t(false, 0, c->ra_start, 0, c->leftmost_rseq, 0, 0, 0, 0, 0, 0);
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
		}

		if (del && -del->svlen() >= config.min_sv_size) {
			del->disc_pairs = c->count;
			del->disc_pairs_maxmapq = c->max_mapq;
			del->disc_pairs_high_mapq = c->confident_count;
			del->source = "DP";
			deletions.push_back(del);
		}
	}
	
	maps_mtx.lock();
	deletions_by_chr[contig_name] = deletions;
	maps_mtx.unlock();

	for (cluster_t* c : lp_clusters) {
		delete c;
	}
	for (cluster_t* c : ow_clusters) {
		delete c;
	}
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
		std::vector<Interval<sv_t*> > iv = del_tree.findContained(del->left_anchor_aln->start, del->right_anchor_aln->end);

		if (iv.size() == 1) {
			sv_t* corr_sr_del = iv[0].value;
			if (corr_sr_del->start-del->left_anchor_aln->start <= stats.max_is && del->right_anchor_aln->end-corr_sr_del->end <= stats.max_is) {
				corr_sr_del->disc_pairs = del->disc_pairs;
				corr_sr_del->disc_pairs_maxmapq = del->disc_pairs_maxmapq;
				corr_sr_del->disc_pairs_high_mapq = del->disc_pairs_high_mapq;
				deletions[i] = NULL;
			}
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

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");

	stats.parse(workdir + "/stats.txt", config.per_contig_stats);

	chr_seqs.read_fasta_into_map(reference_fname);

	ctpl::thread_pool thread_pool1(config.threads);
	std::vector<std::future<void> > futures;
	std::vector<std::unordered_map<std::string, std::pair<std::string, int> > > mateseqs_w_mapq(contig_map.size());
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);

		std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
		std::ifstream fin(fname);
		std::string qname, read_seq, qual; int mapq;
		while (fin >> qname >> read_seq >> qual >> mapq) {
			mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
		}

		std::future<void> future = thread_pool1.push(cluster_dps, contig_id, contig_name, bam_fname, &mateseqs_w_mapq[contig_id]);
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


	std::string sr_vcf_fname = workdir + "/sr.dedup.vcf.gz";
	htsFile* sr_vcf_file = bcf_open(sr_vcf_fname.c_str(), "r");
	bcf_hdr_t* hdr = bcf_hdr_read(sr_vcf_file);

	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(sr_vcf_file, hdr, bcf_entry) == 0) {
		sv_t* sv = bcf_to_sv(hdr, bcf_entry);
		if (sv->rc_consensus && sv->lc_consensus && sv->svtype() == "DEL") { // only use 2SR/2HSR
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

	std::string dp_vcf_fname = workdir + "/dp.vcf.gz";
	htsFile* dp_vcf_file = bcf_open(dp_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(dp_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + dp_vcf_fname + ".");
	}

	int del_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
		auto& deletions = deletions_by_chr[contig_name];
		std::sort(deletions.begin(), deletions.end(), [](const sv_t* sv1, const sv_t* sv2) {return sv1->start < sv2->start;});
		for (sv_t* del : deletions) {
			del->id = "DEL_DP_" + std::to_string(del_id++);
			sv2bcf(hdr, bcf_entry, del, chr_seqs.get_seq(contig_name));
			if (bcf_write(dp_vcf_file, hdr, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + dp_vcf_fname + ".");
			}
		}
    }

	bcf_close(dp_vcf_file);

	tbx_index_build(dp_vcf_fname.c_str(), 0, &tbx_conf_vcf);

	std::string merged_vcf_fname = workdir + "/survindel2.out.vcf.gz";
	htsFile* merged_vcf_file = bcf_open(merged_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(merged_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + merged_vcf_fname + ".");
	}

	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		std::vector<sv_t*> all_entries;
		all_entries.insert(all_entries.end(), deletions_by_chr[contig_name].begin(), deletions_by_chr[contig_name].end());
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
