#ifndef HSR_H_
#define HSR_H_

#include <unordered_set>
#include <htslib/sam.h>

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "utils.h"
#include "sam_utils.h"
#include "remapping.h"
#include "extend_1sr_consensus.h"

extern config_t config;

const int EXTRA_SEQ = 10;


indel_t* remap_rc_cluster(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree, std::string contig_name,
		char* contig_seq, hts_pos_t contig_len, StripedSmithWaterman::Aligner& aligner,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	hts_pos_t left_ext_target_start = consensus->left_ext_target_start(config);
	hts_pos_t left_ext_target_end = consensus->left_ext_target_end(config);

	hts_pos_t right_ext_target_start = consensus->right_ext_target_start(config);
	hts_pos_t right_ext_target_end = consensus->right_ext_target_end(config);

	int old_consensus_len = consensus->consensus.length();
	extend_consensus_to_right(consensus, candidate_reads_itree, right_ext_target_start, right_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int r_ext_len = consensus->consensus.length() - old_consensus_len;

	old_consensus_len = consensus->consensus.length();
	extend_consensus_to_left(consensus, candidate_reads_itree, left_ext_target_start, left_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int l_ext_len = consensus->consensus.length() - old_consensus_len;

	hts_pos_t ref_start = std::max(hts_pos_t(0), consensus->start - EXTRA_SEQ);
	hts_pos_t ref_end = std::min(consensus->end + config.max_is, contig_len);
	hts_pos_t ref_len = ref_end - ref_start;

	hts_pos_t remap_target_end = consensus->remap_boundary;
	hts_pos_t remap_target_start = consensus->remap_boundary - config.max_is;
	if (consensus->remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
		remap_target_start = consensus->breakpoint - config.max_is - 2*consensus->consensus.length();
		remap_target_end = consensus->breakpoint + config.max_is + 2*consensus->consensus.length();
	}
	remap_target_start = std::max(remap_target_start-l_ext_len, hts_pos_t(0));
	remap_target_end = std::min(remap_target_end+r_ext_len, contig_len);
	hts_pos_t remap_target_len = remap_target_end - remap_target_start;

	if (remap_target_start >= remap_target_end ||
		has_Ns(contig_seq, remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

    indel_t* indel = remap_consensus(contig_name, consensus->consensus, contig_seq, contig_len, ref_start, ref_len, remap_target_start,
    		remap_target_len, aligner, NULL, consensus, "1HSR_RC");
    return indel;
}

indel_t* remap_lc_cluster(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree, std::string contig_name,
		char* contig_seq, hts_pos_t contig_len, StripedSmithWaterman::Aligner& aligner,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	hts_pos_t left_ext_target_start = consensus->left_ext_target_start(config);
	hts_pos_t left_ext_target_end = consensus->left_ext_target_end(config);

	hts_pos_t right_ext_target_start = consensus->right_ext_target_start(config);
	hts_pos_t right_ext_target_end = consensus->right_ext_target_end(config);

	int old_consensus_len = consensus->consensus.length();
	extend_consensus_to_left(consensus, candidate_reads_itree, left_ext_target_start, left_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int l_ext_len = consensus->consensus.length() - old_consensus_len;

	old_consensus_len = consensus->consensus.length();
	extend_consensus_to_right(consensus, candidate_reads_itree, right_ext_target_start, right_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int r_ext_len = consensus->consensus.length() - old_consensus_len;

	hts_pos_t ref_end = std::min(consensus->end + EXTRA_SEQ, contig_len);
	hts_pos_t ref_start = std::max(hts_pos_t(0), consensus->start - config.max_is);
	hts_pos_t ref_len = ref_end - ref_start;

	hts_pos_t remap_target_start = consensus->remap_boundary;
	hts_pos_t remap_target_end = consensus->remap_boundary + config.max_is;
	if (consensus->remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
		remap_target_start = consensus->breakpoint - config.max_is - 2*consensus->consensus.length();
		remap_target_end = consensus->breakpoint + config.max_is + 2*consensus->consensus.length();
	}
	remap_target_start = std::max(remap_target_start-l_ext_len, hts_pos_t(0));
	remap_target_end = std::min(remap_target_end+r_ext_len, contig_len);
	hts_pos_t remap_target_len = remap_target_end - remap_target_start;

	// do not attempt if reference region has Ns - this is because of in our aligner, Ns will always match
	if (remap_target_start >= remap_target_end ||
		has_Ns(contig_seq, remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

	indel_t* indel = remap_consensus(contig_name, consensus->consensus, contig_seq, contig_len, remap_target_start, remap_target_len,
			ref_start, ref_len, aligner, consensus, NULL, "1HSR_LC");
    return indel;
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

#endif /* HSR_H_ */
