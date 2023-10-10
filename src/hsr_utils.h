#ifndef HSR_UTILS_H
#define HSR_UTILS_H

#include <algorithm>
#include <htslib/sam.h>

#include "../libs/ssw_cpp.h"
#include "../libs/IntervalTree.h"
#include "utils.h"

int compute_left_half_Ms(bam1_t* r) {
	int border = r->core.l_qseq/2;
	int left_Ms = 0;

	// compute how many Ms (either matches or mismatches are in the first half of the read)
	// in other words, this is readlen/2 - number of insertions in the first half
	int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
	for (int i = 0; i < r->core.n_cigar; i++) {
		char op_char = bam_cigar_opchr(cigar[i]);
		int op_len = bam_cigar_oplen(cigar[i]);
		if (op_char == 'M') {
			left_Ms += std::min(border, qpos+op_len) - qpos;
			qpos += op_len;
		} else if (op_char == 'I') {
			qpos += op_len;
		}

		if (qpos >= border) break;
	}
	return left_Ms;
}

std::pair<int, int> compute_left_and_right_differences_indel_as_1_diff(bam1_t* r) {
	int border = r->core.l_qseq/2;
	int left_Ms = compute_left_half_Ms(r);

	// we computed left_Ms because the MD tag does not include insertions
	// i.e. if I have CIGAR: 50M50I50M and MD: 74A25, the mismatch is not in position 75 in the query,
	// but it is in position 125
	// here we compute the number of deleted bases and mismatches in each half of the read based on the MD tag
	std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
	int m = 0;
	int left_diffs = 0, right_diffs = 0;
	bool del_mode = false;
	for (int i = 0; i < md_tag.length(); i++) {
		char c = md_tag[i];
		if (c >= '0' && c <= '9') {
			m = m*10 + c-'0';
			del_mode = false;
		} else if (c >= 'A' && c <= 'T') {
			left_Ms -= m;
			m = 0;
			if (!del_mode) {
				left_Ms--;
				if (left_Ms > 0) left_diffs++;
				else right_diffs++;
			}
		} else if (c == '^') {
			del_mode = true;
			if (left_Ms > 0) left_diffs++;
			else right_diffs++;
		}
	}

	int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
	for (int i = 0; i < r->core.n_cigar; i++) {
		char op_char = bam_cigar_opchr(cigar[i]);
		int op_len = bam_cigar_oplen(cigar[i]);
		if (op_char == 'M') {
			qpos += op_len;
		} else if (op_char == 'I') {
			if (qpos < border && qpos + op_len > border) {  // this is the case there an insertion is partially in the left
															// half and partially in the right half
				left_diffs++;
				right_diffs++;
			} else if (qpos < border) {
				left_diffs++;
			} else if (qpos >= border) {
				right_diffs++;
			}
			qpos += op_len;
		}
	}

	return {left_diffs, right_diffs};
}

// computes the differences for the left and the right half of the read
std::pair<int, int> compute_left_and_right_differences_indel_as_n_diffs(bam1_t* r) {
    int border = r->core.l_qseq/2;
    int left_Ms = compute_left_half_Ms(r);

    // we computed left_Ms because the MD tag does not include insertions
    // i.e. if I have CIGAR: 50M50I50M and MD: 74A25, the mismatch is not in position 75 in the query,
    // but it is in position 125
    // here we compute the number of deleted bases and mismatches in each half of the read based on the MD tag
    std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
    int m = 0;
    int left_diffs = 0, right_diffs = 0;
    bool del_mode = false;
    for (int i = 0; i < md_tag.length(); i++) {
        char c = md_tag[i];
        if (c >= '0' && c <= '9') {
            m = m*10 + c-'0';
            del_mode = false;
        } else if (c >= 'A' && c <= 'T') {
            left_Ms -= m;
            m = 0;
            if (left_Ms > 0) left_diffs++;
            else right_diffs++;
            if (!del_mode) left_Ms--;
        } else if (c == '^') {
            del_mode = true;
        }
    }

    // We also need to count the number of inserted bases, since the MD tag does not report this information
    int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
    for (int i = 0; i < r->core.n_cigar; i++) {
        char op_char = bam_cigar_opchr(cigar[i]);
        int op_len = bam_cigar_oplen(cigar[i]);
        if (op_char == 'M') {
            qpos += op_len;
        } else if (op_char == 'I') {
            if (qpos < border && qpos + op_len > border) {  // this is the case there an insertion is partially in the left
                                                            // half and partially in the right half
                left_diffs += border - qpos;
                right_diffs += qpos + op_len - border;
            } else if (qpos < border) {
                left_diffs += op_len;
            } else if (qpos >= border) {
                right_diffs += op_len;
            }
            qpos += op_len;
        }
    }

    return {left_diffs, right_diffs};
}

std::pair<int, int> compute_left_and_right_differences(bam1_t* r, bool indel_as_single_diff) {
	if (indel_as_single_diff) {
		return compute_left_and_right_differences_indel_as_1_diff(r);
	} else {
		return compute_left_and_right_differences_indel_as_n_diffs(r);
	}
}

void filter_well_aligned_to_ref(char* contig_seq, hts_pos_t contig_len, std::vector<consensus_t*>& consensuses, config_t config) {
    std::vector<consensus_t*> retained;

	StripedSmithWaterman::Aligner aligner(2, 2, 4, 1, true);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    for (consensus_t* c : consensuses) {
        hts_pos_t ref_start = c->start - 10, ref_end = c->end + 10;
		if (ref_start < 0) ref_start = 0;
		if (ref_end > contig_len) ref_end = contig_len;
        aligner.Align(c->consensus.c_str(), contig_seq+ref_start, ref_end-ref_start, filter, &aln, 0);
        int differences = std::count(aln.cigar_string.begin(), aln.cigar_string.end(), 'D') +
                          std::count(aln.cigar_string.begin(), aln.cigar_string.end(), 'I') + 
						  aln.mismatches;
        if ((aln.sw_score < c->consensus.length()*2 && differences >= config.min_diff_hsr) || is_clipped(aln)) {
            retained.push_back(c);
        } else {
        	delete c;
        }
    }
    consensuses.swap(retained);
}

void select_nonoverlapping_clusters(std::vector<consensus_t*>& consensuses) {
	// for overlapping pairs, keep the ones with higher count
	std::vector<consensus_t*> to_be_deleted;
	std::vector<Interval<consensus_t*> > rc_iv;
	for (consensus_t* c : consensuses) rc_iv.push_back(Interval<consensus_t*>(c->start, c->end, c));
	IntervalTree<consensus_t*> rc_it(rc_iv);
	std::vector<consensus_t*> kept_consensuses;
	for (consensus_t* c : consensuses) {
		std::vector<Interval<consensus_t*> > ov = rc_it.findOverlapping(c->start, c->end);
		bool keep_consensus = true;
		for (Interval<consensus_t*> ov_c : ov) {
			if (ov_c.value->supp_clipped_reads() > c->supp_clipped_reads()) { // a higher count was found
				keep_consensus = false;
				break;
			}
		}

		if (keep_consensus) {
			kept_consensuses.push_back(c);
		} else {
			to_be_deleted.push_back(c);
		}
	}
	for (consensus_t* c : to_be_deleted) delete c;
	consensuses.swap(kept_consensuses);
}

#endif
