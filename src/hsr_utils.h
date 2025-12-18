#ifndef HSR_UTILS_H
#define HSR_UTILS_H

#include <algorithm>
#include <unordered_set>
#include <htslib/sam.h>

#include "../libs/ssw_cpp.h"
#include "../libs/IntervalTree.h"
#include "utils.h"
#include "sw_utils.h"
#include "SegTree.h"

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
		if (op_char == 'M' || op_char == '=' || op_char == 'X') {
			left_Ms += std::min(border, qpos+op_len) - qpos;
			qpos += op_len;
		} else if (op_char == 'I' || op_char == 'S') {
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
			left_Ms -= m;
			m = 0;
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
			left_Ms -= m;
			m = 0;
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

// If the high quality (i.e., supported by >= 3 reads) part of the consensus sequence aligns too well to the reference,
// it probably means that the differences were due to sequencing errors in the individual reads and they corrected each other
// during consensus building. We filter out such consensuses.
void filter_well_aligned_to_ref(char* contig_seq, hts_pos_t contig_len, std::vector<consensus_t*>& consensuses, config_t config) {
    std::vector<consensus_t*> retained;

	StripedSmithWaterman::Aligner aligner(2, 2, 4, 1, true);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    for (consensus_t* c : consensuses) {
        hts_pos_t ref_start = c->start - 10, ref_end = c->end + 10;
		if (ref_start < 0) ref_start = 0;
		if (ref_end > contig_len) ref_end = contig_len;
        aligner.Align(c->sequence.c_str(), contig_seq+ref_start, ref_end-ref_start, filter, &aln, 0);
		int dops = 0, dbases = 0, iops = 0, ibases = 0, xbases = 0;
		for (uint32_t op : aln.cigar) {
			char opchar = cigar_int_to_op(op);
			int oplen = cigar_int_to_len(op);
			if (opchar == 'D') {
				dops++;
				dbases += oplen;
			} else if (opchar == 'I') {
				iops++;
				ibases += oplen;
			} else if (opchar == 'X') xbases += oplen;
		}
        if (dops + iops + xbases >= config.min_diff_hsr || abs(ibases-dbases) >= config.min_sv_size || is_clipped(aln)) {
			retained.push_back(c);
        } else {
        	delete c;
        }
    }
    consensuses.swap(retained);
}

// Remove consensues that are completely contained within another consensues and their sequence is a substring of the other consensus
void filter_fully_contained(std::vector<consensus_t*>& consensuses) {
	std::vector<consensus_t*> sorted = consensuses;
	std::sort(sorted.begin(), sorted.end(), [](const consensus_t* c1, const consensus_t* c2) {
		return c1->start < c2->start || (c1->start == c2->start && c1->end > c2->end);
	});
	
	std::vector<bool> to_delete(sorted.size(), false);
	for (int i = 0; i < sorted.size(); i++) {
		if (to_delete[i]) continue;
		for (int j = i+1; j < sorted.size(); j++) {
			if (to_delete[j]) continue;
			if (sorted[j]->start > sorted[i]->end) break;
			if (sorted[j]->end <= sorted[i]->end) {
				// j is contained within i
				std::string highq_seq = sorted[j]->sequence.substr(sorted[j]->lowq_prefix, sorted[j]->sequence.length()-sorted[j]->lowq_prefix-sorted[j]->lowq_suffix);
				if (sorted[i]->sequence.find(highq_seq) != std::string::npos) {
					to_delete[j] = true;
				}
			}
		}
	}
	consensuses.clear();
	for (int i = 0; i < sorted.size(); i++) {
		if (!to_delete[i]) consensuses.push_back(sorted[i]);
	}
}

// c1 is assumed to be to the left of c2, and neither cluster is not fully contained within the other
void merge_overlapping_pair_of_clusters(consensus_t* c1, consensus_t* c2, consensus_t* target, std::string& c1_seq, std::string& c2_seq, int overlap) {
	target->breakpoint = c1->left_clipped ? c1->breakpoint : c2->breakpoint;
	target->start = c1->start;
	target->end = c2->end;
	target->sequence = c1_seq + c2_seq.substr(overlap);
	target->fwd_reads = c1->fwd_reads + c2->fwd_reads;
	target->rev_reads = c1->rev_reads + c2->rev_reads;
	target->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
	target->max_mapq = std::max(c1->max_mapq, c2->max_mapq);
	if (c1->left_clipped) {
		target->remap_boundary = std::max(c1->remap_boundary, c2->remap_boundary);
	} else {
		target->remap_boundary = std::min(c1->remap_boundary, c2->remap_boundary);
	}
	target->lowq_prefix = c1->lowq_prefix;;
	target->lowq_suffix = c2->lowq_suffix;
}

bool merge_overlapping_pair_of_clusters(consensus_t* c1, consensus_t* c2, consensus_t* target, int min_overlap) {
	// only merge if the overlap is at least half of the length of the shorter consensus
	hts_pos_t overlap = c1->end - c2->start;
	if (overlap < (int) std::min(c1->sequence.length(), c2->sequence.length())/2) return false;
	
	hts_pos_t c1_hq_end = c1->end - c1->lowq_suffix;
	hts_pos_t c2_hq_start = c2->start + c2->lowq_prefix;
	hts_pos_t hq_overlap = c1_hq_end - c2_hq_start;
	std::string c1_hq_seq = c1->sequence.substr(0, c1->sequence.length()-c1->lowq_suffix);
	std::string c2_hq_seq = c2->sequence.substr(c2->lowq_prefix);
	int min_hq_seq_len = std::min(c1_hq_seq.length(), c2_hq_seq.length());
	if (hq_overlap >= min_hq_seq_len/2 && hq_overlap < min_hq_seq_len) { // first, try to see if the high quality parts are compatible
		suffix_prefix_aln_t spa_hq = aln_suffix_prefix_perfect(c1_hq_seq, c2_hq_seq, min_overlap);
		if (spa_hq.overlap) {
			merge_overlapping_pair_of_clusters(c1, c2, target, c1_hq_seq, c2_hq_seq, spa_hq.overlap);
			return true;
		}
	}
	
	suffix_prefix_aln_t spa = aln_suffix_prefix_perfect(c1->sequence, c2->sequence, min_overlap); // then, try the whole sequences
	if (spa.overlap) {
		merge_overlapping_pair_of_clusters(c1, c2, target, c1->sequence, c2->sequence, spa.overlap);
		return true;
	}
	return false;
}

// This only applies to HSR consensuses
void merge_overlapping_clusters(std::vector<consensus_t*>& consensuses, int min_overlap) {

	if (consensuses.empty()) return;

	std::vector<consensus_t*> sorted_by_start = consensuses;
	std::sort(sorted_by_start.begin(), sorted_by_start.end(), [](const consensus_t* c1, const consensus_t* c2) {
		return c1->start < c2->start;
	});

	std::vector<consensus_t*> sorted_by_end = consensuses;
	std::sort(sorted_by_end.begin(), sorted_by_end.end(), [](const consensus_t* c1, const consensus_t* c2) {
		return c1->end < c2->end;
	});

	std::unordered_map<consensus_t*, int> start_pos, end_pos;
	for (int i = 0; i < sorted_by_start.size(); i++) {
		start_pos[sorted_by_start[i]] = i;
		end_pos[sorted_by_end[i]] = i;
	}

	std::vector<consensus_t*> sorted_by_reads = consensuses;
	std::sort(sorted_by_reads.begin(), sorted_by_reads.end(), [](consensus_t* c1, consensus_t* c2) {
		return c1->reads() > c2->reads();
	});

	std::unordered_set<consensus_t*> removed;
	for (int i = 0; i < sorted_by_reads.size(); i++) {
		consensus_t* c = sorted_by_reads[i];
		if (removed.count(c) > 0) continue;

		int e_pos = end_pos[c];
		for (int j = e_pos-1; j >= 0; j--) {
			consensus_t* o = sorted_by_end[j];
			if (removed.count(o) > 0) continue;
			if (o->end < c->start) break;
			if (o->start >= c->start || o->end == c->end) continue; // o is fully contained within c

			if (merge_overlapping_pair_of_clusters(o, c, c, min_overlap)) {
				removed.insert(o);
				std::swap(sorted_by_start[start_pos[o]], sorted_by_start[start_pos[c]]);
				std::swap(start_pos[o], start_pos[c]);
				// end positions do not change
			}
		}

		int s_pos = start_pos[c];
		for (int j = s_pos+1; j < sorted_by_start.size(); j++) {
			consensus_t* o = sorted_by_start[j];
			if (removed.count(o) > 0) continue;
			if (o->start > c->end) break;
			if (o->start == c->start || o->end <= c->end) continue; // o is fully contained within c
			
			if (merge_overlapping_pair_of_clusters(c, o, c, min_overlap)) {
				removed.insert(o);
				std::swap(sorted_by_end[end_pos[o]], sorted_by_end[end_pos[c]]);
				std::swap(end_pos[o], end_pos[c]);
				// start positions do not change
			}
		}
	}

	consensuses.erase(std::remove_if(consensuses.begin(), consensuses.end(),
		[&removed](consensus_t* c){ return removed.count(c) > 0; }), consensuses.end());
}

void enforce_max_ploidy(std::vector<consensus_t*>& consensuses, int max_ploidy) {

	int contig_len = 0;
	for (consensus_t* c : consensuses) {
		if (c->end > contig_len) contig_len = c->end;
	}

	SegTree segtree(contig_len+1);
	std::vector<consensus_t*> retained;
	std::sort(consensuses.begin(), consensuses.end(), [](consensus_t* c1, consensus_t* c2) {
		return c1->reads() > c2->reads();
	});
	for (consensus_t* c : consensuses) {
		if (!segtree.any_ge(c->start, c->end, max_ploidy)) {
			segtree.add(c->start, c->end, 1);
			retained.push_back(c);
		} else {
			delete c;
		}
	}
	consensuses.swap(retained);
}

#endif
