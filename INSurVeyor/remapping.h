#ifndef REMAPPING_H
#define REMAPPING_H

#include <vector>
#include <string>
#include <htslib/sam.h>
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"

struct region_score_t {
    int total_score = 0;
    int remap_start = 0, remap_end = 0;
};
bool operator > (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score > s2.total_score;
}
bool operator >= (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score >= s2.total_score;
}
bool operator < (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score < s2.total_score;
}

struct region_t {
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    region_score_t score;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}
};
bool operator > (const region_t& r1, const region_t& r2) {
	if (r1.score.total_score != r2.score.total_score) {
		return r1.score.total_score > r2.score.total_score;
	}
	return r1.end-r1.start < r2.end-r2.start;
}

struct remap_info_t {
    int start, end;
    int score;
    std::string cigar;
    bool accepted, left_clipped, right_clipped;

    remap_info_t() : start(0), end(0), score(0), cigar(""), accepted(false), left_clipped(false), right_clipped(false) {}

    remap_info_t(StripedSmithWaterman::Alignment& aln, bool accepted, int min_clip_len) : start(aln.ref_begin), end(aln.ref_end),
    score(aln.sw_score), cigar(aln.cigar_string), accepted(accepted) {
        uint32_t f = aln.cigar[0], l = aln.cigar[aln.cigar.size()-1];
        left_clipped = cigar_int_to_op(f) == 'S' && cigar_int_to_len(f) >= min_clip_len;
        right_clipped = cigar_int_to_op(l) == 'S' && cigar_int_to_len(l) >= min_clip_len;
    }
};

region_t get_candidate_region(std::vector<bam1_t*> subcluster, std::string& m_contig_name, int m_contig_id, hts_pos_t m_contig_len, int max_is, int max_trans_size) {
    bam1_t* leftmost_reverse_mate = NULL, * rightmost_reverse_mate = NULL;
    bam1_t* leftmost_forward_mate = NULL, * rightmost_forward_mate = NULL;
    /* gets two regions
     * 1. from reverse mates, max_is to the left of right-most read and max_insertion_size to the right of the left-most
     * 2. from right-most forward read, max_is to the right and max_insertion_size to the left
     */
    for (bam1_t* r : subcluster) { // get leftmost reverse read
        if (bam_is_mrev(r)) {
            if (leftmost_reverse_mate == NULL || leftmost_reverse_mate->core.mpos > r->core.mpos) {
                leftmost_reverse_mate = r;
            }
            if (rightmost_reverse_mate == NULL || rightmost_reverse_mate->core.mpos < r->core.mpos) {
                rightmost_reverse_mate = r;
            }
        }
    }
    for (bam1_t* r : subcluster) { // get rightmost forward read
        if (!bam_is_mrev(r)) {
            if (leftmost_forward_mate == NULL || leftmost_forward_mate->core.mpos > r->core.mpos) {
                leftmost_forward_mate = r;
            }
            if (rightmost_forward_mate == NULL || rightmost_forward_mate->core.mpos < r->core.mpos) {
                rightmost_forward_mate = r;
            }
        }
    }

    hts_pos_t start = INT_MAX;
    hts_pos_t end = 0;
    if (leftmost_reverse_mate != NULL) {
        start = std::min(start, rightmost_reverse_mate->core.mpos-max_is);
        end = std::max(end, leftmost_reverse_mate->core.mpos+max_trans_size);
    }
    if (rightmost_forward_mate != NULL) {
        start = std::min(start, rightmost_forward_mate->core.mpos-max_trans_size);
        end = std::max(end, leftmost_forward_mate->core.mpos+max_is);
    }

    // hts_pos_t contig_len = contigs.get_len(m_contig_name);
    return region_t(m_contig_id, subcluster[0]->core.mtid, std::max(hts_pos_t(0), start), std::min(end, m_contig_len));
}

#endif // REMAPPING_H
