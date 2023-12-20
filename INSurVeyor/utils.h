#ifndef SMALLINSFINDER_UTILS_H
#define SMALLINSFINDER_UTILS_H

#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <emmintrin.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "../src/utils.h"
#include "../src/sw_utils.h"
#include "../src/simd_macros.h"

struct inss_insertion_t {
    std::string id;
    std::string chr;
    hts_pos_t start, end;
    int overlap, r_disc_pairs, l_disc_pairs, rc_fwd_reads = 0, rc_rev_reads = 0, lc_fwd_reads = 0, lc_rev_reads = 0;
    int r_conc_pairs = 0, l_conc_pairs = 0, median_lf_cov = 0, median_rf_cov = 0;
    int mh_len = 0;
    std::string ins_seq;
    std::string left_anchor, right_anchor, left_anchor_cigar, right_anchor_cigar;
    int left_seq_cov = 0, right_seq_cov = 0;
    bool left_bp_precise = false, right_bp_precise = false;
    double rc_avg_nm = 0.0, lc_avg_nm = 0.0;

    inss_insertion_t(std::string chr, hts_pos_t start, hts_pos_t end, int r_disc_pairs, int l_disc_pairs,
    		int rc_fwd_reads, int rc_rev_reads, int lc_fwd_reads, int lc_rev_reads, int overlap, std::string ins_seq) :
	chr(chr), start(start), end(end), r_disc_pairs(r_disc_pairs), l_disc_pairs(l_disc_pairs),
	rc_fwd_reads(rc_fwd_reads), rc_rev_reads(rc_rev_reads), lc_fwd_reads(lc_fwd_reads), lc_rev_reads(lc_rev_reads),
	overlap(overlap), ins_seq(ins_seq) {}

    int rc_reads() { return rc_fwd_reads + rc_rev_reads; }
    int lc_reads() { return lc_fwd_reads + lc_rev_reads; }

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + ins_seq;
    }
};

// Returns a vector scores s.t. scores[N] contains the score of the alignment between reference[aln.ref_begin:aln.ref_begin+N]
// and query[aln.query_begin:aln.query_begin+M]
std::vector<int> inss_ssw_cigar_to_prefix_ref_scores(uint32_t* cigar, int cigar_len,
		int match = 1, int mismatch = -4, int gap_open = -6, int gap_extend = -1) {

	std::vector<int> scores;
	scores.push_back(0);
	int score = 0;

	for (int i = 0; i < cigar_len; i++) {
		uint32_t c = cigar[i];
		char op = cigar_int_to_op(c);
		int len = cigar_int_to_len(c);
		if (op == 'S') {
			continue;
		} else if (op == '=' || op == 'X') {
			for (int i = 0; i < len; i++) {
				score += (op == '=' ? match : mismatch);
				scores.push_back(score);
			}
		} else if (op == 'D') {
			scores.resize(scores.size()+len, score);
			score += gap_open + len*gap_extend;
		} else if (op == 'I') {
			score += gap_open + len*gap_extend;
		}
	}

	return scores;
}

#endif //SMALLINSFINDER_UTILS_H
