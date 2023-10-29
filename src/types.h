#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <htslib/sam.h>

struct consensus_t {
    bool left_clipped;
    hts_pos_t start, breakpoint, end;
    std::string sequence;
    int fwd_clipped, rev_clipped;
    uint8_t max_mapq;
    hts_pos_t remap_boundary;
    int clip_len, lowq_clip_portion;
    int left_ext_reads = 0, right_ext_reads = 0, hq_left_ext_reads = 0, hq_right_ext_reads = 0;
    bool is_hsr = false;
	bool extended_to_left = false, extended_to_right = false;

    static const int LOWER_BOUNDARY_NON_CALCULATED = 0, UPPER_BOUNDARY_NON_CALCULATED = INT32_MAX;
    static const int UNKNOWN_CLIP_LEN = INT16_MAX;

    consensus_t(bool left_clipped, hts_pos_t start, hts_pos_t breakpoint, hts_pos_t end,
                const std::string& sequence, int fwd_clipped, int rev_clipped, int clip_len, uint8_t max_mapq, 
                hts_pos_t remap_boundary, int lowq_clip_portion)
                : left_clipped(left_clipped), start(start), breakpoint(breakpoint), end(end),
                sequence(sequence), fwd_clipped(fwd_clipped), rev_clipped(rev_clipped), clip_len(clip_len), max_mapq(max_mapq), 
                remap_boundary(remap_boundary), lowq_clip_portion(lowq_clip_portion) {}

    std::string to_string() {
        std::stringstream ss;
        ss << start << " " << end << " " << breakpoint << (left_clipped ? " L " : " R ") << sequence << " ";
        ss << fwd_clipped << " " << rev_clipped << " " << (int)max_mapq << " " << remap_boundary << " " << lowq_clip_portion;
        return ss.str();
    }

    int supp_clipped_reads() { return fwd_clipped + rev_clipped; }

    hts_pos_t left_ext_target_start(int max_is, int read_len) {
    	if (!left_clipped) {
    		return start - max_is + read_len;
    	} else {
    		if (remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint - max_is - 2*sequence.length();
			} else {
				return remap_boundary;
			}
    	}
    }
    hts_pos_t left_ext_target_end(int max_is, int read_len) {
		if (!left_clipped) {
			return start;
		} else {
			if (remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint + max_is + 2*sequence.length();
			} else {
				return remap_boundary + max_is;
			}
		}
	}

    hts_pos_t right_ext_target_start(int max_is, int read_len) {
		if (!left_clipped) {
			if (remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint - max_is - 2*sequence.length();
			} else {
				return remap_boundary - max_is;
			}
		} else {
			return end;
		}
	}
    hts_pos_t right_ext_target_end(int max_is, int read_len) {
    	if (!left_clipped) {
			if (remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint + max_is + 2*sequence.length();
			} else {
				return remap_boundary;
			}
		} else {
			return end + max_is - read_len;
		}
    }

    int anchor_len() { return sequence.length() - clip_len; }
};

struct sv_t {
    std::string id;
    std::string chr;
    hts_pos_t start, end;
    std::string ins_seq;
    int rc_fwd_reads = 0, rc_rev_reads = 0, lc_fwd_reads = 0, lc_rev_reads = 0;
    int overlap;
    std::string left_anchor, right_anchor, left_anchor_cigar, right_anchor_cigar;

    sv_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq) : chr(chr), start(start), end(end), ins_seq(ins_seq) {}

    int rc_reads() { return rc_fwd_reads + rc_rev_reads; }
    int lc_reads() { return lc_fwd_reads + lc_rev_reads; }

    virtual std::string unique_key() = 0;

    virtual std::string svtype() = 0;
    virtual hts_pos_t svlen() = 0;

    virtual ~sv_t() {}
};

struct deletion_t : sv_t {
    deletion_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq) : sv_t(chr, start, end, ins_seq) {}

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end);
    }

    std::string svtype() { return "DEL"; }
    hts_pos_t svlen() { return start - end + ins_seq.length(); }
};

struct duplication_t : sv_t {
    duplication_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq) : sv_t(chr, start, end, ins_seq) {}

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end);
    }

    std::string svtype() { return "DUP"; }
    hts_pos_t svlen() { return end - start; }
};

struct insertion_t : sv_t {
    int r_disc_pairs, l_disc_pairs;
    int r_conc_pairs = 0, l_conc_pairs = 0, median_lf_cov = 0, median_rf_cov = 0;
    int left_seq_cov = 0, right_seq_cov = 0;
    bool left_bp_precise = false, right_bp_precise = false;
    double rc_avg_nm = 0.0, lc_avg_nm = 0.0;

    insertion_t(std::string chr, hts_pos_t start, hts_pos_t end, int r_disc_pairs, int l_disc_pairs, std::string ins_seq) :
	sv_t(chr, start, end, ins_seq), r_disc_pairs(r_disc_pairs), l_disc_pairs(l_disc_pairs) {}

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + ins_seq;
    }

    std::string svtype() { return "INS"; }
    hts_pos_t svlen() { return ins_seq.length() - (end-start); }
};


struct clip_consensus_t {
    bool left_clipped;
    hts_pos_t breakpoint, start, end; // we follow the vcf conventions, i.e. this is the base "before" the breakpoint
    std::string seq;
    int fwd_clipped, rev_clipped;
    int clip_len, lowq_clip_portion;
    uint8_t max_mapq;
    hts_pos_t remap_boundary;

    static const int LOWER_BOUNDARY_NON_CALCULATED = 0, UPPER_BOUNDARY_NON_CALCULATED = INT32_MAX;
    static const int UNKNOWN_CLIP_LEN = INT16_MAX;

    clip_consensus_t(bool left_clipped, hts_pos_t start, hts_pos_t end, hts_pos_t breakpoint, std::string& seq, 
    int fwd_clipped, int rev_clipped, int clip_len, int max_mapq, hts_pos_t remap_boundary, int lowq_clip_portion) :
        left_clipped(left_clipped), start(start), end(end), breakpoint(breakpoint), seq(seq), fwd_clipped(fwd_clipped), rev_clipped(rev_clipped),
        clip_len(clip_len), max_mapq(max_mapq), remap_boundary(remap_boundary), lowq_clip_portion(lowq_clip_portion) {}
};


#endif /* TYPES_H */
