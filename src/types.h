#ifndef TYPES_H
#define TYPES_H

#include "htslib/hts.h"
#include "utils.h"
#include <vector>
#include <string>
#include <sstream>
#include "htslib/sam.h"
#include "htslib/vcf.h"

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
    
    consensus_t(std::string& line, bool is_hsr) : is_hsr(is_hsr) {
        std::stringstream ss(line);
        char dir;
        int max_mapq_int;
        ss >> start >> end >> breakpoint >> dir >> sequence >> fwd_clipped >> rev_clipped >> max_mapq_int >> remap_boundary >> lowq_clip_portion;
        left_clipped = dir == 'L';
        max_mapq = (uint8_t) max_mapq_int;
        if (is_hsr) {
            clip_len = UNKNOWN_CLIP_LEN;
        } else {
            clip_len = left_clipped ? breakpoint - start : end - breakpoint;
        }
    }

    std::string name() {
        return std::to_string(start) + "_" + std::to_string(end) + "_" + (left_clipped ? "L" : "R");
    }

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

    hts_pos_t anchor_start() {
        if (left_clipped) {
            return breakpoint;
        } else {
            return start;
        }
    }
    hts_pos_t anchor_end() {
        if (left_clipped) {
            return end;
        } else {
            return breakpoint;
        }
    }

    int anchor_len() { return sequence.length() - clip_len; }

    std::string clip_sequence() {
        if (left_clipped) {
            return sequence.substr(0, clip_len);
        } else {
            return sequence.substr(sequence.length()-clip_len);
        }
    }
};

struct sv_t {
    
    struct anchor_aln_t {
        hts_pos_t start, end;
        int seq_len;
        int best_score, next_best_score;
        std::string cigar;

        anchor_aln_t(hts_pos_t start, hts_pos_t end, int seq_len, int best_score, int next_best_score, std::string cigar) : 
            start(start), end(end), seq_len(seq_len), best_score(best_score), next_best_score(next_best_score), cigar(cigar) {}
    };

    std::string id;
    std::string chr;
    hts_pos_t start, end;
    std::string ins_seq;
    int prefix_mh_len = 0, suffix_mh_len = 0;
    anchor_aln_t* left_anchor_aln,* right_anchor_aln,* full_junction_aln;
    consensus_t* rc_consensus, * lc_consensus;
    int disc_pairs_lf = 0, disc_pairs_rf = 0, disc_pairs_lf_high_mapq = 0, disc_pairs_rf_high_mapq = 0, disc_pairs_lf_maxmapq = 0, 
        disc_pairs_rf_maxmapq = 0, conc_pairs = 0;
    double disc_pairs_lf_avg_nm = 0, disc_pairs_rf_avg_nm = 0;

    int median_left_flanking_cov = 0, median_indel_left_cov = 0, median_indel_right_cov = 0, median_right_flanking_cov = 0;
    int median_left_cluster_cov = 0, median_right_cluster_cov = 0;
    int median_left_flanking_cov_highmq = 0, median_indel_left_cov_highmq = 0;
    int median_indel_right_cov_highmq = 0, median_right_flanking_cov_highmq = 0;
    int median_left_cluster_cov_highmq = 0, median_right_cluster_cov_highmq = 0;
    int l_cluster_region_disc_pairs = 0, r_cluster_region_disc_pairs = 0;

    int overlap = 0;
    double mismatch_rate = 0.0;
    std::string source;
    bool imprecise = false;

    static constexpr const double KS_PVAL_NOT_COMPUTED = -1.0;
    static const int SIZE_NOT_COMPUTED = INT32_MAX;
    
    double ks_pval = KS_PVAL_NOT_COMPUTED, ks_pval_highmq = KS_PVAL_NOT_COMPUTED;
    int min_conf_size = SIZE_NOT_COMPUTED, max_conf_size = SIZE_NOT_COMPUTED, estimated_size = SIZE_NOT_COMPUTED;
    int min_conf_size_highmq = SIZE_NOT_COMPUTED, max_conf_size_highmq = SIZE_NOT_COMPUTED;

    base_frequencies_t left_anchor_base_freqs, right_anchor_base_freqs;
    base_frequencies_t prefix_ref_base_freqs, suffix_ref_base_freqs;
    base_frequencies_t ins_prefix_base_freqs, ins_suffix_base_freqs;

    struct regenotyping_info_t {
        int alt_bp1_better_reads = 0, alt_bp2_better_reads = 0;
        int alt_bp1_better_disc_pairs = 0, alt_bp2_better_disc_pairs = 0;
        int alt_bp1_better_consistent_reads = 0, alt_bp2_better_consistent_reads = 0;
        double alt_bp1_better_consistent_avg_score = 0, alt_bp2_better_consistent_avg_score = 0;
        int ref_bp1_better_reads = 0, ref_bp2_better_reads = 0;
        int ref_bp1_better_consistent_reads = 0, ref_bp2_better_consistent_reads = 0;
        int alt_ref_equal_reads = 0;
        int alt_ext_reads = 0, hq_alt_ext_reads = 0;
        int ext_alt_consensus_length = 0;
        int ext_alt_consensus_to_alt_score = 0, ext_alt_consensus_to_ref_score = 0;
        bool too_deep = false;
    } regenotyping_info;

    int* gt = NULL, ngt;
    bcf1_t* vcf_entry = NULL;

    std::vector<std::string> filters;

    sv_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq, consensus_t* rc_consensus, consensus_t* lc_consensus, 
        anchor_aln_t* left_anchor_aln, anchor_aln_t* right_anchor_aln, anchor_aln_t* full_junction_aln) : 
        chr(chr), start(start), end(end), ins_seq(ins_seq), rc_consensus(rc_consensus), lc_consensus(lc_consensus),
        left_anchor_aln(left_anchor_aln), right_anchor_aln(right_anchor_aln), full_junction_aln(full_junction_aln) {}

    int rc_reads() { return rc_consensus ? rc_consensus->fwd_clipped + rc_consensus->rev_clipped : 0; }
    int lc_reads() { return lc_consensus ? lc_consensus->fwd_clipped + lc_consensus->rev_clipped : 0; }

    int rc_fwd_reads() { return rc_consensus ? rc_consensus->fwd_clipped : 0; }
    int rc_rev_reads() { return rc_consensus ? rc_consensus->rev_clipped : 0; }
    int lc_fwd_reads() { return lc_consensus ? lc_consensus->fwd_clipped : 0; }
    int lc_rev_reads() { return lc_consensus ? lc_consensus->rev_clipped : 0; }

    bool is_pass() { return filters.size() == 1 && filters[0] == "PASS"; }
    bool is_fail() { return !filters.empty() && filters[0] != "PASS"; }

    hts_pos_t remap_boundary_upper() {
        if (rc_consensus == NULL) return consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
        return rc_consensus->remap_boundary;
    }

    hts_pos_t remap_boundary_lower() {
        if (lc_consensus == NULL) return consensus_t::LOWER_BOUNDARY_NON_CALCULATED;
        return lc_consensus->remap_boundary;
    }

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + svtype() + ":" + ins_seq;
    }

    virtual std::string svtype() = 0;
    virtual hts_pos_t svlen() = 0;

    std::string left_anchor_aln_string() {
        if (left_anchor_aln == NULL) return "NA";
        return std::to_string(left_anchor_aln->start+1) + "-" + std::to_string(left_anchor_aln->end+1);
    }
    std::string right_anchor_aln_string() {
        if (right_anchor_aln == NULL) return "NA";
        return std::to_string(right_anchor_aln->start+1) + "-" + std::to_string(right_anchor_aln->end+1);
    }
    std::string full_junction_aln_string() {
        if (full_junction_aln == NULL) return "NA";
        return std::to_string(full_junction_aln->start+1) + "-" + std::to_string(full_junction_aln->end+1);
    }

    bool incomplete_assembly() { return ins_seq.find("-") != std::string::npos; }

    std::string print_gt() {
        std::stringstream ss;
        for (int i = 0; i < ngt; i++) {
            if (i > 0) ss << "/";
            ss << (bcf_gt_is_missing(gt[i]) ? "." : std::to_string(bcf_gt_allele(gt[i])));
        }
        return ss.str();
    }

    int allele_count(int allele) {
        int ac = 0;
        for (int i = 0; i < ngt; i++) {
            if (!bcf_gt_is_missing(gt[i])) ac += (bcf_gt_allele(gt[i]) == allele);
        }
        return ac;
    }

    int missing_alleles() {
        int ac = 0;
        for (int i = 0; i < ngt; i++) {
            ac += bcf_gt_is_missing(gt[i]);
        }
        return ac;
    }

    base_frequencies_t get_left_anchor_base_freqs(char* chr_seq) {
        if (left_anchor_base_freqs.empty()) {
            left_anchor_base_freqs = get_base_frequencies(chr_seq+left_anchor_aln->start, left_anchor_aln->end-left_anchor_aln->start);
        }
        return left_anchor_base_freqs;
    }
    base_frequencies_t get_right_anchor_base_freqs(char* chr_seq) {
        if (right_anchor_base_freqs.empty()) {
            right_anchor_base_freqs = get_base_frequencies(chr_seq+right_anchor_aln->start, right_anchor_aln->end-right_anchor_aln->start);
        }
        return right_anchor_base_freqs;
    }

    base_frequencies_t get_prefix_ref_base_freqs(char* chr_seq) {
        if (prefix_ref_base_freqs.empty()) {
            prefix_ref_base_freqs = get_base_frequencies(chr_seq+start, std::min(end-start, hts_pos_t(5000)));
        }
        return prefix_ref_base_freqs;
    }
    base_frequencies_t get_suffix_ref_base_freqs(char* chr_seq) {
        if (suffix_ref_base_freqs.empty()) {
            suffix_ref_base_freqs = get_base_frequencies(chr_seq+std::max(start, end-5000), std::min(end-start, hts_pos_t(5000)));
        }
        return suffix_ref_base_freqs;
    }

    base_frequencies_t get_ins_prefix_base_freqs() {
        if (ins_prefix_base_freqs.empty()) {
            int d = ins_seq.find("-");
            std::string ins_seq_fh = ins_seq.substr(0, d);
            ins_prefix_base_freqs = get_base_frequencies(ins_seq_fh.c_str(), ins_seq_fh.length());
        }
        return ins_prefix_base_freqs;
    }

    base_frequencies_t get_ins_suffix_base_freqs() {
        if (ins_suffix_base_freqs.empty()) {
            int d = ins_seq.find("-");
            if (d == std::string::npos) {
                ins_suffix_base_freqs = get_ins_prefix_base_freqs();
            } else {
                std::string ins_seq_sh = ins_seq.substr(d+1);
                ins_suffix_base_freqs = get_base_frequencies(ins_seq_sh.c_str(), ins_seq_sh.length());
            }
        }
        return ins_suffix_base_freqs;
    }

    void precompute_base_frequencies(char* chr_seq) {
        get_left_anchor_base_freqs(chr_seq);
        get_right_anchor_base_freqs(chr_seq);
        get_prefix_ref_base_freqs(chr_seq);
        get_suffix_ref_base_freqs(chr_seq);
        get_ins_prefix_base_freqs();
        get_ins_suffix_base_freqs();
    }

    virtual ~sv_t() {}
};

struct deletion_t : sv_t {

    bool remapped = false;
    std::string original_range;

    using sv_t::sv_t;

    std::string svtype() { return "DEL"; }
    hts_pos_t svlen() { return start - end + ins_seq.length(); }
};

struct duplication_t : sv_t {
    using sv_t::sv_t;

    std::string svtype() { return "DUP"; }
    hts_pos_t svlen() { return end - start + ins_seq.length(); }
};

struct insertion_t : sv_t {
    using sv_t::sv_t;

    static const int NOT_COMPUTED = INT32_MAX;

    int prefix_cov_start = NOT_COMPUTED, prefix_cov_end = NOT_COMPUTED, suffix_cov_start = NOT_COMPUTED, suffix_cov_end = NOT_COMPUTED; // start and end of the prefix and suffix of the inserted sequence actually supported by reads. Only applicable to long transpositions, for which the central part is inferred

    std::string svtype() { return "INS"; }
    hts_pos_t svlen() { return ins_seq.length() - (end-start); }
};

#endif /* TYPES_H */
