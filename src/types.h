#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <htslib/sam.h>

// TODO: extract common members of subclasses to here
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

#endif /* TYPES_H */
