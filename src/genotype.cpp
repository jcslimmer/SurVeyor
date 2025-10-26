#include "genotype.h"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <new>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "assemble.h"
#include "htslib/hts_endian.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

#include "sw_utils.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "extend_1sr_consensus.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "vcf_utils.h"
#include "stat_tests.h"
#include "reference_guided_assembly.h"
#include "SegTree.h"
#include "../libs/IntervalTree.h"

#include "genotype_dels.h"
#include "genotype_dups.h"
#include "genotype_inss.h"
#include "genotype_invs.h"

chr_seqs_map_t chr_seqs;
config_t config;
contig_map_t contig_map;
stats_t stats;
std::mutex mtx;

std::string bam_fname, reference_fname, workdir;
bam_pool_t* bam_pool;

std::vector<hts_pos_t> global_isize_dist;

std::vector<double> global_crossing_isize_dist;

StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, false);

std::vector<std::unordered_map<std::string, std::pair<std::string, int> > > mateseqs_w_mapq;
std::vector<int> active_threads_per_chr;
std::vector<std::mutex> mutex_per_chr;

void update_record_bp_reads_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_reads_info_t bp_reads_info, std::string prefix, int bp_number) {
    if (!bp_reads_info.computed) return;

    std::string bp_number_str = std::to_string(bp_number);
    std::string read_fmt_prefix = prefix + "R" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, read_fmt_prefix.c_str(), &(bp_reads_info.reads), 1);

    int consistent_reads = bp_reads_info.consistent_reads();
    if (bp_reads_info.reads) {
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "C").c_str(), &(consistent_reads), 1);
    }
    if (consistent_reads) {
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CF").c_str(), &(bp_reads_info.consistent_fwd), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CR").c_str(), &(bp_reads_info.consistent_rev), 1);
        float rcas = bp_reads_info.consistent_avg_score;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAS").c_str(), &rcas, 1);
        float rcss = bp_reads_info.consistent_stddev_score;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSS").c_str(), &rcss, 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CHQ").c_str(), &(bp_reads_info.consistent_high_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CmQ").c_str(), &(bp_reads_info.consistent_min_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMQ").c_str(), &(bp_reads_info.consistent_max_mq), 1);
        float rcavmq = bp_reads_info.consistent_avg_mq;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAQ").c_str(), &rcavmq, 1);
        float rcstdmq = bp_reads_info.consistent_stddev_mq;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSQ").c_str(), &rcstdmq, 1);
        int mate_cov_bps[] = {bp_reads_info.fwd_mate_cov_bps, bp_reads_info.rev_mate_cov_bps};
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMSPAN").c_str(), mate_cov_bps, 2);
        int hq_mate_cov_bps[] = {bp_reads_info.fwd_hq_mate_cov_bps, bp_reads_info.rev_hq_mate_cov_bps};
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMHQSPAN").c_str(), hq_mate_cov_bps, 2);
    }
}

void reset_record_bp_reads_info(bcf_hdr_t* out_hdr, bcf1_t* b, std::string prefix, int bp_number) {
    std::string bp_number_str = std::to_string(bp_number);
    std::string read_fmt_prefix = prefix + "R" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, read_fmt_prefix.c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "C").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CF").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CR").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAS").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSS").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CHQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CmQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMSPAN").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMHQSPAN").c_str(), NULL, 0);
}

void update_record_bp_pairs_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_pairs_info_t bp_pairs_info, std::string prefix, int bp_number) {
    if (!bp_pairs_info.computed) return;

    std::string bp_number_str = std::to_string(bp_number);
    std::string pairs_fmt_prefix = prefix + "SP" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, pairs_fmt_prefix.c_str(), &(bp_pairs_info.pairs), 1);

    if (bp_pairs_info.pairs) {
        int pairs_hq[] = {bp_pairs_info.pos_high_mapq, bp_pairs_info.neg_high_mapq};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "HQ").c_str(), pairs_hq, 2);
        int pairs_min_mq[] = {bp_pairs_info.pos_min_mq, bp_pairs_info.neg_min_mq};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "mQ").c_str(), pairs_min_mq, 2);
        int pairs_max_mq[] = {bp_pairs_info.pos_max_mq, bp_pairs_info.neg_max_mq};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "MQ").c_str(), pairs_max_mq, 2);
        float pairs_avg_mq[] = {(float) bp_pairs_info.pos_avg_mq, (float) bp_pairs_info.neg_avg_mq};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "AQ").c_str(), pairs_avg_mq, 2);
        float pairs_stddev_mq[] = {(float) bp_pairs_info.pos_stddev_mq, (float) bp_pairs_info.neg_stddev_mq};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "SQ").c_str(), pairs_stddev_mq, 2);
        int pairs_span[] = {bp_pairs_info.lf_span, bp_pairs_info.rf_span};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "SPAN").c_str(), pairs_span, 2);
        float pairs_avg_nm[] = {(float) bp_pairs_info.pos_avg_nm, (float) bp_pairs_info.neg_avg_nm};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMA").c_str(), pairs_avg_nm, 2);
        float pairs_stddev_nm[] = {(float) bp_pairs_info.pos_stddev_nm, (float) bp_pairs_info.neg_stddev_nm};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMS").c_str(), pairs_stddev_nm, 2);
    }
}

void reset_record_bp_pairs_info(bcf_hdr_t* out_hdr, bcf1_t* b, std::string prefix, int bp_number) {
    std::string bp_number_str = std::to_string(bp_number);
    std::string pairs_fmt_prefix = prefix + "SP" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, pairs_fmt_prefix.c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "HQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "mQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "MQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "AQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "SQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "SPAN").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMA").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMS").c_str(), NULL, 0);
}

void update_record_bp_consensus_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_consensus_info_t bp_consensus_info, std::string prefix, int bp_number) {
    update_record_bp_reads_info(out_hdr, b, bp_consensus_info.reads_info, prefix, bp_number);
    update_record_bp_pairs_info(out_hdr, b, bp_consensus_info.pairs_info, prefix, bp_number);
}

void reset_record_bp_consensus_info(bcf_hdr_t* out_hdr, bcf1_t* b, std::string prefix, int bp_number) {
    reset_record_bp_reads_info(out_hdr, b, prefix, bp_number);
    reset_record_bp_pairs_info(out_hdr, b, prefix, bp_number);
}

void update_record(bcf_hdr_t* in_hdr, bcf_hdr_t* out_hdr, sv_t* sv, char* chr_seq, hts_pos_t chr_len, int sample_idx) {
    
    bcf_translate(out_hdr, in_hdr, sv->vcf_entry);

    bcf_subset(out_hdr, sv->vcf_entry, 1, &sample_idx);

    // update INFO fields
    int sv_end = sv->end+1;
    bcf_update_info_int32(out_hdr, sv->vcf_entry, "END", &sv_end, 1);

    if (sv->ins_seq.find("-") == std::string::npos) {
        if (!sv->ins_seq.empty()) {
            int svinslen = sv->ins_seq.length();
            bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVINSLEN", &svinslen, 1);
        }
        int svlen = sv->svlen();
        bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVLEN", &svlen, 1);
    }

    hts_pos_t la_aln_start = bring_within_range(sv->left_anchor_aln->start, hts_pos_t(0), chr_len);
    hts_pos_t la_aln_end = bring_within_range(sv->left_anchor_aln->end, hts_pos_t(0), chr_len);
    hts_pos_t ra_aln_start = bring_within_range(sv->right_anchor_aln->start, hts_pos_t(0), chr_len);
    hts_pos_t ra_aln_end = bring_within_range(sv->right_anchor_aln->end, hts_pos_t(0), chr_len);
    base_frequencies_t left_anchor_base_freqs = get_base_frequencies(chr_seq+la_aln_start, la_aln_end-la_aln_start);
	int labc[] = {left_anchor_base_freqs.a, left_anchor_base_freqs.c, left_anchor_base_freqs.g, left_anchor_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "LEFT_ANCHOR_BASE_COUNT", labc, 4);

	base_frequencies_t right_anchor_base_freqs = get_base_frequencies(chr_seq+ra_aln_start, ra_aln_end-ra_aln_start);
	int rabc[] = {right_anchor_base_freqs.a, right_anchor_base_freqs.c, right_anchor_base_freqs.g, right_anchor_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "RIGHT_ANCHOR_BASE_COUNT", rabc, 4);

	base_frequencies_t prefix_ref_base_freqs = get_base_frequencies(chr_seq+sv->start, std::min(sv->end-sv->start, hts_pos_t(5000)));
	int svrefpbc[] = {prefix_ref_base_freqs.a, prefix_ref_base_freqs.c, prefix_ref_base_freqs.g, prefix_ref_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "SV_REF_PREFIX_BASE_COUNT", svrefpbc, 4);

	base_frequencies_t suffix_ref_base_freqs = get_base_frequencies(chr_seq+sv->end-std::min(sv->end-sv->start, hts_pos_t(5000)), std::min(sv->end-sv->start, hts_pos_t(5000)));
	int svrefsbc[] = {suffix_ref_base_freqs.a, suffix_ref_base_freqs.c, suffix_ref_base_freqs.g, suffix_ref_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "SV_REF_SUFFIX_BASE_COUNT", svrefsbc, 4);

    int d = sv->ins_seq.find("-");
	std::string ins_seq_fh = sv->ins_seq.substr(0, d);
	std::string ins_seq_sh = sv->ins_seq.substr(d+1);
	base_frequencies_t prefix_base_freqs = get_base_frequencies(ins_seq_fh.c_str(), ins_seq_fh.length());
	base_frequencies_t suffix_base_freqs = ins_seq_sh.empty() ? prefix_base_freqs : get_base_frequencies(ins_seq_sh.c_str(), ins_seq_sh.length());
	int pbc[] = {prefix_base_freqs.a, prefix_base_freqs.c, prefix_base_freqs.g, prefix_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "INS_PREFIX_BASE_COUNT", pbc, 4);
	int sbc[] = {suffix_base_freqs.a, suffix_base_freqs.c, suffix_base_freqs.g, suffix_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "INS_SUFFIX_BASE_COUNT", sbc, 4);

    int mh_len = 0;
    if ((sv->svtype() == "DEL" || sv->svtype() == "DUP") && sv->ins_seq.empty()) {
        while (mh_len < abs(sv->svlen()) && sv->end+mh_len+1 < chr_len && toupper(chr_seq[sv->start+mh_len+1]) == toupper(chr_seq[sv->end+mh_len+1])) {
            mh_len++;
        }
    } else if (sv->svtype() == "INS" && sv->start == sv->end) {
        while (mh_len < sv->ins_seq.length() && sv->end+mh_len+1 < chr_len && toupper(sv->ins_seq[mh_len]) == toupper(chr_seq[sv->end+mh_len+1])) {
            mh_len++;
        }
    }
    if (mh_len > 0) {
        bcf_update_info_int32(out_hdr, sv->vcf_entry, "MH_LEN", &mh_len, 1);
    }

    // update FORMAT fields
    bcf_update_genotypes(out_hdr, sv->vcf_entry, sv->sample_info.gt, sv->n_gt);

    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "A", 1);
    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "A", 2);
    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "R", 1);
    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "R", 2);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "N", 1);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "N", 2);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "S", 1);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "S", 2);

    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.alt_bp1, "A", 1);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.alt_bp2, "A", 2);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.ref_bp1, "R", 1);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.ref_bp2, "R", 2);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.neutral_bp1_pairs, "N", 1);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.neutral_bp2_pairs, "N", 2);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.bp1_stray_pairs, "S", 1);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.bp2_stray_pairs, "S", 2);

    int er = sv->sample_info.alt_ref_equal_reads;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ER", &er, 1);

    if (sv->sample_info.too_deep) {
        int td = 1;
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "TD", &td, 1);
    }

    int median_depths[] = {sv->sample_info.left_flanking_cov, sv->sample_info.indel_left_cov, sv->sample_info.indel_right_cov, sv->sample_info.right_flanking_cov};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MD", median_depths, 4);

    int median_depths_highmq[] = {sv->sample_info.left_flanking_cov_highmq, sv->sample_info.indel_left_cov_highmq, sv->sample_info.indel_right_cov_highmq, sv->sample_info.right_flanking_cov_highmq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDHQ", median_depths_highmq, 4);

    int cluster_depths[] = {sv->sample_info.left_anchor_cov, sv->sample_info.right_anchor_cov};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CLMD", cluster_depths, 2);

    int cluster_depths_highmq[] = {sv->sample_info.left_anchor_cov_highmq, sv->sample_info.right_anchor_cov_highmq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CLMDHQ", cluster_depths_highmq, 2);

    if (sv->min_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZE", &(sv->min_conf_size), 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZE", NULL, 0);
    }
    if (sv->max_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZE", &(sv->max_conf_size), 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZE", NULL, 0);
    }
    if (sv->ks_pval != deletion_t::KS_PVAL_NOT_COMPUTED) {
        float ks_pval = sv->ks_pval;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVAL", &ks_pval, 1);
    } else {
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVAL", NULL, 0);
    }

    int ext_reads[] = {sv->sample_info.alt_lext_reads, sv->sample_info.alt_rext_reads};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXR", ext_reads, 2);
    int hq_ext_reads[] = {sv->sample_info.hq_alt_lext_reads, sv->sample_info.hq_alt_rext_reads};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXRHQ", hq_ext_reads, 2);

    if (sv->sample_info.ext_alt_consensus1_length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS", &(sv->sample_info.ext_alt_consensus1_to_alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS", &(sv->sample_info.ext_alt_consensus1_to_ref_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL", &(sv->sample_info.ext_alt_consensus1_length), 1);
        int exss[] = {sv->sample_info.alt_consensus1_split_size1, sv->sample_info.alt_consensus1_split_size2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS", exss, 2);
        int exssc[] = {sv->sample_info.alt_consensus1_split_score1, sv->sample_info.alt_consensus1_split_score2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC", exssc, 2);
        int exsscia[] = {sv->sample_info.alt_consensus1_split_score1_ind_aln, sv->sample_info.alt_consensus1_split_score2_ind_aln};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSCIA", exsscia, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSCIA", NULL, 0);
    }

    if (sv->sample_info.ext_alt_consensus2_length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS2", &(sv->sample_info.ext_alt_consensus2_to_alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS2", &(sv->sample_info.ext_alt_consensus2_to_ref_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL2", &(sv->sample_info.ext_alt_consensus2_length), 1);
        int exss2[] = {sv->sample_info.alt_consensus2_split_size1, sv->sample_info.alt_consensus2_split_size2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS2", exss2, 2);
        int exssc2[] = {sv->sample_info.alt_consensus2_split_score1, sv->sample_info.alt_consensus2_split_score2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC2", exssc2, 2);
        int exssc2ia[] = {sv->sample_info.alt_consensus2_split_score1_ind_aln, sv->sample_info.alt_consensus2_split_score2_ind_aln};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC2IA", exssc2ia, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC2IA", NULL, 0);
    }

    std::string filters;
    for (size_t i = 0; i < sv->sample_info.filters.size(); ++i) {
        filters += sv->sample_info.filters[i];
        if (i + 1 < sv->sample_info.filters.size()) {
            filters += ",";
        }
    }
    if (filters.empty()) {
        filters = "PASS";
    }

    char* filters_cstr = strdup(filters.c_str());
    bcf_update_format_char(out_hdr, sv->vcf_entry, "FT", filters_cstr, strlen(filters_cstr));
}

hts_pos_t get_covered_bps(std::vector<hts_pair_pos_t>& pos_pairs) {
    if (pos_pairs.empty()) {
        return 0;
    }
    
    std::sort(pos_pairs.begin(), pos_pairs.end(), [](const hts_pair_pos_t &a, const hts_pair_pos_t &b) {
        return a.beg < b.beg;
    });
    
    hts_pos_t total = 0;
    hts_pos_t current_start = pos_pairs[0].beg;
    hts_pos_t current_end = pos_pairs[0].end;
    
    for (size_t i = 1; i < pos_pairs.size(); ++i) {
        if (pos_pairs[i].beg <= current_end) {
            current_end = std::max(current_end, pos_pairs[i].end);
        } else {
            total += current_end - current_start;
            current_start = pos_pairs[i].beg;
            current_end = pos_pairs[i].end;
        }
    }
    total += current_end - current_start;
    
    return total;
}

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<std::shared_ptr<bam1_t>>& consistent_reads, 
    double consistent_avg_score, double consistent_stddev_score) {
    
    bp_reads_info.computed = true;
    bp_reads_info.reads = n_reads;

    std::vector<int> mqs;
    std::vector<hts_pair_pos_t> fwd_mate_positions, rev_mate_positions;
    std::vector<hts_pair_pos_t> fwd_hq_mate_positions, rev_hq_mate_positions;

    double sum_mq = 0;
    for (std::shared_ptr<bam1_t> read : consistent_reads) {
        int mq = get_mq(read.get());
        if (bam_is_mrev(read)) {
            bp_reads_info.consistent_fwd++;
            rev_mate_positions.push_back({read->core.mpos, get_mate_endpos(read.get())});
            if (mq >= config.high_confidence_mapq) {
                rev_hq_mate_positions.push_back({read->core.mpos, get_mate_endpos(read.get())});
            }
        } else {
            bp_reads_info.consistent_rev++;
            fwd_mate_positions.push_back({read->core.mpos, get_mate_endpos(read.get())});
            if (mq >= config.high_confidence_mapq) {
                fwd_hq_mate_positions.push_back({read->core.mpos, get_mate_endpos(read.get())});
            }
        }
        bp_reads_info.consistent_min_mq = std::min(bp_reads_info.consistent_min_mq, mq);
        bp_reads_info.consistent_max_mq = std::max(bp_reads_info.consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            bp_reads_info.consistent_high_mq++;
        }

        sum_mq += mq;
        mqs.push_back(mq);
    }
    bp_reads_info.consistent_avg_score = consistent_avg_score;
    bp_reads_info.consistent_stddev_score = consistent_stddev_score;

    bp_reads_info.consistent_avg_mq = sum_mq/consistent_reads.size();
    bp_reads_info.consistent_stddev_mq = stddev(mqs);

    bp_reads_info.fwd_mate_cov_bps = get_covered_bps(fwd_mate_positions);
    bp_reads_info.rev_mate_cov_bps = get_covered_bps(rev_mate_positions);
    bp_reads_info.fwd_hq_mate_cov_bps = get_covered_bps(fwd_hq_mate_positions);
    bp_reads_info.rev_hq_mate_cov_bps = get_covered_bps(rev_hq_mate_positions);
}

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs) {
    std::vector<std::string> temp1, temp2;
    std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
    
    std::vector<std::string> consensus_seqs; 
    
    consensus_seqs = generate_reference_guided_consensus(ref_seq, temp1, seqs, temp2, aligner, harsh_aligner, consensus_contigs_alns, config, stats, false);
    
    std::vector<seq_w_pp_t> seqs_w_pp, temp3, temp4;
    for (std::string& seq : seqs) {
        seqs_w_pp.push_back({seq, true, true});
    }
    std::vector<std::string> consensus_seqs2 = assemble_reads(temp3, seqs_w_pp, temp4, harsh_aligner, config, stats);
    consensus_seqs.insert(consensus_seqs.end(), consensus_seqs2.begin(), consensus_seqs2.end());
    return consensus_seqs;
}

std::vector<std::shared_ptr<bam1_t>> gen_consensus_and_find_consistent_seqs_subset(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score) {

    if (reads.empty()) {
        avg_score = 0;
        stddev_score = 0;
        consensus_seq = "";
        return reads;
    }

    if (revcomp_read.empty()) {
        revcomp_read.resize(reads.size(), false);
    }

    std::vector<std::string> seqs;
    for (int i = 0; i < reads.size(); i++) {
        std::shared_ptr<bam1_t> read = reads[i];
        std::string seq = get_sequence(read.get());
        if (revcomp_read[i]) rc(seq);
        seqs.push_back(seq);
    }

    avg_score = 0;
    std::vector<std::string> consensus_seqs = gen_consensus_seqs(ref_seq, seqs);

    std::vector<std::shared_ptr<bam1_t>> consistent_reads;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    std::vector<int> start_positions, end_positions;
    std::vector<StripedSmithWaterman::Alignment> alns;
    double cum_score = 0;
    std::vector<double> aln_scores;
    for (std::string cseq : consensus_seqs) {
        std::vector<std::shared_ptr<bam1_t>> curr_consistent_reads;
        std::vector<int> curr_start_positions, curr_end_positions;
        std::vector<StripedSmithWaterman::Alignment> curr_alns;
        double curr_cum_score = 0;
        std::vector<double> curr_aln_scores;
        for (int i = 0; i < reads.size(); i++) {
            std::shared_ptr<bam1_t> read = reads[i];
            std::string seq = get_sequence(read.get());
            if (revcomp_read[i]) rc(seq);

            harsh_aligner.Align(seq.c_str(), cseq.c_str(), cseq.length(), filter, &aln, 0);
            curr_alns.push_back(aln);

            double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
            if (mismatch_rate <= config.max_seq_error && !is_left_clipped(aln, config.min_clip_len) && !is_right_clipped(aln, config.min_clip_len)) {
                curr_consistent_reads.push_back(read);
                curr_cum_score += double(aln.sw_score)/seq.length();
                curr_aln_scores.push_back(double(aln.sw_score)/seq.length());
                curr_start_positions.push_back(aln.ref_begin);
                curr_end_positions.push_back(aln.ref_end);
            }
        }

        if (curr_cum_score > cum_score) {
            consistent_reads = curr_consistent_reads;
            cum_score = curr_cum_score;
            aln_scores = curr_aln_scores;
            start_positions = curr_start_positions;
            end_positions = curr_end_positions;
            alns = curr_alns;
            consensus_seq = cseq;
        }
    }
    std::sort(start_positions.begin(), start_positions.end());
    std::sort(end_positions.begin(), end_positions.end(), std::greater<int>());
    if (!consistent_reads.empty()) avg_score = cum_score/consistent_reads.size();
    else avg_score = 0;
    stddev_score = stddev(aln_scores);

    if (start_positions.size() >= 2) {
        correct_contig(consensus_seq, seqs, alns, config, true);

        consensus_seq = consensus_seq.substr(start_positions[1], end_positions[1]-start_positions[1]);

        int p;
        int window_len = 2*stats.read_len/3;
        while ((p = find_char_in_str(consensus_seq, 'N', 0, window_len)) != -1) {
            consensus_seq = consensus_seq.substr(p+1);
        }

        while ((p = find_char_in_str(consensus_seq, 'N', consensus_seq.length()-window_len, consensus_seq.length())) != -1) {
            consensus_seq = consensus_seq.substr(0, p);
        }
    } else {
        consensus_seq = "";
    }

    return consistent_reads;
}

std::vector<std::shared_ptr<bam1_t>> find_seqs_consistent_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, double& avg_score, double& stddev_score) {

    if (reads.empty()) {
        avg_score = 0;
        stddev_score = 0;
        return reads;
    }

    std::vector<std::shared_ptr<bam1_t>> consistent_reads;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    std::vector<StripedSmithWaterman::Alignment> alns;
    double cum_score = 0;
    std::vector<double> aln_scores;
    for (std::shared_ptr<bam1_t> read : reads) {
        std::string seq = get_sequence(read.get(), true);
        if (!bam_is_mrev(read)) rc(seq);
        
        harsh_aligner.Align(seq.c_str(), ref_seq.c_str(), ref_seq.length(), filter, &aln, 0);
        alns.push_back(aln);

        double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
        if (mismatch_rate <= config.max_seq_error && !is_left_clipped(aln, config.min_clip_len) && !is_right_clipped(aln, config.min_clip_len)) {
            consistent_reads.push_back(read);
            cum_score += double(aln.sw_score)/seq.length();
            aln_scores.push_back(double(aln.sw_score)/seq.length());
        }
    }

    avg_score = mean(aln_scores);
    stddev_score = stddev(aln_scores);

    return consistent_reads;
}

void read_mates(int contig_id) {
    mutex_per_chr[contig_id].lock();
    if (active_threads_per_chr[contig_id] == 0) {
		std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
		std::ifstream fin(fname);
		std::string qname, read_seq, qual; int mapq;
		while (fin >> qname >> read_seq >> qual >> mapq) {
			mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
		}
	}
	active_threads_per_chr[contig_id]++;
	mutex_per_chr[contig_id].unlock();
}

void release_mates(int contig_id) {
    mutex_per_chr[contig_id].lock();
	active_threads_per_chr[contig_id]--;
	if (active_threads_per_chr[contig_id] == 0) {
		mateseqs_w_mapq[contig_id].clear();
	}
	mutex_per_chr[contig_id].unlock();
}

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension) {
    int contig_id = contig_map.get_id(contig_name);
    candidate_reads_for_extension = get_extension_reads(contig_name, target_ivals, contig_len, config, stats, bam_file);
    std::vector<Interval<ext_read_t*>> it_ivals;
    for (ext_read_t* ext_read : candidate_reads_for_extension) {
        Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
        it_ivals.push_back(it_ival);
    }
    return IntervalTree<ext_read_t*>(it_ivals);
}

std::unordered_map<std::string, std::string> assign_reads(std::string in_vcf_fname) {
    
    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (in_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
    }

    bcf_hdr_t* in_vcf_header = bcf_hdr_read(in_vcf_file);
    if (in_vcf_header == NULL) {
        throw std::runtime_error("Failed to read the VCF header.");
    }

    std::unordered_map<std::string, float> sv_epr_map;

    bcf1_t* vcf_record = bcf_init();
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        bcf_unpack(vcf_record, BCF_UN_ALL);

        std::string id = vcf_record->d.id;
        float epr = get_sv_epr(in_vcf_header, vcf_record);
        sv_epr_map[id] = epr;
    }
    hts_close(in_vcf_file);
    bcf_hdr_destroy(in_vcf_header);

    std::string alt_reads_association_fname = workdir + "/alt_reads_to_sv_associations.txt";
    std::ifstream alt_reads_association_fin(alt_reads_association_fname);
    std::string sv_id, read_name, temp;
    std::unordered_map<std::string, std::string> read_to_sv_map;
    std::unordered_map<std::string, float> read_to_epr_map;
    while (alt_reads_association_fin >> sv_id >> read_name) {
        float epr = sv_epr_map[sv_id];
        if (epr > read_to_epr_map[read_name]) {
            read_to_epr_map[read_name] = epr; // Store the highest EPR for the read
            read_to_sv_map[read_name] = sv_id;
        }
    }

    return read_to_sv_map;
}

void clear_invalid_stat_tests(bcf_hdr_t* hdr, std::vector<std::shared_ptr<deletion_t>>& dels) {
    std::vector<std::pair<std::shared_ptr<deletion_t>, float>> dels_w_epr;
    hts_pos_t chr_len = 0;
    for (const auto& del : dels) {
        dels_w_epr.push_back({del, get_sv_epr(hdr, del->vcf_entry)});
        chr_len = std::max(chr_len, del->end);
    }

    std::sort(dels_w_epr.begin(), dels_w_epr.end(), 
    [](const std::pair<std::shared_ptr<deletion_t>, float>& a, const std::pair<std::shared_ptr<deletion_t>, float>& b) {
        return a.second > b.second; // Sort by EPR in descending order
    });

    SegTree seg_tree(chr_len+1);
    for (const auto& del_epr : dels_w_epr) {
        const auto& del = del_epr.first;
        hts_pos_t midpoint = (del->start + del->end) / 2;
        if (seg_tree.any_ge(midpoint-stats.max_is, midpoint, 1)) {
            del->ks_pval = deletion_t::KS_PVAL_NOT_COMPUTED;
            del->min_conf_size = deletion_t::SIZE_NOT_COMPUTED;
            del->max_conf_size = deletion_t::SIZE_NOT_COMPUTED;
        } else {
            seg_tree.add(midpoint-stats.max_is, midpoint, 1);
        }
    }
}

void rebalance_depth_cov(hts_pos_t strong_start, hts_pos_t strong_end, hts_pos_t weak_start, hts_pos_t weak_end, 
    int per_base_cov_delta, int& cov_to_update) {
    int ov = overlap(strong_start, strong_end, weak_start, weak_end);
    if (ov <= 0 || weak_end-weak_start <= 0) {
        return;
    }

    uint64_t weak_total_cov = (weak_end - weak_start) * cov_to_update;
    weak_total_cov += per_base_cov_delta * ov;
    cov_to_update = weak_total_cov / (weak_end - weak_start);
}

void rebalance_covs(bcf_hdr_t* hdr, std::vector<std::shared_ptr<deletion_t>>& dels) {
    std::vector<std::pair<std::shared_ptr<deletion_t>, float>> dels_w_epr;
    for (const auto& del : dels) {
        dels_w_epr.push_back({del, get_sv_epr(hdr, del->vcf_entry)});
    }

    std::vector<Interval<std::pair<std::shared_ptr<deletion_t>, float>>> it_ivals;
    for (const auto& del_epr : dels_w_epr) {
        const auto& del = del_epr.first;
        float epr = del_epr.second;
        Interval<std::pair<std::shared_ptr<deletion_t>, float>> it_ival(del->start, del->end, {del, epr});
        it_ivals.push_back(it_ival);
    }
    IntervalTree<std::pair<std::shared_ptr<deletion_t>, float>> it_tree(it_ivals);

    for (auto& curr_del : dels_w_epr) {
        auto ov_dels = it_tree.findOverlapping(curr_del.first->start, curr_del.first->end);
        for (const auto& it_del : ov_dels) {
            auto& ov_del = it_del.value;
            int ov_del_alt_ac = count_alt_alleles(hdr, ov_del.first->vcf_entry);
            if (ov_del.second <= curr_del.second || ov_del.first == curr_del.first || ov_del_alt_ac == 0) continue;

            deletion_t* strong_del = ov_del.first.get();
            deletion_t* weak_del = curr_del.first.get();
            int strong_del_alt_ac = ov_del_alt_ac;

            int strong_del_avg_flanking_cov = (strong_del->sample_info.left_flanking_cov + strong_del->sample_info.right_flanking_cov)/2;
            int strong_del_avg_indel_cov = (strong_del->sample_info.indel_left_cov + strong_del->sample_info.indel_right_cov)/2;
            int avg_depth_delta = std::max(0, (strong_del_avg_flanking_cov-strong_del_avg_indel_cov)/2 * strong_del_alt_ac);

            // int strong_del_avg_flanking_cov_hq = (strong_del->sample_info.left_flanking_cov_highmq + strong_del->sample_info.right_flanking_cov_highmq)/2;
            // int strong_del_avg_indel_cov_hq = (strong_del->sample_info.indel_left_cov_highmq + strong_del->sample_info.indel_right_cov_highmq)/2;
            // int avg_depth_delta_hq = std::max(0, (strong_del_avg_flanking_cov_hq-strong_del_avg_indel_cov_hq)/2 * strong_del_alt_ac);

            if (weak_del->end - weak_del->start <= config.indel_tested_region_size) {
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_left_cov);
                weak_del->sample_info.indel_right_cov = weak_del->sample_info.indel_left_cov;
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_left_cov_highmq);
                weak_del->sample_info.indel_right_cov_highmq = weak_del->sample_info.indel_left_cov_highmq;
            } else {
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->start+config.indel_tested_region_size, avg_depth_delta, weak_del->sample_info.indel_left_cov);
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->end-config.indel_tested_region_size, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_right_cov);
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->start+config.indel_tested_region_size, avg_depth_delta, weak_del->sample_info.indel_left_cov_highmq);
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->end-config.indel_tested_region_size, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_right_cov_highmq);
                
            }
            rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start-config.flanking_size, weak_del->start, avg_depth_delta, weak_del->sample_info.left_flanking_cov);
            rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->end, weak_del->end+config.flanking_size, avg_depth_delta, weak_del->sample_info.right_flanking_cov);
        }
    }
}

int main(int argc, char* argv[]) {

    std::string in_vcf_fname = argv[1];
    std::string out_vcf_fname = argv[2];
    bam_fname = argv[3];
    reference_fname = argv[4];
    workdir = argv[5];
    std::string sample_name = argv[6];

    std::unordered_map<std::string, std::string> reads_to_sv_map;
    bool reassign_evidence = false;
    if (argc > 7 && std::string(argv[7]) == "--reassign-evidence") {
        reassign_evidence = true;
        reads_to_sv_map = assign_reads(in_vcf_fname);
    }

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    chr_seqs.read_fasta_into_map(reference_fname);
    bam_pool = new bam_pool_t(config.threads, bam_fname, reference_fname);

    // read crossing isize distribution
    std::ifstream crossing_isizes_dist_fin(workdir + "/crossing_isizes.txt");
	int isize, count;
	while (crossing_isizes_dist_fin >> isize >> count) {
		for (int i = 0; i < count; i++) global_crossing_isize_dist.push_back(isize);
	}
	std::random_shuffle(global_crossing_isize_dist.begin(), global_crossing_isize_dist.end());
	global_crossing_isize_dist.resize(100000);
	crossing_isizes_dist_fin.close();

    std::string full_cmd_fname = workdir + "/genotype_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
    std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    mateseqs_w_mapq.resize(contig_map.size());
    active_threads_per_chr = std::vector<int>(contig_map.size());
	mutex_per_chr = std::vector<std::mutex>(contig_map.size());

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (in_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
    }

    bcf_hdr_t* in_vcf_header = bcf_hdr_read(in_vcf_file);
    if (in_vcf_header == NULL) {
        throw std::runtime_error("Failed to read the VCF header.");
    }

    bcf1_t* vcf_record = bcf_init();
    std::unordered_map<std::string, std::vector<std::shared_ptr<deletion_t>>> dels_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<duplication_t>>> dups_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<insertion_t>>> inss_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<inversion_t>>> invs_by_chr;
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(in_vcf_header, vcf_record);
        if (sv == nullptr) {
            std::cout << "Ignoring SV of unsupported type: " << vcf_record->d.id << std::endl; 
            continue;
        }

        sv->vcf_entry = bcf_dup(vcf_record);
        if (sv->svtype() == "DEL") {
            dels_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<deletion_t>(sv));
        } else if (sv->svtype() == "DUP") {
            if (sv->end-sv->start <= 0) {
                std::cout << "Discarding SV with invalid coordinates: " << sv->id << std::endl;
                continue;
            }
            dups_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<duplication_t>(sv));
        } else if (sv->svtype() == "INS") {
        	inss_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<insertion_t>(sv));
        } else if (sv->svtype() == "INV") {
            invs_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<inversion_t>(sv));
        }
    }

    int* imap = NULL;
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    bcf_hdr_t* out_vcf_header = bcf_subset_header(in_vcf_header, sample_name, imap);
    add_fmt_tags(out_vcf_header);
    if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
    	throw std::runtime_error("Failed to read the VCF header.");
    }

    evidence_logger_t* evidence_logger = NULL;
    if (!reassign_evidence) {
        evidence_logger = new evidence_logger_t(workdir);
    }

    // genotype chrs in descending order of svs
    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    const int BLOCK_SIZE = 20;

    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
    	std::string contig_name = contig_map.get_name(contig_id);
        std::vector<std::shared_ptr<deletion_t>>& dels = dels_by_chr[contig_name];
        for (int i = 0; i < dels.size(); i += BLOCK_SIZE) {
            std::vector<deletion_t*> block_dels;
            for (int j = i; j < std::min(i + BLOCK_SIZE, (int)dels.size()); j++) {
                block_dels.push_back(dels[j].get());
            }
            std::future<void> future = thread_pool.push(genotype_dels, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_dels, in_vcf_header, out_vcf_header, stats, config, 
                    contig_map, bam_pool, &mateseqs_w_mapq[contig_id], workdir, &global_crossing_isize_dist, evidence_logger, reassign_evidence, &reads_to_sv_map);
            futures.push_back(std::move(future));
        }

        std::vector<std::shared_ptr<duplication_t>>& dups = dups_by_chr[contig_name];
        for (int i = 0; i < dups.size(); i += BLOCK_SIZE) {
            std::vector<duplication_t*> block_dups;
            for (int j = i; j < std::min(i + BLOCK_SIZE, (int)dups.size()); j++) {
                block_dups.push_back(dups[j].get());
            }
            std::future<void> future = thread_pool.push(genotype_dups, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_dups, in_vcf_header, out_vcf_header, stats, config,
                    contig_map, bam_pool, &mateseqs_w_mapq[contig_id], workdir, &global_crossing_isize_dist, evidence_logger, reassign_evidence, &reads_to_sv_map);
            futures.push_back(std::move(future));
        }

        std::vector<std::shared_ptr<insertion_t>>& inss = inss_by_chr[contig_name];
        for (int i = 0; i < inss.size(); i += BLOCK_SIZE) {
            std::vector<insertion_t*> block_inss;
            for (int j = i; j < std::min(i + BLOCK_SIZE, (int)inss.size()); j++) {
                block_inss.push_back(inss[j].get());
            }
            std::future<void> future = thread_pool.push(genotype_inss, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_inss, in_vcf_header, out_vcf_header, stats, config,
                    contig_map, bam_pool, &mateseqs_w_mapq[contig_id], &global_crossing_isize_dist, evidence_logger, reassign_evidence, &reads_to_sv_map);
            futures.push_back(std::move(future));
        }

        std::vector<std::shared_ptr<inversion_t>>& invs = invs_by_chr[contig_name];
        for (int i = 0; i < invs.size(); i += BLOCK_SIZE) {
            std::vector<inversion_t*> block_invs;
            for (int j = i; j < std::min(i + BLOCK_SIZE, (int)invs.size()); j++) {
                block_invs.push_back(invs[j].get());
            }
            std::future<void> future = thread_pool.push(genotype_invs, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_invs, in_vcf_header, out_vcf_header, stats, config,
                    contig_map, bam_pool, &mateseqs_w_mapq[contig_id]);
            futures.push_back(std::move(future));
        }
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
        } catch (std::string s) {
            std::cerr << s << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
    futures.clear();

    if (reassign_evidence) {
        for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
            std::string contig_name = contig_map.get_name(contig_id);
            std::vector<std::shared_ptr<deletion_t>>& dels = dels_by_chr[contig_name];
            clear_invalid_stat_tests(in_vcf_header, dels);
            rebalance_covs(in_vcf_header, dels);
        }
    }

    // print contigs in vcf order
    int n_seqs;
    const char** seqnames = bcf_hdr_seqnames(in_vcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; i++) {
    	std::string contig_name = seqnames[i];
    	std::vector<std::shared_ptr<sv_t>> contig_svs;
    	if (dels_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dels_by_chr[contig_name].begin(), dels_by_chr[contig_name].end());
    	if (dups_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dups_by_chr[contig_name].begin(), dups_by_chr[contig_name].end());
    	if (inss_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), inss_by_chr[contig_name].begin(), inss_by_chr[contig_name].end());
        if (invs_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), invs_by_chr[contig_name].begin(), invs_by_chr[contig_name].end());
    	std::sort(contig_svs.begin(), contig_svs.end(), [](const std::shared_ptr<sv_t>& sv1, const std::shared_ptr<sv_t>& sv2) {return sv1->start < sv2->start;});

		for (auto& sv : contig_svs) {
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AC", NULL, 0);
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AN", NULL, 0);
            update_record(in_vcf_header, out_vcf_header, sv.get(), chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), imap[0]);
			if (bcf_write(out_vcf_file, out_vcf_header, sv->vcf_entry) != 0) {
				throw std::runtime_error("Failed to write VCF record to " + out_vcf_fname);
			}
		}
    }
    delete[] imap;
    free(seqnames);

    bcf_destroy(vcf_record);
    bcf_hdr_destroy(out_vcf_header);
    bcf_hdr_destroy(in_vcf_header);
    bcf_close(in_vcf_file);
    bcf_close(out_vcf_file);
    delete bam_pool;
    delete evidence_logger;
}
