#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

#include "sw_utils.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "assemble.h"
#include "extend_1sr_consensus.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "vcf_utils.h"
#include "stat_tests.h"

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

    base_frequencies_t left_anchor_base_freqs = get_base_frequencies(chr_seq+sv->left_anchor_aln->start, sv->left_anchor_aln->end-sv->left_anchor_aln->start);
	int labc[] = {left_anchor_base_freqs.a, left_anchor_base_freqs.c, left_anchor_base_freqs.g, left_anchor_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "LEFT_ANCHOR_BASE_COUNT", labc, 4);

	base_frequencies_t right_anchor_base_freqs = get_base_frequencies(chr_seq+sv->right_anchor_aln->start, sv->right_anchor_aln->end-sv->right_anchor_aln->start);
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
    bcf_update_genotypes(out_hdr, sv->vcf_entry, sv->regenotyping_info.gt, sv->regenotyping_info.n_gt);

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1", &(sv->regenotyping_info.alt_bp1_better_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARC1", &(sv->regenotyping_info.alt_bp1_better_consistent_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARCF1", &(sv->regenotyping_info.alt_bp1_better_consistent_reads_fwd), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARCR1", &(sv->regenotyping_info.alt_bp1_better_consistent_reads_rev), 1);
    float arcas1 = sv->regenotyping_info.alt_bp1_better_consistent_avg_score;
    bcf_update_format_float(out_hdr, sv->vcf_entry, "ARCAS1", &arcas1, 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARC1HQ", &(sv->regenotyping_info.alt_bp1_better_consistent_high_mq), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARD1", &(sv->regenotyping_info.alt_bp1_better_disc_pairs), 1);
    
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR2", &(sv->regenotyping_info.alt_bp2_better_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARC2", &(sv->regenotyping_info.alt_bp2_better_consistent_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARCF2", &(sv->regenotyping_info.alt_bp2_better_consistent_reads_fwd), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARCR2", &(sv->regenotyping_info.alt_bp2_better_consistent_reads_rev), 1);
    float arcas2 = sv->regenotyping_info.alt_bp2_better_consistent_avg_score;
    bcf_update_format_float(out_hdr, sv->vcf_entry, "ARCAS2", &arcas2, 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARC2HQ", &(sv->regenotyping_info.alt_bp2_better_consistent_high_mq), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARD2", &(sv->regenotyping_info.alt_bp2_better_disc_pairs), 1);

    int arcmq[] = {sv->regenotyping_info.alt_bp1_better_consistent_max_mq, sv->regenotyping_info.alt_bp2_better_consistent_max_mq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARCMQ", arcmq, 2);

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1", &(sv->regenotyping_info.ref_bp1_better_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR2", &(sv->regenotyping_info.ref_bp2_better_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "RRC1", &(sv->regenotyping_info.ref_bp1_better_consistent_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "RRC2", &(sv->regenotyping_info.ref_bp2_better_consistent_reads), 1);

    int er = sv->regenotyping_info.alt_ref_equal_reads;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ER", &er, 1);

    if (sv->regenotyping_info.too_deep) {
        int td = 1;
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "TD", &td, 1);
    }

    int median_depths[] = {sv->median_left_flanking_cov, sv->median_indel_left_cov, sv->median_indel_right_cov, sv->median_right_flanking_cov};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MD", median_depths, 4);

    int median_depths_highmq[] = {sv->median_left_flanking_cov_highmq, sv->median_indel_left_cov_highmq, sv->median_indel_right_cov_highmq, sv->median_right_flanking_cov_highmq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDHQ", median_depths_highmq, 4);

    int cluster_depths[] = {sv->median_left_cluster_cov, sv->median_right_cluster_cov};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CLMD", cluster_depths, 2);

    int cluster_depths_highmq[] = {sv->median_left_cluster_cov_highmq, sv->median_right_cluster_cov_highmq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CLMDHQ", cluster_depths_highmq, 2);

    if (sv->min_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZE", &(sv->min_conf_size), 1);
    }
    if (sv->max_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZE", &(sv->max_conf_size), 1);
    }
    if (sv->ks_pval != deletion_t::KS_PVAL_NOT_COMPUTED) {
        float ks_pval = sv->ks_pval;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVAL", &ks_pval, 1);
    }

    int dp[] = {sv->disc_pairs_lf, sv->disc_pairs_rf};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP", dp, 2);

    int disc_pairs_surr[] = {sv->l_cluster_region_disc_pairs, sv->r_cluster_region_disc_pairs};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPS", disc_pairs_surr, 2);
    
    int disc_pairs_surr_hq[] = {sv->l_cluster_region_disc_pairs_high_mapq, sv->r_cluster_region_disc_pairs_high_mapq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPSHQ", disc_pairs_surr_hq, 2);
    
    if (sv->disc_pairs_lf + sv->disc_pairs_rf > 0) {
        int dphq[] = {sv->disc_pairs_lf_high_mapq, sv->disc_pairs_rf_high_mapq};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPHQ", dphq, 2);
        
        int dpmq[] = {sv->disc_pairs_lf_maxmapq, sv->disc_pairs_rf_maxmapq};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPMQ", dpmq, 2);

        float dpnm[] = {(float) sv->disc_pairs_lf_avg_nm, (float) sv->disc_pairs_rf_avg_nm};
        bcf_update_format_float(out_hdr, sv->vcf_entry, "DPNM", &dpnm, 2);

        int dpsp[] = {sv->disc_pairs_lf_span, sv->disc_pairs_rf_span};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPSP", dpsp, 2);
    }

    int cp[] = {sv->conc_pairs_lbp, sv->conc_pairs_midp, sv->conc_pairs_rbp};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CP", cp, 3);
    int cphq[] = {sv->conc_pairs_lbp_high_mapq, sv->conc_pairs_midp_high_mapq, sv->conc_pairs_rbp_high_mapq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CPHQ", cphq, 3);

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXR", &(sv->regenotyping_info.alt_ext_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXRHQ", &(sv->regenotyping_info.hq_alt_ext_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL", &(sv->regenotyping_info.ext_alt_consensus_length), 1);
    if (sv->regenotyping_info.ext_alt_consensus_length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS", &(sv->regenotyping_info.ext_alt_consensus_to_alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS", &(sv->regenotyping_info.ext_alt_consensus_to_ref_score), 1);
        
        int exss[] = {sv->regenotyping_info.alt_consensus_split_size1, sv->regenotyping_info.alt_consensus_split_size2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS", exss, 2);
    }
}

void reset_stats(sv_t* sv) {
    sv->median_left_flanking_cov = 0;
    sv->median_indel_left_cov = 0;
    sv->median_indel_right_cov = 0;
    sv->median_right_flanking_cov = 0;
    sv->median_left_flanking_cov_highmq = 0;
    sv->median_indel_left_cov_highmq = 0;
    sv->median_indel_right_cov_highmq = 0;
    sv->median_right_flanking_cov_highmq = 0;
    sv->median_left_cluster_cov = 0;
    sv->median_right_cluster_cov = 0;
    sv->disc_pairs_lf = 0;
    sv->disc_pairs_rf = 0;
    sv->disc_pairs_lf_high_mapq = 0;
    sv->disc_pairs_rf_high_mapq = 0;
    sv->disc_pairs_lf_maxmapq = 0;
    sv->disc_pairs_rf_maxmapq = 0;
    sv->disc_pairs_lf_avg_nm = 0;
    sv->disc_pairs_rf_avg_nm = 0;
    sv->conc_pairs_lbp = 0;
    sv->conc_pairs_midp = 0;
    sv->conc_pairs_rbp = 0;
    sv->l_cluster_region_disc_pairs = 0;
    sv->r_cluster_region_disc_pairs = 0;
}

std::vector<bam1_t*> find_consistent_seqs_subset(std::string ref_seq, std::vector<bam1_t*>& reads, std::string& consensus_seq, double& avg_score) {
    std::vector<std::string> seqs;
    for (bam1_t* read : reads) {
        std::string seq = get_sequence(read, true);
        if (!bam_is_mrev(read)) rc(seq);
        seqs.push_back(seq);
    }

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    std::vector<bam1_t*> consistent_reads;
    avg_score = 0;
    if (!reads.empty()) {
        std::vector<seq_w_pp_t> seqs_w_pp, temp1, temp2;
        for (std::string& seq : seqs) {
            seqs_w_pp.push_back({seq, true, true});
        }
        std::vector<std::string> consensuses = assemble_reads(temp1, seqs_w_pp, temp2, harsh_aligner, config, stats);

        if (consensuses.empty()) {
            return consistent_reads;
        }

        consensus_seq = consensuses[0];
        std::vector<int> start_positions, end_positions;
        for (bam1_t* read : reads) {
            std::string seq = get_sequence(read, true);
            if (!bam_is_mrev(read)) rc(seq);
            aligner.Align(seq.c_str(), consensus_seq.c_str(), consensus_seq.length(), filter, &aln, 0);
            if (!is_left_clipped(aln, config.min_clip_len) && !is_right_clipped(aln, config.min_clip_len)) {
                consistent_reads.push_back(read);
                avg_score += double(aln.sw_score)/seq.length();
                start_positions.push_back(aln.ref_begin);
                end_positions.push_back(aln.ref_end);
            }
        }
        std::sort(start_positions.begin(), start_positions.end());
        std::sort(end_positions.begin(), end_positions.end());
        if (!consistent_reads.empty()) avg_score /= consistent_reads.size();

        if (start_positions.size() > 2) {
            int start = start_positions[2];
            int end = end_positions[start_positions.size()-3];
            consensus_seq = consensus_seq.substr(start, end-start);
        } else {
            consensus_seq = "";
        }
    }
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
    candidate_reads_for_extension = get_extension_reads(contig_name, target_ivals, contig_len, stats, bam_file);
    std::vector<Interval<ext_read_t*>> it_ivals;
    for (ext_read_t* ext_read : candidate_reads_for_extension) {
        Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
        it_ivals.push_back(it_ival);
    }
    return IntervalTree<ext_read_t*>(it_ivals);
}

void genotype_del(deletion_t* del, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
    int del_start = del->start, del_end = del->end;

    hts_pos_t extend = stats.read_len + 20;

    // build alt allele
    /* POS in VCF is the base BEFORE the deletion - i.e., the first deleted base is POS+1.
     * Therefore, we want the ALT allele to *include* base POS
     * (note that POS is 1-based in the VCF file, but htslib kindly returns the 0-based coordinate here).
     * As for the END coordinate, my current understanding (which may change) is that it represents the last base deleted.
     * Therefore, the ALT allele should NOT include base END, i.e. it should start at END+1.
     * Here we shift both coordinates by 1, to make them the base immediately AFTER the breakpoints, which is a bit more intuitive for me. */
    del_start++; del_end++;

    char* contig_seq = chr_seqs.get_seq(del->chr);
    hts_pos_t contig_len = chr_seqs.get_len(del->chr);

    // all ranges will be start-inclusive and end-exclusive, i.e. [a,b)
    hts_pos_t alt_start = std::max(hts_pos_t(0), del_start-extend);
    hts_pos_t alt_end = std::min(del_end+extend, contig_len);
    int alt_lh_len = del_start-alt_start, alt_rh_len = alt_end-del_end;
    int alt_len = alt_lh_len + del->ins_seq.length() + alt_rh_len;
    char* alt_seq = new char[alt_len + 1];
    strncpy(alt_seq, contig_seq+alt_start, alt_lh_len);
    strncpy(alt_seq+alt_lh_len, del->ins_seq.c_str(), del->ins_seq.length());
    strncpy(alt_seq+alt_lh_len+del->ins_seq.length(), contig_seq+del_end, alt_rh_len);
    alt_seq[alt_len] = 0;

    // extract ref alleles - will be useful for consensus generation
    int ref_bp1_start = alt_start, ref_bp1_end = std::min(del_start+extend, contig_len);
    int ref_bp1_len = ref_bp1_end - ref_bp1_start;
    char* ref_bp1_seq = new char[ref_bp1_len + 1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;

    int ref_bp2_start = std::max(hts_pos_t(0), del_end-extend), ref_bp2_end = alt_end;
    int ref_bp2_len = ref_bp2_end - ref_bp2_start;
    char* ref_bp2_seq = new char[ref_bp2_len + 1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;

    std::stringstream l_region, r_region;
    l_region << del->chr << ":" << alt_start << "-" << ref_bp1_end;
    r_region << del->chr << ":" << ref_bp2_start << "-" << alt_end;
    
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int same = 0;
    std::vector<bam1_t*> alt_better_reads, ref_bp1_better_seqs, ref_bp2_better_seqs;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < del_start || del_end < get_unclipped_start(read)) continue;
        if (del_start < get_unclipped_start(read) && get_unclipped_end(read) < del_end) continue;

        std::string seq;

        if (!is_samechr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < del_start-stats.max_is) continue; // positive mate and after deletion, potentially discordant...
            // however, there is an exception: both reads in a pair can be left-clipped on the end-side of the deletion (we allow 5bp tolerance)
            if (read->core.mpos > del_start && (abs(read->core.pos-read->core.mpos) > 5 || !is_left_clipped(read, config.min_clip_len))) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > del_end+stats.max_is) continue; // negative mate and before deletion, potentially discordant...
            // however, there is an exception: both reads in a pair can be right-clipped on the start-side of the deletion (we allow 5bp tolerance)
            if (mate_endpos < del_end && (abs(mate_endpos-bam_endpos(read)) > 5 || !is_right_clipped(read, config.min_clip_len))) continue;
            seq = get_sequence(read, true);
        }

        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

        // align to REF (two breakpoints)
        uint16_t ref_aln_score = 0;
        bool increase_ref_bp1_better = false, increase_ref_bp2_better = false;
        if (is_perfectly_aligned(read)) {
            ref_aln_score = read->core.l_qseq;
            if (read->core.pos < del_start && bam_endpos(read) > del_start) {
                increase_ref_bp1_better = true;
            }
            if (read->core.pos < del_end && bam_endpos(read) > del_end) {
                increase_ref_bp2_better = true;
            }
        } else {
            aligner.Align(seq.c_str(), ref_bp1_seq, ref_bp1_len, filter, &ref1_aln, 0);
            aligner.Align(seq.c_str(), ref_bp2_seq, ref_bp2_len, filter, &ref2_aln, 0);
            ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                increase_ref_bp1_better = true;
            }
            if (ref2_aln.sw_score >= ref1_aln.sw_score) {
                increase_ref_bp2_better = true;
            }
        }

        if (alt_aln.sw_score > ref_aln_score) {
            alt_better_reads.push_back(bam_dup1(read));
        } else if (alt_aln.sw_score < ref_aln_score) {
            if (increase_ref_bp1_better) {
                ref_bp1_better_seqs.push_back(bam_dup1(read));
            }
            if (increase_ref_bp2_better) {
                ref_bp2_better_seqs.push_back(bam_dup1(read));
            }
        } else {
            same++;
        }

        if (alt_better_reads.size() + ref_bp1_better_seqs.size() + ref_bp2_better_seqs.size() + same > 4 * stats.get_max_depth(del->chr)) {
            alt_better_reads.clear();
            ref_bp1_better_seqs.clear();
            ref_bp2_better_seqs.clear();
            same = 0;
            del->regenotyping_info.too_deep = true;
            break;
        }
    }

    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seq, alt_better_reads, alt_consensus_seq, alt_avg_score);
    std::vector<bam1_t*> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score);
    std::vector<bam1_t*> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score);

    if (!alt_consensus_seq.empty()) {
       // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, del->start-stats.max_is, del->start, del->chr, chr_seqs.get_len(del->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, del->end, del->end+stats.max_is, del->chr, chr_seqs.get_len(del->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        del->regenotyping_info.alt_ext_reads = alt_consensus->left_ext_reads + alt_consensus->right_ext_reads;
        del->regenotyping_info.hq_alt_ext_reads = alt_consensus->hq_left_ext_reads + alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        delete alt_consensus;

        hts_pos_t lh_start = del->start-alt_consensus_seq.length();
        if (lh_start < 0) lh_start = 0;
        hts_pos_t lh_len = del->start-lh_start;
        hts_pos_t rh_end = del->end+alt_consensus_seq.length();
        if (rh_end > contig_len) rh_end = contig_len;
        hts_pos_t rh_len = rh_end-del->end;
    
        delete[] alt_seq;
        alt_seq = new char[2*alt_consensus_seq.length() + 1];
        strncpy(alt_seq, contig_seq+lh_start, lh_len);
        strncpy(alt_seq+lh_len, contig_seq+del_end, rh_len);
        alt_seq[lh_len+rh_len] = 0;

        // align to ref+SV
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, lh_len+rh_len, filter, &alt_aln, 0);

        // length of the left and right flanking regions of the deletion covered by alt_consensus_seq
        int lf_aln_rlen = std::max(hts_pos_t(0), lh_len - alt_aln.ref_begin);
        int rf_aln_rlen = std::max(hts_pos_t(0), alt_aln.ref_end - lh_len);

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        int temp;
        int lf_aln_qlen = query_prefix_len(alt_aln.cigar, lf_aln_rlen, temp);
        int rf_aln_qlen = query_suffix_len(alt_aln.cigar, rf_aln_rlen, temp);
        del->regenotyping_info.alt_consensus_split_size1 = lf_aln_qlen;
        del->regenotyping_info.alt_consensus_split_size2 = rf_aln_qlen;

        // align to ref
        hts_pos_t lbp_start = lh_start, lbp_end = del->start + alt_consensus_seq.length();
        hts_pos_t rbp_start = del->end - alt_consensus_seq.length(), rbp_end = rh_end;
        if (lbp_end > contig_len) lbp_end = contig_len;
        if (rbp_start < 0) rbp_start = 0;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+lbp_start, lbp_end-lbp_start, filter, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+rbp_start, rbp_end-rbp_start, filter, &ref2_aln, 0);
    
        del->regenotyping_info.ext_alt_consensus_length = alt_consensus_seq.length();
        del->regenotyping_info.ext_alt_consensus_to_alt_score = alt_aln.sw_score;
        del->regenotyping_info.ext_alt_consensus_to_ref_score = std::max(ref1_aln.sw_score, ref2_aln.sw_score);
    }

    del->regenotyping_info.alt_bp1_better_reads = alt_better_reads.size();
    del->regenotyping_info.alt_bp2_better_reads = alt_better_reads.size();
    del->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_reads_consistent.size();
    del->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_reads_consistent.size();
    for (bam1_t* read : alt_better_reads_consistent) {
        if (bam_is_mrev(read)) {
            del->regenotyping_info.alt_bp1_better_consistent_reads_fwd++;
        } else {
            del->regenotyping_info.alt_bp1_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        del->regenotyping_info.alt_bp1_better_consistent_max_mq = std::max(del->regenotyping_info.alt_bp1_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            del->regenotyping_info.alt_bp1_better_consistent_high_mq++;
        }
    }
    del->regenotyping_info.alt_bp2_better_consistent_reads_fwd = del->regenotyping_info.alt_bp1_better_consistent_reads_fwd;
    del->regenotyping_info.alt_bp2_better_consistent_reads_rev = del->regenotyping_info.alt_bp1_better_consistent_reads_rev;
    del->regenotyping_info.alt_bp2_better_consistent_max_mq = del->regenotyping_info.alt_bp1_better_consistent_max_mq;
    del->regenotyping_info.alt_bp2_better_consistent_high_mq = del->regenotyping_info.alt_bp1_better_consistent_high_mq;
    del->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_avg_score;
    del->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_avg_score;
    del->regenotyping_info.ref_bp1_better_reads = ref_bp1_better_seqs.size();
    del->regenotyping_info.ref_bp2_better_reads = ref_bp2_better_seqs.size();
    del->regenotyping_info.ref_bp1_better_consistent_reads = ref_bp1_better_seqs_consistent.size();
    del->regenotyping_info.ref_bp2_better_consistent_reads = ref_bp2_better_seqs_consistent.size();
    del->regenotyping_info.alt_ref_equal_reads = same;

    delete[] alt_seq;
    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    free(regions[0]);
    free(regions[1]);

    for (bam1_t* read : alt_better_reads) bam_destroy1(read);
    for (bam1_t* read : ref_bp1_better_seqs) bam_destroy1(read);
    for (bam1_t* read : ref_bp2_better_seqs) bam_destroy1(read);
    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_dels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<deletion_t*> dels,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (deletion_t* del : dels) {
        target_ivals.push_back({del->start-stats.max_is, del->start+stats.max_is});
        target_ivals.push_back({del->end-stats.max_is, del->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_file, candidate_reads_for_extension);

    std::vector<deletion_t*> small_deletions, large_deletions;       
    std::vector<sv_t*> small_svs;         
    for (deletion_t* del : dels) {
        genotype_del(del, bam_file, candidate_reads_for_extension_itree, mateseqs_w_mapq[contig_id]);
        if (-del->svlen() >= stats.max_is) {
            large_deletions.push_back(del);
        } else {
            small_deletions.push_back(del);
            small_svs.push_back(del);
        }
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    release_mates(contig_id);

    depth_filter_del(contig_name, dels, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, global_crossing_isize_dist, small_svs, bam_file, config, stats, config.min_sv_size, true);
    std::string mates_nms_file = workdir + "/workspace/long-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dels, bam_file, config, stats, true, mates_nms_file);
    calculate_cluster_region_disc(contig_name, dels, bam_file, config);
}

void genotype_small_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {

	hts_pos_t dup_start = dup->start, dup_end = dup->end;
	hts_pos_t contig_len = chr_seqs.get_len(dup->chr);

	hts_pos_t extend = stats.read_len + 20;
	hts_pos_t svlen = dup->svlen();

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;

    hts_pos_t ref_start = std::max(hts_pos_t(0), dup_start-extend), ref_end = std::min(dup_end+extend, contig_len);
	hts_pos_t ref_len = ref_end - ref_start;
	char* ref_seq = new char[ref_len + 1];
	strncpy(ref_seq, chr_seqs.get_seq(dup->chr)+ref_start, ref_len);
	ref_seq[ref_len] = 0;

    std::vector<char*> alt_seqs;
	for (int copies = 1; copies*svlen < stats.read_len; copies++) {
		int alt_len = ref_len + copies*svlen;

		char* alt_seq = new char[alt_len+1];
		int pos = 0;
		strncpy(alt_seq, chr_seqs.get_seq(dup->chr)+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < copies; i++) {
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
			strncpy(alt_seq+pos, chr_seqs.get_seq(dup->chr)+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
		}
		strncpy(alt_seq+pos, chr_seqs.get_seq(dup->chr)+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
		alt_seqs.push_back(alt_seq);
	}

    std::stringstream region;
    region << dup->chr << ":" << ref_start << "-" << ref_end;

	hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region.str().c_str());
	
    bam1_t* read = bam_init1();

    std::vector<bam1_t*> ref_better_reads;
    std::vector<std::vector<bam1_t*>> alt_better_reads(alt_seqs.size());
    int same = 0;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;

        std::string seq;

        if (!is_samechr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < dup_start-stats.max_is) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            if (get_mate_endpos(read) > dup_end+stats.max_is) continue;
            seq = get_sequence(read, true);
        }

        aligner.Align(seq.c_str(), ref_seq, ref_len, filter, &ref_aln, 0);

        uint16_t best_aln_score = 0;
        std::vector<uint16_t> alt_aln_scores(alt_seqs.size());
        for (int i = 0; i < alt_seqs.size(); i++) {
            aligner.Align(seq.c_str(), alt_seqs[i], strlen(alt_seqs[i]), filter, &alt_aln, 0);
            alt_aln_scores[i] = alt_aln.sw_score;
            if (alt_aln.sw_score > best_aln_score) {
                best_aln_score = alt_aln.sw_score;
            }
        }

        if (best_aln_score > ref_aln.sw_score) {
            for (int i = 0; i < alt_seqs.size(); i++) {
                if (alt_aln_scores[i] == best_aln_score) {
                    alt_better_reads[i].push_back(bam_dup1(read));
                }
            }
        } else if (best_aln_score < ref_aln.sw_score) {
            ref_better_reads.push_back(bam_dup1(read));
        } else {
            same++;
        }
    }

    int alt_with_most_reads = 0;
    for (int i = 1; i < alt_better_reads.size(); i++) {
        if (alt_better_reads[i].size() > alt_better_reads[alt_with_most_reads].size()) {
            alt_with_most_reads = i;
        }
    }

    if (alt_better_reads[alt_with_most_reads].size() > 20*stats.get_max_depth(dup->chr) || ref_better_reads.size() > 20*stats.get_max_depth(dup->chr)) {
        alt_better_reads[alt_with_most_reads].clear();
        ref_better_reads.clear();
        same = 0;
        dup->regenotyping_info.too_deep = true;
    }

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seqs[alt_with_most_reads], alt_better_reads[alt_with_most_reads], alt_consensus_seq, alt_avg_score);
    std::vector<bam1_t*> ref_better_reads_consistent = find_consistent_seqs_subset(ref_seq, ref_better_reads, ref_consensus_seq, ref_avg_score);

    if (!alt_consensus_seq.empty()) {
       // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, dup->start-stats.max_is, dup->start, dup->chr, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->end, dup->end+stats.max_is, dup->chr, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        dup->regenotyping_info.alt_ext_reads = alt_consensus->left_ext_reads + alt_consensus->right_ext_reads;
        dup->regenotyping_info.hq_alt_ext_reads = alt_consensus->hq_left_ext_reads + alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        delete alt_consensus;

        hts_pos_t ref_start = dup->start-alt_consensus_seq.length();
        if (ref_start < 0) ref_start = 0;
        hts_pos_t ref_end = dup->end+alt_consensus_seq.length();
        if (ref_end > contig_len) ref_end = contig_len;
        aligner.Align(alt_consensus_seq.c_str(), chr_seqs.get_seq(dup->chr)+ref_start, ref_end-ref_start, filter, &ref_aln, 0);
        dup->regenotyping_info.ext_alt_consensus_to_ref_score = ref_aln.sw_score;

        int n_extra_copies = alt_with_most_reads+1;
        int alt_len = ref_end - ref_start + n_extra_copies*svlen;
		char* alt_seq = new char[alt_len+1];
		int pos = 0;
		strncpy(alt_seq, chr_seqs.get_seq(dup->chr)+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < n_extra_copies; i++) {
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
			strncpy(alt_seq+pos, chr_seqs.get_seq(dup->chr)+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
		}
		strncpy(alt_seq+pos, chr_seqs.get_seq(dup->chr)+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);
        dup->regenotyping_info.ext_alt_consensus_to_alt_score = alt_aln.sw_score;

        delete[] alt_seq;

        dup->regenotyping_info.ext_alt_consensus_length = alt_consensus_seq.length();
    }

    dup->regenotyping_info.alt_bp1_better_reads = alt_better_reads[alt_with_most_reads].size();
    dup->regenotyping_info.alt_bp2_better_reads = alt_better_reads[alt_with_most_reads].size();
    dup->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_reads_consistent.size();
    dup->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_reads_consistent.size();
    for (bam1_t* read : alt_better_reads_consistent) {
        if (bam_is_mrev(read)) {
            dup->regenotyping_info.alt_bp1_better_consistent_reads_fwd++;
        } else {
            dup->regenotyping_info.alt_bp1_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        dup->regenotyping_info.alt_bp1_better_consistent_max_mq = std::max(dup->regenotyping_info.alt_bp1_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            dup->regenotyping_info.alt_bp1_better_consistent_high_mq++;
        }
    }
    dup->regenotyping_info.alt_bp2_better_consistent_reads_fwd = dup->regenotyping_info.alt_bp1_better_consistent_reads_fwd;
    dup->regenotyping_info.alt_bp2_better_consistent_reads_rev = dup->regenotyping_info.alt_bp1_better_consistent_reads_rev;
    dup->regenotyping_info.alt_bp2_better_consistent_max_mq = dup->regenotyping_info.alt_bp1_better_consistent_max_mq;
    dup->regenotyping_info.alt_bp2_better_consistent_high_mq = dup->regenotyping_info.alt_bp1_better_consistent_high_mq;
    dup->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.ref_bp1_better_reads = ref_better_reads.size();
    dup->regenotyping_info.ref_bp2_better_reads = ref_better_reads.size();
    dup->regenotyping_info.ref_bp1_better_consistent_reads = ref_better_reads_consistent.size();
    dup->regenotyping_info.ref_bp2_better_consistent_reads = ref_better_reads_consistent.size();
    dup->regenotyping_info.alt_ref_equal_reads = same;

    delete[] ref_seq;
    for (char* alt_seq : alt_seqs) {
        delete[] alt_seq;
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_large_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
    
    hts_pos_t dup_start = dup->start, dup_end = dup->end;

    char* contig_seq = chr_seqs.get_seq(dup->chr);
	hts_pos_t contig_len = chr_seqs.get_len(dup->chr);

	hts_pos_t extend = stats.read_len + 20;

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;

	// all ranges will be start-inclusive and end-exclusive, i.e. [a,b)

	hts_pos_t ref_bp1_start = std::max(hts_pos_t(0), dup_start-extend), ref_bp1_end = std::min(dup_start+extend, contig_len);
	hts_pos_t ref_bp1_len = ref_bp1_end - ref_bp1_start;
	hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), dup_end-extend), ref_bp2_end = std::min(dup_end+extend, contig_len);
	hts_pos_t ref_bp2_len = ref_bp2_end - ref_bp2_start;

	// build alt allele
	hts_pos_t alt_lh_len = dup_end - ref_bp2_start;
	hts_pos_t alt_rh_len = ref_bp1_end - dup_start;
	hts_pos_t alt_len = alt_lh_len + dup->ins_seq.length() + alt_rh_len;
	char* alt_seq = new char[alt_len + 1];
	strncpy(alt_seq, contig_seq+ref_bp2_start, alt_lh_len);
    strncpy(alt_seq+alt_lh_len, dup->ins_seq.c_str(), dup->ins_seq.length());
    strncpy(alt_seq+alt_lh_len+dup->ins_seq.length(), contig_seq+dup_start, alt_rh_len);
	alt_seq[alt_len] = 0;

    std::stringstream l_region, r_region;
    l_region << dup->chr << ":" << ref_bp1_start << "-" << ref_bp1_end;
    r_region << dup->chr << ":" << ref_bp2_start << "-" << ref_bp2_end;
    
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    std::vector<bam1_t*> alt_better_reads, ref_bp1_better_reads, ref_bp2_better_reads;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;

        std::string seq;
        
        if (!is_samechr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < dup_start-stats.max_is) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            if (get_mate_endpos(read) > dup_end+stats.max_is) continue;
            seq = get_sequence(read, true);
        }
        
        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

        // align to REF (two breakpoints)
        uint16_t ref_aln_score = 0;
        bool increase_ref_bp1_better = false, increase_ref_bp2_better = false;
        if (is_perfectly_aligned(read)) {
            ref_aln_score = read->core.l_qseq;
            if (read->core.pos < dup_start && bam_endpos(read) > dup_start) {
                increase_ref_bp1_better = true;
            }
            if (read->core.pos < dup_end && bam_endpos(read) > dup_end) {
                increase_ref_bp2_better = true;
            }
        } else {
            aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln, 0);
            aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln, 0);
            ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                increase_ref_bp1_better = true;
            }
            if (ref2_aln.sw_score >= ref1_aln.sw_score) {
                increase_ref_bp2_better = true;
            }
        }

        if (alt_aln.sw_score > ref_aln_score) {
            alt_better_reads.push_back(bam_dup1(read));
        } else {
            if (increase_ref_bp1_better) {
                ref_bp1_better_reads.push_back(bam_dup1(read));
            }
            if (increase_ref_bp2_better) {
                ref_bp2_better_reads.push_back(bam_dup1(read));
            }
        }

        if (ref_bp1_better_reads.size() + ref_bp2_better_reads.size() > 4*stats.get_max_depth(dup->chr)) {
            alt_better_reads.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            dup->regenotyping_info.too_deep = true;
            break;
        }
    }

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;

    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seq, alt_better_reads, alt_consensus_seq, alt_avg_score);
    std::vector<bam1_t*> ref_bp1_better_reads_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_reads, ref_bp1_consensus_seq, ref_bp1_avg_score);
    std::vector<bam1_t*> ref_bp2_better_reads_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_reads, ref_bp2_consensus_seq, ref_bp2_avg_score);

    if (!alt_consensus_seq.empty()) {
       // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, dup->end-stats.max_is, dup->end, dup->chr, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->start, dup->start+stats.max_is, dup->chr, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        dup->regenotyping_info.alt_ext_reads = alt_consensus->left_ext_reads + alt_consensus->right_ext_reads;
        dup->regenotyping_info.hq_alt_ext_reads = alt_consensus->hq_left_ext_reads + alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        delete alt_consensus;

        hts_pos_t lh_start = dup->end-alt_consensus_seq.length();
        if (lh_start < 0) lh_start = 0;
        hts_pos_t lh_len = dup->end-lh_start;
        hts_pos_t rh_start = dup->start;
        hts_pos_t rh_end = dup->start+alt_consensus_seq.length();
        if (rh_end > contig_len) rh_end = contig_len;
        hts_pos_t rh_len = rh_end-dup->start;
    
        delete[] alt_seq;
        alt_seq = new char[lh_len + dup->ins_seq.length() + rh_len + 1];
        strncpy(alt_seq, contig_seq+lh_start, lh_len);
        strncpy(alt_seq+lh_len, dup->ins_seq.c_str(), dup->ins_seq.length());
        strncpy(alt_seq+lh_len+dup->ins_seq.length(), contig_seq+rh_start, rh_len);
        alt_seq[lh_len+dup->ins_seq.length()+rh_len] = 0;

        // align to ref+SV
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, lh_len+dup->ins_seq.length()+rh_len, filter, &alt_aln, 0);

        // align to ref
        hts_pos_t ref_bp1_start = dup->start - alt_consensus_seq.length(), ref_bp1_end = dup->start + alt_consensus_seq.length();
        if (ref_bp1_start < 0) ref_bp1_start = 0;
        if (ref_bp1_end > contig_len) ref_bp1_end = contig_len;
        hts_pos_t ref_bp2_start = dup->end - alt_consensus_seq.length(), ref_bp2_end = dup->end + alt_consensus_seq.length();
        if (ref_bp2_start < 0) ref_bp2_start = 0;
        if (ref_bp2_end > contig_len) ref_bp2_end = contig_len;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &ref2_aln, 0);

        dup->regenotyping_info.ext_alt_consensus_length = alt_consensus_seq.length();
        dup->regenotyping_info.ext_alt_consensus_to_alt_score = alt_aln.sw_score;
        dup->regenotyping_info.ext_alt_consensus_to_ref_score = std::max(ref1_aln.sw_score, ref2_aln.sw_score);
    }

    dup->regenotyping_info.alt_bp1_better_reads = alt_better_reads.size();
    dup->regenotyping_info.alt_bp2_better_reads = alt_better_reads.size();
    dup->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_reads_consistent.size();
    dup->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_reads_consistent.size();
    for (bam1_t* read : alt_better_reads_consistent) {
        if (bam_is_mrev(read)) {
            dup->regenotyping_info.alt_bp1_better_consistent_reads_fwd++;
        } else {
            dup->regenotyping_info.alt_bp1_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        dup->regenotyping_info.alt_bp1_better_consistent_max_mq = std::max(dup->regenotyping_info.alt_bp1_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            dup->regenotyping_info.alt_bp1_better_consistent_high_mq++;
        }
    }
    dup->regenotyping_info.alt_bp2_better_consistent_reads_fwd = dup->regenotyping_info.alt_bp1_better_consistent_reads_fwd;
    dup->regenotyping_info.alt_bp2_better_consistent_reads_rev = dup->regenotyping_info.alt_bp1_better_consistent_reads_rev;
    dup->regenotyping_info.alt_bp2_better_consistent_max_mq = dup->regenotyping_info.alt_bp1_better_consistent_max_mq;
    dup->regenotyping_info.alt_bp2_better_consistent_high_mq = dup->regenotyping_info.alt_bp1_better_consistent_high_mq;
    dup->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.ref_bp1_better_reads = ref_bp1_better_reads.size();
    dup->regenotyping_info.ref_bp2_better_reads = ref_bp2_better_reads.size();
    dup->regenotyping_info.ref_bp1_better_consistent_reads = ref_bp1_better_reads_consistent.size();
    dup->regenotyping_info.ref_bp2_better_consistent_reads = ref_bp2_better_reads_consistent.size();

    delete[] alt_seq;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_dups(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<duplication_t*> dups,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (duplication_t* dup : dups) {
        target_ivals.push_back({dup->start-stats.max_is, dup->start+stats.max_is});
        target_ivals.push_back({dup->end-stats.max_is, dup->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_file, candidate_reads_for_extension);

    std::vector<sv_t*> small_dups;
    for (duplication_t* dup : dups) {
        if (dup->svlen() <= stats.read_len-2*config.min_clip_len) {
			genotype_small_dup(dup, bam_file, candidate_reads_for_extension_itree, mateseqs_w_mapq[contig_id]);
            small_dups.push_back(dup);
		} else {
			genotype_large_dup(dup, bam_file, candidate_reads_for_extension_itree, mateseqs_w_mapq[contig_id]);
		}
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    release_mates(contig_id);

    depth_filter_dup(contig_name, dups, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, global_crossing_isize_dist, small_dups, bam_file, config, stats, config.min_sv_size, true);
    std::string mates_nms_file = workdir + "/workspace/outward-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dups, bam_file, config, stats, true, mates_nms_file);
    calculate_cluster_region_disc(contig_name, dups, bam_file, config, stats);
}

void genotype_ins(insertion_t* ins, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
    hts_pos_t ins_start = ins->start, ins_end = ins->end;

	hts_pos_t extend = stats.read_len + 20;

	// build alt allele
	/*
	 * POS in VCF is the base BEFORE the insertion
	 * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" insertion (no deletion),POS == END, otherwise the last base deleted
	 * As usual, in order to make intervals [ ), we increase the coordinates by 1
	 */
	ins_start++; ins_end++;

    char* contig_seq = chr_seqs.get_seq(ins->chr);
    hts_pos_t contig_len = chr_seqs.get_len(ins->chr);

	hts_pos_t alt_start = std::max(hts_pos_t(0), ins_start-extend);
	hts_pos_t alt_end = std::min(ins_end+extend, contig_len);
	int alt_lf_len = ins_start-alt_start, alt_rf_len = alt_end-ins_end;
    hts_pos_t ins_lh_len = std::min(extend, (hts_pos_t) ins->ins_seq.length());
	hts_pos_t ins_rh_len = ins_lh_len;

	int alt_bp1_len = alt_lf_len + ins_lh_len;
	char* alt_bp1_seq = new char[alt_bp1_len+1];
	strncpy(alt_bp1_seq, contig_seq+alt_start, alt_lf_len);
	strncpy(alt_bp1_seq+alt_lf_len, ins->ins_seq.c_str(), ins_lh_len);
	alt_bp1_seq[alt_bp1_len] = 0;

    int alt_bp2_len = ins_rh_len + alt_rf_len;
	char* alt_bp2_seq = new char[alt_bp2_len+1];
	strncpy(alt_bp2_seq, ins->ins_seq.c_str()+(ins->ins_seq.length()-ins_rh_len), ins_rh_len);
	strncpy(alt_bp2_seq+ins_rh_len, contig_seq+ins_end, alt_rf_len);
	alt_bp2_seq[alt_bp2_len] = 0;

    hts_pos_t ref_bp1_start = alt_start, ref_bp1_end = std::min(ins_start+extend, contig_len);
    hts_pos_t ref_bp1_len = ref_bp1_end - ref_bp1_start;
    hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), ins_end-extend), ref_bp2_end = alt_end;
    hts_pos_t ref_bp2_len = ref_bp2_end - ref_bp2_start;

    std::stringstream l_region, r_region;
    l_region << ins->chr << ":" << ref_bp1_start << "-" << ref_bp1_end;
    r_region << ins->chr << ":" << ref_bp2_start << "-" << ref_bp2_end;

    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int same = 0;
    std::vector<bam1_t*> alt_bp1_better_seqs, alt_bp2_better_seqs, ref_bp1_better_seqs, ref_bp2_better_seqs;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < ins_start || ins_end < get_unclipped_start(read)) continue;
        if (ins_start < get_unclipped_start(read) && get_unclipped_end(read) < ins_end) continue;

        std::string seq = get_sequence(read);
        
        // align to ALT
        aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter, &alt1_aln, 0);
        aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter, &alt2_aln, 0);

        // align to REF (two breakpoints)
        aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln, 0);
        aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln, 0);

        StripedSmithWaterman::Alignment& alt_aln = alt1_aln.sw_score >= alt2_aln.sw_score ? alt1_aln : alt2_aln;
        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            if (alt1_aln.sw_score >= alt2_aln.sw_score) {
                alt_bp1_better_seqs.push_back(bam_dup1(read));
            } 
            if (alt1_aln.sw_score <= alt2_aln.sw_score) {
                alt_bp2_better_seqs.push_back(bam_dup1(read));
            }
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                ref_bp1_better_seqs.push_back(bam_dup1(read));
            } 
            if (ref1_aln.sw_score <= ref2_aln.sw_score) {
                ref_bp2_better_seqs.push_back(bam_dup1(read));
            }
        } else {
            same++;
        }

        if (alt_bp1_better_seqs.size() + alt_bp2_better_seqs.size() + ref_bp1_better_seqs.size() + ref_bp2_better_seqs.size() + same > 4 * stats.get_max_depth(ins->chr)) {
            alt_bp1_better_seqs.clear();
            alt_bp2_better_seqs.clear();
            ref_bp1_better_seqs.clear();
            ref_bp2_better_seqs.clear();
            same = 0;
            ins->regenotyping_info.too_deep = true;
            break;
        }
    }

    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    std::vector<bam1_t*> alt_bp1_better_seqs_consistent = find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_seqs, alt_bp1_consensus_seq, alt_bp1_avg_score);
    delete[] alt_bp1_seq;
    std::vector<bam1_t*> alt_bp2_better_seqs_consistent = find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_seqs, alt_bp2_consensus_seq, alt_bp2_avg_score);
    delete[] alt_bp2_seq;

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    std::vector<bam1_t*> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    std::vector<bam1_t*> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score);
    delete[] ref_bp2_seq;

    if (!alt_bp1_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_bp1_consensus = new consensus_t(false, 0, 0, 0, alt_bp1_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, ins->start-stats.max_is, ins->start, ins->chr, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, ins->start, ins->start+stats.max_is, ins->chr, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->regenotyping_info.alt_ext_reads = alt_bp1_consensus->left_ext_reads + alt_bp1_consensus->right_ext_reads;
        ins->regenotyping_info.hq_alt_ext_reads = alt_bp1_consensus->hq_left_ext_reads + alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;
        delete alt_bp1_consensus;
    }
    if (!alt_bp2_consensus_seq.empty()) {
        consensus_t* alt_bp2_consensus = new consensus_t(false, 0, 0, 0, alt_bp2_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, ins->end-stats.max_is, ins->end, ins->chr, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, ins->end, ins->end+stats.max_is, ins->chr, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->regenotyping_info.alt_ext_reads += alt_bp2_consensus->left_ext_reads + alt_bp2_consensus->right_ext_reads;
        ins->regenotyping_info.hq_alt_ext_reads += alt_bp2_consensus->hq_left_ext_reads + alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;
        delete alt_bp2_consensus;
    }

    ins->regenotyping_info.ext_alt_consensus_length = alt_bp1_consensus_seq.length() + alt_bp2_consensus_seq.length();

    ref1_aln.Clear();
    alt1_aln.Clear();
    if (!alt_bp1_consensus_seq.empty()) {
        hts_pos_t ref_bp1_start = ins->start-alt_bp1_consensus_seq.length();
        if (ref_bp1_start < 0) ref_bp1_start = 0;
        hts_pos_t ref_bp1_end = ins->start+alt_bp1_consensus_seq.length();
        if (ref_bp1_end > contig_len) ref_bp1_end = contig_len;
        aligner.Align(alt_bp1_consensus_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &ref1_aln, 0);

        alt_bp1_seq = new char[2*alt_bp1_consensus_seq.length()+1];
        strncpy(alt_bp1_seq, contig_seq+ref_bp1_start, ins->start-ref_bp1_start);
        int ins_seq_portion_len = std::min(ins->ins_seq.length(), alt_bp1_consensus_seq.length());
        strncpy(alt_bp1_seq+ins->start-ref_bp1_start, ins->ins_seq.c_str(), ins_seq_portion_len);
        int extra_len = alt_bp1_consensus_seq.length() - ins->ins_seq.length();
        if (extra_len > 0) {
            if (ins->end+extra_len > contig_len) extra_len = contig_len-ins->end;
            strncpy(alt_bp1_seq+ins->start-ref_bp1_start+ins_seq_portion_len, contig_seq+ins->end, extra_len);
        } else {
            extra_len = 0;
        }
        int alt_bp1_seq_len = ins->start-ref_bp1_start + ins_seq_portion_len + extra_len;
        alt_bp1_seq[alt_bp1_seq_len] = 0;

        aligner.Align(alt_bp1_consensus_seq.c_str(), alt_bp1_seq, alt_bp1_seq_len, filter, &alt1_aln, 0);
        delete[] alt_bp1_seq;
    }

    ref2_aln.Clear();
    alt2_aln.Clear();
    if (!alt_bp2_consensus_seq.empty()) {
        hts_pos_t ref_bp2_start = ins->end-alt_bp2_consensus_seq.length();
        if (ref_bp2_start < 0) ref_bp2_start = 0;
        hts_pos_t ref_bp2_end = ins->end+alt_bp2_consensus_seq.length();
        if (ref_bp2_end > contig_len) ref_bp2_end = contig_len; 
        aligner.Align(alt_bp2_consensus_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &ref2_aln, 0);

        alt_bp2_seq = new char[2*alt_bp2_consensus_seq.length()+1];
        int extra_len = alt_bp2_consensus_seq.length() - ins->ins_seq.length();
        if (extra_len > 0) {
            if (ins->start-extra_len < 0) extra_len = ins->start;
            strncpy(alt_bp2_seq, contig_seq+ins->start-extra_len, extra_len);
        } else {
            extra_len = 0;
        }
        int ins_seq_portion_len = std::min(ins->ins_seq.length(), alt_bp2_consensus_seq.length());
        strncpy(alt_bp2_seq+extra_len, ins->ins_seq.c_str()+(ins->ins_seq.length()-ins_seq_portion_len), ins_seq_portion_len);
        strncpy(alt_bp2_seq+extra_len+ins_seq_portion_len, contig_seq+ins->end, ref_bp2_end-ins->end);
        int alt_bp2_seq_len = extra_len + ins_seq_portion_len + ref_bp2_end-ins->end;
        alt_bp2_seq[alt_bp2_seq_len] = 0;

        aligner.Align(alt_bp2_consensus_seq.c_str(), alt_bp2_seq, alt_bp2_seq_len, filter, &alt2_aln, 0);
        delete[] alt_bp2_seq;
    }

    ins->regenotyping_info.ext_alt_consensus_to_ref_score = ref1_aln.sw_score + ref2_aln.sw_score;
    ins->regenotyping_info.ext_alt_consensus_to_alt_score = alt1_aln.sw_score + alt2_aln.sw_score;

    ins->regenotyping_info.alt_bp1_better_reads = alt_bp1_better_seqs.size();
    ins->regenotyping_info.alt_bp2_better_reads = alt_bp2_better_seqs.size();
    ins->regenotyping_info.alt_bp1_better_consistent_reads = alt_bp1_better_seqs_consistent.size();
    ins->regenotyping_info.alt_bp2_better_consistent_reads = alt_bp2_better_seqs_consistent.size();
    for (bam1_t* read : alt_bp1_better_seqs_consistent) {
        if (!bam_is_rev(read)) {
            ins->regenotyping_info.alt_bp1_better_consistent_reads_fwd++;
        } else {
            ins->regenotyping_info.alt_bp1_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        ins->regenotyping_info.alt_bp1_better_consistent_max_mq = std::max(ins->regenotyping_info.alt_bp1_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            ins->regenotyping_info.alt_bp1_better_consistent_high_mq++;
        }
    }
    for (bam1_t* read : alt_bp2_better_seqs_consistent) {
        if (!bam_is_rev(read)) {
            ins->regenotyping_info.alt_bp2_better_consistent_reads_fwd++;
        } else {
            ins->regenotyping_info.alt_bp2_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        ins->regenotyping_info.alt_bp2_better_consistent_max_mq = std::max(ins->regenotyping_info.alt_bp2_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            ins->regenotyping_info.alt_bp2_better_consistent_high_mq++;
        }
    }
    ins->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_bp1_avg_score;
    ins->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_bp2_avg_score;
    ins->regenotyping_info.ref_bp1_better_reads = ref_bp1_better_seqs.size();
    ins->regenotyping_info.ref_bp2_better_reads = ref_bp2_better_seqs.size();
    ins->regenotyping_info.ref_bp1_better_consistent_reads = ref_bp1_better_seqs_consistent.size();
    ins->regenotyping_info.ref_bp2_better_consistent_reads = ref_bp2_better_seqs_consistent.size();
    ins->regenotyping_info.alt_ref_equal_reads = same;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void find_discordant_pairs(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, stats_t& stats,
                           std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
    
    std::vector<char*> regions;
    for (insertion_t* ins : insertions) {
        std::stringstream ss;
        ss << ins->chr << ":" << std::max(hts_pos_t(1), ins->start-stats.max_is) << "-" << ins->start;
        regions.push_back(strdup(ss.str().c_str()));
    }

    std::sort(insertions.begin(), insertions.end(), [](insertion_t* a, insertion_t* b) { return a->start < b->start; });

    std::vector<hts_pos_t> dp1_start(insertions.size(), INT32_MAX), dp1_end(insertions.size(), 0);
    std::vector<hts_pos_t> dp2_start(insertions.size(), INT32_MAX), dp2_end(insertions.size(), 0);

    int curr_pos = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        while (curr_pos < insertions.size() && insertions[curr_pos]->start < read->core.pos) curr_pos++;

        if (bam_is_rev(read)) continue;

        std::string qname = bam_get_qname(read);
        if (is_samechr(read)) {
            if (read->core.flag & BAM_FREAD1) {
                qname += "_2";
            } else {
                qname += "_1";
            }
        }
        if (mateseqs_w_mapq_chr.count(qname) == 0) continue;

        std::string mate_seq = mateseqs_w_mapq_chr[qname].first;
        rc(mate_seq);

        StripedSmithWaterman::Filter filter;
        StripedSmithWaterman::Alignment aln;
        for (int i = curr_pos; i < insertions.size() && insertions[i]->start-stats.max_is < read->core.pos; i++) {
            harsh_aligner.Align(mate_seq.c_str(), insertions[i]->ins_seq.c_str(), insertions[i]->ins_seq.length(), filter, &aln, 0);

            double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
            int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);
            
            if (mismatch_rate <= config.max_seq_error && (lc_size < config.min_clip_len || aln.ref_begin == 0) && 
                (rc_size < config.min_clip_len || aln.ref_end >= insertions[i]->ins_seq.length()-1)) {
                insertions[i]->disc_pairs_lf++;

                if (read->core.pos < dp1_start[i]) dp1_start[i] = read->core.pos;
                if (bam_endpos(read) > dp1_end[i]) dp1_end[i] = bam_endpos(read);

                if (read->core.qual >= config.high_confidence_mapq) insertions[i]->disc_pairs_lf_high_mapq++;
                insertions[i]->disc_pairs_lf_maxmapq = std::max(insertions[i]->disc_pairs_lf_maxmapq, mateseqs_w_mapq_chr[qname].second);
                insertions[i]->disc_pairs_lf_avg_nm += get_nm(read);
            } else {
                insertions[i]->l_cluster_region_disc_pairs++;
                if (read->core.qual >= config.high_confidence_mapq) insertions[i]->l_cluster_region_disc_pairs_high_mapq++;
            }
        }
    }
    for (int i = 0; i < insertions.size(); i++) {
        insertion_t* ins = insertions[i];
        if (ins->disc_pairs_lf > 0) ins->disc_pairs_lf_avg_nm /= ins->disc_pairs_lf;
        ins->disc_pairs_lf_span = std::max(hts_pos_t(0), dp1_end[i] - dp1_start[i]);
    }

    for (char* region : regions) free(region);
    hts_itr_destroy(iter);

    regions.clear();
    for (insertion_t* ins : insertions) {
        std::stringstream ss;
        ss << ins->chr << ":" << ins->end << "-" << ins->end+stats.max_is;
        regions.push_back(strdup(ss.str().c_str()));
    }

    std::sort(insertions.begin(), insertions.end(), [](insertion_t* a, insertion_t* b) { return a->end < b->end; });

    curr_pos = 0;
    iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        while (curr_pos < insertions.size() && insertions[curr_pos]->end+stats.max_is < read->core.pos) curr_pos++;

        if (!bam_is_rev(read)) continue;

        std::string qname = bam_get_qname(read);
        if (is_samechr(read)) {
            if (read->core.flag & BAM_FREAD1) {
                qname += "_2";
            } else {
                qname += "_1";
            }
        }
        if (mateseqs_w_mapq_chr.count(qname) == 0) continue;

        std::string mate_seq = mateseqs_w_mapq_chr[qname].first;

        StripedSmithWaterman::Filter filter;
        StripedSmithWaterman::Alignment aln;
        for (int i = curr_pos; i < insertions.size() && insertions[i]->end < read->core.pos; i++) {
            harsh_aligner.Align(mate_seq.c_str(), insertions[i]->ins_seq.c_str(), insertions[i]->ins_seq.length(), filter, &aln, 0);

            double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
            int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);

            if (mismatch_rate <= config.max_seq_error && (lc_size < config.min_clip_len || aln.ref_begin == 0) && 
                (rc_size < config.min_clip_len || aln.ref_end >= insertions[i]->ins_seq.length()-1)) {
                insertions[i]->disc_pairs_rf++;

                if (read->core.pos < dp2_start[i]) dp2_start[i] = read->core.pos;
                if (bam_endpos(read) > dp2_end[i]) dp2_end[i] = bam_endpos(read);

                if (read->core.qual >= config.high_confidence_mapq) insertions[i]->disc_pairs_rf_high_mapq++;
                insertions[i]->disc_pairs_rf_maxmapq = std::max(insertions[i]->disc_pairs_rf_maxmapq, mateseqs_w_mapq_chr[qname].second);
                insertions[i]->disc_pairs_rf_avg_nm += get_nm(read);
            } else {
                insertions[i]->r_cluster_region_disc_pairs++;
                if (read->core.qual >= config.high_confidence_mapq) insertions[i]->r_cluster_region_disc_pairs_high_mapq++;
            }
        }
    }
    for (int i = 0; i < insertions.size(); i++) {
        insertion_t* ins = insertions[i];
        if (ins->disc_pairs_rf > 0) ins->disc_pairs_rf_avg_nm /= ins->disc_pairs_rf;
        ins->disc_pairs_rf_span = std::max(hts_pos_t(0), dp2_end[i] - dp2_start[i]);
    }

    for (char* region : regions) free(region);
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void genotype_inss(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<insertion_t*> inss,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (insertion_t* ins : inss) {
        target_ivals.push_back({ins->start-stats.max_is, ins->start+stats.max_is});
        target_ivals.push_back({ins->end-stats.max_is, ins->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_pool->get_bam_reader(id), candidate_reads_for_extension);
                    
    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);
    for (insertion_t* ins : inss) { 
        genotype_ins(ins, bam_file, candidate_reads_for_extension_itree, mateseqs_w_mapq[contig_id]);
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    depth_filter_ins(contig_name, inss, bam_file, config, stats);
    calculate_ptn_ratio(contig_name, inss, bam_file, config, stats);
    find_discordant_pairs(contig_name, inss, bam_file, stats, mateseqs_w_mapq[contig_id]);

    release_mates(contig_id);
}

void genotype_small_inv(inversion_t* inv, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
    hts_pos_t inv_start = inv->start, inv_end = inv->end;

    hts_pos_t extend = stats.read_len + 20;

    // build alt allele
    /*
     * POS in VCF is the base BEFORE the inversion
     * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" inversion (no deletion),POS == END, otherwise the last base deleted
     * As usual, in order to make intervals [ ), we increase the coordinates by 1
     */
    inv_start++; inv_end++;

    char* contig_seq = chr_seqs.get_seq(inv->chr);
    hts_pos_t contig_len = chr_seqs.get_len(inv->chr);

    hts_pos_t alt_start = std::max(hts_pos_t(0), inv_start-extend);
    hts_pos_t alt_end = std::min(inv_end+extend, contig_len);
    int alt_lf_len = inv_start-alt_start, alt_rf_len = alt_end-inv_end;
    int alt_len = alt_end - alt_start;
    char* alt_seq = new char[alt_len+1];
    strncpy(alt_seq, contig_seq+alt_start, alt_lf_len);
    char* inv_seq = new char[inv->end-inv->start+1];
    strncpy(inv_seq, contig_seq+inv->start, inv->end-inv->start);
    inv_seq[inv->end-inv->start] = 0;
    rc(inv_seq);
    strncpy(alt_seq+alt_lf_len, inv_seq, inv->end-inv->start);
    strncpy(alt_seq+alt_lf_len+(inv->end-inv->start), contig_seq+inv->end, alt_rf_len);
    alt_seq[alt_len] = 0;

    hts_pos_t ref_start = std::max(hts_pos_t(0), inv_start-extend);
    hts_pos_t ref_end = std::min(inv_end+extend, contig_len);
    int ref_len = ref_end - ref_start;
    char* ref_seq = new char[ref_len+1];
    strncpy(ref_seq, contig_seq+ref_start, ref_len);
    ref_seq[ref_len] = 0;

    std::stringstream region;
    region << inv->chr << ":" << ref_start << "-" << ref_end;
    char* regions[1];
    regions[0] = strdup(region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 1);

    bam1_t* read = bam_init1();

    int same = 0;
    std::vector<bam1_t*> alt_better_seqs, ref_better_seqs;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < inv_start || inv_end < get_unclipped_start(read)) continue;
        if (inv_start < get_unclipped_start(read) && get_unclipped_end(read) < inv_end) continue;

        std::string seq = get_sequence(read);
        if (bam_is_rev(read) && bam_is_mrev(read) && inv->end+stats.read_len/2 <= get_mate_endpos(read)) rc(seq);
        if (!bam_is_rev(read) && !bam_is_mrev(read) && read->core.mpos <= inv->start-stats.read_len/2) rc(seq);
        
        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

        // // align to REF
        aligner.Align(seq.c_str(), ref_seq, ref_len, filter, &ref_aln, 0);

        if (alt_aln.sw_score > ref_aln.sw_score) {
            alt_better_seqs.push_back(bam_dup1(read));
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            ref_better_seqs.push_back(bam_dup1(read));
        } else {
            same++;
        }

        if (alt_better_seqs.size() + ref_better_seqs.size() + same > 4 * stats.get_max_depth(inv->chr)) {
            alt_better_seqs.clear();
            ref_better_seqs.clear();
            same = 0;
            inv->regenotyping_info.too_deep = true;
            break;
        }
    }

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seq, alt_better_seqs, alt_consensus_seq, alt_avg_score);
    std::vector<bam1_t*> ref_better_reads_consistent = find_consistent_seqs_subset(ref_seq, ref_better_seqs, ref_consensus_seq, ref_avg_score);

    if (!alt_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, inv->start-stats.max_is, inv->start, inv->chr, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, inv->start, inv->start+stats.max_is, inv->chr, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        inv->regenotyping_info.alt_ext_reads = alt_consensus->left_ext_reads + alt_consensus->right_ext_reads;
        inv->regenotyping_info.hq_alt_ext_reads = alt_consensus->hq_left_ext_reads + alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        delete alt_consensus;

        delete[] alt_seq;
        alt_start = std::max(hts_pos_t(0), inv->start-stats.max_is);
        alt_end = std::min(inv->end+stats.max_is, contig_len);
        alt_seq = new char[alt_end-alt_start+1];
        strncpy(alt_seq, contig_seq+alt_start, inv->start-alt_start);
        strncpy(alt_seq+(inv->start-alt_start), inv_seq, inv->end-inv->start);
        strncpy(alt_seq+(inv->start-alt_start)+(inv->end-inv->start), contig_seq+inv->end, alt_end-inv->end);
        alt_seq[alt_end-alt_start] = 0;

        // align to ref+SV
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, alt_end-alt_start, filter, &alt_aln, 0);

        delete[] ref_seq;
        ref_start = std::max(hts_pos_t(0), inv->start-stats.max_is);
        ref_end = std::min(inv->end+stats.max_is, contig_len);
        ref_len = ref_end-ref_start;
        ref_seq = new char[ref_len+1];
        strncpy(ref_seq, contig_seq+ref_start, ref_len);
        ref_seq[ref_len] = 0;

        // align to ref
        aligner.Align(alt_consensus_seq.c_str(), ref_seq, ref_len, filter, &ref_aln, 0);

        inv->regenotyping_info.ext_alt_consensus_length = alt_consensus_seq.length();
        inv->regenotyping_info.ext_alt_consensus_to_ref_score = ref_aln.sw_score;
        inv->regenotyping_info.ext_alt_consensus_to_alt_score = alt_aln.sw_score;
    }

    inv->regenotyping_info.alt_bp1_better_reads = alt_better_seqs.size();
    inv->regenotyping_info.alt_bp2_better_reads = alt_better_seqs.size();
    inv->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_reads_consistent.size();
    inv->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_reads_consistent.size();
    for (bam1_t* read : alt_better_reads_consistent) {
        if (!bam_is_rev(read)) {
            inv->regenotyping_info.alt_bp1_better_consistent_reads_fwd++;
        } else {
            inv->regenotyping_info.alt_bp1_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        inv->regenotyping_info.alt_bp1_better_consistent_max_mq = std::max(inv->regenotyping_info.alt_bp1_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            inv->regenotyping_info.alt_bp1_better_consistent_high_mq++;
        }
    }
    inv->regenotyping_info.alt_bp2_better_consistent_reads_fwd = inv->regenotyping_info.alt_bp1_better_consistent_reads_fwd;
    inv->regenotyping_info.alt_bp2_better_consistent_reads_rev = inv->regenotyping_info.alt_bp1_better_consistent_reads_rev;
    inv->regenotyping_info.alt_bp2_better_consistent_max_mq = inv->regenotyping_info.alt_bp1_better_consistent_max_mq;
    inv->regenotyping_info.alt_bp2_better_consistent_high_mq = inv->regenotyping_info.alt_bp1_better_consistent_high_mq;
    inv->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_avg_score;
    inv->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_avg_score;
    inv->regenotyping_info.ref_bp1_better_reads = ref_better_seqs.size();
    inv->regenotyping_info.ref_bp2_better_reads = ref_better_seqs.size();
    inv->regenotyping_info.ref_bp1_better_consistent_reads = ref_better_reads_consistent.size();
    inv->regenotyping_info.ref_bp2_better_consistent_reads = ref_better_reads_consistent.size();
    inv->regenotyping_info.alt_ref_equal_reads = same;

    delete[] alt_seq;
    delete[] ref_seq;
    delete[] inv_seq;
    
    free(regions[0]);
    
    for (bam1_t* read : alt_better_seqs) bam_destroy1(read);
    for (bam1_t* read : ref_better_seqs) bam_destroy1(read);
    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_large_inv(inversion_t* inv, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
    hts_pos_t inv_start = inv->start, inv_end = inv->end;

    hts_pos_t extend = stats.read_len + 20;

    // build alt allele
    /*
     * POS in VCF is the base BEFORE the inversion
     * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" inversion (no deletion),POS == END, otherwise the last base deleted
     * As usual, in order to make intervals [ ), we increase the coordinates by 1
     */
    inv_start++; inv_end++;

    char* contig_seq = chr_seqs.get_seq(inv->chr);
    hts_pos_t contig_len = chr_seqs.get_len(inv->chr);

    hts_pos_t alt_start = std::max(hts_pos_t(0), inv_start-extend);
    hts_pos_t alt_end = std::min(inv_end+extend, contig_len);
    int alt_lf_len = inv_start-alt_start, alt_rf_len = alt_end-inv_end;
    hts_pos_t inv_border_len = std::min(extend, inv->end-inv->start);

    char* inv_prefix = new char[inv_border_len+1];
    strncpy(inv_prefix, contig_seq+inv->start, inv_border_len);
    inv_prefix[inv_border_len] = 0;

    char* inv_prefix_rc = strdup(inv_prefix);
    rc(inv_prefix_rc);

    char* inv_suffix = new char[inv_border_len+1];
    strncpy(inv_suffix, contig_seq+inv->end-inv_border_len, inv_border_len);
    inv_suffix[inv_border_len] = 0;

    char* inv_suffix_rc = strdup(inv_suffix);
    rc(inv_suffix_rc);

    int alt_bp1_len = alt_lf_len + inv_border_len;
    char* alt_bp1_seq = new char[alt_bp1_len+1];
    strncpy(alt_bp1_seq, contig_seq+alt_start, alt_lf_len);
    strncpy(alt_bp1_seq+alt_lf_len, inv_suffix_rc, inv_border_len);
    alt_bp1_seq[alt_bp1_len] = 0;

    int alt_bp2_len = inv_border_len + alt_rf_len;
    char* alt_bp2_seq = new char[alt_bp2_len+1];
    strncpy(alt_bp2_seq, inv_prefix_rc, inv_border_len);
    strncpy(alt_bp2_seq+inv_border_len, contig_seq+inv->end, alt_rf_len);
    alt_bp2_seq[alt_bp2_len] = 0;

    hts_pos_t ref_bp1_start = std::max(hts_pos_t(0), inv_start-extend);
    hts_pos_t ref_bp1_end = std::min(inv->start+extend, contig_len);
    hts_pos_t ref_bp1_len = ref_bp1_end-ref_bp1_start;
    hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), inv->end-extend);
    hts_pos_t ref_bp2_end = std::min(inv_end+extend, contig_len);
    hts_pos_t ref_bp2_len = ref_bp2_end-ref_bp2_start;

    std::stringstream l_region, r_region;
    l_region << inv->chr << ":" << ref_bp1_start << "-" << ref_bp1_end;
    r_region << inv->chr << ":" << ref_bp2_start << "-" << ref_bp2_end;

    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int same = 0;
    std::vector<bam1_t*> alt_bp1_better_reads, alt_bp2_better_reads, ref_bp1_better_reads, ref_bp2_better_reads;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < inv_start || inv_end < get_unclipped_start(read)) continue;
        if (inv_start < get_unclipped_start(read) && get_unclipped_end(read) < inv_end) continue;

        std::string seq = get_sequence(read);
        if (bam_is_rev(read) && bam_is_mrev(read) && inv->end+stats.read_len/2 <= get_mate_endpos(read)) rc(seq);
        if (!bam_is_rev(read) && !bam_is_mrev(read) && read->core.mpos <= inv->start-stats.read_len/2) rc(seq);

        // if mate is inside the inversion, we can't infer the orientation and we have to test both
        if (inv->start-10 <= read->core.mpos && get_mate_endpos(read) <= inv->end+10) {
            // align to ALT
            StripedSmithWaterman::Alignment alt1_aln_fwd, alt1_aln_rev, alt2_aln_fwd, alt2_aln_rev;
            aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter, &alt1_aln_fwd, 0);
            aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter, &alt2_aln_fwd, 0);

            rc(seq);
            aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter, &alt1_aln_rev, 0);
            aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter, &alt2_aln_rev, 0);

            alt1_aln = alt1_aln_fwd.sw_score >= alt1_aln_rev.sw_score ? alt1_aln_fwd : alt1_aln_rev;
            alt2_aln = alt2_aln_fwd.sw_score >= alt2_aln_rev.sw_score ? alt2_aln_fwd : alt2_aln_rev;

            // align to REF
            StripedSmithWaterman::Alignment ref1_aln_fwd, ref1_aln_rev, ref2_aln_fwd, ref2_aln_rev;
            aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln_fwd, 0);
            aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln_fwd, 0);

            rc(seq);
            aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln_rev, 0);
            aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln_rev, 0);

            ref1_aln = ref1_aln_fwd.sw_score >= ref1_aln_rev.sw_score ? ref1_aln_fwd : ref1_aln_rev;
            ref2_aln = ref2_aln_fwd.sw_score >= ref2_aln_rev.sw_score ? ref2_aln_fwd : ref2_aln_rev;
        } else {
            // align to ALT
            aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter, &alt1_aln, 0);
            aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter, &alt2_aln, 0);

            // align to REF
            aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln, 0);
            aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln, 0);
        }

        StripedSmithWaterman::Alignment& alt_aln = alt1_aln.sw_score >= alt2_aln.sw_score ? alt1_aln : alt2_aln;
        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            if (alt1_aln.sw_score >= alt2_aln.sw_score) {
                alt_bp1_better_reads.push_back(bam_dup1(read));
            } 
            if (alt1_aln.sw_score <= alt2_aln.sw_score) {
                alt_bp2_better_reads.push_back(bam_dup1(read));
            }
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                ref_bp1_better_reads.push_back(bam_dup1(read));
            } 
            if (ref1_aln.sw_score <= ref2_aln.sw_score) {
                ref_bp2_better_reads.push_back(bam_dup1(read));
            }
        } else {
            same++;
        }

        if (alt_bp1_better_reads.size() + alt_bp2_better_reads.size() + ref_bp1_better_reads.size() + ref_bp2_better_reads.size() + same > 4 * stats.get_max_depth(inv->chr)) {
            alt_bp1_better_reads.clear();
            alt_bp2_better_reads.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            same = 0;
            inv->regenotyping_info.too_deep = true;
            break;
        }
    }

    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    std::vector<bam1_t*> alt_bp1_better_reads_consistent = find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_reads, alt_bp1_consensus_seq, alt_bp1_avg_score);
    delete[] alt_bp1_seq;
    std::vector<bam1_t*> alt_bp2_better_reads_consistent = find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_reads, alt_bp2_consensus_seq, alt_bp2_avg_score);
    delete[] alt_bp2_seq;

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    std::vector<bam1_t*> ref_bp1_better_reads_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_reads, ref_bp1_consensus_seq, ref_bp1_avg_score);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    std::vector<bam1_t*> ref_bp2_better_reads_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_reads, ref_bp2_consensus_seq, ref_bp2_avg_score);
    delete[] ref_bp2_seq;

    if (!alt_bp1_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_bp1_consensus = new consensus_t(false, 0, 0, 0, alt_bp1_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, inv->start-stats.max_is, inv->start, inv->chr, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, inv->start, inv->start+stats.max_is, inv->chr, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        inv->regenotyping_info.alt_ext_reads = alt_bp1_consensus->left_ext_reads + alt_bp1_consensus->right_ext_reads;
        inv->regenotyping_info.hq_alt_ext_reads = alt_bp1_consensus->hq_left_ext_reads + alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;
        delete alt_bp1_consensus;
    }

    if (!alt_bp2_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_bp2_consensus = new consensus_t(false, 0, 0, 0, alt_bp2_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, inv->end-stats.max_is, inv->end, inv->chr, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, inv->end, inv->end+stats.max_is, inv->chr, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        inv->regenotyping_info.alt_ext_reads += alt_bp2_consensus->left_ext_reads + alt_bp2_consensus->right_ext_reads;
        inv->regenotyping_info.hq_alt_ext_reads += alt_bp2_consensus->hq_left_ext_reads + alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;
        delete alt_bp2_consensus;
    }
    inv->regenotyping_info.ext_alt_consensus_length = alt_bp1_consensus_seq.length() + alt_bp2_consensus_seq.length();

    ref1_aln.Clear();
    alt1_aln.Clear();
    if (!alt_bp1_consensus_seq.empty()) {
        hts_pos_t ref_bp1_start = inv->start-alt_bp1_consensus_seq.length();
        if (ref_bp1_start < 0) ref_bp1_start = 0;
        hts_pos_t ref_bp1_end = inv->start+alt_bp1_consensus_seq.length();
        if (ref_bp1_end > contig_len) ref_bp1_end = contig_len;
        aligner.Align(alt_bp1_consensus_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &ref1_aln, 0);

        hts_pos_t alt_bp1_start = std::max(hts_pos_t(0), hts_pos_t(inv->start-alt_bp1_consensus_seq.length()));
        hts_pos_t alt_bp1_end = std::min(hts_pos_t(inv->start+alt_bp1_consensus_seq.length()), contig_len);
        hts_pos_t alt_bp1_len = alt_bp1_end-alt_bp1_start;
        alt_bp1_seq = new char[alt_bp1_len+1];
        strncpy(alt_bp1_seq, contig_seq+alt_bp1_start, inv->start-alt_bp1_start);
        
        hts_pos_t ext_inv_border_len = std::min((hts_pos_t) alt_bp1_consensus_seq.length(), inv->end-inv->start);
        char* ext_inv_suffix_rc = new char[ext_inv_border_len+1];
        strncpy(ext_inv_suffix_rc, contig_seq+inv->end-ext_inv_border_len, ext_inv_border_len);
        ext_inv_suffix_rc[ext_inv_border_len] = 0;
        rc(ext_inv_suffix_rc);
        strncpy(alt_bp1_seq+(inv->start-alt_bp1_start), ext_inv_suffix_rc, ext_inv_border_len);
        int extra_len = alt_bp1_consensus_seq.length() - ext_inv_border_len;
        if (extra_len > 0) {
            strncpy(alt_bp1_seq+(inv->start-alt_bp1_start)+ext_inv_border_len, contig_seq+inv->end, extra_len);
        }
        alt_bp1_seq[alt_bp1_len] = 0;

        aligner.Align(alt_bp1_consensus_seq.c_str(), alt_bp1_seq, alt_bp1_len, filter, &alt1_aln, 0);
        delete[] ext_inv_suffix_rc;
        delete[] alt_bp1_seq;
    }

    ref2_aln.Clear();
    alt2_aln.Clear();
    if (!alt_bp2_consensus_seq.empty()) {
        hts_pos_t ref_bp2_start = inv->end-alt_bp2_consensus_seq.length();
        if (ref_bp2_start < 0) ref_bp2_start = 0;
        hts_pos_t ref_bp2_end = inv->end+alt_bp2_consensus_seq.length();
        if (ref_bp2_end > contig_len) ref_bp2_end = contig_len;
        aligner.Align(alt_bp2_consensus_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &ref2_aln, 0);

        hts_pos_t alt_bp2_start = std::max(hts_pos_t(0), hts_pos_t(inv->end-alt_bp2_consensus_seq.length()));
        hts_pos_t alt_bp2_end = std::min(hts_pos_t(inv->end+alt_bp2_consensus_seq.length()), contig_len);
        hts_pos_t alt_bp2_len = alt_bp2_end-alt_bp2_start;
        alt_bp2_seq = new char[alt_bp2_len+1];
        
        hts_pos_t ext_inv_border_len = std::min((hts_pos_t) alt_bp2_consensus_seq.length(), inv->end-inv->start);
        int extra_len = alt_bp2_consensus_seq.length() - ext_inv_border_len;
        if (extra_len > 0) {
            strncpy(alt_bp2_seq, contig_seq+inv->start-extra_len, extra_len);
        } else extra_len = 0;
        char* ext_inv_prefix_rc = new char[ext_inv_border_len+1];
        strncpy(ext_inv_prefix_rc, contig_seq+inv->start, ext_inv_border_len);
        ext_inv_prefix_rc[ext_inv_border_len] = 0;
        rc(ext_inv_prefix_rc);
        strncpy(alt_bp2_seq+extra_len, ext_inv_prefix_rc, ext_inv_border_len);
        strncpy(alt_bp2_seq+extra_len+ext_inv_border_len, contig_seq+inv->end, alt_bp2_end-inv->end);
        alt_bp2_seq[alt_bp2_len] = 0;

        aligner.Align(alt_bp2_consensus_seq.c_str(), alt_bp2_seq, alt_bp2_len, filter, &alt2_aln, 0);
        delete[] ext_inv_prefix_rc;
        delete[] alt_bp2_seq;
    }

    inv->regenotyping_info.ext_alt_consensus_to_ref_score = ref1_aln.sw_score + ref2_aln.sw_score;
    inv->regenotyping_info.ext_alt_consensus_to_alt_score = alt1_aln.sw_score + alt2_aln.sw_score;

    inv->regenotyping_info.alt_bp1_better_reads = alt_bp1_better_reads.size();
    inv->regenotyping_info.alt_bp2_better_reads = alt_bp2_better_reads.size();
    inv->regenotyping_info.alt_bp1_better_consistent_reads = alt_bp1_better_reads_consistent.size();
    inv->regenotyping_info.alt_bp2_better_consistent_reads = alt_bp2_better_reads_consistent.size();
    for (bam1_t* read : alt_bp1_better_reads_consistent) {
        if (!bam_is_rev(read)) {
            inv->regenotyping_info.alt_bp1_better_consistent_reads_fwd++;
        } else {
            inv->regenotyping_info.alt_bp1_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        inv->regenotyping_info.alt_bp1_better_consistent_max_mq = std::max(inv->regenotyping_info.alt_bp1_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            inv->regenotyping_info.alt_bp1_better_consistent_high_mq++;
        }
    }
    for (bam1_t* read : alt_bp2_better_reads_consistent) {
        if (!bam_is_rev(read)) {
            inv->regenotyping_info.alt_bp2_better_consistent_reads_fwd++;
        } else {
            inv->regenotyping_info.alt_bp2_better_consistent_reads_rev++;
        }
        int mq = get_mq(read);
        inv->regenotyping_info.alt_bp2_better_consistent_max_mq = std::max(inv->regenotyping_info.alt_bp2_better_consistent_max_mq, mq);
        if (mq >= config.high_confidence_mapq) {
            inv->regenotyping_info.alt_bp2_better_consistent_high_mq++;
        }
    }
    inv->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_bp1_avg_score;
    inv->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_bp2_avg_score;
    inv->regenotyping_info.ref_bp1_better_reads = ref_bp1_better_reads.size();
    inv->regenotyping_info.ref_bp2_better_reads = ref_bp2_better_reads.size();
    inv->regenotyping_info.ref_bp1_better_consistent_reads = ref_bp1_better_reads_consistent.size();
    inv->regenotyping_info.ref_bp2_better_consistent_reads = ref_bp2_better_reads_consistent.size();
    inv->regenotyping_info.alt_ref_equal_reads = same;

    free(regions[0]);
    free(regions[1]);

    for (bam1_t* read : alt_bp1_better_reads) bam_destroy1(read);
    for (bam1_t* read : alt_bp2_better_reads) bam_destroy1(read);
    for (bam1_t* read : ref_bp1_better_reads) bam_destroy1(read);
    for (bam1_t* read : ref_bp2_better_reads) bam_destroy1(read);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_invs(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<inversion_t*> invs,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (inversion_t* inv : invs) {
        target_ivals.push_back({inv->start-stats.max_is, inv->start+stats.max_is});
        target_ivals.push_back({inv->end-stats.max_is, inv->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_pool->get_bam_reader(id), candidate_reads_for_extension);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);
    for (inversion_t* inv : invs) {
        if (inv->svlen() < stats.read_len-2*config.min_clip_len) {
            genotype_small_inv(inv, bam_file, candidate_reads_for_extension_itree, mateseqs_w_mapq[contig_id]);
        } else {
            genotype_large_inv(inv, bam_file, candidate_reads_for_extension_itree, mateseqs_w_mapq[contig_id]);
        }
    }

    calculate_ptn_ratio(contig_name, invs, bam_file, config, stats, true);
    depth_filter_inv(contig_name, invs, bam_file, config, stats);

    for (inversion_t* inv : invs) {
        double ptn1 = inv->disc_pairs_lf/double(inv->conc_pairs_lbp);
        double ptn2 = inv->disc_pairs_rf/double(inv->conc_pairs_rbp);
        double ptn = std::min(ptn1, ptn2);
        inv->regenotyping_info.n_gt = 2;
        inv->regenotyping_info.gt = new int[2];
        if (ptn >= 0.75) {
            inv->regenotyping_info.gt[0] = bcf_gt_unphased(1);
            inv->regenotyping_info.gt[1] = bcf_gt_unphased(1);
        } else if (ptn <= 0.25) {
            inv->regenotyping_info.gt[0] = bcf_gt_unphased(0);
            inv->regenotyping_info.gt[1] = bcf_gt_unphased(0);
        } else {
            inv->regenotyping_info.gt[0] = bcf_gt_unphased(0);
            inv->regenotyping_info.gt[1] = bcf_gt_unphased(1);
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

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    open_samFile_t* bam_file = open_samFile(bam_fname);
	if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + reference_fname;
	}

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
    std::unordered_map<std::string, std::vector<deletion_t*> > dels_by_chr;
    std::unordered_map<std::string, std::vector<duplication_t*> > dups_by_chr;
    std::unordered_map<std::string, std::vector<insertion_t*> > inss_by_chr;
    std::unordered_map<std::string, std::vector<inversion_t*> > invs_by_chr;
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        sv_t* sv = bcf_to_sv(in_vcf_header, vcf_record);
        if (sv == NULL) {
            std::cout << "Ignoring SV of unsupported type: " << vcf_record->d.id << std::endl; 
            continue;
        }

        sv->vcf_entry = bcf_dup(vcf_record);
        reset_stats(sv);
        if (sv->svtype() == "DEL") {
            dels_by_chr[sv->chr].push_back((deletion_t*) sv);
        } else if (sv->svtype() == "DUP") {
            dups_by_chr[sv->chr].push_back((duplication_t*) sv);
        } else if (sv->svtype() == "INS") {
        	inss_by_chr[sv->chr].push_back((insertion_t*) sv);
        } else if (sv->svtype() == "INV") {
            invs_by_chr[sv->chr].push_back((inversion_t*) sv);
        }
    }

    int* imap = new int[1];
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    bcf_hdr_t* out_vcf_header;
    int sample_idx = find_sample_index(in_vcf_header, sample_name);
    if (sample_idx >= 0) {
        char** samples = new char*[1];
        samples[0] = strdup(sample_name.c_str());
        out_vcf_header = bcf_hdr_subset(in_vcf_header, 1, samples, imap);
    } else {
        out_vcf_header = bcf_hdr_subset(in_vcf_header, 0, NULL, NULL);
        bcf_hdr_add_sample(out_vcf_header, sample_name.c_str());
        imap[0] = -1;
    }
    add_fmt_tags(out_vcf_header);
    if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
    	throw std::runtime_error("Failed to read the VCF header.");
    }

    // genotype chrs in descending order of svs
    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    const int BLOCK_SIZE = 20;

    for (int i = 0; i < contig_map.size(); i++) {
    	std::string contig_name = contig_map.get_name(i);
        std::vector<deletion_t*>& dels = dels_by_chr[contig_name];
        for (int i = 0; i < dels.size(); i += BLOCK_SIZE) {
            std::vector<deletion_t*> block_dels(dels.begin() + i, dels.begin() + std::min(i + BLOCK_SIZE, static_cast<int>(dels.size())));
            std::future<void> future = thread_pool.push(genotype_dels, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_dels, in_vcf_header, out_vcf_header, stats, config);
            futures.push_back(std::move(future));
        }

        std::vector<duplication_t*>& dups = dups_by_chr[contig_name];
        for (int i = 0; i < dups.size(); i += BLOCK_SIZE) {
            std::vector<duplication_t*> block_dups(dups.begin() + i, dups.begin() + std::min(i + BLOCK_SIZE, static_cast<int>(dups.size())));
            std::future<void> future = thread_pool.push(genotype_dups, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_dups, in_vcf_header, out_vcf_header, stats, config);
            futures.push_back(std::move(future));
        }

        std::vector<insertion_t*>& inss = inss_by_chr[contig_name];
        for (int i = 0; i < inss.size(); i += BLOCK_SIZE) {
            std::vector<insertion_t*> block_inss(inss.begin() + i, inss.begin() + std::min(i + BLOCK_SIZE, static_cast<int>(inss.size())));
            std::future<void> future = thread_pool.push(genotype_inss, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_inss, in_vcf_header, out_vcf_header, stats, config);
            futures.push_back(std::move(future));
        }

        std::vector<inversion_t*>& invs = invs_by_chr[contig_name];
        for (int i = 0; i < invs.size(); i += BLOCK_SIZE) {
            std::vector<inversion_t*> block_invs(invs.begin() + i, invs.begin() + std::min(i + BLOCK_SIZE, static_cast<int>(invs.size())));
            std::future<void> future = thread_pool.push(genotype_invs, contig_name, chr_seqs.get_seq(contig_name),
                    chr_seqs.get_len(contig_name), block_invs, in_vcf_header, out_vcf_header, stats, config);
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

    // print contigs in vcf order
    int n_seqs;
    const char** seqnames = bcf_hdr_seqnames(in_vcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; i++) {
    	std::string contig_name = seqnames[i];
    	std::vector<sv_t*> contig_svs;
    	if (dels_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dels_by_chr[contig_name].begin(), dels_by_chr[contig_name].end());
    	if (dups_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dups_by_chr[contig_name].begin(), dups_by_chr[contig_name].end());
    	if (inss_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), inss_by_chr[contig_name].begin(), inss_by_chr[contig_name].end());
        if (invs_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), invs_by_chr[contig_name].begin(), invs_by_chr[contig_name].end());
    	std::sort(contig_svs.begin(), contig_svs.end(), [](const sv_t* sv1, const sv_t* sv2) {return sv1->start < sv2->start;});

		for (auto& sv : contig_svs) {
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AC", NULL, 0);
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AN", NULL, 0);
            update_record(in_vcf_header, out_vcf_header, sv, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), imap[0]);
			if (bcf_write(out_vcf_file, out_vcf_header, sv->vcf_entry) != 0) {
				throw std::runtime_error("Failed to write VCF record to " + out_vcf_fname);
			}
		}
    }
    delete[] seqnames;

    bcf_destroy(vcf_record);
    bcf_hdr_destroy(out_vcf_header);
    bcf_hdr_destroy(in_vcf_header);
    bcf_close(in_vcf_file);
    bcf_close(out_vcf_file);
}
