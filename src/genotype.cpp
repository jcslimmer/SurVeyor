#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "assemble.h"
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
    }
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

void update_record_bp_consensus_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_consensus_info_t bp_consensus_info, std::string prefix, int bp_number) {
    update_record_bp_reads_info(out_hdr, b, bp_consensus_info.reads_info, prefix, bp_number);
    update_record_bp_pairs_info(out_hdr, b, bp_consensus_info.pairs_info, prefix, bp_number);
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
    bcf_update_genotypes(out_hdr, sv->vcf_entry, sv->sample_info.gt, sv->n_gt);

    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.alt_bp1, "A", 1);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.alt_bp2, "A", 2);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.ref_bp1, "R", 1);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.ref_bp2, "R", 2);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.pairs_crossing_midpoint, "M", 1);
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
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL", &(sv->sample_info.ext_alt_consensus1_length), 1);

    if (sv->sample_info.ext_alt_consensus1_length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS", &(sv->sample_info.ext_alt_consensus1_to_alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS", &(sv->sample_info.ext_alt_consensus1_to_ref_score), 1);
        int exss[] = {sv->sample_info.alt_consensus1_split_size1, sv->sample_info.alt_consensus1_split_size2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS", exss, 2);
        int exssc[] = {sv->sample_info.alt_consensus1_split_score1, sv->sample_info.alt_consensus1_split_score2};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC", exssc, 2);
        int exsscia[] = {sv->sample_info.alt_consensus1_split_score1_ind_aln, sv->sample_info.alt_consensus1_split_score2_ind_aln};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSCIA", exsscia, 2);

        if (sv->svtype() == "INS") {
            bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS2", &(sv->sample_info.ext_alt_consensus2_to_alt_score), 1);
            bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS2", &(sv->sample_info.ext_alt_consensus2_to_ref_score), 1);
            bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL2", &(sv->sample_info.ext_alt_consensus2_length), 1);
            int exss2[] = {sv->sample_info.alt_consensus2_split_size1, sv->sample_info.alt_consensus2_split_size2};
            bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSS2", exss2, 2);
            int exssc2[] = {sv->sample_info.alt_consensus2_split_score1, sv->sample_info.alt_consensus2_split_score2};
            bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC2", exssc2, 2);
            int exssc2ia[] = {sv->sample_info.alt_consensus2_split_score1_ind_aln, sv->sample_info.alt_consensus2_split_score2_ind_aln};
            bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXSSC2IA", exssc2ia, 2);
        }
    }
}

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<bam1_t*>& consistent_reads, 
    double consistent_avg_score, double consistent_stddev_score) {
    
    bp_reads_info.computed = true;
    bp_reads_info.reads = n_reads;

    std::vector<int> mqs;

    double sum_mq = 0;
    for (bam1_t* read : consistent_reads) {
        if (bam_is_mrev(read)) {
            bp_reads_info.consistent_fwd++;
        } else {
            bp_reads_info.consistent_rev++;
        }
        int mq = get_mq(read);
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
}

std::vector<bam1_t*> find_consistent_seqs_subset(std::string ref_seq, std::vector<bam1_t*>& reads, std::string& consensus_seq, double& avg_score, double& stddev_score) {

    if (reads.empty()) {
        avg_score = 0;
        stddev_score = 0;
        consensus_seq = "";
        return reads;
    }

    std::vector<std::string> seqs;
    for (bam1_t* read : reads) {
        std::string seq = get_sequence(read, true);
        if (!bam_is_mrev(read)) rc(seq);
        seqs.push_back(seq);
    }

    std::vector<bam1_t*> consistent_reads;
    avg_score = 0;
    std::vector<std::string> temp1, temp2;
    std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
    
    std::vector<std::string> consensus_seqs = generate_reference_guided_consensus(ref_seq, temp1, seqs, temp2, aligner, harsh_aligner, consensus_contigs_alns, config, stats, false);
    
    std::vector<seq_w_pp_t> seqs_w_pp, temp3, temp4;
    for (std::string& seq : seqs) {
        seqs_w_pp.push_back({seq, true, true});
    }
    std::vector<std::string> consensus_seqs2 = assemble_reads(temp3, seqs_w_pp, temp4, harsh_aligner, config, stats);
    consensus_seqs.insert(consensus_seqs.end(), consensus_seqs2.begin(), consensus_seqs2.end());

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    std::vector<int> start_positions, end_positions;
    std::vector<StripedSmithWaterman::Alignment> alns;
    double cum_score = 0;
    std::vector<double> aln_scores;
    for (std::string cseq : consensus_seqs) {
        std::vector<bam1_t*> curr_consistent_reads;
        std::vector<int> curr_start_positions, curr_end_positions;
        std::vector<StripedSmithWaterman::Alignment> curr_alns;
        double curr_cum_score = 0;
        std::vector<double> curr_aln_scores;
        for (bam1_t* read : reads) {
            std::string seq = get_sequence(read, true);
            if (!bam_is_mrev(read)) rc(seq);
            
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
            del->sample_info.too_deep = true;
            break;
        }
    }

    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seq, alt_better_reads, alt_consensus_seq, alt_avg_score, alt_stddev_score);
    std::vector<bam1_t*> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score);
    std::vector<bam1_t*> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score);

    if (!alt_consensus_seq.empty()) {
       // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, del->start-stats.max_is, del->start, chr_seqs.get_len(del->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, del->end, del->end+stats.max_is, chr_seqs.get_len(del->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);

        del->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        del->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        del->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        del->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
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
        auto query_lh_aln_score = find_aln_prefix_score(alt_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        del->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second - get_left_clip_size(alt_aln);
        del->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second - get_right_clip_size(alt_aln);
        del->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        del->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        del->left_anchor_aln->start = del->start - lf_aln_rlen;
        del->left_anchor_aln->end = del->start;
        del->left_anchor_aln->seq_len = lf_aln_rlen;
        del->right_anchor_aln->start = del_end;
        del->right_anchor_aln->end = del_end + rf_aln_rlen;
        del->right_anchor_aln->seq_len = rf_aln_rlen;

        // align to ref
        hts_pos_t lbp_start = lh_start, lbp_end = del->start + alt_consensus_seq.length();
        hts_pos_t rbp_start = del->end - alt_consensus_seq.length(), rbp_end = rh_end;
        if (lbp_end > contig_len) lbp_end = contig_len;
        if (rbp_start < 0) rbp_start = 0;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+lbp_start, lbp_end-lbp_start, filter, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+rbp_start, rbp_end-rbp_start, filter, &ref2_aln, 0);
    
        del->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();
        del->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.sw_score;
        del->sample_info.ext_alt_consensus1_to_ref_score = std::max(ref1_aln.sw_score, ref2_aln.sw_score);

        ref1_aln.Clear();
        std::string lh_query = alt_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+lbp_start, lbp_end-lbp_start, filter, &ref1_aln, 0);
        del->sample_info.alt_consensus1_split_score1_ind_aln = ref1_aln.sw_score;

        ref2_aln.Clear();
        std::string rh_query = alt_consensus_seq.substr(alt_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+rbp_start, rbp_end-rbp_start, filter, &ref2_aln, 0);
        del->sample_info.alt_consensus1_split_score2_ind_aln = ref2_aln.sw_score;
    }

    set_bp_consensus_info(del->sample_info.alt_bp1.reads_info, alt_better_reads.size(), alt_better_reads_consistent, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(del->sample_info.ref_bp1.reads_info, ref_bp1_better_seqs.size(), ref_bp1_better_seqs_consistent, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(del->sample_info.ref_bp2.reads_info, ref_bp2_better_seqs.size(), ref_bp2_better_seqs_consistent, ref_bp2_avg_score, ref_bp2_stddev_score);
    
    del->sample_info.alt_ref_equal_reads = same;

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
    count_stray_pairs(contig_name, dels, bam_file, config, stats);
}

void genotype_small_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {

	hts_pos_t dup_start = dup->start, dup_end = dup->end;
    char* contig_seq = chr_seqs.get_seq(dup->chr);
	hts_pos_t contig_len = chr_seqs.get_len(dup->chr);

	hts_pos_t extend = stats.read_len + 20;
	hts_pos_t svlen = dup->svlen();

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;

    hts_pos_t ref_start = std::max(hts_pos_t(0), dup_start-extend), ref_end = std::min(dup_end+extend, contig_len);
	hts_pos_t ref_len = ref_end - ref_start;
	char* ref_seq = new char[ref_len + 1];
	strncpy(ref_seq, contig_seq+ref_start, ref_len);
	ref_seq[ref_len] = 0;

    std::vector<char*> alt_seqs;
	for (int copies = 1; copies*svlen < stats.read_len; copies++) {
		int alt_len = ref_len + copies*svlen;

		char* alt_seq = new char[alt_len+1];
		int pos = 0;
		strncpy(alt_seq, contig_seq+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < copies; i++) {
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
			strncpy(alt_seq+pos, contig_seq+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
		}
		strncpy(alt_seq+pos, contig_seq+dup_end, ref_end-dup_end);
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
        dup->sample_info.too_deep = true;
    }

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    double alt_stddev_score, ref_stddev_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seqs[alt_with_most_reads], alt_better_reads[alt_with_most_reads], alt_consensus_seq, alt_avg_score, alt_stddev_score);
    std::vector<bam1_t*> ref_better_reads_consistent = find_consistent_seqs_subset(ref_seq, ref_better_reads, ref_consensus_seq, ref_avg_score, ref_stddev_score);

    if (!alt_consensus_seq.empty()) {
       // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, dup->start-stats.max_is, dup->start, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->end, dup->end+stats.max_is, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        dup->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        dup->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        dup->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        dup->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        delete alt_consensus;

        hts_pos_t ref_start = dup->start-alt_consensus_seq.length();
        if (ref_start < 0) ref_start = 0;
        hts_pos_t ref_end = dup->end+alt_consensus_seq.length();
        if (ref_end > contig_len) ref_end = contig_len;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_start, ref_end-ref_start, filter, &ref_aln, 0);
        dup->sample_info.ext_alt_consensus1_to_ref_score = ref_aln.sw_score;

        int n_extra_copies = alt_with_most_reads+1;
        int alt_len = ref_end - ref_start + n_extra_copies*svlen;
		char* alt_seq = new char[alt_len+1];
		int pos = 0;
		strncpy(alt_seq, contig_seq+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < n_extra_copies; i++) {
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
			strncpy(alt_seq+pos, contig_seq+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
		}
		strncpy(alt_seq+pos, contig_seq+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);
        dup->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.sw_score;

        // delete[] alt_seq;

        int lf_seq_end = dup_start - ref_start;
        int lf_aln_rlen = std::max(0, lf_seq_end-alt_aln.ref_begin);
        int rf_seq_len = ref_end - dup_end;
        int rf_seq_start = pos - rf_seq_len;
        int rf_aln_rlen = std::max(0, alt_aln.ref_end-rf_seq_start);

        int temp;
        auto query_lh_aln_score = find_aln_prefix_score(alt_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        dup->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second - get_left_clip_size(alt_aln);
        dup->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second - get_right_clip_size(alt_aln);
        dup->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        dup->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        dup->left_anchor_aln->start = dup_end - lf_aln_rlen;
        dup->left_anchor_aln->end = dup_end;
        dup->left_anchor_aln->seq_len = lf_aln_rlen;
        dup->right_anchor_aln->start = dup_start;
        dup->right_anchor_aln->end = dup_start + rf_aln_rlen;
        dup->right_anchor_aln->seq_len = rf_aln_rlen;

        dup->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();

        ref_aln.Clear();
        std::string lh_query = alt_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+ref_start, ref_end-ref_start, filter, &ref_aln, 0);
        dup->sample_info.alt_consensus1_split_score1_ind_aln = ref_aln.sw_score;

        ref_aln.Clear();
        std::string rh_query = alt_consensus_seq.substr(alt_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+ref_start, ref_end-ref_start, filter, &ref_aln, 0);
        dup->sample_info.alt_consensus1_split_score2_ind_aln = ref_aln.sw_score;
    }

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads[alt_with_most_reads].size(), alt_better_reads_consistent, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_better_reads.size(), ref_better_reads_consistent, ref_avg_score, ref_stddev_score);
    
    dup->sample_info.alt_ref_equal_reads = same;

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
            dup->sample_info.too_deep = true;
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
    double alt_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seq, alt_better_reads, alt_consensus_seq, alt_avg_score, alt_stddev_score);
    std::vector<bam1_t*> ref_bp1_better_reads_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_reads, ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score);
    std::vector<bam1_t*> ref_bp2_better_reads_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_reads, ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score);

    if (!alt_consensus_seq.empty()) {
       // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, dup->end-stats.max_is, dup->end, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->start, dup->start+stats.max_is, chr_seqs.get_len(dup->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        dup->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        dup->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        dup->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        dup->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
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

        // length of the left and right flanking regions of the deletion covered by alt_consensus_seq
        int lf_aln_rlen = std::max(hts_pos_t(0), lh_len - alt_aln.ref_begin);
        int rf_aln_rlen = std::max(hts_pos_t(0), alt_aln.ref_end - lh_len);

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        int temp;
        auto query_lh_aln_score = find_aln_prefix_score(alt_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        dup->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second - get_left_clip_size(alt_aln);
        dup->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second - get_right_clip_size(alt_aln);
        dup->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        dup->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        dup->left_anchor_aln->start = dup->end - lf_aln_rlen;
        dup->left_anchor_aln->end = dup->end;
        dup->left_anchor_aln->seq_len = lf_aln_rlen;
        dup->right_anchor_aln->start = dup->start;
        dup->right_anchor_aln->end = dup->start + rf_aln_rlen;
        dup->right_anchor_aln->seq_len = rf_aln_rlen;

        // align to ref
        hts_pos_t ref_bp1_start = dup->start - alt_consensus_seq.length(), ref_bp1_end = dup->start + alt_consensus_seq.length();
        if (ref_bp1_start < 0) ref_bp1_start = 0;
        if (ref_bp1_end > contig_len) ref_bp1_end = contig_len;
        hts_pos_t ref_bp2_start = dup->end - alt_consensus_seq.length(), ref_bp2_end = dup->end + alt_consensus_seq.length();
        if (ref_bp2_start < 0) ref_bp2_start = 0;
        if (ref_bp2_end > contig_len) ref_bp2_end = contig_len;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &ref2_aln, 0);

        dup->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();
        dup->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.sw_score;
        dup->sample_info.ext_alt_consensus1_to_ref_score = std::max(ref1_aln.sw_score, ref2_aln.sw_score);

        ref1_aln.Clear();
        std::string lh_query = alt_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &ref1_aln, 0);
        dup->sample_info.alt_consensus1_split_score1_ind_aln = ref1_aln.sw_score;

        ref2_aln.Clear();
        std::string rh_query = alt_consensus_seq.substr(alt_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &ref2_aln, 0);
        dup->sample_info.alt_consensus1_split_score2_ind_aln = ref2_aln.sw_score;
    }

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads.size(), alt_better_reads_consistent, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_bp1_better_reads.size(), ref_bp1_better_reads_consistent, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp2.reads_info, ref_bp2_better_reads.size(), ref_bp2_better_reads_consistent, ref_bp2_avg_score, ref_bp2_stddev_score);
    
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
    count_stray_pairs(contig_name, dups, bam_file, config, stats);
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
            ins->sample_info.too_deep = true;
            break;
        }
    }

    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_bp1_stddev_score, alt_bp2_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bam1_t*> alt_bp1_better_seqs_consistent = find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_seqs, alt_bp1_consensus_seq, alt_bp1_avg_score, alt_bp1_stddev_score);
    delete[] alt_bp1_seq;
    std::vector<bam1_t*> alt_bp2_better_seqs_consistent = find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_seqs, alt_bp2_consensus_seq, alt_bp2_avg_score, alt_bp2_stddev_score);
    delete[] alt_bp2_seq;

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    std::vector<bam1_t*> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    std::vector<bam1_t*> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score);
    delete[] ref_bp2_seq;

    if (!alt_bp1_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_bp1_consensus = new consensus_t(false, 0, 0, 0, alt_bp1_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, ins->start-stats.max_is, ins->start, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, ins->start, ins->start+stats.max_is, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->sample_info.alt_lext_reads = alt_bp1_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads = alt_bp1_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads = alt_bp1_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads = alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;
        delete alt_bp1_consensus;
    }
    if (!alt_bp2_consensus_seq.empty()) {
        consensus_t* alt_bp2_consensus = new consensus_t(false, 0, 0, 0, alt_bp2_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, ins->end-stats.max_is, ins->end, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, ins->end, ins->end+stats.max_is, chr_seqs.get_len(ins->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->sample_info.alt_lext_reads += alt_bp2_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads += alt_bp2_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads += alt_bp2_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads += alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;
        delete alt_bp2_consensus;
    }

    ins->sample_info.ext_alt_consensus1_length = alt_bp1_consensus_seq.length();
    ins->sample_info.ext_alt_consensus2_length = alt_bp2_consensus_seq.length();

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

        // length of the left and right flanking regions of the deletion covered by alt_consensus_seq
        hts_pos_t lh_len = ins->start - ref_bp1_start;
        int lf_aln_rlen = std::max(hts_pos_t(0), lh_len - alt1_aln.ref_begin);
        int rf_aln_rlen = std::max(hts_pos_t(0), alt1_aln.ref_end - lh_len);

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        int temp;
        auto query_lh_aln_score = find_aln_prefix_score(alt1_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt1_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        ins->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second - get_left_clip_size(alt1_aln);
        ins->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second - get_right_clip_size(alt1_aln);
        ins->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        ins->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        ins->left_anchor_aln->start = ins->start - lf_aln_rlen;
        ins->left_anchor_aln->end = ins->start;
        ins->left_anchor_aln->seq_len = lf_aln_rlen;

        StripedSmithWaterman::Alignment aln;

        aln.Clear();
        std::string lh_query = alt_bp1_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &aln, 0);
        ins->sample_info.alt_consensus1_split_score1_ind_aln = aln.sw_score;

        aln.Clear();
        std::string rh_query = alt_bp1_consensus_seq.substr(alt_bp1_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter, &aln, 0);
        ins->sample_info.alt_consensus1_split_score2_ind_aln = aln.sw_score;
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

        // length of the left and right flanking regions of the deletion covered by alt_consensus_seq
        hts_pos_t lh_len = extra_len + ins_seq_portion_len;
        int lf_aln_rlen = std::max(hts_pos_t(0), lh_len - alt2_aln.ref_begin);
        int rf_aln_rlen = std::max(hts_pos_t(0), alt2_aln.ref_end - lh_len);

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        int temp;
        auto query_lh_aln_score = find_aln_prefix_score(alt2_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt2_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        ins->sample_info.alt_consensus2_split_size1 = query_lh_aln_score.second - get_left_clip_size(alt2_aln);
        ins->sample_info.alt_consensus2_split_size2 = query_rh_aln_score.second - get_right_clip_size(alt2_aln);
        ins->sample_info.alt_consensus2_split_score1 = query_lh_aln_score.first;
        ins->sample_info.alt_consensus2_split_score2 = query_rh_aln_score.first;

        ins->right_anchor_aln->start = ins->end;
        ins->right_anchor_aln->end = ins->end + rf_aln_rlen;
        ins->right_anchor_aln->seq_len = rf_aln_rlen;

        StripedSmithWaterman::Alignment aln;

        aln.Clear();
        std::string lh_query = alt_bp2_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &aln, 0);
        ins->sample_info.alt_consensus2_split_score1_ind_aln = aln.sw_score;

        aln.Clear();
        std::string rh_query = alt_bp2_consensus_seq.substr(alt_bp2_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter, &aln, 0);
        ins->sample_info.alt_consensus2_split_score2_ind_aln = aln.sw_score;
    }

    ins->sample_info.ext_alt_consensus1_to_ref_score = ref1_aln.sw_score;
    ins->sample_info.ext_alt_consensus2_to_ref_score = ref2_aln.sw_score;
    ins->sample_info.ext_alt_consensus1_to_alt_score = alt1_aln.sw_score;
    ins->sample_info.ext_alt_consensus2_to_alt_score = alt2_aln.sw_score;

    set_bp_consensus_info(ins->sample_info.alt_bp1.reads_info, alt_bp1_better_seqs.size(), alt_bp1_better_seqs_consistent, alt_bp1_avg_score, alt_bp1_stddev_score);
    set_bp_consensus_info(ins->sample_info.alt_bp2.reads_info, alt_bp2_better_seqs.size(), alt_bp2_better_seqs_consistent, alt_bp2_avg_score, alt_bp2_stddev_score);
    set_bp_consensus_info(ins->sample_info.ref_bp1.reads_info, ref_bp1_better_seqs.size(), ref_bp1_better_seqs_consistent, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(ins->sample_info.ref_bp2.reads_info, ref_bp2_better_seqs.size(), ref_bp2_better_seqs_consistent, ref_bp2_avg_score, ref_bp2_stddev_score);
    
    ins->sample_info.alt_ref_equal_reads = same;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
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
    find_discordant_pairs(contig_name, inss, bam_file, stats, mateseqs_w_mapq[contig_id], harsh_aligner, config);

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
            inv->sample_info.too_deep = true;
            break;
        }
    }

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    double alt_stddev_score, ref_stddev_score;
    std::vector<bam1_t*> alt_better_reads_consistent = find_consistent_seqs_subset(alt_seq, alt_better_seqs, alt_consensus_seq, alt_avg_score, alt_stddev_score);
    std::vector<bam1_t*> ref_better_reads_consistent = find_consistent_seqs_subset(ref_seq, ref_better_seqs, ref_consensus_seq, ref_avg_score, ref_stddev_score);

    if (!alt_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_consensus = new consensus_t(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, inv->start-stats.max_is, inv->start, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, inv->start, inv->start+stats.max_is, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        inv->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        inv->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        inv->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        inv->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
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

        inv->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();
        inv->sample_info.ext_alt_consensus1_to_ref_score = ref_aln.sw_score;
        inv->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.sw_score;
    }

    set_bp_consensus_info(inv->sample_info.alt_bp1.reads_info, alt_better_seqs.size(), alt_better_reads_consistent, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(inv->sample_info.ref_bp1.reads_info, ref_better_seqs.size(), ref_better_reads_consistent, ref_avg_score, ref_stddev_score);
    
    inv->sample_info.alt_ref_equal_reads = same;

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
            inv->sample_info.too_deep = true;
            break;
        }
    }

    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_bp1_stddev_score, alt_bp2_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bam1_t*> alt_bp1_better_reads_consistent = find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_reads, alt_bp1_consensus_seq, alt_bp1_avg_score, alt_bp1_stddev_score);
    delete[] alt_bp1_seq;
    std::vector<bam1_t*> alt_bp2_better_reads_consistent = find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_reads, alt_bp2_consensus_seq, alt_bp2_avg_score, alt_bp2_stddev_score);
    delete[] alt_bp2_seq;

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    std::vector<bam1_t*> ref_bp1_better_reads_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_reads, ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    std::vector<bam1_t*> ref_bp2_better_reads_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_reads, ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score);
    delete[] ref_bp2_seq;

    if (!alt_bp1_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_bp1_consensus = new consensus_t(false, 0, 0, 0, alt_bp1_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, inv->start-stats.max_is, inv->start, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, inv->start, inv->start+stats.max_is, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        inv->sample_info.alt_lext_reads = alt_bp1_consensus->left_ext_reads;
        inv->sample_info.alt_rext_reads = alt_bp1_consensus->right_ext_reads;
        inv->sample_info.hq_alt_lext_reads = alt_bp1_consensus->hq_left_ext_reads;
        inv->sample_info.hq_alt_rext_reads = alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;
        delete alt_bp1_consensus;
    }
    if (!alt_bp2_consensus_seq.empty()) {
        // all we care about is the consensus sequence
        consensus_t* alt_bp2_consensus = new consensus_t(false, 0, 0, 0, alt_bp2_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, inv->end-stats.max_is, inv->end, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, inv->end, inv->end+stats.max_is, chr_seqs.get_len(inv->chr), config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        inv->sample_info.alt_lext_reads += alt_bp2_consensus->left_ext_reads;
        inv->sample_info.alt_rext_reads += alt_bp2_consensus->right_ext_reads;
        inv->sample_info.hq_alt_lext_reads += alt_bp2_consensus->hq_left_ext_reads;
        inv->sample_info.hq_alt_rext_reads += alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;
        delete alt_bp2_consensus;
    }
    inv->sample_info.ext_alt_consensus1_length = alt_bp1_consensus_seq.length() + alt_bp2_consensus_seq.length();

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

    inv->sample_info.ext_alt_consensus1_to_ref_score = ref1_aln.sw_score + ref2_aln.sw_score;
    inv->sample_info.ext_alt_consensus1_to_alt_score = alt1_aln.sw_score + alt2_aln.sw_score;

    set_bp_consensus_info(inv->sample_info.alt_bp1.reads_info, alt_bp1_better_reads.size(), alt_bp1_better_reads_consistent, alt_bp1_avg_score, alt_bp1_stddev_score);
    set_bp_consensus_info(inv->sample_info.alt_bp2.reads_info, alt_bp2_better_reads.size(), alt_bp2_better_reads_consistent, alt_bp2_avg_score, alt_bp2_stddev_score);
    set_bp_consensus_info(inv->sample_info.ref_bp1.reads_info, ref_bp1_better_reads.size(), ref_bp1_better_reads_consistent, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(inv->sample_info.ref_bp2.reads_info, ref_bp2_better_reads.size(), ref_bp2_better_reads_consistent, ref_bp2_avg_score, ref_bp2_stddev_score);

    inv->sample_info.alt_ref_equal_reads = same;

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

    calculate_ptn_ratio(contig_name, invs, bam_file, config, stats);
    depth_filter_inv(contig_name, invs, bam_file, config, stats);

    for (inversion_t* inv : invs) {
        double ptn1 = inv->sample_info.alt_bp1.pairs_info.pairs/double(inv->sample_info.ref_bp1.pairs_info.pairs);
        double ptn2 = inv->sample_info.alt_bp2.pairs_info.pairs/double(inv->sample_info.ref_bp2.pairs_info.pairs);
        double ptn = std::min(ptn1, ptn2);
        inv->n_gt = 2;
        inv->sample_info.gt = new int[2];
        if (ptn >= 0.75) {
            inv->sample_info.gt[0] = bcf_gt_unphased(1);
            inv->sample_info.gt[1] = bcf_gt_unphased(1);
        } else if (ptn <= 0.25) {
            inv->sample_info.gt[0] = bcf_gt_unphased(0);
            inv->sample_info.gt[1] = bcf_gt_unphased(0);
        } else {
            inv->sample_info.gt[0] = bcf_gt_unphased(0);
            inv->sample_info.gt[1] = bcf_gt_unphased(1);
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

    // apply hard filters to inversions
    int n_seqs;
    const char** seqnames = bcf_hdr_seqnames(in_vcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; i++) {
        std::string contig_name = seqnames[i];
        for (inversion_t* inv : invs_by_chr[contig_name]) {
            std::vector<std::string> filters;
            if (inv->sample_info.left_flanking_cov > stats.get_max_depth(inv->chr) || inv->sample_info.right_flanking_cov > stats.get_max_depth(inv->chr) ||
                inv->sample_info.left_flanking_cov < stats.get_min_depth(inv->chr) || inv->sample_info.right_flanking_cov < stats.get_min_depth(inv->chr) ||
                inv->sample_info.left_anchor_cov > stats.get_max_depth(inv->chr) || inv->sample_info.right_anchor_cov > stats.get_max_depth(inv->chr)) {
                inv->sample_info.filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (inv->sample_info.alt_bp1.reads_info.consistent_reads() > stats.get_max_depth(inv->chr) 
             || inv->sample_info.alt_bp2.reads_info.consistent_reads() > stats.get_max_depth(inv->chr)) {
                inv->sample_info.filters.push_back("ANOMALOUS_SC_NUMBER");
            }
            if (inv->sample_info.alt_bp1.pairs_info.pos_max_mq < config.high_confidence_mapq || inv->sample_info.alt_bp1.pairs_info.neg_max_mq < config.high_confidence_mapq
             || inv->sample_info.alt_bp2.pairs_info.pos_max_mq < config.high_confidence_mapq || inv->sample_info.alt_bp2.pairs_info.neg_max_mq < config.high_confidence_mapq) {
                inv->sample_info.filters.push_back("LOW_MAPQ_DISC_PAIRS");
            }

            if (inv->sample_info.alt_bp1.pairs_info.pairs < stats.get_min_disc_pairs_by_insertion_size(inv->svlen())/2) {
                inv->sample_info.filters.push_back("NOT_ENOUGH_DISC_PAIRS");
            }

            double ptn_ratio_bp1 = double(inv->sample_info.alt_bp1.pairs_info.pairs)/(inv->sample_info.alt_bp1.pairs_info.pairs+inv->sample_info.ref_bp1.pairs_info.pairs);
            double ptn_ratio_bp2 = double(inv->sample_info.alt_bp2.pairs_info.pairs)/(inv->sample_info.alt_bp2.pairs_info.pairs+inv->sample_info.ref_bp2.pairs_info.pairs);
            if (ptn_ratio_bp1 < 0.25 || ptn_ratio_bp2 < 0.25) {
                inv->sample_info.filters.push_back("LOW_PTN_RATIO");
            }

            if (inv->right_anchor_aln->end-inv->left_anchor_aln->start < stats.max_is/2 || inv->rbp_right_anchor_aln->end-inv->rbp_left_anchor_aln->start < stats.max_is/2) {
                inv->sample_info.filters.push_back("SHORT_ANCHOR");
            }

            if (inv->sample_info.filters.empty()) {
                inv->sample_info.filters.push_back("PASS");
            }
        }
    }

    // print contigs in vcf order
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
