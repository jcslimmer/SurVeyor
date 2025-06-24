#ifndef GENOTYPE_INSS_H
#define GENOTYPE_INSS_H

#include "stat_tests.h"

#include "genotype.h"

void genotype_ins(insertion_t* ins, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t* evidence_logger,
                bool reassign_evidence, std::unordered_map<std::string, std::string>& reads_to_sv_map) {
    
    hts_pos_t ins_start = ins->start, ins_end = ins->end;

	hts_pos_t extend = stats.read_len + 20;

	// build alt allele
	/*
	 * POS in VCF is the base BEFORE the insertion
	 * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" insertion (no deletion),POS == END, otherwise the last base deleted
	 * As usual, in order to make intervals [ ), we increase the coordinates by 1
	 */
	ins_start++; ins_end++;

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
    std::vector<std::shared_ptr<bam1_t>> alt_bp1_better_seqs, alt_bp2_better_seqs, ref_bp1_better_seqs, ref_bp2_better_seqs;

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
            if (reassign_evidence && reads_to_sv_map[bam_get_qname(read)] != ins->id) continue;
            if (alt1_aln.sw_score >= alt2_aln.sw_score) {
                alt_bp1_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            } 
            if (alt1_aln.sw_score <= alt2_aln.sw_score) {
                alt_bp2_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            }
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                ref_bp1_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            } 
            if (ref1_aln.sw_score <= ref2_aln.sw_score) {
                ref_bp2_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
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
    auto alt_bp1_better_seqs_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_seqs, std::vector<bool>(), alt_bp1_consensus_seq, alt_bp1_avg_score, alt_bp1_stddev_score);
    delete[] alt_bp1_seq;
    auto alt_bp2_better_seqs_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_seqs, std::vector<bool>(), alt_bp2_consensus_seq, alt_bp2_avg_score, alt_bp2_stddev_score);
    delete[] alt_bp2_seq;

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    auto ref_bp1_better_seqs_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    auto ref_bp2_better_seqs_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score);
    delete[] ref_bp2_seq;

    if (evidence_logger) evidence_logger->log_reads_associations(ins->id, alt_bp1_better_seqs);
    if (evidence_logger) evidence_logger->log_reads_associations(ins->id, alt_bp2_better_seqs);

    if (alt_bp1_consensus_seq.length() >= 2*config.min_clip_len) {
        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_bp1_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp1_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, ins->start-stats.max_is, ins->start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, ins->start, ins->start+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->sample_info.alt_lext_reads = alt_bp1_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads = alt_bp1_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads = alt_bp1_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads = alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;
    }
    if (alt_bp2_consensus_seq.length() >= 2*config.min_clip_len) {
        std::shared_ptr<consensus_t> alt_bp2_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp2_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, ins->end-stats.max_is, ins->end, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, ins->end, ins->end+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->sample_info.alt_lext_reads += alt_bp2_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads += alt_bp2_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads += alt_bp2_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads += alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;
    }

    ins->sample_info.ext_alt_consensus1_length = alt_bp1_consensus_seq.length();
    ins->sample_info.ext_alt_consensus2_length = alt_bp2_consensus_seq.length();

    ref1_aln.Clear();
    alt1_aln.Clear();
    if (alt_bp1_consensus_seq.length() >= 2*config.min_clip_len) {
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
    if (alt_bp2_consensus_seq.length() >= 2*config.min_clip_len) {
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
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq_chr,
    std::vector<double>* global_crossing_isize_dist, evidence_logger_t* evidence_logger,
    bool reassign_evidence, std::unordered_map<std::string, std::string>& reads_to_sv_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

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
        genotype_ins(ins, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_logger, reassign_evidence, reads_to_sv_map);
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    depth_filter_ins(contig_name, inss, bam_file, config, stats);
    calculate_ptn_ratio(contig_name, inss, bam_file, config, stats, evidence_logger, reassign_evidence, reads_to_sv_map, *mateseqs_w_mapq_chr);
    std::vector<sv_t*> inss_sv(inss.begin(), inss.end());
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, inss_sv, bam_file, config, stats, config.min_sv_size, true);

    release_mates(contig_id);
}

#endif // GENOTYPE_INSS_H