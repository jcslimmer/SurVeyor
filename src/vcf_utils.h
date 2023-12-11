#ifndef VCF_UTILS_H
#define VCF_UTILS_H

#include <chrono>
#include <ctime>
#include <htslib/vcf.h>
#include "types.h"
#include "utils.h"

bcf_hrec_t* generate_contig_hrec() {
	bcf_hrec_t* contig_hrec = new bcf_hrec_t;
	contig_hrec->type = BCF_HL_CTG;
	contig_hrec->key = strdup("contig");
	contig_hrec->value = NULL;
	contig_hrec->keys = contig_hrec->vals = NULL;
	contig_hrec->nkeys = 0;
	int r1 = bcf_hrec_add_key(contig_hrec, "ID", 2);
	int r2 = bcf_hrec_add_key(contig_hrec, "length", 6);
	if (r1 || r2) {
		throw std::runtime_error("Failed to create contig to VCF header.");
	}
	return contig_hrec;
}
bcf_hdr_t* generate_vcf_header_base(chr_seqs_map_t& contigs, std::string sample_name, config_t config, std::string command) {
	bcf_hdr_t* header = bcf_hdr_init("w");

	// add contigs
	for (std::string contig_name : contigs.ordered_contigs) {
		bcf_hrec_t* hrec = generate_contig_hrec();
		int r1 = bcf_hrec_set_val(hrec, 0, contig_name.c_str(), contig_name.length(), false);
		std::string len_str = std::to_string(contigs.get_len(contig_name));
		int r2 = bcf_hrec_set_val(hrec, 1, len_str.c_str(), len_str.length(), false);
		if (r1 || r2) {
			throw std::runtime_error("Failed to create contig to VCF header.");
		}
		bcf_hdr_add_hrec(header, hrec);
	}

	int len;

	// add FILTER tags
	char size_flt_tag[1000];
	sprintf(size_flt_tag, "##FILTER=<ID=SMALL,Description=\"SV smaller than %d bp.\">", config.min_sv_size); // TODO: remove this filter
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	const char* remap_boundary_flt_tag = "##FILTER=<ID=REMAP_BOUNDARY_FILTER,Description=\"One of the breakpoints is incompatible with mate locations.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_boundary_flt_tag, &len));

	const char* failed_to_ext_flt_tag = "##FILTER=<ID=FAILED_TO_EXTEND,Description=\"No reads can extend the consensus.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, failed_to_ext_flt_tag, &len));

	const char* anomalous_fl_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_FLANKING_DEPTH,Description=\"The insertion region has anomalous depth.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_fl_depth_flt_tag, &len));

	const char* homopolymer_flt_tag = "##FILTER=<ID=HOMOPOLYMER_INSSEQ,Description=\"Inserted sequence is a homopolymer.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, homopolymer_flt_tag, &len));

	const char* anomalous_sc_flt_tag = "##FILTER=<ID=ANOMALOUS_SC_NUMBER,Description=\"The number of soft-clipped reads supporting this call is too large.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_sc_flt_tag, &len));

	// add INFO tags
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, end_tag, &len));

	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svlen_tag, &len));

	const char* svinslen_tag = "##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinslen_tag, &len));

	const char* svinsseq_tag = "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinsseq_tag, &len));

	const char* source_tag = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, source_tag, &len));

	const char* sr_tag = "##INFO=<ID=SPLIT_READS,Number=2,Type=Integer,Description=\"Split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sr_tag, &len));

	const char* fwd_sr_tag = "##INFO=<ID=FWD_SPLIT_READS,Number=2,Type=Integer,Description=\"Forward split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fwd_sr_tag, &len));

	const char* rev_sr_tag = "##INFO=<ID=REV_SPLIT_READS,Number=2,Type=Integer,Description=\"Reverse split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rev_sr_tag, &len));
	
	const char* rcc_ext_1sr_reads_tag = "##INFO=<ID=RCC_EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads extending a the right-clipped consensus to the left and to the right, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rcc_ext_1sr_reads_tag, &len));

	const char* lcc_ext_1sr_reads_tag = "##INFO=<ID=LCC_EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads extending a the left-clipped consensus to the left and to the right, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, lcc_ext_1sr_reads_tag, &len));

	const char* rcc_hq_ext_1sr_reads_tag = "##INFO=<ID=RCC_HQ_EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads with high MAPQ extending a the right-clipped consensus to the left and to the right, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rcc_hq_ext_1sr_reads_tag, &len));

	const char* lcc_hq_ext_1sr_reads_tag = "##INFO=<ID=LCC_HQ_EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads with high MAPQ extending a the left-clipped consensus to the left and to the right, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, lcc_hq_ext_1sr_reads_tag, &len));

	const char* fullj_score_tag = "##INFO=<ID=FULL_JUNCTION_SCORE,Number=1,Type=Integer,Description=\"Score of the best alignment of the full junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fullj_score_tag, &len));

	const char* fullj_cigar_tag = "##INFO=<ID=FULL_JUNCTION_CIGAR,Number=1,Type=String,Description=\"CIGAR of the best alignment of the full junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fullj_cigar_tag, &len));

	const char* splitj_score_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE,Number=2,Type=Integer,Description=\"Score of the best alignment of the left-half and right-half of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score_tag, &len));

	const char* splitj_score2_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE2,Number=2,Type=Integer,Description=\"Score of the second best alignment of the left-half and right-half of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score2_tag, &len));

	const char* splitj_size_tag = "##INFO=<ID=SPLIT_JUNCTION_SIZE,Number=2,Type=Integer,Description=\"Size of the the left-half and right-half of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_size_tag, &len));

	const char* splitj_cigar_tag = "##INFO=<ID=SPLIT_JUNCTION_CIGAR,Number=2,Type=String,Description=\"CIGAR of the best alignment of the left-half and right-half of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_cigar_tag, &len));

	const char* splitj_maprange_tag = "##INFO=<ID=SPLIT_JUNCTION_MAPPING_RANGE,Number=1,Type=String,Description=\"Mapping locations of the left-half and right-half of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_maprange_tag, &len));

	const char* max_mapq_tag = "##INFO=<ID=MAX_MAPQ,Number=2,Type=Integer,Description=\"Maximum MAPQ of clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_mapq_tag, &len));

	const char* remap_lb_tag = "##INFO=<ID=REMAP_LB,Number=1,Type=Integer,Description=\"Minimum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_lb_tag, &len));

	const char* remap_ub_tag = "##INFO=<ID=REMAP_UB,Number=1,Type=Integer,Description=\"Maximum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_ub_tag, &len));

	const char* disc_pairs_tag = "##INFO=<ID=DISC_PAIRS,Number=1,Type=Integer,Description=\"Discordant pairs supporting the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_tag, &len));

	const char* dp_max_mapq_tag = "##INFO=<ID=DISC_PAIRS_MAXMAPQ,Number=1,Type=Integer,Description=\"Maximum MAPQ of supporting discordant pairs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dp_max_mapq_tag, &len));

	const char* disc_pairs_hmapq_tag = "##INFO=<ID=DISC_PAIRS_HIGHMAPQ,Number=1,Type=Integer,Description=\"HDiscordant pairs with high MAPQ supporting the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_hmapq_tag, &len));

	// add ALT
	const char* del_alt_tag = "##ALT=<ID=DEL,Description=\"Deletion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, del_alt_tag, &len));

	const char* dup_alt_tag = "##ALT=<ID=DUP,Description=\"Tandem Duplication\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dup_alt_tag, &len));

	const char* ins_alt_tag = "##ALT=<ID=INS,Description=\"Insertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ins_alt_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, gt_tag, &len));

	const char* ft_tag = "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS indicates a reliable call. Any other value means the call is not reliable.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ft_tag, &len));

	std::string cmd_tag = "##SurVeyorCommand=" + command;
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cmd_tag.c_str(), &len));

	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::string version_tag = "##SurVeyorVersion=" + config.version + "; Date=" + std::ctime(&now_time);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, version_tag.c_str(), &len));

	std::stringstream called_by_ss;
	called_by_ss << "##calledBy=SurVeyor " << config.version << "; ";
	called_by_ss << "seed: " << config.seed << "; ";
	called_by_ss << "max-clipped-pos-dist: " << config.max_clipped_pos_dist << "; ";
	called_by_ss << "min-sv-size: " << config.min_sv_size << "; ";
	called_by_ss << "min-clip-len: " << config.min_clip_len << "; ";
	called_by_ss << "max-seq-error: " << config.max_seq_error << "; ";
	called_by_ss << "max-clipped-pos-dist: " << config.max_clipped_pos_dist << "; ";
	called_by_ss << "max-trans-size: " << config.max_trans_size << "; ";
	called_by_ss << "min-stable-mapq: " << config.min_stable_mapq << "; ";
	called_by_ss << "min-size-for-depth-filtering: " << config.min_size_for_depth_filtering << "; ";
	called_by_ss << "min-diff-hsr: " << config.min_diff_hsr << "; ";
	called_by_ss << "sampling-regions: " << (config.sampling_regions.empty() ? "no" : config.sampling_regions) << "; ";
	called_by_ss << "per-contig-stats: " << (config.per_contig_stats ? "true" : "false") << "; ";
	std::string called_by = called_by_ss.str();
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, called_by.c_str(), &len));

	return header;
}
bcf_hdr_t* generate_vcf_header(chr_seqs_map_t& contigs, std::string sample_name, config_t config, std::string command) {
	bcf_hdr_t* header = generate_vcf_header_base(contigs, sample_name, config, command);

	int len;

	// add INFO tags
	const char* overlap_tag = "##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Overlap (in bp) between the left and right contigs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, overlap_tag, &len));

	const char* mismatch_rate_tag = "##INFO=<ID=MISMATCH_RATE,Number=1,Type=Float,Description=\"Mismatch rate of overlap between the left and right contigs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mismatch_rate_tag, &len));

	const char* sd_tag = "##INFO=<ID=STABLE_DEPTHS,Number=2,Type=Integer,Description=\"Depths of the stable regions (in practice, the regions left and right of the insertion site).\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sd_tag, &len));

	const char* algo_tag = "##INFO=<ID=ALGORITHM,Number=1,Type=String,Description=\"Algorithm used to report the call.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, algo_tag, &len));

	// add samples
	bcf_hdr_add_sample(header, sample_name.c_str());

	return header;
}

void sv2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, sv_t* sv, char* chr_seq) {
	bcf_clear(bcf_entry);
	
	bcf_entry->rid = bcf_hdr_name2id(hdr, sv->chr.c_str());
	bcf_entry->pos = sv->start;
	bcf_update_id(hdr, bcf_entry, sv->id.c_str());
	
	std::string alleles = std::string(1, chr_seq[sv->start]) + ",<" + sv->svtype() + ">";
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	// add FILTERs
	for (std::string& filter : sv->filters) {
		int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter.c_str());
		bcf_add_filter(hdr, bcf_entry, filter_id);
	}
	
	// add INFO
	int int_conv = sv->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", sv->svtype().c_str());
	if (sv->ins_seq.find("-") == std::string::npos) {
		int_conv = sv->svlen();
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
		if (!sv->ins_seq.empty()) {
			bcf_update_info_int32(hdr, bcf_entry, "SVINSLEN", &int_conv, 1);
		}
	}
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", sv->source.c_str());
	if (!sv->ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", sv->ins_seq.c_str());
	}

	int int2_conv[2];
	int2_conv[0] = sv->rc_reads(), int2_conv[1] = sv->lc_reads();
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_READS", int2_conv, 2);
	int2_conv[0] = sv->rc_fwd_reads(), int2_conv[1] = sv->lc_fwd_reads();
	bcf_update_info_int32(hdr, bcf_entry, "FWD_SPLIT_READS", int2_conv, 2);
	int2_conv[0] = sv->rc_rev_reads(), int2_conv[1] = sv->lc_rev_reads();
	bcf_update_info_int32(hdr, bcf_entry, "REV_SPLIT_READS", int2_conv, 2);

	if (sv->full_junction_aln != NULL) {
		bcf_update_info_int32(hdr, bcf_entry, "FULL_JUNCTION_SCORE", &sv->full_junction_aln->best_score, 1);
		bcf_update_info_string(hdr, bcf_entry, "FULL_JUNCTION_CIGAR", sv->full_junction_aln->cigar.c_str());
	}
	int2_conv[0] = sv->left_anchor_aln->best_score, int2_conv[1] = sv->right_anchor_aln->best_score;
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE", int2_conv, 2);
	int2_conv[0] = sv->left_anchor_aln->next_best_score, int2_conv[1] = sv->right_anchor_aln->next_best_score;
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE2", int2_conv, 2);
	int2_conv[0] = sv->left_anchor_aln->seq_len, int2_conv[1] = sv->right_anchor_aln->seq_len;
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE", int2_conv, 2);
	char* split_junction_cigar = (char*) malloc(sv->left_anchor_aln->cigar.length() + sv->right_anchor_aln->cigar.length() + 2);
	sprintf(split_junction_cigar, "%s,%s", sv->left_anchor_aln->cigar.c_str(), sv->right_anchor_aln->cigar.c_str());
	bcf_update_info_string(hdr, bcf_entry, "SPLIT_JUNCTION_CIGAR", split_junction_cigar);
	std::string split_junction_mapping_range = sv->left_anchor_aln_string() + "," + sv->right_anchor_aln_string();
	bcf_update_info_string(hdr, bcf_entry, "SPLIT_JUNCTION_MAPPING_RANGE", split_junction_mapping_range.c_str());

	int max_mapq[] = {sv->rc_consensus ? (int) sv->rc_consensus->max_mapq : 0, sv->lc_consensus ? (int) sv->lc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	bcf_update_info_int32(hdr, bcf_entry, "OVERLAP", &sv->overlap, 1);
	bcf_update_info_float(hdr, bcf_entry, "MISMATCH_RATE", &sv->mismatch_rate, 1);

	if (sv->rc_consensus) {
		int ext_1sr_reads[] = { sv->rc_consensus->left_ext_reads, sv->rc_consensus->right_ext_reads };
		bcf_update_info_int32(hdr, bcf_entry, "RCC_EXT_1SR_READS", ext_1sr_reads, 2);
		int hq_ext_1sr_reads[] = { sv->rc_consensus->hq_left_ext_reads, sv->rc_consensus->hq_right_ext_reads };
		bcf_update_info_int32(hdr, bcf_entry, "RCC_HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
		if (sv->rc_consensus->remap_boundary != consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
			bcf_update_info_int32(hdr, bcf_entry, "REMAP_UB", &sv->rc_consensus->remap_boundary, 1);
		}
	}
	if (sv->lc_consensus) {
		int ext_1sr_reads[] = { sv->lc_consensus->left_ext_reads, sv->lc_consensus->right_ext_reads };
		bcf_update_info_int32(hdr, bcf_entry, "LCC_EXT_1SR_READS", ext_1sr_reads, 2);
		int hq_ext_1sr_reads[] = { sv->lc_consensus->hq_left_ext_reads, sv->lc_consensus->hq_right_ext_reads };
		bcf_update_info_int32(hdr, bcf_entry, "LCC_HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
		if (sv->lc_consensus->remap_boundary != consensus_t::LOWER_BOUNDARY_NON_CALCULATED) {
			bcf_update_info_int32(hdr, bcf_entry, "REMAP_LB", &sv->lc_consensus->remap_boundary, 1);
		}
	}
	
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &sv->disc_pairs, 1);
	if (sv->disc_pairs > 0) {
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_MAXMAPQ", &sv->disc_pairs_maxmapq, 1);
	}
	if (sv->source == "DP") {
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_HIGHMAPQ", &sv->disc_pairs_high_mapq, 1);
	}

	// add GT info
	int gt[1];
	gt[0] = bcf_gt_unphased(1);
	bcf_update_genotypes(hdr, bcf_entry, gt, 1);

	const char* ft_val = sv->is_pass() ? "PASS" : "FAIL";
	bcf_update_format_string(hdr, bcf_entry, "FT", &ft_val, 1);

	if (sv->svtype() == "INS") {
		if (sv->ins_seq.find("-") != std::string::npos) {
			bcf_update_info_flag(hdr, bcf_entry, "INCOMPLETE_ASSEMBLY", "", 1);
		}
	}
}

std::string get_sv_type(bcf_hdr_t* hdr, bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
        throw std::runtime_error("Failed to determine SVTYPE for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

int get_sv_end(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        delete[] data;
        return end-1; // return 0-based
    }

    bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
    if (size > 0) {
        int svlen = data[0];
        delete[] data;
        return sv->pos + abs(svlen);
    }

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_sv_info_str(bcf_hdr_t* hdr, bcf1_t* sv, std::string info) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, info.c_str(), &data, &len) < 0) {
        throw std::runtime_error("Failed to fetch " + info + " for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

std::string get_ins_seq(bcf_hdr_t* hdr, bcf1_t* sv) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	bcf_unpack(sv, BCF_UN_INFO);
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		return sv->d.allele[1];
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "SVINSSEQ", (void**) &data, &size);
	if (data) return data;

	return "";
}

sv_t* bcf_to_sv(bcf_hdr_t* hdr, bcf1_t* b) {

	int* data = NULL;
	int len = 0;
	bcf_get_info_int32(hdr, b, "MAX_MAPQ", &data, &len);
	int max_rc_mapq = 0, max_lc_mapq = 0;
	if (len > 0) {
		max_rc_mapq = data[0];
		max_lc_mapq = data[1];
	}

	int rc_fwd_reads = 0, rc_rev_reads = 0, lc_fwd_reads = 0, lc_rev_reads = 0;
	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "FWD_SPLIT_READS", &data, &len);
	if (len > 0) {
		rc_fwd_reads = data[0];
		lc_fwd_reads = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "REV_SPLIT_READS", &data, &len);
	if (len > 0) {
		rc_rev_reads = data[0];
		lc_rev_reads = data[1];
	}


	consensus_t* rc_consensus = NULL;
	if (rc_fwd_reads + rc_rev_reads > 0) {
		rc_consensus = new consensus_t(false, 0, 0, 0, "", rc_fwd_reads, rc_rev_reads, 0, max_rc_mapq, consensus_t::UPPER_BOUNDARY_NON_CALCULATED, 0);

		data = NULL;
		len = 0;
		bcf_get_info_int32(hdr, b, "RCC_EXT_1SR_READS", &data, &len);
		int rcc_left_ext_reads = 0, rcc_right_ext_reads = 0;
		if (len > 0) {
			rc_consensus->left_ext_reads = data[0];
			rc_consensus->right_ext_reads = data[1];
		}

		data = NULL;
		len = 0;
		bcf_get_info_int32(hdr, b, "RCC_HQ_EXT_1SR_READS", &data, &len);
		if (len > 0) {
			rc_consensus->hq_left_ext_reads = data[0];
			rc_consensus->hq_right_ext_reads = data[1];
		}

		data = NULL;
		len = 0;
		bcf_get_info_int32(hdr, b, "REMAP_UB", &data, &len);
		if (len > 0) rc_consensus->remap_boundary = data[0];
	} 

	consensus_t* lc_consensus = NULL;
	if (lc_fwd_reads + lc_rev_reads > 0) {
		lc_consensus = new consensus_t(true, 0, 0, 0, "", lc_fwd_reads, lc_rev_reads, 0, max_lc_mapq, consensus_t::LOWER_BOUNDARY_NON_CALCULATED, 0);

		data = NULL;
		len = 0;
		bcf_get_info_int32(hdr, b, "LCC_EXT_1SR_READS", &data, &len);
		if (len > 0) {
			lc_consensus->left_ext_reads = data[0];
			lc_consensus->right_ext_reads = data[1];
		}

		data = NULL;
		len = 0;
		bcf_get_info_int32(hdr, b, "LCC_HQ_EXT_1SR_READS", &data, &len);
		if (len > 0) {
			lc_consensus->hq_left_ext_reads = data[0];
			lc_consensus->hq_right_ext_reads = data[1];
		}

		data = NULL;
		len = 0;
		bcf_get_info_int32(hdr, b, "REMAP_LB", &data, &len);
		if (len > 0) lc_consensus->remap_boundary = data[0];
	}

	char* data2 = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "SPLIT_JUNCTION_MAPPING_RANGE", (void**) &data2, &len);
	std::string left_split_mapping_range = "", right_split_mapping_range = "";
	if (data2) {
		std::string mapping_range = data2;
		size_t comma_pos = mapping_range.find(",");
		left_split_mapping_range = mapping_range.substr(0, comma_pos);
		right_split_mapping_range = mapping_range.substr(comma_pos+1);
	}

	int left_split_mapping_start = std::stoi(left_split_mapping_range.substr(0, left_split_mapping_range.find("-")))-1;
	int left_split_mapping_end = std::stoi(left_split_mapping_range.substr(left_split_mapping_range.find("-")+1))-1;
	int right_split_mapping_start = std::stoi(right_split_mapping_range.substr(0, right_split_mapping_range.find("-")))-1;
	int right_split_mapping_end = std::stoi(right_split_mapping_range.substr(right_split_mapping_range.find("-")+1))-1;

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SCORE", &data, &len);
	int left_split_score = data[0], right_split_score = data[1];

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SCORE2", &data, &len);
	int left_split_score2 = data[0], right_split_score2 = data[1];

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SIZE", &data, &len);
	int left_split_size = data[0], right_split_size = data[1];

	data2 = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "SPLIT_JUNCTION_CIGAR", (void**) &data2, &len);
	std::string left_split_cigar = "", right_split_cigar = "";
	if (data2) {
		std::string cigar = data2;
		size_t comma_pos = cigar.find(",");
		left_split_cigar = cigar.substr(0, comma_pos);
		right_split_cigar = cigar.substr(comma_pos+1);
	}

	int full_junction_score = 0;
	len = 0;
	bcf_get_info_int32(hdr, b, "FULL_JUNCTION_SCORE", &data, &len);
	if (len > 0) full_junction_score = data[0];

	std::string full_junction_cigar = "";
	data2 = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "FULL_JUNCTION_CIGAR", (void**) &data2, &len);
	if (data2) full_junction_cigar = data2;

	sv_t::anchor_aln_t* left_anchor_aln = new sv_t::anchor_aln_t(left_split_mapping_start, left_split_mapping_end, left_split_size, left_split_score, left_split_score2, left_split_cigar);
	sv_t::anchor_aln_t* right_anchor_aln = new sv_t::anchor_aln_t(right_split_mapping_start, right_split_mapping_end, right_split_size, right_split_score, right_split_score2, right_split_cigar);
	sv_t::anchor_aln_t* full_junction_aln = full_junction_score > 0 ? new sv_t::anchor_aln_t(0, 0, 0, full_junction_score, 0, full_junction_cigar) : NULL;

	std::string svtype = get_sv_type(hdr, b);
	sv_t* sv;
	if (svtype == "DEL") {
		sv = new deletion_t(bcf_seqname_safe(hdr, b), b->pos, get_sv_end(hdr, b), get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, full_junction_aln);
	} else if (svtype == "DUP") {
		sv = new duplication_t(bcf_seqname_safe(hdr, b), b->pos, get_sv_end(hdr, b), get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, full_junction_aln);
	} else if (svtype == "INS") {
		sv = new insertion_t(bcf_seqname_safe(hdr, b), b->pos, get_sv_end(hdr, b), get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, full_junction_aln);
	} else {
		throw std::runtime_error("Unsupported SV type: " + svtype);
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "DISC_PAIRS", &data, &len);
	if (len > 0) sv->disc_pairs = data[0];

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "DISC_PAIRS_MAXMAPQ", &data, &len);
	if (len > 0) sv->disc_pairs_maxmapq = data[0];

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "DISC_PAIRS_HIGHMAPQ", &data, &len);
	if (len > 0) sv->disc_pairs_high_mapq = data[0];

	sv->id = b->d.id;
	sv->source = get_sv_info_str(hdr, b, "SOURCE");

	for (int i = 0; i < b->d.n_flt; i++) {
		sv->filters.push_back(bcf_hdr_int2id(hdr, BCF_DT_ID, b->d.flt[i]));
	}

	return sv;
}

#endif /* VCF_UTILS_H */
