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
bcf_hdr_t* generate_vcf_header(chr_seqs_map_t& contigs, std::string sample_name, config_t config, std::string command) {
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
	sprintf(size_flt_tag, "##FILTER=<ID=SMALL,Description=\"Insertion smaller than %d bp.\">", config.min_sv_size);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	// TODO: remove this filter
	const char* alt_short_flt_tag = "##FILTER=<ID=ALT_SHORTER_THAN_REF,Description=\"If this insertion/replacement was real, alternative"
			"allele would be shorter than reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, alt_short_flt_tag, &len));

	const char* anomalous_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_DEPTH,Description=\"The insertion region has anomalous depth.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_depth_flt_tag, &len));

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

	const char* overlap_tag = "##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Overlap (in bp) between the left and right contigs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, overlap_tag, &len));

	const char* mismatch_rate_tag = "##INFO=<ID=MISMATCH_RATE,Number=1,Type=Float,Description=\"Mismatch rate of overlap between the left and right contigs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mismatch_rate_tag, &len));

	const char* sr_tag = "##INFO=<ID=SPLIT_READS,Number=2,Type=Integer,Description=\"Split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sr_tag, &len));

	const char* fwd_sr_tag = "##INFO=<ID=FWD_SPLIT_READS,Number=2,Type=Integer,Description=\"Forward split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fwd_sr_tag, &len));

	const char* rev_sr_tag = "##INFO=<ID=REV_SPLIT_READS,Number=2,Type=Integer,Description=\"Reverse split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rev_sr_tag, &len));

	const char* sd_tag = "##INFO=<ID=STABLE_DEPTHS,Number=2,Type=Integer,Description=\"Depths of the stable regions (in practice, the regions left and right of the insertion site).\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sd_tag, &len));

	const char* consensuses_id_tag = "##INFO=<ID=CONSENSUSES_ID,Number=1,Type=String,Description=\"ID of the consensuses used to generate this call (for development/debug).\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, consensuses_id_tag, &len));

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

	const char* algo_tag = "##INFO=<ID=ALGORITHM,Number=1,Type=String,Description=\"Algorithm used to report the call.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, algo_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, gt_tag, &len));

	// add ALT
	const char* ins_alt_tag = "##ALT=<ID=INS,Description=\"Insertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ins_alt_tag, &len));

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
	called_by_ss << "max-trans-size: " << config.max_trans_size << "; ";
	called_by_ss << "min-stable-mapq: " << config.min_stable_mapq << "; ";
	called_by_ss << "min-clip-len: " << config.min_clip_len << "; ";
	called_by_ss << "max-seq-error: " << config.max_seq_error << "; ";
	called_by_ss << "sampling-regions: " << (config.sampling_regions.empty() ? "no" : config.sampling_regions) << "; ";
	called_by_ss << "per-contig-stats: " << (config.per_contig_stats ? "true" : "false") << "; ";
	std::string called_by = called_by_ss.str();
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, called_by.c_str(), &len));

	// add samples
	bcf_hdr_add_sample(header, sample_name.c_str());

	return header;
}

void sv2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, sv_t* sv, char* chr_seq, std::vector<std::string>& filters) {
	bcf_clear(bcf_entry);
	
	bcf_entry->rid = bcf_hdr_name2id(hdr, sv->chr.c_str());
	bcf_entry->pos = sv->start;
	bcf_update_id(hdr, bcf_entry, sv->id.c_str());
	
	std::string alleles = std::string(1, chr_seq[sv->start]) + ",<" + sv->svtype() + ">";
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	// add FILTERs
	for (std::string& filter : filters) {
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
	int2_conv[0] = sv->rc_fwd_reads, int2_conv[1] = sv->lc_fwd_reads;
	bcf_update_info_int32(hdr, bcf_entry, "FWD_SPLIT_READS", int2_conv, 2);
	int2_conv[0] = sv->rc_rev_reads, int2_conv[1] = sv->lc_rev_reads;
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
	}
	if (sv->lc_consensus) {
		int ext_1sr_reads[] = { sv->lc_consensus->left_ext_reads, sv->lc_consensus->right_ext_reads };
		bcf_update_info_int32(hdr, bcf_entry, "LCC_EXT_1SR_READS", ext_1sr_reads, 2);
		int hq_ext_1sr_reads[] = { sv->lc_consensus->hq_left_ext_reads, sv->lc_consensus->hq_right_ext_reads };
		bcf_update_info_int32(hdr, bcf_entry, "LCC_HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
	}

	// add GT info
	int gt[1];
	gt[0] = bcf_gt_unphased(1);
	bcf_update_genotypes(hdr, bcf_entry, gt, 1);

	if (!filters.empty()) {
		const char* ft_val = (filters[0] == "PASS") ? "PASS" : "FAIL";
		bcf_update_format_string(hdr, bcf_entry, "FT", &ft_val, 1);
	}
}

void del2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, deletion_t* del, char* chr_seq, std::vector<std::string>& filters) {
	sv2bcf(hdr, bcf_entry, del, chr_seq, filters);
}

void dup2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, duplication_t* dup, char* chr_seq, std::vector<std::string>& filters) {
	sv2bcf(hdr, bcf_entry, dup, chr_seq, filters);
}

void ins2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, insertion_t* ins, char* chr_seq, std::vector<std::string>& filters) {
	sv2bcf(hdr, bcf_entry, ins, chr_seq, filters);

	if (ins->ins_seq.find("-") != std::string::npos) {
		bcf_update_info_flag(hdr, bcf_entry, "INCOMPLETE_ASSEMBLY", "", 1);
	}
}

#endif /* VCF_UTILS_H */
