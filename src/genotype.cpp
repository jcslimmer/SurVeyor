#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <chrono>

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"

#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
#include "vcf_utils.h"

chr_seqs_map_t chr_seqs;
config_t config;
stats_t stats;
std::mutex mtx;

std::string bam_fname, reference_fname;
bam_pool_t* bam_pool;

StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

void add_tags(bcf_hdr_t* hdr) {
    int len = 0;

    bcf_hdr_remove(hdr, BCF_HL_INFO, "END");
    const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, end_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "SVLEN");
    const char* svlen_tag = "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, svlen_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "SVINSLEN");
    const char* svinslen_tag = "##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, svinslen_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "LEFT_ANCHOR_BASE_COUNT");
    const char* labc_tag = "##INFO=<ID=LEFT_ANCHOR_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the left anchor region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,labc_tag, &len));
	
    bcf_hdr_remove(hdr, BCF_HL_INFO, "RIGHT_ANCHOR_BASE_COUNT");
	const char* rabc_tag = "##INFO=<ID=RIGHT_ANCHOR_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the right anchor region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rabc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "SV_REF_PREFIX_BASE_COUNT");
	const char* svrefpbc_tag = "##INFO=<ID=SV_REF_PREFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 5000 bp of the reference region of the SV (if the SV affects less than 5000 bp, the whole SV is considered).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,svrefpbc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "SV_REF_SUFFIX_BASE_COUNT");
	const char* svrefsbc_tag = "##INFO=<ID=SV_REF_SUFFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the last 5000 bp of the reference region of the SV (if the SV affects less than 5000 bp, the whole SV is considered).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,svrefsbc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "INS_PREFIX_BASE_COUNT");
	const char* pbc_tag = "##INFO=<ID=INS_PREFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the prefix of the inserted sequence (for incomplete assemblies) or in the full inserted sequence.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,pbc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_INFO, "INS_SUFFIX_BASE_COUNT");
	const char* sbc_tag = "##INFO=<ID=INS_SUFFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the suffix of the inserted sequence (for incomplete assemblies) or in the full inserted sequence. For insertion not marked with INCOMPLETE_ASSEMBLY, this will be identical to PREFIX_BASE_COUNT.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,sbc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AR");
    const char* ar_tag = "##FORMAT=<ID=AR,Number=1,Type=Integer,Description=\"Number of reads supporting the alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,ar_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "RR");
    const char* rr_tag = "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Number of reads supporting the reference allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rr_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ER");
    const char* er_tag = "##FORMAT=<ID=ER,Number=1,Type=Integer,Description=\"Number of reads supporting equally well reference and alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,er_tag, &len));
}

void update_record(bcf_hdr_t* in_hdr, bcf_hdr_t* out_hdr, sv_t* sv, char* chr_seq) {
    
    bcf_translate(out_hdr, in_hdr, sv->vcf_entry);

    // update INFO fields
    bcf_update_info_int32(out_hdr, sv->vcf_entry, "END", &sv->end, 1);

    int svlen = sv->svlen();
    bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVLEN", &svlen, 1);

    if (!sv->ins_seq.empty()) {
        int svinslen = sv->ins_seq.length();
        bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVINSLEN", &svinslen, 1);
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

    // update FORMAT fields
    int ar = sv->regenotyping_info.alt_better_reads;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR", &ar, 1);

    int rr = sv->regenotyping_info.ref_better_reads;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR", &rr, 1);

    int er = sv->regenotyping_info.alt_ref_equal_reads;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ER", &er, 1);
}

void genotype_del(deletion_t* del) {
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

    open_samFile_t* bam_file = bam_pool->get_bam_reader();
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int ref_better = 0, alt_better = 0, same = 0;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        std::string seq = get_sequence(read);
        
        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

        // align to REF (two breakpoints)
        aligner.Align(seq.c_str(), ref_bp1_seq, ref_bp1_len, filter, &ref1_aln, 0);
        aligner.Align(seq.c_str(), ref_bp2_seq, ref_bp2_len, filter, &ref2_aln, 0);

        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            alt_better++;
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            ref_better++;
        } else {
            same++;
        }
    }

    del->regenotyping_info.alt_better_reads = alt_better;
    del->regenotyping_info.ref_better_reads = ref_better;
    del->regenotyping_info.alt_ref_equal_reads = same;

    bam_pool->release_bam_reader(bam_file);
    delete[] alt_seq;
    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_dels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<deletion_t*> dels,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {
                    
    for (deletion_t* del : dels) {
        genotype_del(del);
    }
}

void genotype_small_dup(duplication_t* dup) {
    std::stringstream log_ss;

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
			strncpy(alt_seq+pos, chr_seqs.get_seq(dup->chr)+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
		}
		strncpy(alt_seq+pos, chr_seqs.get_seq(dup->chr)+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
		alt_seqs.push_back(alt_seq);
	}

    std::stringstream region;
    region << dup->chr << ":" << ref_start << "-" << ref_end;

    open_samFile_t* bam_file = bam_pool->get_bam_reader();
	hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region.str().c_str());
	
    bam1_t* read = bam_init1();

    int ref_better = 0, alt_better = 0, same = 0;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_best_aln, alt_aln, ref_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

        std::string seq = get_sequence(read);
        aligner.Align(seq.c_str(), ref_seq, ref_len, filter, &ref_aln, 0);

        for (int i = 0; i < alt_seqs.size(); i++) {
            aligner.Align(seq.c_str(), alt_seqs[i], strlen(alt_seqs[i]), filter, &alt_aln, 0);
            if (alt_aln.sw_score > alt_best_aln.sw_score) {
                alt_best_aln = alt_aln;
            }
        }

        if (alt_best_aln.sw_score > ref_aln.sw_score) {
            alt_better++;
        } else if (alt_best_aln.sw_score < ref_aln.sw_score) {
            ref_better++;
        } else {
            same++;
        }
    }

    dup->regenotyping_info.alt_better_reads = alt_better;
    dup->regenotyping_info.ref_better_reads = ref_better;
    dup->regenotyping_info.alt_ref_equal_reads = same;

    bam_pool->release_bam_reader(bam_file);
    delete[] ref_seq;
    for (char* alt_seq : alt_seqs) {
        delete[] alt_seq;
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_large_dup(duplication_t* dup) {
    
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

    open_samFile_t* bam_file = bam_pool->get_bam_reader();
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int alt_better = 0, same = 0;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        std::string seq = get_sequence(read);
        
        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

        // align to REF (two breakpoints)
        aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln, 0);
        aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln, 0);

        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            alt_better++;
        } else {
            same++;
        }
    }

    dup->regenotyping_info.alt_better_reads = alt_better;
    dup->regenotyping_info.ref_better_reads = 0;
    dup->regenotyping_info.alt_ref_equal_reads = same;

    bam_pool->release_bam_reader(bam_file);
    delete[] alt_seq;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_dups(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<duplication_t*> dups,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {
                    
    for (duplication_t* dup : dups) {
        if (dup->svlen() <= stats.read_len-2*config.min_clip_len) {
			genotype_small_dup(dup);
		} else {
			genotype_large_dup(dup);
		}
    }
}

void genotype_ins(insertion_t* ins) {
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
	int alt_len = alt_lf_len + ins->ins_seq.length() + alt_rf_len;
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

    open_samFile_t* bam_file = bam_pool->get_bam_reader();
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int ref_better = 0, alt_better = 0, same = 0;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

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
            alt_better++;
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            ref_better++;
        } else {
            same++;
        }
    }

    ins->regenotyping_info.alt_better_reads = alt_better;
    ins->regenotyping_info.ref_better_reads = ref_better;
    ins->regenotyping_info.alt_ref_equal_reads = same;

    bam_pool->release_bam_reader(bam_file);
    delete[] alt_bp1_seq;
    delete[] alt_bp2_seq;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_inss(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<insertion_t*> inss,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {
                    
    for (insertion_t* ins : inss) { 
        genotype_ins(ins);
    }
}

int main(int argc, char* argv[]) {

    std::string in_vcf_fname = argv[1];
    std::string out_vcf_fname = argv[2];
    bam_fname = argv[3];
    reference_fname = argv[4];
    std::string workdir = argv[5];
    std::string sample_name = argv[6];

    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    open_samFile_t* bam_file = open_samFile(bam_fname);
	if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + reference_fname;
	}

    chr_seqs.read_fasta_into_map(reference_fname);
    bam_pool = new bam_pool_t(bam_fname, reference_fname);

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
    std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

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
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        sv_t* sv = bcf_to_sv(in_vcf_header, vcf_record);
        sv->vcf_entry = bcf_dup(vcf_record);
        if (sv->svtype() == "DEL") {
            dels_by_chr[sv->chr].push_back((deletion_t*) sv);
        } else if (sv->svtype() == "DUP") {
            dups_by_chr[sv->chr].push_back((duplication_t*) sv);
        } else if (sv->svtype() == "INS") {
        	inss_by_chr[sv->chr].push_back((insertion_t*) sv);
        }
    }

    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    bcf_hdr_t* out_vcf_header = bcf_hdr_dup(in_vcf_header);
    add_tags(out_vcf_header);
    if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
    	throw std::runtime_error("Failed to read the VCF header.");
    }

    // genotype chrs in descending order of svs
    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (auto& p : dels_by_chr) {
    	std::string contig_name = p.first;
        std::vector<deletion_t*>& dels = p.second;
        std::future<void> future = thread_pool.push(genotype_dels, contig_name, chr_seqs.get_seq(contig_name),
        		chr_seqs.get_len(contig_name), dels_by_chr[contig_name], in_vcf_header, out_vcf_header, stats, config);
        futures.push_back(std::move(future));
    }
    for (auto& p : dups_by_chr) {
    	std::string contig_name = p.first;
        std::vector<duplication_t*>& dups = p.second;
        std::future<void> future = thread_pool.push(genotype_dups, contig_name, chr_seqs.get_seq(contig_name),
        		chr_seqs.get_len(contig_name), dups_by_chr[contig_name], in_vcf_header, out_vcf_header, stats, config);
        futures.push_back(std::move(future));
    }

    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
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
    	std::sort(contig_svs.begin(), contig_svs.end(), [](const sv_t* sv1, const sv_t* sv2) {return sv1->start < sv2->start;});

		for (auto& sv : contig_svs) {
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AC", NULL, 0);
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AN", NULL, 0);
            update_record(in_vcf_header, out_vcf_header, sv, chr_seqs.get_seq(contig_name));
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

    // tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);

}
