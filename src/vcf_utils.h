#ifndef VCF_UTILS_H
#define VCF_UTILS_H

#include <cstddef>
#include <iostream>
#include <chrono>
#include <ctime>
#include <htslib/vcf.h>
#include <sstream>
#include "htslib/hts.h"
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

// Define a struct to hold all properties for each format tag
typedef struct {
    const char* suffix;       // Tag suffix (e.g., "C", "CF", etc.)
	int n;					  // Number of values	
    const char* type;         // Data type ("Integer" or "Float")
    const char* desc_format;  // Description format string
} format_tag_def_t;

void add_headers(bcf_hdr_t* hdr, const char prefix, int pos, const char* allele_type, const format_tag_def_t formats[], int n_formats) {
	char tag_id[20];
    char tag_str[1024];
    char desc_str[200];
    int len;

    for (size_t i = 0; i < n_formats; i++) {
        const format_tag_def_t* fmt = &formats[i];

        snprintf(tag_id, sizeof(tag_id), fmt->suffix, prefix, pos, fmt->suffix);
        snprintf(desc_str, sizeof(desc_str), fmt->desc_format, pos, allele_type);
        snprintf(tag_str, sizeof(tag_str), 
                "##FORMAT=<ID=%s,Number=%d,Type=%s,Description=\"%s\">",
                tag_id, fmt->n, fmt->type, desc_str);

        bcf_hdr_remove(hdr, BCF_HL_FMT, tag_id);
        bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, tag_str, &len));
    }
}

void add_read_support_headers(bcf_hdr_t* hdr, const char prefix, int pos, const char* allele_type) {
    // Array of tag definitions
    const format_tag_def_t formats[] = {
        {
            "%cR%d", 1, "Integer",
            "Number of reads supporting breakpoint %d in the %s allele."
        },
        {
            "%cR%dC", 1, "Integer", 
            "Number of consistent reads supporting breakpoint %d in the %s allele."
        },
        {
            "%cR%dCF", 1, "Integer",
            "Number of consistent forward reads supporting breakpoint %d in the %s allele."
        },
        {
            "%cR%dCR", 1, "Integer",
            "Number of consistent reverse reads supporting breakpoint %d in the %s allele."
        },
        {
            "%cR%dCAS", 1, "Float",
            "Average aln score of consistent reads supporting breakpoint %d of the SV to the %s allele consensus."
        },
		{
			"%cR%dCSS", 1, "Float",
			"Standard deviation of aln score of consistent reads supporting breakpoint %d of the SV to the %s allele consensus."
		},
        {
            "%cR%dCHQ", 1, "Integer",
            "Number of high-quality consistent reads supporting breakpoint %d in the %s allele."
        },
		{
            "%cR%dCmQ", 1, "Integer",
            "Minimum mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
        },
        {
            "%cR%dCMQ", 1, "Integer",
            "Maximum mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
        },
		{
            "%cR%dCAQ", 1, "Float",
            "Average mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
        },
		{
			"%cR%dCSQ", 1, "Float",
			"Standard deviation of mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
		}
    };

    add_headers(hdr, prefix, pos, allele_type, formats, sizeof(formats)/sizeof(formats[0]));
}

void add_pairs_support_headers(bcf_hdr_t* hdr, const char prefix, int pos, const char* allele_type) {
	const format_tag_def_t formats[] = {
		{
			"%cSP%d", 1, "Integer",
			"Number of pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dHQ", 2, "Integer",
			"Number of high-quality positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dmQ", 2, "Integer",
			"Minimum mapping quality of positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dMQ", 2, "Integer",
			"Maximum mapping quality of positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dAQ", 2, "Float",
			"Average mapping quality of positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dSQ", 2, "Float",
			"Standard deviation of mapping quality of positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dSPAN", 2, "Integer",
			"Number of base pairs to the left and right, respectively, of the breakpoint %d in the %s allele covered by reads within its supporting pairs."
		},
		{
			"%cSP%dNMA", 2, "Float",
			"Average NM (non-matches) of positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		},
		{
			"%cSP%dNMS", 2, "Float",
			"Standard deviation of NM (non-matches) of positive and negative reads, respectively, within pairs supporting breakpoint %d in the %s allele."
		}
	};

	add_headers(hdr, prefix, pos, allele_type, formats, sizeof(formats)/sizeof(formats[0]));
}

void add_stray_pairs_headers(bcf_hdr_t* hdr, int bp_number) {
	const format_tag_def_t formats[] = {
		{
			"%cSP%d", 1, "Integer",
			"Number of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dHQ", 1, "Integer",
			"Number of high-quality pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dmQ", 1, "Integer",
			"Minimum mapping quality of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dMQ", 1, "Integer",
			"Maximum mapping quality of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dAQ", 1, "Float",
			"Average mapping quality of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dSQ", 1, "Float",
			"Standard deviation of mapping quality of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dNMA", 1, "Float",
			"Average NM (non-matches) of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dNMS", 1, "Float",
			"Standard deviation of NM (non-matches) of pairs around the breakpoint %d in the %s allele that are discordant and yet do not support the SV."
		}
	};

	add_headers(hdr, 'S', bp_number, "reference", formats, sizeof(formats)/sizeof(formats[0]));
}

void add_fmt_tags(bcf_hdr_t* hdr) {
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

	bcf_hdr_remove(hdr, BCF_HL_INFO, "MH_LEN");
	const char* mh_len = "##INFO=<ID=MH_LEN,Number=1,Type=Integer,Description=\"Length of the microhomology.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,mh_len, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "GT");
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, gt_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "FT");
	const char* ft_tag = "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS indicates a reliable call. Any other value means the call is not reliable.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ft_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "MD");
	const char* md_tag = "##FORMAT=<ID=MD,Number=4,Type=Integer,Description=\"Median depth in the left flanking, prefix, suffix and right flanking regions of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, md_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "MDHQ");
	const char* mdhq_tag = "##FORMAT=<ID=MDHQ,Number=4,Type=Integer,Description=\"Median depth in the left flanking, prefix, suffix and right flanking regions of the SV. Only high MAPQ reads are considered.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdhq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "CLMD");
	const char* clmd_tag = "##FORMAT=<ID=CLMD,Number=2,Type=Integer,Description=\"Median depth in the left and right flanking regions, considering only the regions where the SV evidence (DP and SR) lies.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, clmd_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "CLMDHQ");
	const char* clmdhq_tag = "##FORMAT=<ID=CLMDHQ,Number=2,Type=Integer,Description=\"Median depth in the left and right flanking regions, considering only the regions where the SV evidence (DP and SR) lies. Only high MAPQ reads are considered.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, clmdhq_tag, &len));
	
	add_read_support_headers(hdr, 'A', 1, "alternate");
	add_pairs_support_headers(hdr, 'A', 1, "alternate");
	
	add_read_support_headers(hdr, 'A', 2, "alternate");
	add_pairs_support_headers(hdr, 'A', 2, "alternate");

	add_read_support_headers(hdr, 'R', 1, "reference");
	add_pairs_support_headers(hdr, 'R', 1, "reference");

	add_read_support_headers(hdr, 'R', 2, "reference");
	add_pairs_support_headers(hdr, 'R', 2, "reference");

	add_stray_pairs_headers(hdr, 1);
	add_stray_pairs_headers(hdr, 2);

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ER");
    const char* er_tag = "##FORMAT=<ID=ER,Number=1,Type=Integer,Description=\"Number of reads supporting equally well reference and alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,er_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "TD");
    const char* nc_tag = "##FORMAT=<ID=TD,Number=1,Type=Integer,Description=\"The variant region is too deep to be genotyped reliably.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,nc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MINSIZE");
    const char* minsize_tag = "##FORMAT=<ID=MINSIZE,Number=1,Type=Integer,Description=\"Minimum size of the event calculated based on insert size distribution."
            "Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, minsize_tag, &len));
    
    bcf_hdr_remove(hdr, BCF_HL_FMT, "MAXSIZE");
    const char* maxsize_tag = "##FORMAT=<ID=MAXSIZE,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, maxsize_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "KSPVAL");
    const char* kspval_tag = "##FORMAT=<ID=KSPVAL,Number=1,Type=Float,Description=\"p-value of the KS test.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, kspval_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AXR");
    const char* axr_tag = "##FORMAT=<ID=AXR,Number=2,Type=Integer,Description=\"Number of reads used to extend the alternative allele consensus to the left and the right, respectively.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, axr_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AXRHQ");
    const char* axrhq_tag = "##FORMAT=<ID=AXRHQ,Number=2,Type=Integer,Description=\"Number of high-quality reads used to extend the alternative allele consensus to the left and the right, respectively.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, axrhq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "EXL");
    const char* exl_tag = "##FORMAT=<ID=EXL,Number=1,Type=Integer,Description=\"Length of the extended alternative allele consensus.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exl_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXL2");
	const char* exl2_tag = "##FORMAT=<ID=EXL2,Number=1,Type=Integer,Description=\"Length of the extended alternative allele consensus for the second breakpoint (only for insertions).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exl2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "EXAS");
    const char* exas_tag = "##FORMAT=<ID=EXAS,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus and the original alternative allele (the reference with the SV applied).\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exas_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXAS2");
	const char* exas2_tag = "##FORMAT=<ID=EXAS2,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus for the second breakpoint and the original alternative allele (the reference with the SV applied).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exas2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "EXRS");
    const char* exrs_tag = "##FORMAT=<ID=EXRS,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus and the reference.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exrs_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXRS2");
	const char* exrs2_tag = "##FORMAT=<ID=EXRS2,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus for the second breakpoint and the reference.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exrs2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSS");
	const char* exss_tag = "##FORMAT=<ID=EXSS,Number=2,Type=Integer,Description=\"How much of the extended alternative allele aligns to the left and right, respectively, of the breakpoint.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exss_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSS2");
	const char* exss2_tag = "##FORMAT=<ID=EXSS2,Number=2,Type=Integer,Description=\"How much of the extended alternative allele aligns to the left and right, respectively, of the second breakpoint (only for insertions).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exss2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSC");
	const char* exssc_tag = "##FORMAT=<ID=EXSSC,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exssc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSCIA");
	const char* exsscia_tag = "##FORMAT=<ID=EXSSCIA,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus for the first breakpoint, when allowed to map to the reference independently.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exsscia_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSC2");
	const char* exssc2_tag = "##FORMAT=<ID=EXSSC2,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus for the second breakpoint (only for insertions).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exssc2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSC2IA");
	const char* exssc2ia_tag = "##FORMAT=<ID=EXSSC2IA,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus for the second breakpoint, when allowed to map to the reference independently.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exssc2ia_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EPR");
	const char* epr_tag = "##FORMAT=<ID=EPR,Number=1,Type=Float,Description=\"Probability of the SV existing in the sample, according to the ML model.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, epr_tag, &len));
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
	const char* remap_boundary_flt_tag = "##FILTER=<ID=REMAP_BOUNDARY_FILTER,Description=\"One of the breakpoints is incompatible with mate locations.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_boundary_flt_tag, &len));

	const char* failed_to_ext_flt_tag = "##FILTER=<ID=FAILED_TO_EXTEND,Description=\"No reads can extend the consensus.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, failed_to_ext_flt_tag, &len));

	const char* anomalous_fl_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_FLANKING_DEPTH,Description=\"The insertion region has anomalous depth.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_fl_depth_flt_tag, &len));

	const char* anomalous_sc_flt_tag = "##FILTER=<ID=ANOMALOUS_SC_NUMBER,Description=\"The number of soft-clipped reads supporting this call is too large.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_sc_flt_tag, &len));

	const char* not_enough_disc_pairs_flt_tag = "##FILTER=<ID=NOT_ENOUGH_DISC_PAIRS,Description=\"Not enough discordant pairs supporting the event.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, not_enough_disc_pairs_flt_tag, &len));

	const char* low_mapq_disc_pairs_flt_tag = "##FILTER=<ID=LOW_MAPQ_DISC_PAIRS,Description=\"Discordant pairs supporting the SV have low MAPQ.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_mapq_disc_pairs_flt_tag, &len));

	const char* weak_support_flt_tag = "##FILTER=<ID=WEAK_SUPPORT,Description=\"Remapped breakpoint has low support from local reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_support_flt_tag, &len));

	const char* low_ptn_ratio_flt_tag = "##FILTER=<ID=LOW_PTN_RATIO,Description=\"Low positive-to-negative ratio.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_ptn_ratio_flt_tag, &len));

	const char* short_anchor_flt_tag = "##FILTER=<ID=SHORT_ANCHOR,Description=\"Read pairs supporting the SV are limited to a very small portion of the genome.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, short_anchor_flt_tag, &len));

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

	const char* infsvinsseq_tag = "##INFO=<ID=SVINSSEQ_INFERRED,Number=1,Type=String,Description=\"Inferred insertion sequence. When the inserted sequence is too long "
		"to be fully assembled but SurVeyor suspects it to be a transposition, it uses the reference to infer the content of the insertion. Not guaranteed to be accurate. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, infsvinsseq_tag, &len));

	const char* source_tag = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, source_tag, &len));

	const char* remap_lb_tag = "##INFO=<ID=REMAP_LB,Number=1,Type=Integer,Description=\"Minimum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_lb_tag, &len));

	const char* remap_ub_tag = "##INFO=<ID=REMAP_UB,Number=1,Type=Integer,Description=\"Maximum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_ub_tag, &len));

	const char* est_size_tag = "##INFO=<ID=EST_SIZE,Number=1,Type=Integer,Description=\"Estimated size of the imprecise event. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, est_size_tag, &len));

	const char* og_range_tag = "##INFO=<ID=ORIGINAL_RANGE,Number=1,Type=String,Description=\"Unadjusted imprecise range predicted by discordant pairs. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, og_range_tag, &len));

	const char* imprecise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The reported boundaries are not precise.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, imprecise_tag, &len));

	const char* ia_tag = "##INFO=<ID=INCOMPLETE_ASSEMBLY,Number=0,Type=Flag,Description=\"The inserted sequence is too long, and it could not be fully assembled using short reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ia_tag, &len));

	// add ALT
	const char* del_alt_tag = "##ALT=<ID=DEL,Description=\"Deletion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, del_alt_tag, &len));

	const char* dup_alt_tag = "##ALT=<ID=DUP,Description=\"Tandem Duplication\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dup_alt_tag, &len));

	const char* ins_alt_tag = "##ALT=<ID=INS,Description=\"Insertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ins_alt_tag, &len));

	const char* inv_alt_tag = "##ALT=<ID=INV,Description=\"Invertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, inv_alt_tag, &len));

	add_fmt_tags(header);

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

	if (!sample_name.empty()) {
		bcf_hdr_add_sample(header, sample_name.c_str());
	}
	
	return header;
}

void sv2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, sv_t* sv, char* chr_seq, bool for_gt = false) {
	bcf_clear(bcf_entry);
	
	bcf_entry->rid = bcf_hdr_name2id(hdr, sv->chr.c_str());
	bcf_entry->pos = sv->start;
	bcf_update_id(hdr, bcf_entry, sv->id.c_str());
	
	std::string ref_allele = std::string(1, chr_seq[sv->start]);
	std::string alleles = ref_allele + ",";
	if (sv->svtype() == "BND") {
		breakend_t* bnd = (breakend_t*) sv;
		if (bnd->direction == '+') {
			alleles += std::string("]") + bcf_seqname(hdr, bcf_entry) + ":" + std::to_string(bnd->end+1) + "]" + ref_allele; 
		} else {
			alleles += std::string("[") + bcf_seqname(hdr, bcf_entry) + ":" + std::to_string(bnd->end+1) + "[" + ref_allele;
		}
	} else {
		alleles += "<" + sv->svtype() + ">";
	}
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	int n_samples = bcf_hdr_nsamples(hdr);

	// add FILTERs
	if (n_samples) {
		std::string ft;
		for (std::string& filter : sv->sample_info.filters) {
			if (!ft.empty()) ft += ";";
			ft += filter;
		}
		if (ft.empty()) {
			ft = "PASS";
		}
		const char* ft_cstr[1];
		ft_cstr[0] = (char*) ft.c_str();
		bcf_update_format_string(hdr, bcf_entry, "FT", ft_cstr, 1);
	}

	// add INFO
	int int_conv = sv->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", sv->svtype().c_str());
	if (sv->ins_seq.find("-") == std::string::npos) {
		int_conv = sv->svlen();
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
		if (!sv->ins_seq.empty()) {
			int_conv = sv->ins_seq.length();
			bcf_update_info_int32(hdr, bcf_entry, "SVINSLEN", &int_conv, 1);
		}
	}
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", sv->source.c_str());
	if (!sv->ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", sv->ins_seq.c_str());
	}
	if (!sv->inferred_ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ_INFERRED", sv->inferred_ins_seq.c_str());
	}

	if (!for_gt) {
		if (sv->rc_consensus && sv->rc_consensus->remap_boundary != consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
			bcf_update_info_int32(hdr, bcf_entry, "REMAP_UB", &sv->rc_consensus->remap_boundary, 1);
		}
		if (sv->lc_consensus && sv->lc_consensus->remap_boundary != consensus_t::LOWER_BOUNDARY_NON_CALCULATED) {
			bcf_update_info_int32(hdr, bcf_entry, "REMAP_LB", &sv->lc_consensus->remap_boundary, 1);
		}
	}

	bcf_update_info_flag(hdr, bcf_entry, "IMPRECISE", "", sv->imprecise);
	bcf_update_genotypes(hdr, bcf_entry, sv->sample_info.gt, sv->n_gt);

	if (sv->incomplete_ins_seq()) {
		bcf_update_info_flag(hdr, bcf_entry, "INCOMPLETE_ASSEMBLY", "", 1);
	}

	if (sv->svtype() == "DEL") {
		if (!for_gt) {
			deletion_t* del = (deletion_t*) sv;
			if (!del->original_range.empty()) {
				bcf_update_info_string(hdr, bcf_entry, "ORIGINAL_RANGE", del->original_range.c_str());
			}
			if (del->imprecise && del->estimated_size != deletion_t::SIZE_NOT_COMPUTED) {
				bcf_update_info_int32(hdr, bcf_entry, "EST_SIZE", &(del->estimated_size), 1);
			}
		}
	}
}

std::string get_sv_type(bcf_hdr_t* hdr, bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
		// try to determine svtype from ALT
		bcf_unpack(sv, BCF_UN_STR);
		if (sv->d.allele[1][0] == '<') {
			std::string alt = sv->d.allele[1];
			return alt.substr(1, alt.length()-2);
		}
		std::cerr << "Failed to determine SVTYPE for sv " << std::string(sv->d.id) << std::endl;
		return "";
    }
    std::string svtype = data;
    free(data);
    return svtype;
}

int get_sv_len(bcf_hdr_t* hdr, bcf1_t* sv) {
	int* data = NULL;
	int size = 0;
	bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
	if (size > 0) {
		int len = data[0];
		free(data);
		return len;
	}

	if (sv->d.allele[1][0] != '<') {
		return strlen(sv->d.allele[1])-strlen(sv->d.allele[0]);
	}

	return bcf_int32_missing;
}

int get_sv_end(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        free(data);
        return end-1; // return 0-based
    }

	if (get_sv_type(hdr, sv) == "INS") return sv->pos;

	int svlen = get_sv_len(hdr, sv);
	if (svlen != bcf_int32_missing) {
		return sv->pos + abs(svlen);
	}

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_sv_info_str(bcf_hdr_t* hdr, bcf1_t* sv, std::string info) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, info.c_str(), &data, &len) < 0) {
		return "";
    }
    std::string svtype = data;
    free(data);
    return svtype;
}

std::string get_ins_seq(bcf_hdr_t* hdr, bcf1_t* sv) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	bcf_unpack(sv, BCF_UN_INFO);
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		// return everything after the padding base
		return sv->d.allele[1] + 1;
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "SVINSSEQ", (void**) &data, &size);
	if (data) {
		std::string ins_seq = data;
		free(data);
		return ins_seq;
	}

	// LEFT_SVINSSEQ + "-" + RIGHT_SVINSSEQ, if they exist
	std::string left_ins_seq = get_sv_info_str(hdr, sv, "LEFT_SVINSSEQ");
	std::string right_ins_seq = get_sv_info_str(hdr, sv, "RIGHT_SVINSSEQ");
	if (!left_ins_seq.empty() && !right_ins_seq.empty()) {
		return left_ins_seq + "-" + right_ins_seq;
	}

	return "";
}

sv_t* bcf_to_sv(bcf_hdr_t* hdr, bcf1_t* b) {

	int* data = NULL;
	float* f_data = NULL;
	int len = 0;

	consensus_t* rc_consensus = NULL;
	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "REMAP_UB", &data, &len);
	if (len > 0) {
		rc_consensus = new consensus_t(false, 0, 0, 0, "", 0, 0, 0, 0, data[0], 0);
	}

	consensus_t* lc_consensus = NULL;
	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "REMAP_LB", &data, &len);
	if (len > 0) {
		lc_consensus = new consensus_t(true, 0, 0, 0, "", 0, 0, 0, 0, data[0], 0);
	}

	data = NULL;
	len = 0;
	int imprecise = bcf_get_info_flag(hdr, b, "IMPRECISE", &data, &len);

	std::string svtype = get_sv_type(hdr, b);
	hts_pos_t start = b->pos, end = get_sv_end(hdr, b);
	hts_pos_t left_anchor_mapping_start, left_anchor_mapping_end, right_anchor_mapping_start, right_anchor_mapping_end;
	if (svtype == "DUP") {
		left_anchor_mapping_start = start, left_anchor_mapping_end = start + 150;
		right_anchor_mapping_start = end - 150, right_anchor_mapping_end = end;
	} else if (svtype == "INV") {
		left_anchor_mapping_start = start - 150, left_anchor_mapping_end = start;
		right_anchor_mapping_start = start, right_anchor_mapping_end = start + 150;
	} else {
		left_anchor_mapping_start = start - 150, left_anchor_mapping_end = start;
		right_anchor_mapping_start = end, right_anchor_mapping_end = end + 150;
	}

	sv_t::anchor_aln_t* left_anchor_aln = new sv_t::anchor_aln_t(std::max(hts_pos_t(0), left_anchor_mapping_start), left_anchor_mapping_end, 0, 0);
	sv_t::anchor_aln_t* right_anchor_aln = new sv_t::anchor_aln_t(std::max(hts_pos_t(0), right_anchor_mapping_start), right_anchor_mapping_end, 0, 0);

	sv_t* sv;
	if (svtype == "DEL") {
		sv = new deletion_t(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln);
	} else if (svtype == "DUP") {
		sv = new duplication_t(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln);
	} else if (svtype == "INS") {
		sv = new insertion_t(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln);
	} else if (svtype == "INV") {
		sv_t::anchor_aln_t* lbp_left_anchor_aln = left_anchor_aln, *lbp_right_anchor_aln = right_anchor_aln;
		hts_pos_t rbp_left_split_mapping_start, rbp_left_split_mapping_end, rbp_right_split_mapping_start, rbp_right_split_mapping_end;
		rbp_left_split_mapping_start = end - 150, rbp_left_split_mapping_end = end;
		rbp_right_split_mapping_start = end, rbp_right_split_mapping_end = end + 150;
		sv_t::anchor_aln_t* rbp_left_anchor_aln = new sv_t::anchor_aln_t(std::max(hts_pos_t(0), rbp_left_split_mapping_start), rbp_left_split_mapping_end, 0, 0);
		sv_t::anchor_aln_t* rbp_right_anchor_aln = new sv_t::anchor_aln_t(std::max(hts_pos_t(0), rbp_right_split_mapping_start), rbp_right_split_mapping_end, 0, 0);
		sv = new inversion_t(bcf_seqname_safe(hdr, b), b->pos, get_sv_end(hdr, b), get_ins_seq(hdr, b), rc_consensus, lc_consensus, lbp_left_anchor_aln, lbp_right_anchor_aln, rbp_left_anchor_aln, rbp_right_anchor_aln);
	} else if (svtype == "BND") {
		std::string alt = b->d.allele[1];
		size_t colon_pos = alt.find(':');
		size_t pos_start = colon_pos + 1, pos_end = alt.find_last_of("[]");
    	hts_pos_t pos2 = std::stoi(alt.substr(pos_start, pos_end-pos_start));
		char dir;
		if (alt[pos_end] == ']') {
			dir = '+';
		} else if (alt[pos_end] == '[') {
			dir = '-';
		}
		sv = new breakend_t(bcf_seqname_safe(hdr, b), b->pos, pos2-1, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, dir);
	} else {
		return NULL;
	}
	sv->imprecise = (imprecise == 1);

	char* s_data = NULL;
	len = 0;
	int res = bcf_get_info_string(hdr, b, "SVINSSEQ_INFERRED", (void**) &s_data, &len);
	if (len > 0) {
		sv->inferred_ins_seq = s_data;
	}

	sv->id = b->d.id;
	sv->source = get_sv_info_str(hdr, b, "SOURCE");

	char** ss_data = NULL;
	len = 0;
	bcf_get_format_string(hdr, b, "FT", &ss_data, &len);
	std::string ft = "PASS";
	if (len > 0) {
		ft = ss_data[0];
	}

	std::istringstream tokenStream(ft);
	std::string token;
	while (std::getline(tokenStream, token, ';')) {
		std::string filter = token;
		if (filter != "PASS") sv->sample_info.filters.push_back(token);
	}

	int n = 0;
	sv->n_gt = bcf_get_genotypes(hdr, b, &(sv->sample_info.gt), &n);
	if (sv->n_gt < 0) sv->n_gt = 0;
	std::sort(sv->sample_info.gt, sv->sample_info.gt+sv->n_gt);

	return sv;
}

int find_sample_index(bcf_hdr_t* hdr, std::string sample_name) {
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        if (sample_name == hdr->samples[i]) {
            return i;
        }
    }
    return -1;  // Sample not found
}

#endif /* VCF_UTILS_H */
