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
    const char* type;         // Data type ("Integer" or "Float")
    const char* desc_format;  // Description format string
} format_tag_def_t;

void add_read_support_headers(bcf_hdr_t* hdr, char letter, int pos, const char* allele_type) {
    // Array of tag definitions
    const format_tag_def_t formats[] = {
        {
            "", 
            "Integer",
            "Number of reads supporting breakpoint %d in the %s allele."
        },
        {
            "C",
            "Integer", 
            "Number of consistent reads supporting breakpoint %d in the %s allele."
        },
        {
            "CF",
            "Integer",
            "Number of consistent forward reads supporting breakpoint %d in the %s allele."
        },
        {
            "CR",
            "Integer",
            "Number of consistent reverse reads supporting breakpoint %d in the %s allele."
        },
        {
            "CAS",
            "Float",
            "Average aln score of consistent reads supporting breakpoint %d of the SV to the %s allele consensus."
        },
        {
            "CHQ",
            "Integer",
            "Number of high-quality consistent reads supporting breakpoint %d in the %s allele."
        },
		{
            "CmQ",
            "Integer",
            "Minimum mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
        },
        {
            "CMQ",
            "Integer",
            "Maximum mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
        },
		{
            "CAQ",
            "Float",
            "Average mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
        },
		{
			"CSQ",
			"Float",
			"Standard deviation of mate mapping quality of consistent reads supporting breakpoint %d in the %s allele."
		}
    };

    char tag_id[10];
    char tag_str[256];
    char desc_str[200];
    int len;

    for (size_t i = 0; i < sizeof(formats) / sizeof(formats[0]); i++) {
        const format_tag_def_t* fmt = &formats[i];

        snprintf(tag_id, sizeof(tag_id), "%cR%d%s", letter, pos, fmt->suffix);
        snprintf(desc_str, sizeof(desc_str), fmt->desc_format, pos, allele_type);
        snprintf(tag_str, sizeof(tag_str), 
                "##FORMAT=<ID=%s,Number=1,Type=%s,Description=\"%s\">",
                tag_id, fmt->type, desc_str);

        bcf_hdr_remove(hdr, BCF_HL_FMT, tag_id);
        bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, tag_str, &len));
    }
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

	bcf_hdr_remove(hdr, BCF_HL_INFO, "HARD_FILTERS");
	const char* hf_tag = "##INFO=<ID=HARD_FILTERS,Number=.,Type=String,Description=\"PASS or not according to hard filters.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,hf_tag, &len));

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

	bcf_hdr_remove(hdr, BCF_HL_FMT, "DPS");
    const char* dps_tag = "##FORMAT=<ID=DPS,Number=2,Type=Integer,Description=\"Number of discordant pairs surrounding but not supporting the left and the right breakpoint of the SV, respectively.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dps_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "DPSHQ");
	const char* dpshq_tag = "##FORMAT=<ID=DPSHQ,Number=1,Type=Integer,Description=\"Number of high-quality discordant pairs surrounding but not supporting the left and the right breakpoint of the SV, respectively.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpshq_tag, &len));
	
	bcf_hdr_remove(hdr, BCF_HL_FMT, "DPSP");
	const char* dpsp_tag = "##FORMAT=<ID=DPSP,Number=2,Type=Integer,Description=\"Size of the region covered by discordant pairs supporting the left and right breakpoints of the SV, respectively.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpsp_tag, &len));

	add_read_support_headers(hdr, 'A', 1, "alternate");
	add_read_support_headers(hdr, 'A', 2, "alternate");

	add_read_support_headers(hdr, 'R', 1, "reference");
	add_read_support_headers(hdr, 'R', 2, "reference");

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

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP");
    const char* dp_tag = "##FORMAT=<ID=DP,Number=2,Type=Integer,Description=\"Number of discordant pairs supporting the first and second breakpoints of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DPHQ");
    const char* dphq_tag = "##FORMAT=<ID=DPHQ,Number=2,Type=Integer,Description=\"Number of high quality discordant pairs supporting the first and second breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dphq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DPMQ");
    const char* dpmq_tag = "##FORMAT=<ID=DPMQ,Number=2,Type=Integer,Description=\"Maximum mapping quality of discordant pairs supporting the first and the second breakpoints of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "DPNM");
	const char* dpnm_tag = "##FORMAT=<ID=DPNM,Number=2,Type=Float,Description=\"Average NM value of the left and right reads in the discordant pairs supporting this SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpnm_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "CP");
    const char* cp_tag = "##FORMAT=<ID=CP,Number=3,Type=Integer,Description=\"Number of concordant pairs crossing the left breakpoint, the mid point and the right breakpoint.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, cp_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "CPHQ");
	const char* cphq_tag = "##FORMAT=<ID=CPHQ,Number=3,Type=Integer,Description=\"Number of high-quality concordant pairs crossing the left breakpoint, the mid point and the right breakpoint.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, cphq_tag, &len));

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

	const char* homopolymer_flt_tag = "##FILTER=<ID=HOMOPOLYMER_INSSEQ,Description=\"Inserted sequence is a homopolymer.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, homopolymer_flt_tag, &len));

	const char* anomalous_sc_flt_tag = "##FILTER=<ID=ANOMALOUS_SC_NUMBER,Description=\"The number of soft-clipped reads supporting this call is too large.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_sc_flt_tag, &len));

	const char* size_flt_tag = "##FILTER=<ID=SIZE_FILTER,Description=\"Size of the event is outside the predicted confidence interval.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	const char* depth_flt_tag = "##FILTER=<ID=DEPTH_FILTER,Description=\"Depth of the region is incompatible with type of the event. Only applicable to long events.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, depth_flt_tag, &len));

	const char* anom_del_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_DEL_DEPTH,Description=\"Depth of the deleted region is anomalous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anom_del_depth_flt_tag, &len));

	const char* not_enough_disc_pairs_flt_tag = "##FILTER=<ID=NOT_ENOUGH_DISC_PAIRS,Description=\"Not enough discordant pairs supporting the event.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, not_enough_disc_pairs_flt_tag, &len));

	const char* weak_split_aln_flt_tag = "##FILTER=<ID=WEAK_SPLIT_ALIGNMENT,Description=\"Split alignment not significantly better than full junction alignment.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_split_aln_flt_tag, &len));

	const char* low_mapq_disc_pairs_flt_tag = "##FILTER=<ID=LOW_MAPQ_DISC_PAIRS,Description=\"Discordant pairs supporting the SV have low MAPQ.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_mapq_disc_pairs_flt_tag, &len));

	const char* low_mapq_cons_flt_tag = "##FILTER=<ID=LOW_MAPQ_CONSENSUSES,Description=\"No high MAPQ read supports the consensus(es) used to call this SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_mapq_cons_flt_tag, &len));

	const char* weak_support_flt_tag = "##FILTER=<ID=WEAK_SUPPORT,Description=\"Remapped breakpoint has low support from local reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_support_flt_tag, &len));

	const char* low_ptn_ratio_flt_tag = "##FILTER=<ID=LOW_PTN_RATIO,Description=\"Low positive-to-negative ratio.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_ptn_ratio_flt_tag, &len));

	const char* ks_filter_flt_tag = "##FILTER=<ID=KS_FILTER,Description=\"According to KS test, the local IS distribution is not "
			"sufficiently different from the global IS distribution.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ks_filter_flt_tag, &len));

	const char* ambiguous_flt_tag = "##FILTER=<ID=AMBIGUOUS_REGION,Description=\"Region containing the deletion is ambiguous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ambiguous_flt_tag, &len));

	const char* low_support_flt_tag = "##FILTER=<ID=LOW_SUPPORT,Description=\"Insertion has low support.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_support_flt_tag, &len));

	const char* weak_anchor_flt_tag = "##FILTER=<ID=WEAK_ANCHOR,Description=\"Left or right split junction alignment score is too low.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_anchor_flt_tag, &len));

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

	const char* mh_len = "##INFO=<ID=MH_LEN,Number=1,Type=Integer,Description=\"Length of the microhomology.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header,mh_len, &len));

	const char* source_tag = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, source_tag, &len));

	const char* fwd_sr_info_tag = "##INFO=<ID=FWD_SPLIT_READS,Number=2,Type=Integer,Description=\"Forward split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fwd_sr_info_tag, &len));

	const char* rev_sr_info_tag = "##INFO=<ID=REV_SPLIT_READS,Number=2,Type=Integer,Description=\"Reverse split reads supporting the left and right breakpoints of this ins.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rev_sr_info_tag, &len));

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

	const char* splitj_score_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE,Number=2,Type=Integer,Description=\"Score of the best alignment of the prefix and suffix of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score_tag, &len));

	const char* splitj_score2_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE2,Number=2,Type=Integer,Description=\"Score of the second best alignment of the prefix and suffix of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score2_tag, &len));

	const char* splitj_size_tag = "##INFO=<ID=SPLIT_JUNCTION_SIZE,Number=2,Type=Integer,Description=\"Size of the the prefix and suffix of the junction sequenc alignment to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_size_tag, &len));

	const char* splitj_size_lbp_tag = "##INFO=<ID=SPLIT_JUNCTION_SIZE_LBP,Number=2,Type=Integer,Description=\"Size of the the prefix and suffix of the alignment of the junction sequence supporting the left breakpoint to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_size_lbp_tag, &len));

	const char* splitj_size_rbp_tag = "##INFO=<ID=SPLIT_JUNCTION_SIZE_RBP,Number=2,Type=Integer,Description=\"Size of the the prefix and suffix of the alignment of the junction sequence supporting the right breakpoint to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_size_rbp_tag, &len));

	const char* splitj_cigar_tag = "##INFO=<ID=SPLIT_JUNCTION_CIGAR,Number=2,Type=String,Description=\"CIGAR of the best alignment of the prefix and suffix of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_cigar_tag, &len));

	const char* splitj_maprange_tag = "##INFO=<ID=SPLIT_JUNCTION_MAPPING_RANGE,Number=1,Type=String,Description=\"Original mapping locations of the prefix and suffix of the junction sequence to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_maprange_tag, &len));

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
			int_conv = sv->ins_seq.length();
			bcf_update_info_int32(hdr, bcf_entry, "SVINSLEN", &int_conv, 1);
		}
	}
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", sv->source.c_str());
	if (!sv->ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", sv->ins_seq.c_str());
	}
	if (sv->mh_len > 0) {
		bcf_update_info_int32(hdr, bcf_entry, "MH_LEN", &sv->mh_len, 1);
	}
	if (!sv->inferred_ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ_INFERRED", sv->inferred_ins_seq.c_str());
	}

	if (!for_gt) {
		int int2_conv[2];
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
		if (sv->svtype() == "INV") {
			inversion_t* inv = (inversion_t*) sv;
			int2_conv[0] = inv->left_anchor_aln->seq_len, int2_conv[1] = inv->right_anchor_aln->seq_len;
			bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE_LBP", int2_conv, 2);
			int2_conv[0] = inv->rbp_left_anchor_aln->seq_len, int2_conv[1] = inv->rbp_right_anchor_aln->seq_len;
			bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE_RBP", int2_conv, 2);
		} else {
			int2_conv[0] = sv->left_anchor_aln->seq_len, int2_conv[1] = sv->right_anchor_aln->seq_len;
			bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE", int2_conv, 2);
		}
		char* split_junction_cigar = (char*) malloc(sv->left_anchor_aln->cigar.length() + sv->right_anchor_aln->cigar.length() + 2);
		std::stringstream ss;
		ss << sv->left_anchor_aln->cigar << "," << sv->right_anchor_aln->cigar;
		bcf_update_info_string(hdr, bcf_entry, "SPLIT_JUNCTION_CIGAR", ss.str().c_str());
	}

	std::string split_junction_mapping_range = sv->left_anchor_aln_string() + "," + sv->right_anchor_aln_string();
	if (sv->svtype() == "INV") {
		inversion_t* inv = (inversion_t*) sv;
		split_junction_mapping_range += "," + inv->rbp_left_anchor_aln->to_string() + "," + inv->rbp_right_anchor_aln->to_string();
	}
	bcf_update_info_string(hdr, bcf_entry, "SPLIT_JUNCTION_MAPPING_RANGE", split_junction_mapping_range.c_str());

	if (!for_gt) {
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

		int max_mapq1 = sv->rc_consensus ? (int) sv->rc_consensus->max_mapq : 0;
		int max_mapq2 = sv->lc_consensus ? (int) sv->lc_consensus->max_mapq : 0;
		bcf_update_info_int32(hdr, bcf_entry, "ARC1MQ", &max_mapq1, 1);
		bcf_update_info_int32(hdr, bcf_entry, "ARC2MQ", &max_mapq2, 1);

		int median_depths[] = {sv->median_left_flanking_cov, sv->median_indel_left_cov, sv->median_indel_right_cov, sv->median_right_flanking_cov};
		bcf_update_format_int32(hdr, bcf_entry, "MD", median_depths, 4);

		int median_depths_highmq[] = {sv->median_left_flanking_cov_highmq, sv->median_indel_left_cov_highmq, sv->median_indel_right_cov_highmq, sv->median_right_flanking_cov_highmq};
		bcf_update_format_int32(hdr, bcf_entry, "MDHQ", median_depths_highmq, 4);

		int cluster_depths[] = {sv->median_left_cluster_cov, sv->median_right_cluster_cov};
		bcf_update_format_int32(hdr, bcf_entry, "CLMD", cluster_depths, 2);

		int cluster_depths_highmq[] = {sv->median_left_cluster_cov_highmq, sv->median_right_cluster_cov_highmq};
		bcf_update_format_int32(hdr, bcf_entry, "CLMDHQ", cluster_depths_highmq, 2);

		int disc_pairs_surr[] = {sv->l_cluster_region_disc_pairs, sv->r_cluster_region_disc_pairs};
		bcf_update_format_int32(hdr, bcf_entry, "DPS", disc_pairs_surr, 2);
		
		int dp[] = {sv->disc_pairs_lf, sv->disc_pairs_rf};
		bcf_update_format_int32(hdr, bcf_entry, "DP", dp, 2);
		if (sv->disc_pairs_lf + sv->disc_pairs_rf > 0) {
			int disc_pairs_surr_hq[] = {sv->l_cluster_region_disc_pairs_high_mapq, sv->r_cluster_region_disc_pairs_high_mapq};
			bcf_update_format_int32(hdr, bcf_entry, "DPSHQ", disc_pairs_surr_hq, 2);

			int dphq[] = {sv->disc_pairs_lf_high_mapq, sv->disc_pairs_rf_high_mapq};
			bcf_update_format_int32(hdr, bcf_entry, "DPHQ", dphq, 2);

			int dpmq[] = {sv->disc_pairs_lf_high_mapq, sv->disc_pairs_rf_high_mapq};
			bcf_update_format_int32(hdr, bcf_entry, "DPMQ", dpmq, 2);
		}

		float dpnm[] = {(float) sv->disc_pairs_lf_avg_nm, (float) sv->disc_pairs_rf_avg_nm};
    	bcf_update_format_float(hdr, bcf_entry, "DPNM", &dpnm, 2);

		int conc_pairs[] = {sv->conc_pairs_lbp, sv->conc_pairs_midp, sv->conc_pairs_rbp};
		bcf_update_format_int32(hdr, bcf_entry, "CP", conc_pairs, 3);

		int conc_pairs_hq[] = {sv->conc_pairs_lbp_high_mapq, sv->conc_pairs_midp_high_mapq, sv->conc_pairs_rbp_high_mapq};
		bcf_update_format_int32(hdr, bcf_entry, "CPHQ", conc_pairs_hq, 3);

		base_frequencies_t left_anchor_base_freqs = sv->get_left_anchor_base_freqs(chr_seq);
		int labc[] = {left_anchor_base_freqs.a, left_anchor_base_freqs.c, left_anchor_base_freqs.g, left_anchor_base_freqs.t};
		bcf_update_info_int32(hdr, bcf_entry, "LEFT_ANCHOR_BASE_COUNT", labc, 4);

		base_frequencies_t right_anchor_base_freqs = sv->get_right_anchor_base_freqs(chr_seq);
		int rabc[] = {right_anchor_base_freqs.a, right_anchor_base_freqs.c, right_anchor_base_freqs.g, right_anchor_base_freqs.t};
		bcf_update_info_int32(hdr, bcf_entry, "RIGHT_ANCHOR_BASE_COUNT", rabc, 4);

		base_frequencies_t prefix_ref_base_freqs = sv->get_prefix_ref_base_freqs(chr_seq);
		int svrefpbc[] = {prefix_ref_base_freqs.a, prefix_ref_base_freqs.c, prefix_ref_base_freqs.g, prefix_ref_base_freqs.t};
		bcf_update_info_int32(hdr, bcf_entry, "SV_REF_PREFIX_BASE_COUNT", svrefpbc, 4);
		
		base_frequencies_t suffix_ref_base_freqs = sv->get_suffix_ref_base_freqs(chr_seq);
		int svrefsbc[] = {suffix_ref_base_freqs.a, suffix_ref_base_freqs.c, suffix_ref_base_freqs.g, suffix_ref_base_freqs.t};
		bcf_update_info_int32(hdr, bcf_entry, "SV_REF_SUFFIX_BASE_COUNT", svrefsbc, 4);

		base_frequencies_t ins_prefix_base_freqs = sv->get_ins_prefix_base_freqs();
		int pbc[] = {ins_prefix_base_freqs.a, ins_prefix_base_freqs.c, ins_prefix_base_freqs.g, ins_prefix_base_freqs.t};
		bcf_update_info_int32(hdr, bcf_entry, "INS_PREFIX_BASE_COUNT", pbc, 4);
		
		base_frequencies_t ins_suffix_base_freqs = sv->get_ins_suffix_base_freqs();
		int sbc[] = {ins_suffix_base_freqs.a, ins_suffix_base_freqs.c, ins_suffix_base_freqs.g, ins_suffix_base_freqs.t};
		bcf_update_info_int32(hdr, bcf_entry, "INS_SUFFIX_BASE_COUNT", sbc, 4);
	}

	bcf_update_info_flag(hdr, bcf_entry, "IMPRECISE", "", sv->imprecise);
	bcf_update_genotypes(hdr, bcf_entry, sv->sample_info.gt, sv->n_gt);

	if (sv->incomplete_ins_seq()) {
		bcf_update_info_flag(hdr, bcf_entry, "INCOMPLETE_ASSEMBLY", "", 1);
	}

	const char* ft_val = sv->is_pass() ? "PASS" : "FAIL";
	bcf_update_format_string(hdr, bcf_entry, "FT", &ft_val, sv->n_gt);

	if (sv->svtype() == "DEL") {
		if (!for_gt) {
			deletion_t* del = (deletion_t*) sv;
			if (!del->original_range.empty()) {
				bcf_update_info_string(hdr, bcf_entry, "ORIGINAL_RANGE", del->original_range.c_str());
			}
			if (del->min_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
				bcf_update_format_int32(hdr, bcf_entry, "MINSIZE", &(del->min_conf_size), 1);
			}
			if (del->max_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
				bcf_update_format_int32(hdr, bcf_entry, "MAXSIZE", &(del->max_conf_size), 1);
			}
			if (del->imprecise && del->estimated_size != deletion_t::SIZE_NOT_COMPUTED) {
				bcf_update_info_int32(hdr, bcf_entry, "EST_SIZE", &(del->estimated_size), 1);
			}
			if (del->ks_pval != deletion_t::KS_PVAL_NOT_COMPUTED) {
				float ks_pval = del->ks_pval;
				bcf_update_format_float(hdr, bcf_entry, "KSPVAL", &ks_pval, 1);
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

	int rc_fwd_reads = 0, rc_rev_reads = 0, lc_fwd_reads = 0, lc_rev_reads = 0;
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

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "ARC1MQ", &data, &len);
	int max_rc_mapq = 0, max_lc_mapq = 0;
	if (len > 0) {
		max_rc_mapq = data[0];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "ARC2MQ", &data, &len);
	if (len > 0) {
		max_lc_mapq = data[0];
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

	int full_junction_score = 0;
	len = 0;
	bcf_get_info_int32(hdr, b, "FULL_JUNCTION_SCORE", &data, &len);
	if (len > 0) full_junction_score = data[0];

	std::string full_junction_cigar = "";
	char* s_data = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "FULL_JUNCTION_CIGAR", (void**) &s_data, &len);
	if (s_data) full_junction_cigar = s_data;

	std::string svtype = get_sv_type(hdr, b);

	s_data = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "SPLIT_JUNCTION_MAPPING_RANGE", (void**) &s_data, &len);

	std::vector<std::string> split_junction_mapping_ranges;
	std::string mapping_range_str;
	if (s_data != NULL) {
		mapping_range_str = s_data;
		std::string mapping_range;
		std::istringstream tokenStream(mapping_range_str);
		while (std::getline(tokenStream, mapping_range, ',')) {
			split_junction_mapping_ranges.push_back(mapping_range);
		}
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SCORE", &data, &len);
	int left_split_score = 0, right_split_score = 0;
	if (len > 0) {
		left_split_score = data[0], right_split_score = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SCORE2", &data, &len);
	int left_split_score2 = 0, right_split_score2 = 0;
	if (len > 0) {
		left_split_score2 = data[0], right_split_score2 = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SIZE", &data, &len);
	int left_split_size = 0, right_split_size = 0;
	if (len > 0) {
		left_split_size = data[0], right_split_size = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SIZE_LBP", &data, &len);
	if (len > 0) {
		left_split_size = data[0], right_split_size = data[1];
	}

	data = NULL;
	len = 0;
	int rbp_left_split_size = 0, rbp_right_split_size = 0;
	bcf_get_info_int32(hdr, b, "SPLIT_JUNCTION_SIZE_RBP", &data, &len);
	if (len > 0) {
		rbp_left_split_size = data[0], rbp_right_split_size = data[1];
	}

	s_data = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "SPLIT_JUNCTION_CIGAR", (void**) &s_data, &len);
	std::string left_split_cigar = "", right_split_cigar = "";
	if (s_data) {
		std::string cigar = s_data;
		size_t comma_pos = cigar.find(",");
		left_split_cigar = cigar.substr(0, comma_pos);
		right_split_cigar = cigar.substr(comma_pos+1);
	}

	data = NULL;
	len = 0;
	int imprecise = bcf_get_info_flag(hdr, b, "IMPRECISE", &data, &len);

	// if the SV does not have a SPLIT_JUNCTION_MAPPING_RANGE, set one using a 150 bp window around the SV
	hts_pos_t left_split_mapping_start, left_split_mapping_end, right_split_mapping_start, right_split_mapping_end;
	if (split_junction_mapping_ranges.empty()) {
		if (svtype == "DUP") {
			left_split_mapping_start = b->pos, left_split_mapping_end = b->pos + 150;
			right_split_mapping_start = std::max(0, get_sv_end(hdr, b) - 150), right_split_mapping_end = get_sv_end(hdr, b);
		} else {
			left_split_mapping_start = std::max(hts_pos_t(0), b->pos - 150), left_split_mapping_end = b->pos;
			right_split_mapping_start = get_sv_end(hdr, b), right_split_mapping_end = get_sv_end(hdr, b) + 150;
		}
	} else {
		left_split_mapping_start = std::stoi(split_junction_mapping_ranges[0].substr(0, split_junction_mapping_ranges[0].find("-")))-1;
		left_split_mapping_end = std::stoi(split_junction_mapping_ranges[0].substr(split_junction_mapping_ranges[0].find("-")+1))-1;
		right_split_mapping_start = std::stoi(split_junction_mapping_ranges[1].substr(0, split_junction_mapping_ranges[1].find("-")))-1;
		right_split_mapping_end = std::stoi(split_junction_mapping_ranges[1].substr(split_junction_mapping_ranges[1].find("-")+1))-1;
	}
	sv_t::anchor_aln_t* left_anchor_aln = new sv_t::anchor_aln_t(left_split_mapping_start, left_split_mapping_end, left_split_size, left_split_score, left_split_score2, left_split_cigar);
	sv_t::anchor_aln_t* right_anchor_aln = new sv_t::anchor_aln_t(right_split_mapping_start, right_split_mapping_end, right_split_size, right_split_score, right_split_score2, right_split_cigar);
	sv_t::anchor_aln_t* full_junction_aln = full_junction_score > 0 ? new sv_t::anchor_aln_t(0, 0, 0, full_junction_score, 0, full_junction_cigar) : NULL;

	sv_t* sv;
	hts_pos_t end = get_sv_end(hdr, b);
	if (svtype == "DEL") {
		sv = new deletion_t(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, full_junction_aln);
	} else if (svtype == "DUP") {
		sv = new duplication_t(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, full_junction_aln);
	} else if (svtype == "INS") {
		sv = new insertion_t(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln, full_junction_aln);
	} else if (svtype == "INV") {
		sv_t::anchor_aln_t* lbp_left_anchor_aln = left_anchor_aln, *lbp_right_anchor_aln = right_anchor_aln;
		hts_pos_t rbp_left_split_mapping_start, rbp_left_split_mapping_end, rbp_right_split_mapping_start, rbp_right_split_mapping_end;
		if (split_junction_mapping_ranges.size() > 2){
			rbp_left_split_mapping_start = std::stoi(split_junction_mapping_ranges[2].substr(0, split_junction_mapping_ranges[2].find("-")))-1;
			rbp_left_split_mapping_end = std::stoi(split_junction_mapping_ranges[2].substr(split_junction_mapping_ranges[2].find("-")+1))-1;
			rbp_right_split_mapping_start = std::stoi(split_junction_mapping_ranges[3].substr(0, split_junction_mapping_ranges[3].find("-")))-1;
			rbp_right_split_mapping_end = std::stoi(split_junction_mapping_ranges[3].substr(split_junction_mapping_ranges[3].find("-")+1))-1;
		} else {
			left_split_mapping_end = b->pos, right_split_mapping_start = b->pos, right_split_mapping_end = b->pos + 150;
			rbp_left_split_mapping_start = std::max(hts_pos_t(0), end - 150), rbp_left_split_mapping_end = end;
			rbp_right_split_mapping_start = end, rbp_right_split_mapping_end = end + 150;
		}
		sv_t::anchor_aln_t* rbp_left_anchor_aln = new sv_t::anchor_aln_t(rbp_left_split_mapping_start, rbp_left_split_mapping_end, rbp_left_split_size, 0, 0, "");
		sv_t::anchor_aln_t* rbp_right_anchor_aln = new sv_t::anchor_aln_t(rbp_right_split_mapping_start, rbp_right_split_mapping_end, rbp_right_split_size, 0, 0, "");
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

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "MH_LEN", &data, &len);
	if (len > 0) {
		sv->mh_len = data[0];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "DP", &data, &len);
	if (len > 0) {
		sv->disc_pairs_lf = data[0];
		sv->disc_pairs_rf = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "DPHQ", &data, &len);
	if (len > 0) {
		sv->disc_pairs_lf_high_mapq = data[0];
		sv->disc_pairs_rf_high_mapq = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "DPMQ", &data, &len);
	if (len > 0) {
		sv->disc_pairs_lf_maxmapq = data[0];
		sv->disc_pairs_rf_maxmapq = data[1];
	}

	f_data = NULL;
	len = 0;
	bcf_get_format_float(hdr, b, "DPNM", &f_data, &len);
	if (len > 0) {
		sv->disc_pairs_lf_avg_nm = f_data[0];
		sv->disc_pairs_rf_avg_nm = f_data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "CP", &data, &len);
	if (len > 0) {
		sv->conc_pairs_lbp = data[0];
		sv->conc_pairs_midp = data[1];
		sv->conc_pairs_rbp = data[2];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "CPHQ", &data, &len);
	if (len > 0) {
		sv->conc_pairs_lbp_high_mapq = data[0];
		sv->conc_pairs_midp_high_mapq = data[1];
		sv->conc_pairs_rbp_high_mapq = data[2];
	}

	f_data = NULL;
	len = 0;
	bcf_get_format_float(hdr, b, "KSPVAL", &f_data, &len);
	if (len > 0) {
		sv->ks_pval = f_data[0];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "MINSIZE", &data, &len);
	if (len > 0) {
		sv->min_conf_size = data[0];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "MAXSIZE", &data, &len);
	if (len > 0) {
		sv->max_conf_size = data[0];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "MD", &data, &len);
	if (len > 0) {
		sv->median_left_flanking_cov = data[0];
		sv->median_indel_left_cov = data[1];
		sv->median_indel_right_cov = data[2];
		sv->median_right_flanking_cov = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "MDHQ", &data, &len);
	if (len > 0) {
		sv->median_left_flanking_cov_highmq = data[0];
		sv->median_indel_left_cov_highmq = data[1];
		sv->median_indel_right_cov_highmq = data[2];
		sv->median_right_flanking_cov_highmq = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "CLMD", &data, &len);
	if (len > 0) {
		sv->median_left_cluster_cov = data[0];
		sv->median_right_cluster_cov = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "CLMDHQ", &data, &len);
	if (len > 0) {
		sv->median_left_cluster_cov_highmq = data[0];
		sv->median_right_cluster_cov_highmq = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "DPS", &data, &len);
	if (len > 0) {
		sv->l_cluster_region_disc_pairs = data[0];
		sv->r_cluster_region_disc_pairs = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_format_int32(hdr, b, "DPSHQ", &data, &len);
	if (len > 0) {
		sv->l_cluster_region_disc_pairs_high_mapq = data[0];
		sv->r_cluster_region_disc_pairs_high_mapq = data[1];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "LEFT_ANCHOR_BASE_COUNT", &data, &len);
	if (len > 0) {
		sv->left_anchor_base_freqs.a = data[0];
		sv->left_anchor_base_freqs.c = data[1];
		sv->left_anchor_base_freqs.g = data[2];
		sv->left_anchor_base_freqs.t = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "RIGHT_ANCHOR_BASE_COUNT", &data, &len);
	if (len > 0) {
		sv->right_anchor_base_freqs.a = data[0];
		sv->right_anchor_base_freqs.c = data[1];
		sv->right_anchor_base_freqs.g = data[2];
		sv->right_anchor_base_freqs.t = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SV_REF_PREFIX_BASE_COUNT", &data, &len);
	if (len > 0) {
		sv->prefix_ref_base_freqs.a = data[0];
		sv->prefix_ref_base_freqs.c = data[1];
		sv->prefix_ref_base_freqs.g = data[2];
		sv->prefix_ref_base_freqs.t = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "SV_REF_SUFFIX_BASE_COUNT", &data, &len);
	if (len > 0) {
		sv->suffix_ref_base_freqs.a = data[0];
		sv->suffix_ref_base_freqs.c = data[1];
		sv->suffix_ref_base_freqs.g = data[2];
		sv->suffix_ref_base_freqs.t = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "INS_PREFIX_BASE_COUNT", &data, &len);
	if (len > 0) {
		sv->ins_prefix_base_freqs.a = data[0];
		sv->ins_prefix_base_freqs.c = data[1];
		sv->ins_prefix_base_freqs.g = data[2];
		sv->ins_prefix_base_freqs.t = data[3];
	}

	data = NULL;
	len = 0;
	bcf_get_info_int32(hdr, b, "INS_SUFFIX_BASE_COUNT", &data, &len);
	if (len > 0) {
		sv->ins_suffix_base_freqs.a = data[0];
		sv->ins_suffix_base_freqs.c = data[1];
		sv->ins_suffix_base_freqs.g = data[2];
		sv->ins_suffix_base_freqs.t = data[3];
	}

	s_data = NULL;
	len = 0;
	int res = bcf_get_info_string(hdr, b, "SVINSSEQ_INFERRED", (void**) &s_data, &len);
	if (len > 0) {
		sv->inferred_ins_seq = s_data;
	}

	sv->id = b->d.id;
	sv->source = get_sv_info_str(hdr, b, "SOURCE");

	for (int i = 0; i < b->d.n_flt; i++) {
		sv->filters.push_back(bcf_hdr_int2id(hdr, BCF_DT_ID, b->d.flt[i]));
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
