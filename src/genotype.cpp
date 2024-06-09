#include <cstdlib>
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

#include "sw_utils.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "reference_guided_assembly.h"
#include "extend_1sr_consensus.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
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

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1");
    const char* ar1_tag = "##FORMAT=<ID=AR1,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 in the alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,ar1_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ARC1");
    const char* arc1_tag = "##FORMAT=<ID=ARC1,Number=1,Type=Integer,Description=\"Number of consistent reads supporting the breakpoint 1 in the alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,arc1_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AR2");
    const char* ar2_tag = "##FORMAT=<ID=AR2,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 in the alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,ar2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ARC2");
    const char* arc2_tag = "##FORMAT=<ID=ARC2,Number=1,Type=Integer,Description=\"Number of consistent reads supporting the breakpoint 2 in the alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,arc2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ARCAS1");
    const char* arcas1_tag = "##FORMAT=<ID=ARCAS1,Number=1,Type=Float,Description=\"Average aln score of consistent reads supporting the first breakpoint of the SV to the alternate allele consensus.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,arcas1_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ARCAS2");
    const char* arcas2_tag = "##FORMAT=<ID=ARCAS2,Number=1,Type=Float,Description=\"Average aln score of consistent reads supporting the second breakpoint of the SV to the alternate allele consensus.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,arcas2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1");
    const char* rr1_tag = "##FORMAT=<ID=RR1,Number=1,Type=Integer,Description=\"Number of reads supporting the breakpoint 1 reference allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rr1_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "RRC1");
    const char* rrc1_tag = "##FORMAT=<ID=RRC1,Number=1,Type=Integer,Description=\"Number of consistent reads supporting the breakpoint 1 reference allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rrc1_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "RR2");
    const char* rr2_tag = "##FORMAT=<ID=RR2,Number=1,Type=Integer,Description=\"Number of reads supporting the breakpoint 2 reference allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rr2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "RRC2");
    const char* rrc2_tag = "##FORMAT=<ID=RRC2,Number=1,Type=Integer,Description=\"Number of consistent reads supporting the breakpoint 2 reference allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rrc2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "ER");
    const char* er_tag = "##FORMAT=<ID=ER,Number=1,Type=Integer,Description=\"Number of reads supporting equally well reference and alternate allele.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,er_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "TD");
    const char* nc_tag = "##FORMAT=<ID=TD,Number=1,Type=Integer,Description=\"The variant region is too deep to be genotyped reliably.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,nc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDLF");
    const char* mdlf_tag = "##FORMAT=<ID=MDLF,Number=1,Type=Integer,Description=\"Median depth of coverage in the left flanking region of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdlf_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDSP");
    const char* mdsp_tag = "##FORMAT=<ID=MDSP,Number=1,Type=Integer,Description=\"Median depth of coverage in the SV prefix region of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdsp_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDSF");
    const char* mdsf_tag = "##FORMAT=<ID=MDSF,Number=1,Type=Integer,Description=\"Median depth of coverage in the SV suffix region of the SV (for short SVs, it will be the same as MDSP).\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdsf_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDRF");
    const char* mdrf_tag = "##FORMAT=<ID=MDRF,Number=1,Type=Integer,Description=\"Median depth of coverage in the right flanking region of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdrf_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDLFHQ");
    const char* mdlfhq_tag = "##FORMAT=<ID=MDLFHQ,Number=1,Type=Integer,Description=\"Median depth of coverage in the left flanking region of the SV for high-quality reads.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdlfhq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDSPHQ");
    const char* mdsp_hq_tag = "##FORMAT=<ID=MDSPHQ,Number=1,Type=Integer,Description=\"Median depth of coverage in the SV prefix region of the SV for high-quality reads.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdsp_hq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDSFHQ");
    const char* mdsf_hq_tag = "##FORMAT=<ID=MDSFHQ,Number=1,Type=Integer,Description=\"Median depth of coverage in the SV suffix region of the SV for high-quality reads (for short SVs, it will be the same as MDSPHQ).\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdsf_hq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDRFHQ");
    const char* mdrf_hq_tag = "##FORMAT=<ID=MDRFHQ,Number=1,Type=Integer,Description=\"Median depth of coverage in the right flanking region of the SV for high-quality reads.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdrf_hq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDLC");
    const char* mdlc_tag = "##FORMAT=<ID=MDLC,Number=1,Type=Integer,Description=\"Median depth of coverage in the left cluster region of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdlc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MDRC");
    const char* mdrc_tag = "##FORMAT=<ID=MDRC,Number=1,Type=Integer,Description=\"Median depth of coverage in the right cluster region of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdrc_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MINSIZE");
    const char* minsize_tag = "##FORMAT=<ID=MINSIZE,Number=1,Type=Integer,Description=\"Minimum size of the event calculated based on insert size distribution."
            "Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, minsize_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MINSIZEHQ");
    const char* minsizehq_tag = "##FORMAT=<ID=MINSIZEHQ,Number=1,Type=Integer,Description=\"Minimum size of the event calculated based on insert size distribution for high-quality reads."
            "Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, minsizehq_tag, &len));
    
    bcf_hdr_remove(hdr, BCF_HL_FMT, "MAXSIZE");
    const char* maxsize_tag = "##FORMAT=<ID=MAXSIZE,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, maxsize_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "MAXSIZEHQ");
    const char* maxsizehq_tag = "##FORMAT=<ID=MAXSIZEHQ,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution for high-quality reads."
        	"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, maxsizehq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "KSPVAL");
    const char* kspval_tag = "##FORMAT=<ID=KSPVAL,Number=1,Type=Float,Description=\"p-value of the KS test.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, kspval_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "KSPVALHQ");
    const char* kspvalhq_tag = "##FORMAT=<ID=KSPVALHQ,Number=1,Type=Float,Description=\"p-value of the KS test for high-quality reads.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, kspvalhq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP1");
    const char* dp1_tag = "##FORMAT=<ID=DP1,Number=1,Type=Integer,Description=\"Number of discordant pairs supporting the first breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp1_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP2");
    const char* dp2_tag = "##FORMAT=<ID=DP2,Number=1,Type=Integer,Description=\"Number of discordant pairs supporting the second breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp2_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP1HQ");
    const char* dp1hq_tag = "##FORMAT=<ID=DP1HQ,Number=1,Type=Integer,Description=\"Number of high quality discordant pairs supporting the first breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp1hq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP2HQ");
    const char* dp2hq_tag = "##FORMAT=<ID=DP2HQ,Number=1,Type=Integer,Description=\"Number of high quality discordant pairs supporting the second breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp2hq_tag, &len));
    
    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP1MQ");
    const char* dp1mq_tag = "##FORMAT=<ID=DP1MQ,Number=1,Type=Integer,Description=\"Maximum mapping quality of discordant pairs supporting the first breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp1mq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DP2MQ");
    const char* dp2mq_tag = "##FORMAT=<ID=DP2MQ,Number=1,Type=Integer,Description=\"Maximum mapping quality of discordant pairs supporting the second breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dp2mq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DPLANM");
    const char* dplanm_tag = "##FORMAT=<ID=DPLANM,Number=1,Type=Float,Description=\"Average NM value of the left-most reads in the discordant pairs supporting this SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dplanm_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DPRANM");
    const char* dpranm_tag = "##FORMAT=<ID=DPRANM,Number=1,Type=Float,Description=\"Average NM value of the right-most reads in the discordant pairs supporting this SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpranm_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DPSL");
    const char* dpsl_tag = "##FORMAT=<ID=DPSL,Number=1,Type=Integer,Description=\"Number of discordant pairs surrounding but not supporting the left breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpsl_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "DPSR");
    const char* dpsr_tag = "##FORMAT=<ID=DPSR,Number=1,Type=Integer,Description=\"Number of discordant pairs surrounding but not supporting the right breakpoint of the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dpsr_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "CP");
    const char* cp_tag = "##FORMAT=<ID=CP,Number=1,Type=Integer,Description=\"Number of concordant pairs disproving the SV.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, cp_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AXR");
    const char* axr_tag = "##FORMAT=<ID=AXR,Number=1,Type=Integer,Description=\"Number of reads used to extend the alternative allele consensus.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, axr_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "AXRHQ");
    const char* axrhq_tag = "##FORMAT=<ID=AXRHQ,Number=1,Type=Integer,Description=\"Number of high-quality reads used to extend the alternative allele consensus.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, axrhq_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "EXL");
    const char* exl_tag = "##FORMAT=<ID=EXL,Number=1,Type=Integer,Description=\"Length of the extended alternative allele consensus.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exl_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "EXAS");
    const char* exas_tag = "##FORMAT=<ID=EXAS,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus and the original alternative allele (the reference with the SV applied).\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exas_tag, &len));

    bcf_hdr_remove(hdr, BCF_HL_FMT, "EXRS");
    const char* exrs_tag = "##FORMAT=<ID=EXRS,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus and the reference.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exrs_tag, &len));
}

void update_record(bcf_hdr_t* in_hdr, bcf_hdr_t* out_hdr, sv_t* sv, char* chr_seq, int sample_idx) {
    
    bcf_translate(out_hdr, in_hdr, sv->vcf_entry);

    bcf_subset(out_hdr, sv->vcf_entry, 1, &sample_idx);

    // update INFO fields
    int sv_end = sv->end+1;
    bcf_update_info_int32(out_hdr, sv->vcf_entry, "END", &sv_end, 1);

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
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1", &(sv->regenotyping_info.alt_bp1_better_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR2", &(sv->regenotyping_info.alt_bp2_better_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARC1", &(sv->regenotyping_info.alt_bp1_better_consistent_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARC2", &(sv->regenotyping_info.alt_bp2_better_consistent_reads), 1);
    float arcas1 = sv->regenotyping_info.alt_bp1_better_consistent_avg_score;
    bcf_update_format_float(out_hdr, sv->vcf_entry, "ARCAS1", &arcas1, 1);
    float arcas2 = sv->regenotyping_info.alt_bp2_better_consistent_avg_score;
    bcf_update_format_float(out_hdr, sv->vcf_entry, "ARCAS2", &arcas2, 1);

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

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDLF", &sv->median_left_flanking_cov, 1);
    if (sv->start != sv->end) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDSP", &sv->median_indel_left_cov, 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDSF", &sv->median_indel_right_cov, 1);
    }
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDRF", &sv->median_right_flanking_cov, 1);
    
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDLFHQ", &sv->median_left_flanking_cov_highmq, 1);
    if (sv->start != sv->end) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDSPHQ", &sv->median_indel_left_cov_highmq, 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDSFHQ", &sv->median_indel_right_cov_highmq, 1);
    }
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDRFHQ", &sv->median_right_flanking_cov_highmq, 1);
    
    if (sv->svtype() == "DEL") {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDLC", &sv->median_left_cluster_cov, 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDRC", &sv->median_right_cluster_cov, 1);
    }

    if (sv->min_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZE", &(sv->min_conf_size), 1);
    }
    if (sv->min_conf_size_highmq != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZEHQ", &(sv->min_conf_size_highmq), 1);
    }
    if (sv->max_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZE", &(sv->max_conf_size), 1);
    }
    if (sv->max_conf_size_highmq != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZEHQ", &(sv->max_conf_size_highmq), 1);
    }
    if (sv->ks_pval != deletion_t::KS_PVAL_NOT_COMPUTED) {
        float ks_pvalue = sv->ks_pval;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVAL", &ks_pvalue, 1);
    }
    if (sv->ks_pval_highmq != deletion_t::KS_PVAL_NOT_COMPUTED) {
        float ks_pvalue_highmq = sv->ks_pval_highmq;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVALHQ", &ks_pvalue_highmq, 1);
    }

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP1", &(sv->disc_pairs_lf), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP1HQ", &(sv->disc_pairs_lf_high_mapq), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP1MQ", &(sv->disc_pairs_lf_maxmapq), 1);
    float dplanm = sv->disc_pairs_lf_avg_nm;
    bcf_update_format_float(out_hdr, sv->vcf_entry, "DPLANM", &dplanm, 1);

    if (sv->svtype() == "INS") {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP2", &(sv->disc_pairs_rf), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP2HQ", &(sv->disc_pairs_rf_high_mapq), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "DP2MQ", &(sv->disc_pairs_rf_maxmapq), 1);
        float dpranm = sv->disc_pairs_rf_avg_nm;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "DPRANM", &dpranm, 1);
    }

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPSL", &(sv->l_cluster_region_disc_pairs), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "DPSR", &(sv->r_cluster_region_disc_pairs), 1);
    if (sv->svtype() != "DUP") {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "CP", &(sv->conc_pairs), 1);
    }

    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXR", &(sv->regenotyping_info.alt_ext_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXRHQ", &(sv->regenotyping_info.hq_alt_ext_reads), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXL", &(sv->regenotyping_info.ext_alt_consensus_length), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXAS", &(sv->regenotyping_info.ext_alt_consensus_to_alt_score), 1);
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "EXRS", &(sv->regenotyping_info.ext_alt_consensus_to_ref_score), 1);
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
    sv->l_cluster_region_disc_pairs = 0;
    sv->r_cluster_region_disc_pairs = 0;
}

std::vector<std::string> find_consistent_seqs_subset(std::string ref_seq, std::vector<std::string>& seqs, std::string& consensus_seq, double& avg_score) {
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln;
    std::vector<std::string> consistent_seqs;
    avg_score = 0;
    if (!seqs.empty()) {
        std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
        std::vector<std::string> v1, v2;
        std::vector<std::string> consensuses = generate_reference_guided_consensus(ref_seq, v1, seqs, v2, aligner, harsh_aligner, consensus_contigs_alns, config, stats);

        if (consensuses.empty()) {
            return consistent_seqs;
        }

        consensus_seq = consensuses[0];
        std::vector<int> start_positions, end_positions;
        for (std::string& seq : seqs) {
            aligner.Align(seq.c_str(), consensus_seq.c_str(), consensus_seq.length(), filter, &aln, 0);
            if (!is_left_clipped(aln, config.min_clip_len) && !is_right_clipped(aln, config.min_clip_len)) {
                consistent_seqs.push_back(seq);
                avg_score += double(aln.sw_score)/seq.length();
                start_positions.push_back(aln.ref_begin);
                end_positions.push_back(aln.ref_end);
            }
        }
        std::sort(start_positions.begin(), start_positions.end());
        std::sort(end_positions.begin(), end_positions.end());
        if (!consistent_seqs.empty()) avg_score /= consistent_seqs.size();

        if (start_positions.size() > 2) {
            int start = start_positions[2];
            int end = end_positions[start_positions.size()-3];
            consensus_seq = consensus_seq.substr(start, end-start);
        } else {
            consensus_seq = "";
        }
    }
    return consistent_seqs;
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
    candidate_reads_for_extension = get_extension_reads(contig_name, target_ivals, contig_len, mateseqs_w_mapq[contig_id], stats, bam_file);
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

    std::vector<std::string> alt_better_seqs, ref_bp1_better_seqs, ref_bp2_better_seqs;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < del_start || del_end < get_unclipped_start(read)) continue;
        if (del_start < get_unclipped_start(read) && get_unclipped_end(read) < del_end) continue;

        std::string seq = get_sequence(read);

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
            alt_better_seqs.push_back(seq);
        } else if (alt_aln.sw_score < ref_aln_score) {
            if (increase_ref_bp1_better) {
                ref_bp1_better_seqs.push_back(seq);
            }
            if (increase_ref_bp2_better) {
                ref_bp2_better_seqs.push_back(seq);
            }
        } else {
            same++;
        }

        if (alt_better_seqs.size() + ref_bp1_better_seqs.size() + ref_bp2_better_seqs.size() + same > 4 * stats.get_max_depth(del->chr)) {
            alt_better_seqs.clear();
            ref_bp1_better_seqs.clear();
            ref_bp2_better_seqs.clear();
            same = 0;
            del->regenotyping_info.too_deep = true;
            break;
        }
    }

    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    std::vector<std::string> alt_better_seqs_consistent = find_consistent_seqs_subset(alt_seq, alt_better_seqs, alt_consensus_seq, alt_avg_score);
    std::vector<std::string> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score);
    std::vector<std::string> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score);

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

        // align to ref
        hts_pos_t lbp_start = lh_start, lbp_end = del->start + alt_consensus_seq.length();
        hts_pos_t rbp_start = del->end - alt_consensus_seq.length(), rbp_end = rh_end;
        if (lbp_end > contig_len) lbp_end = contig_len;
        if (rbp_start < 0) rbp_start = 0;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+lbp_start, lbp_end-lbp_start, filter, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+rbp_start, rbp_end-rbp_start, filter, &ref2_aln, 0);
    
        del->regenotyping_info.ext_alt_consensus_length = alt_consensus->sequence.length();
        del->regenotyping_info.ext_alt_consensus_to_alt_score = alt_aln.sw_score;
        del->regenotyping_info.ext_alt_consensus_to_ref_score = std::max(ref1_aln.sw_score, ref2_aln.sw_score);
    }

    del->regenotyping_info.alt_bp1_better_reads = alt_better_seqs.size();
    del->regenotyping_info.alt_bp2_better_reads = alt_better_seqs.size();
    del->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_seqs_consistent.size();
    del->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_seqs_consistent.size();
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
    calculate_ptn_ratio(contig_name, large_deletions, bam_file, config, stats, true);
    calculate_cluster_region_disc(contig_name, dels, bam_file);
}

void genotype_small_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr) {
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

    std::vector<std::string> ref_better_seqs;
    std::vector<std::vector<std::string>> alt_better_seqs(alt_seqs.size());
    int same = 0;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;

        std::string seq = get_sequence(read);
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
                    alt_better_seqs[i].push_back(seq);
                }
            }
        } else if (best_aln_score < ref_aln.sw_score) {
            ref_better_seqs.push_back(seq);
        } else {
            same++;
        }
    }

    int alt_with_most_reads = 0;
    for (int i = 1; i < alt_better_seqs.size(); i++) {
        if (alt_better_seqs[i].size() > alt_better_seqs[alt_with_most_reads].size()) {
            alt_with_most_reads = i;
        }
    }

    if (alt_better_seqs[alt_with_most_reads].size() > 20*stats.get_max_depth(dup->chr) || ref_better_seqs.size() > 20*stats.get_max_depth(dup->chr)) {
        alt_better_seqs[alt_with_most_reads].clear();
        ref_better_seqs.clear();
        same = 0;
        dup->regenotyping_info.too_deep = true;
    }

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    std::vector<std::string> alt_better_seqs_consistent = find_consistent_seqs_subset(alt_seqs[alt_with_most_reads], alt_better_seqs[alt_with_most_reads], alt_consensus_seq, alt_avg_score);
    std::vector<std::string> ref_better_seqs_consistent = find_consistent_seqs_subset(ref_seq, ref_better_seqs, ref_consensus_seq, ref_avg_score);

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

    dup->regenotyping_info.alt_bp1_better_reads = alt_better_seqs[alt_with_most_reads].size();
    dup->regenotyping_info.alt_bp2_better_reads = alt_better_seqs[alt_with_most_reads].size();
    dup->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_seqs_consistent.size();
    dup->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_seqs_consistent.size();
    dup->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.ref_bp1_better_reads = ref_better_seqs.size();
    dup->regenotyping_info.ref_bp2_better_reads = ref_better_seqs.size();
    dup->regenotyping_info.ref_bp1_better_consistent_reads = ref_better_seqs_consistent.size();
    dup->regenotyping_info.ref_bp2_better_consistent_reads = ref_better_seqs_consistent.size();
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

    std::vector<std::string> alt_better_seqs, ref_bp1_better_seqs, ref_bp2_better_seqs;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;

        std::string seq = get_sequence(read);
        
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
            alt_better_seqs.push_back(seq);
        } else {
            if (increase_ref_bp1_better) {
                ref_bp1_better_seqs.push_back(seq);
            }
            if (increase_ref_bp2_better) {
                ref_bp2_better_seqs.push_back(seq);
            }
        }

        if (ref_bp1_better_seqs.size() + ref_bp2_better_seqs.size() > 4 * stats.get_max_depth(dup->chr)) {
            alt_better_seqs.clear();
            ref_bp1_better_seqs.clear();
            ref_bp2_better_seqs.clear();
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
    std::vector<std::string> alt_better_seqs_consistent = find_consistent_seqs_subset(alt_seq, alt_better_seqs, alt_consensus_seq, alt_avg_score);
    std::vector<std::string> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score);
    std::vector<std::string> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score);

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

        dup->regenotyping_info.ext_alt_consensus_length = alt_consensus->sequence.length();
        dup->regenotyping_info.ext_alt_consensus_to_alt_score = alt_aln.sw_score;
        dup->regenotyping_info.ext_alt_consensus_to_ref_score = std::max(ref1_aln.sw_score, ref2_aln.sw_score);
    }

    dup->regenotyping_info.alt_bp1_better_reads = alt_better_seqs.size();
    dup->regenotyping_info.alt_bp2_better_reads = alt_better_seqs.size();
    dup->regenotyping_info.alt_bp1_better_consistent_reads = alt_better_seqs_consistent.size();
    dup->regenotyping_info.alt_bp2_better_consistent_reads = alt_better_seqs_consistent.size();
    dup->regenotyping_info.alt_bp1_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.alt_bp2_better_consistent_avg_score = alt_avg_score;
    dup->regenotyping_info.ref_bp1_better_reads = ref_bp1_better_seqs.size();
    dup->regenotyping_info.ref_bp2_better_reads = ref_bp2_better_seqs.size();
    dup->regenotyping_info.ref_bp1_better_consistent_reads = ref_bp1_better_seqs_consistent.size();
    dup->regenotyping_info.ref_bp2_better_consistent_reads = ref_bp2_better_seqs_consistent.size();

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
    calculate_cluster_region_disc(contig_name, dups, bam_file, stats);
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

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    int same = 0;
    std::vector<std::string> alt_bp1_better_seqs, alt_bp2_better_seqs, ref_bp1_better_seqs, ref_bp2_better_seqs;

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
                alt_bp1_better_seqs.push_back(seq);
            } 
            if (alt1_aln.sw_score <= alt2_aln.sw_score) {
                alt_bp2_better_seqs.push_back(seq);
            }
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                ref_bp1_better_seqs.push_back(seq);
            } 
            if (ref1_aln.sw_score <= ref2_aln.sw_score) {
                ref_bp2_better_seqs.push_back(seq);
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
    std::vector<std::string> alt_bp1_better_seqs_consistent = find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_seqs, alt_bp1_consensus_seq, alt_bp1_avg_score);
    delete[] alt_bp1_seq;
    std::vector<std::string> alt_bp2_better_seqs_consistent = find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_seqs, alt_bp2_consensus_seq, alt_bp2_avg_score);
    delete[] alt_bp2_seq;
    
    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    std::vector<std::string> ref_bp1_better_seqs_consistent = find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, ref_bp1_consensus_seq, ref_bp1_avg_score);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    std::vector<std::string> ref_bp2_better_seqs_consistent = find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, ref_bp2_consensus_seq, ref_bp2_avg_score);
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

    int curr_pos = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        while (curr_pos < insertions.size() && insertions[curr_pos]->start < read->core.pos) curr_pos++;

        if (bam_is_rev(read)) continue;

        std::string qname = bam_get_qname(read);
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
                if (read->core.qual >= config.high_confidence_mapq) insertions[i]->disc_pairs_lf_high_mapq++;
                insertions[i]->disc_pairs_lf_maxmapq = std::max(insertions[i]->disc_pairs_lf_maxmapq, mateseqs_w_mapq_chr[qname].second);
                insertions[i]->disc_pairs_lf_avg_nm += get_nm(read);
            } else {
                insertions[i]->l_cluster_region_disc_pairs++;
            }
        }
    }
    for (insertion_t* ins : insertions) {
        if (ins->disc_pairs_lf > 0) ins->disc_pairs_lf_avg_nm /= ins->disc_pairs_lf;
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
                if (read->core.qual >= config.high_confidence_mapq) insertions[i]->disc_pairs_rf_high_mapq++;
                insertions[i]->disc_pairs_rf_maxmapq = std::max(insertions[i]->disc_pairs_rf_maxmapq, mateseqs_w_mapq_chr[qname].second);
                insertions[i]->disc_pairs_rf_avg_nm += get_nm(read);
            } else {
                insertions[i]->r_cluster_region_disc_pairs++;
            }
        }
    }
    for (insertion_t* ins : insertions) {
        if (ins->disc_pairs_rf > 0) ins->disc_pairs_rf_avg_nm /= ins->disc_pairs_rf;
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
    calculate_ptn_ratio(contig_name, inss, bam_file, stats);
    find_discordant_pairs(contig_name, inss, bam_file, stats, mateseqs_w_mapq[contig_id]);
    
    release_mates(contig_id);
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
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        sv_t* sv = bcf_to_sv(in_vcf_header, vcf_record);
        sv->vcf_entry = bcf_dup(vcf_record);
        reset_stats(sv);
        if (sv->svtype() == "DEL") {
            dels_by_chr[sv->chr].push_back((deletion_t*) sv);
        } else if (sv->svtype() == "DUP") {
            dups_by_chr[sv->chr].push_back((duplication_t*) sv);
        } else if (sv->svtype() == "INS") {
        	inss_by_chr[sv->chr].push_back((insertion_t*) sv);
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
    add_tags(out_vcf_header);
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
    	std::sort(contig_svs.begin(), contig_svs.end(), [](const sv_t* sv1, const sv_t* sv2) {return sv1->start < sv2->start;});

		for (auto& sv : contig_svs) {
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AC", NULL, 0);
			// bcf_update_info_int32(out_vcf_header, vcf_record, "AN", NULL, 0);
            update_record(in_vcf_header, out_vcf_header, sv, chr_seqs.get_seq(contig_name), imap[0]);
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
