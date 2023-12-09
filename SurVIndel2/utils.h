#ifndef SURVINDEL2_UTILS_H
#define SURVINDEL2_UTILS_H

#include <atomic>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unistd.h>
#include <numeric>
#include <chrono>
#include <ctime>

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "../libs/IntervalTree.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
#include "../src/types.h"
#include "../src/vcf_utils.h"


struct indel_t {
	sv_t* sv;
    int med_left_flanking_cov = 0, med_indel_left_cov = 0, med_indel_right_cov = 0, med_right_flanking_cov = 0;
    int med_left_cluster_cov = 0, med_right_cluster_cov = 0;
    bool remapped = false;

    indel_t(sv_t* sv) : sv(sv) {}

    bool imprecise() { return sv->lc_consensus == NULL && sv->rc_consensus == NULL && remapped == false; }
};

struct sv2_deletion_t : indel_t {
	static const int SIZE_NOT_COMPUTED = INT32_MAX;
	static const int REMAP_LB_NOT_COMPUTED = 0, REMAP_UB_NOT_COMPUTED = INT32_MAX;

    int max_conf_size = SIZE_NOT_COMPUTED, estimated_size = SIZE_NOT_COMPUTED, conc_pairs = 0;
    double ks_pval = -1.0;
    hts_pos_t remap_boundary_lower = REMAP_LB_NOT_COMPUTED, remap_boundary_upper = REMAP_UB_NOT_COMPUTED;
    int l_cluster_region_disc_pairs = 0, r_cluster_region_disc_pairs = 0;

    std::string original_range;
    std::string genotype;

    sv2_deletion_t(sv_t* sv) :
        indel_t(sv) {
		if (sv->rc_consensus != NULL) this->remap_boundary_upper = sv->rc_consensus->remap_boundary;
		if (sv->lc_consensus != NULL) this->remap_boundary_lower = sv->lc_consensus->remap_boundary;
	}
};

struct sv2_duplication_t : indel_t {
    hts_pos_t original_start, original_end;
    int lc_cluster_region_disc_pairs = 0, rc_cluster_region_disc_pairs = 0;

    sv2_duplication_t(sv_t* sv) :
		indel_t(sv), original_start(sv->start), original_end(sv->end) {
		if (sv->rc_consensus != NULL) this->original_end = sv->rc_consensus->breakpoint;
		if (sv->lc_consensus != NULL) this->original_start = sv->lc_consensus->breakpoint; 
	}
};

bcf_hdr_t* sv2_generate_vcf_header(chr_seqs_map_t& contigs, std::string& sample_name, config_t config, std::string command) {
	bcf_hdr_t* header = generate_vcf_header_base(contigs, sample_name, config, command);

	int len;

	// add FILTER tags
	const char* size_flt_tag = "##FILTER=<ID=SIZE_FILTER,Description=\"Size of the event is outside the predicted confidence interval.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	const char* remap_boundary_flt_tag = "##FILTER=<ID=REMAP_BOUNDARY_FILTER,Description=\"One of the breakpoints is incompatible with mate locations.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_boundary_flt_tag, &len));

	const char* depth_flt_tag = "##FILTER=<ID=DEPTH_FILTER,Description=\"Depth of the region is incompatible with type of the event. Only applicable to long events.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, depth_flt_tag, &len));

	const char* anom_del_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_DEL_DEPTH,Description=\"Depth of the deleted region is anomalous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anom_del_depth_flt_tag, &len));

	const char* not_enough_ow_pairs_flt_tag = "##FILTER=<ID=NOT_ENOUGH_OW_PAIRS,Description=\"Not enough outward oriented pairs supporting a large duplication.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, not_enough_ow_pairs_flt_tag, &len));

	const char* weak_split_aln_flt_tag = "##FILTER=<ID=WEAK_SPLIT_ALIGNMENT,Description=\"Split alignment not significantly better than full junction alignment.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_split_aln_flt_tag, &len));

	const char* low_mapq_cons_flt_tag = "##FILTER=<ID=LOW_MAPQ_CONSENSUSES,Description=\"No high MAPQ read supports the consensus(es) used to call this SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_mapq_cons_flt_tag, &len));

	const char* weak_support_flt_tag = "##FILTER=<ID=WEAK_SUPPORT,Description=\"Remapped breakpoint has low support from local reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_support_flt_tag, &len));

	const char* failed_to_ext_flt_tag = "##FILTER=<ID=FAILED_TO_EXTEND,Description=\"No reads can extend the consensus.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, failed_to_ext_flt_tag, &len));

	const char* low_ptn_ratio_flt_tag = "##FILTER=<ID=LOW_PTN_RATIO,Description=\"Low positive-to-negative ratio.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_ptn_ratio_flt_tag, &len));

	const char* ks_filter_flt_tag = "##FILTER=<ID=KS_FILTER,Description=\"According to KS test, the local IS distribution is not "
			"sufficiently different from the global IS distribution.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ks_filter_flt_tag, &len));

	const char* ambiguous_flt_tag = "##FILTER=<ID=AMBIGUOUS_REGION,Description=\"Region containing the deletion is ambiguous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ambiguous_flt_tag, &len));

	// add INFO tags
	const char* max_size_tag = "##INFO=<ID=MAX_SIZE,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_size_tag, &len));

	const char* ks_pval_tag = "##INFO=<ID=KS_PVAL,Number=1,Type=Float,Description=\"p-value of the KS test. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ks_pval_tag, &len));

	const char* est_size_tag = "##INFO=<ID=EST_SIZE,Number=1,Type=Integer,Description=\"Estimated size of the imprecise event. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, est_size_tag, &len));

	const char* median_depths_tag = "##INFO=<ID=MEDIAN_DEPTHS,Number=4,Type=Integer,Description=\"Depths of, respectively, the region flanking the indel to the left,"
			"the left portion of the indel, the right portion of the indel, the region flanking the indel to the right. Numbers 2 and 3 will be identical for short indels.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, median_depths_tag, &len));

	const char* cluster_depths_tag = "##INFO=<ID=CLUSTER_DEPTHS,Number=2,Type=Integer,Description=\"Depths of the left and right cluster regions.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cluster_depths_tag, &len));

	const char* disc_pairs_surr_tag = "##INFO=<ID=DISC_PAIRS_SURROUNDING,Number=2,Type=Integer,Description=\"Discordant pairs around the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_surr_tag, &len));

	const char* conc_pairs_tag = "##INFO=<ID=CONC_PAIRS,Number=1,Type=Integer,Description=\"Concordant pairs supporting the absence of a SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, conc_pairs_tag, &len));

	const char* imprecise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The reported boundaries are not precise.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, imprecise_tag, &len));

	const char* og_range_tag = "##INFO=<ID=ORIGINAL_RANGE,Number=1,Type=String,Description=\"Unadjusted imprecise range predicted by discordant pairs. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, og_range_tag, &len));

	// add samples
	bcf_hdr_add_sample(header, sample_name.c_str());

	return header;
}

void sv2_del2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, std::string& contig_name, sv2_deletion_t* del, std::vector<std::string>& filters) {
	sv2bcf(hdr, bcf_entry, del->sv, chr_seq, filters);
	
	int int_conv;

	if (del->max_conf_size != sv2_deletion_t::SIZE_NOT_COMPUTED) {
		bcf_update_info_int32(hdr, bcf_entry, "MAX_SIZE", &(del->max_conf_size), 1);
	}
	if (del->remap_boundary_lower != sv2_deletion_t::REMAP_LB_NOT_COMPUTED) {
		int_conv = del->remap_boundary_lower;
		bcf_update_info_int32(hdr, bcf_entry, "REMAP_LB", &int_conv, 1);
	}
	if (del->remap_boundary_upper != sv2_deletion_t::REMAP_UB_NOT_COMPUTED) {
		int_conv = del->remap_boundary_upper;
		bcf_update_info_int32(hdr, bcf_entry, "REMAP_UB", &int_conv, 1);
	}
	int median_depths[] = {del->med_left_flanking_cov, del->med_indel_left_cov, del->med_indel_right_cov, del->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	int cluster_depths[] = {del->med_left_cluster_cov, del->med_right_cluster_cov};
	bcf_update_info_int32(hdr, bcf_entry, "CLUSTER_DEPTHS", cluster_depths, 2);
	int disc_pairs_surr[] = {del->l_cluster_region_disc_pairs, del->r_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	bcf_update_info_int32(hdr, bcf_entry, "CONC_PAIRS", &del->conc_pairs, 1);

	bcf_update_info_flag(hdr, bcf_entry, "IMPRECISE", "", del->imprecise());
}

void sv2_dup2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, std::string& contig_name, sv2_duplication_t* dup, std::vector<std::string>& filters) {
	sv2bcf(hdr, bcf_entry, dup->sv, chr_seq, filters);

	// add info
	int int_conv;

	int median_depths[] = {dup->med_left_flanking_cov, dup->med_indel_left_cov, dup->med_indel_right_cov, dup->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	int disc_pairs_surr[] = {dup->rc_cluster_region_disc_pairs, dup->lc_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
}

#endif //SURVINDEL2_UTILS_H
