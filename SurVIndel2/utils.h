#ifndef SURVINDEL2_UTILS_H
#define SURVINDEL2_UTILS_H

#include <atomic>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "htslib/kseq.h"
#include <unistd.h>
#include <numeric>
#include <chrono>
#include <ctime>
KSEQ_INIT(int, read)

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "../libs/IntervalTree.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
#include "../src/types.h"

struct config_t {

    int threads, seed;
    int min_is, max_is; // find a way to move this to stats_t
    int min_sv_size;
    int read_len; // this is not exactly "config", but it is more convenient to place it here
    int min_clip_len;
    double max_seq_error;
    int max_clipped_pos_dist;
    int min_size_for_depth_filtering;
    int min_diff_hsr;
    std::string sampling_regions, version;
    bool log;

    int clip_penalty = 7;
    int min_score_diff = 15;
    int high_confidence_mapq = 60;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        threads = std::stoi(config_params["threads"]);
        seed = std::stoi(config_params["seed"]);
        min_is = std::stoi(config_params["min_is"]);
        max_is = std::stoi(config_params["max_is"]);
        min_sv_size = std::stoi(config_params["min_sv_size"]);
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        read_len = std::stod(config_params["read_len"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        max_clipped_pos_dist = std::stoi(config_params["max_clipped_pos_dist"]);
        min_size_for_depth_filtering = std::stoi(config_params["min_size_for_depth_filtering"]);
        min_diff_hsr = std::stoi(config_params["min_diff_hsr"]);
        sampling_regions = config_params["sampling_regions"];
        version = config_params["version"];
        log = std::stoi(config_params["log"]);
    }

    std::string clip_penalty_padding() { return std::string(this->clip_penalty, 'N'); }
};

struct stats_t {
    double median_depth;
    int min_depth, max_depth;
    int pop_avg_crossing_is = 0;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        median_depth = std::stod(config_params["median_depth"]);
        pop_avg_crossing_is = std::stoi(config_params["pop_avg_crossing_is"]);
        min_depth = std::stoi(config_params["min_depth"]);
        max_depth = std::stoi(config_params["max_depth"]);
    }
};

struct indel_t {
	sv_t* sv;
    int disc_pairs = 0, disc_pairs_high_mapq = 0, disc_pairs_maxmapq = 0;
    consensus_t* lc_consensus = NULL,* rc_consensus = NULL;
    int med_left_flanking_cov = 0, med_indel_left_cov = 0, med_indel_right_cov = 0, med_right_flanking_cov = 0;
    int med_left_cluster_cov = 0, med_right_cluster_cov = 0;
    bool remapped = false;
    std::string rightmost_rightfacing_seq, leftmost_leftfacing_seq;

    indel_t(sv_t* sv, consensus_t* lc_consensus, consensus_t* rc_consensus)
    : sv(sv), lc_consensus(lc_consensus), rc_consensus(rc_consensus) { }

	bool is_single_consensus() { return (lc_consensus == NULL || rc_consensus == NULL) && lc_consensus != rc_consensus; }
    bool imprecise() { return lc_consensus == NULL && rc_consensus == NULL && remapped == false; }
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

    sv2_deletion_t(sv_t* sv, consensus_t* lc_consensus, consensus_t* rc_consensus) :
        indel_t(sv, lc_consensus, rc_consensus) {}
};

struct sv2_duplication_t : indel_t {
    hts_pos_t original_start, original_end;
    int lc_cluster_region_disc_pairs = 0, rc_cluster_region_disc_pairs = 0;

    sv2_duplication_t(sv_t* sv, consensus_t* lc_consensus, consensus_t* rc_consensus) :
		indel_t(sv, lc_consensus, rc_consensus), original_start(sv->start), original_end(sv->end) {}
};

indel_t* sv_to_indel(sv_t* sv, consensus_t* rc_consensus, consensus_t* lc_consensus) {
	if (sv == NULL) return NULL;

	indel_t* indel;
	if (sv->svtype() == "DEL") {
		sv2_deletion_t* del = new sv2_deletion_t(sv, lc_consensus, rc_consensus);
		if (rc_consensus != NULL) del->remap_boundary_upper = rc_consensus->remap_boundary;
		if (lc_consensus != NULL) del->remap_boundary_lower = lc_consensus->remap_boundary;
		indel = del;
	} else {
		sv2_duplication_t* dup = new sv2_duplication_t(sv, lc_consensus, rc_consensus);
		if (rc_consensus != NULL) dup->original_end = rc_consensus->breakpoint;
		if (lc_consensus != NULL) dup->original_start = lc_consensus->breakpoint; 
		indel = dup;
	}
	return indel;
}

struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    contig_map_t() {}
    contig_map_t(std::string workdir) { load(workdir); }

    void load(std::string workdir) {
        std::ifstream fin(workdir + "/contig_map");
        std::string name;
        int id = 0;
        while (fin >> name) {
            name_to_id[name] = id;
            id_to_name.push_back(name);
            id++;
        }
    }

    size_t size() {return id_to_name.size();}
    std::string get_name(size_t id) {return id_to_name[id];};
};


struct chr_seq_t {
    char* seq;
    hts_pos_t len;

    chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
    ~chr_seq_t() {delete[] seq;}
};
struct chr_seqs_map_t {
    std::unordered_map<std::string, chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        return seqs[seq_name]->seq;
    }

    hts_pos_t get_len(std::string seq_name) {
        return seqs[seq_name]->len;
    }

    void clear() {
        for (auto& e : seqs) {
            delete e.second;
            e.second = NULL;
        }
    }

    ~chr_seqs_map_t() {
        clear();
    }
};

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

bool is_homopolymer(const char* seq, int len) {
	int a = 0, c = 0, g = 0, t = 0;
	for (int i = 0; i < len; i++) {
		char b = std::toupper(seq[i]);
		if (b == 'A') a++;
		else if (b == 'C') c++;
		else if (b == 'G') g++;
		else if (b == 'T') t++;
	}
	return max(a, c, g, t)/double(a+c+g+t) >= 0.8;
}

template<typename T>
T mean(std::vector<T>& v) {
    return std::accumulate(v.begin(), v.end(), (T)0.0)/v.size();
}

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
bcf_hdr_t* generate_vcf_header(chr_seqs_map_t& contigs, std::string& sample_name, config_t config, std::string command) {
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
	const char* small_flt_tag = "##FILTER=<ID=SMALL,Description=\"Event is smaller than what required by the user.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, small_flt_tag, &len));

	const char* size_flt_tag = "##FILTER=<ID=SIZE_FILTER,Description=\"Size of the event is outside the predicted confidence interval.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	const char* remap_boundary_flt_tag = "##FILTER=<ID=REMAP_BOUNDARY_FILTER,Description=\"One of the breakpoints is incompatible with mate locations.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_boundary_flt_tag, &len));

	const char* depth_flt_tag = "##FILTER=<ID=DEPTH_FILTER,Description=\"Depth of the region is incompatible with type of the event. Only applicable to long events.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, depth_flt_tag, &len));

	const char* anom_flanking_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_FLANKING_DEPTH,Description=\"Depth of region(s) flanking this event is anomalous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anom_flanking_depth_flt_tag, &len));

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
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the indel (DEL or DUP).\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, end_tag, &len));

	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svlen_tag, &len));

	const char* max_size_tag = "##INFO=<ID=MAX_SIZE,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_size_tag, &len));

	const char* ks_pval_tag = "##INFO=<ID=KS_PVAL,Number=1,Type=Float,Description=\"p-value of the KS test. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ks_pval_tag, &len));

	const char* est_size_tag = "##INFO=<ID=EST_SIZE,Number=1,Type=Integer,Description=\"Estimated size of the imprecise event. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, est_size_tag, &len));

	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_size_tag, &len));
	const char* remap_lb_tag = "##INFO=<ID=REMAP_LB,Number=1,Type=Integer,Description=\"Minimum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_lb_tag, &len));

	const char* remap_ub_tag = "##INFO=<ID=REMAP_UB,Number=1,Type=Integer,Description=\"Maximum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_ub_tag, &len));

//	const char* depths_tag = "##INFO=<ID=DEPTHS,Number=4,Type=Integer,Description=\"Depths of, respectively, the region flanking the indel to the left,"
//			"the left portion of the indel, the right portion of the indel, the region flanking the indel to the right. Numbers 2 and 3 will be identical for short indels.\">";
//	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, depths_tag, &len));

	const char* median_depths_tag = "##INFO=<ID=MEDIAN_DEPTHS,Number=4,Type=Integer,Description=\"Depths of, respectively, the region flanking the indel to the left,"
			"the left portion of the indel, the right portion of the indel, the region flanking the indel to the right. Numbers 2 and 3 will be identical for short indels.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, median_depths_tag, &len));

	const char* cluster_depths_tag = "##INFO=<ID=CLUSTER_DEPTHS,Number=2,Type=Integer,Description=\"Depths of the left and right cluster regions.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cluster_depths_tag, &len));

	const char* disc_pairs_tag = "##INFO=<ID=DISC_PAIRS,Number=1,Type=Integer,Description=\"Discordant pairs supporting the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_tag, &len));

	const char* disc_pairs_hmapq_tag = "##INFO=<ID=DISC_PAIRS_HIGHMAPQ,Number=1,Type=Integer,Description=\"HDiscordant pairs with high MAPQ supporting the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_hmapq_tag, &len));

	const char* disc_pairs_surr_tag = "##INFO=<ID=DISC_PAIRS_SURROUNDING,Number=2,Type=Integer,Description=\"Discordant pairs around the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_surr_tag, &len));

	const char* conc_pairs_tag = "##INFO=<ID=CONC_PAIRS,Number=1,Type=Integer,Description=\"Concordant pairs supporting the absence of a SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, conc_pairs_tag, &len));

	const char* clipped_reads_tag = "##INFO=<ID=CLIPPED_READS,Number=2,Type=Integer,Description=\"Reads supporting the right and the left breakpoints, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, clipped_reads_tag, &len));

	const char* max_mapq_tag = "##INFO=<ID=MAX_MAPQ,Number=2,Type=Integer,Description=\"Maximum MAPQ of clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_mapq_tag, &len));

	const char* dp_max_mapq_tag = "##INFO=<ID=DISC_PAIRS_MAXMAPQ,Number=1,Type=Integer,Description=\"Maximum MAPQ of supporting discordant pairs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dp_max_mapq_tag, &len));

	const char* ext_1sr_reads_tag = "##INFO=<ID=EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads extending a 1SR consensus to the left and to the right.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ext_1sr_reads_tag, &len));

	const char* hq_ext_1sr_reads_tag = "##INFO=<ID=HQ_EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads with high MAPQ extending a 1SR consensus to the left and to the right.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, hq_ext_1sr_reads_tag, &len));

	const char* fullj_score_tag = "##INFO=<ID=FULL_JUNCTION_SCORE,Number=1,Type=Integer,Description=\"Full junction score.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fullj_score_tag, &len));

	const char* splitj_score_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE,Number=2,Type=Integer,Description=\"Score of the best alignment of the left-half and right-half of the junction.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score_tag, &len));

	const char* splitj_score2_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE2,Number=2,Type=Integer,Description=\"Score of the second best alignment of the left-half and right-half of the junction.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score2_tag, &len));

	const char* splitj_size_tag = "##INFO=<ID=SPLIT_JUNCTION_SIZE,Number=2,Type=Integer,Description=\"Size of the the left-half and right-half of the junction.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_size_tag, &len));

	const char* source_tag = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source algorithm of the indel.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, source_tag, &len));

	const char* svinsseq_tag = "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinsseq_tag, &len));

	const char* imprecise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The reported boundaries are not precise.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, imprecise_tag, &len));

	const char* og_range_tag = "##INFO=<ID=ORIGINAL_RANGE,Number=1,Type=String,Description=\"Unadjusted imprecise range predicted by discordant pairs. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, og_range_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, gt_tag, &len));

	const char* ft_tag = "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS indicates a reliable call. Any other value means the call is not reliable.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ft_tag, &len));

	// add ALT
	const char* del_alt_tag = "##ALT=<ID=DEL,Description=\"Deletion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, del_alt_tag, &len));

	const char* dup_alt_tag = "##ALT=<ID=DUP,Description=\"Tandem Duplication\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dup_alt_tag, &len));

	std::string cmd_tag = "##SurVeyorCommand=" + command;
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cmd_tag.c_str(), &len));

	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::string version_tag = "##SurVeyor=" + config.version + "; Date=" + std::ctime(&now_time);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, version_tag.c_str(), &len));

	std::stringstream called_by_ss;
	called_by_ss << "##calledBy=SurVeyor " << config.version << "; ";
	called_by_ss << "seed: " << config.seed << "; ";
	called_by_ss << "min_sv_size: " << config.min_sv_size << "; ";
	called_by_ss << "min_clip_len: " << config.min_clip_len << "; ";
	called_by_ss << "max_seq_error: " << config.max_seq_error << "; ";
	called_by_ss << "max_clipped_pos_dist: " << config.max_clipped_pos_dist << "; ";
	called_by_ss << "min_size_for_depth_filtering: " << config.min_size_for_depth_filtering << "; ";
	called_by_ss << "sampling-regions: " << (config.sampling_regions.empty() ? "no" : config.sampling_regions) << "; ";
	std::string called_by = called_by_ss.str();
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, called_by.c_str(), &len));

	// add samples
	bcf_hdr_add_sample(header, sample_name.c_str());

	return header;
}

void indel2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, std::string& contig_name, char* chr_seq, sv_t* sv, std::vector<std::string>& filters) {
	bcf_entry->rid = bcf_hdr_name2id(hdr, contig_name.c_str());
	bcf_entry->pos = sv->start;
	bcf_update_id(hdr, bcf_entry, sv->id.c_str());
	
	std::string alleles = std::string(1, chr_seq[sv->start]) + ",<" + sv->svtype() + ">";
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	// add filters
	for (std::string& filter : filters) {
		int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter.c_str());
		bcf_add_filter(hdr, bcf_entry, filter_id);
	}

	// add info
	int int_conv; // using hts_pos_t + bcf_update_info_int64 throws [E::bcf_update_info] The type 257 not implemented yet
	
	int_conv = sv->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	int_conv = sv->svlen();
	bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", sv->svtype().c_str());
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", sv->source.c_str());
	if (!sv->ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", sv->ins_seq.c_str());
	}

	int split_junction_score[] = {sv->left_anchor_aln->best_score, sv->right_anchor_aln->best_score};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE", split_junction_score, 2);
	int split_junction_score2[] = {sv->left_anchor_aln->next_best_score, sv->right_anchor_aln->next_best_score};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE2", split_junction_score2, 2);
	int split_junction_size[] = {sv->left_anchor_aln->seq_len, sv->right_anchor_aln->seq_len};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE", split_junction_size, 2);
	if (sv->full_junction_aln != NULL) {
		bcf_update_info_int32(hdr, bcf_entry, "FULL_JUNCTION_SCORE", &sv->full_junction_aln->best_score, 1);
	}

	// add GT info
	int gt[1];
	gt[0] = bcf_gt_unphased(1);
	bcf_update_genotypes(hdr, bcf_entry, gt, 1);

	const char* ft_val = (filters[0] == "PASS") ? "PASS" : "FAIL";
	bcf_update_format_string(hdr, bcf_entry, "FT", &ft_val, 1);
}

void del2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, std::string& contig_name, sv2_deletion_t* del, std::vector<std::string>& filters) {
	bcf_clear(bcf_entry);

	indel2bcf(hdr, bcf_entry, contig_name, chr_seq, del->sv, filters);
	
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
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &del->disc_pairs, 1);
	if (del->disc_pairs > 0) {
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_MAXMAPQ", &del->disc_pairs_maxmapq, 1);
	}
	if (del->sv->source == "DP") {
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_HIGHMAPQ", &del->disc_pairs_high_mapq, 1);
	}
	int disc_pairs_surr[] = {del->l_cluster_region_disc_pairs, del->r_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	bcf_update_info_int32(hdr, bcf_entry, "CONC_PAIRS", &del->conc_pairs, 1);
	int clipped_reads[] = {del->rc_consensus ? (int) del->rc_consensus->supp_clipped_reads() : 0, del->lc_consensus ? (int) del->lc_consensus->supp_clipped_reads() : 0};
	bcf_update_info_int32(hdr, bcf_entry, "CLIPPED_READS", clipped_reads, 2);
	int max_mapq[] = {del->rc_consensus ? (int) del->rc_consensus->max_mapq : 0, del->lc_consensus ? (int) del->lc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	if (del->is_single_consensus()) {
		if (del->lc_consensus) {
			int ext_1sr_reads[] = { del->lc_consensus->left_ext_reads, del->lc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { del->lc_consensus->hq_left_ext_reads, del->lc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
		} else if (del->rc_consensus) {
			int ext_1sr_reads[] = { del->rc_consensus->left_ext_reads, del->rc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { del->rc_consensus->hq_left_ext_reads, del->rc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
		}
	}

	bcf_update_info_flag(hdr, bcf_entry, "IMPRECISE", "", del->imprecise());
}

void dup2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, std::string& contig_name, sv2_duplication_t* dup, std::vector<std::string>& filters) {
	bcf_clear(bcf_entry);
	
	indel2bcf(hdr, bcf_entry, contig_name, chr_seq, dup->sv, filters);

	// add info
	int int_conv;

	int median_depths[] = {dup->med_left_flanking_cov, dup->med_indel_left_cov, dup->med_indel_right_cov, dup->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &dup->disc_pairs, 1);
	int disc_pairs_surr[] = {dup->rc_cluster_region_disc_pairs, dup->lc_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	int clipped_reads[] = {dup->lc_consensus ? (int) dup->lc_consensus->supp_clipped_reads() : 0, dup->rc_consensus ? (int) dup->rc_consensus->supp_clipped_reads() : 0};
	bcf_update_info_int32(hdr, bcf_entry, "CLIPPED_READS", clipped_reads, 2);
	int max_mapq[] = {dup->lc_consensus ? (int) dup->lc_consensus->max_mapq : 0, dup->rc_consensus ? (int) dup->rc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	if (dup->is_single_consensus()) {
		if (dup->lc_consensus) {
			int ext_1sr_reads[] = { dup->lc_consensus->left_ext_reads, dup->lc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { dup->lc_consensus->hq_left_ext_reads, dup->lc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
		} else if (dup->rc_consensus) {
			int ext_1sr_reads[] = { dup->rc_consensus->left_ext_reads, dup->rc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { dup->rc_consensus->hq_left_ext_reads, dup->rc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
		}
	}
}

void remove_marked_consensuses(std::vector<consensus_t*>& consensuses, std::vector<bool>& used) {
	for (int i = 0; i < consensuses.size(); i++) if (used[i]) consensuses[i] = NULL;
	consensuses.erase(std::remove(consensuses.begin(), consensuses.end(), (consensus_t*) NULL), consensuses.end());
}

bool file_exists(std::string& fname) {
	return std::ifstream(fname).good();
}

#endif //SURVINDEL2_UTILS_H
