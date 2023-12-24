#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>

#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "utils.h"
#include "vcf_utils.h"

config_t config;
stats_t stats;
double min_ptn_ratio;

contig_map_t contig_map;

std::vector<uint32_t> min_disc_pairs_by_size;

std::pair<int, int> support(inss_insertion_t* insertion) {
	return {insertion->r_disc_pairs+insertion->rc_reads(), insertion->l_disc_pairs+insertion->lc_reads()};
}

std::pair<double, double> ptn_score(inss_insertion_t* insertion) {
	int r_positive = insertion->r_disc_pairs + insertion->rc_reads();
	int r_negative = insertion->r_conc_pairs;
	int l_positive = insertion->l_disc_pairs + insertion->lc_reads();
	int l_negative = insertion->l_conc_pairs;
	return {double(r_positive)/(r_positive+r_negative), double(l_positive)/(l_positive+l_negative)};
}

// filters for all categories
void add_AST_filters(inss_insertion_t* insertion, std::vector<std::string>& filters) {
	if (insertion->ins->ins_seq.length()-(insertion->ins->end-insertion->ins->start) < config.min_sv_size) filters.push_back("SMALL");
	if (insertion->median_lf_cov > stats.get_max_depth(insertion->ins->chr) || insertion->median_rf_cov > stats.get_max_depth(insertion->ins->chr) ||
		insertion->median_lf_cov < stats.get_min_depth(insertion->ins->chr) || insertion->median_rf_cov < stats.get_min_depth(insertion->ins->chr))
		filters.push_back("ANOMALOUS_FLANKING_DEPTH");
	int d = insertion->ins->ins_seq.find("-");
	if (is_homopolymer(insertion->ins->ins_seq.substr(0, d))) filters.push_back("HOMOPOLYMER_INSSEQ");
	else if (d != std::string::npos && is_homopolymer(insertion->ins->ins_seq.substr(d+1))) filters.push_back("HOMOPOLYMER_INSSEQ");
	if (insertion->rc_reads() > stats.get_max_depth(insertion->ins->chr) || insertion->lc_reads() > stats.get_max_depth(insertion->ins->chr)) filters.push_back("ANOMALOUS_SC_NUMBER");
}

void add_AT_filters(inss_insertion_t* insertion, std::vector<std::string>& filters) {
	auto supports = support(insertion);
	if (insertion->l_disc_pairs + insertion->r_disc_pairs == 0) filters.push_back("NO_DISC_SUPPORT");
	else if (supports.first < stats.get_median_depth(insertion->ins->chr)/5 && supports.second < stats.get_median_depth(insertion->ins->chr)/5)
		filters.push_back("LOW_SUPPORT");
}

std::vector<std::string> get_transurveyor_insertions_filterlist(inss_insertion_t* insertion) {
	std::vector<std::string> filters;
	add_AST_filters(insertion, filters);
	add_AT_filters(insertion, filters);
	auto scores = ptn_score(insertion);
	if (scores.first < min_ptn_ratio && scores.second < min_ptn_ratio) filters.push_back("LOW_SCORE");
	int svlen = insertion->ins->ins_seq.length();
	if (svlen >= min_disc_pairs_by_size.size()) svlen = min_disc_pairs_by_size.size()-1;
	if (insertion->r_disc_pairs + insertion->l_disc_pairs < min_disc_pairs_by_size[svlen]) filters.push_back("NOT_ENOUGH_DISC_PAIRS");
	if (filters.empty()) filters.push_back("PASS");
	return filters;
}

std::vector<std::string> get_assembled_insertions_filterlist(inss_insertion_t* insertion) {
	std::vector<std::string> filters;
	add_AST_filters(insertion, filters);
	add_AT_filters(insertion, filters);
	size_t divide = insertion->ins->ins_seq.find("-");
	bool incomplete_assembly = (divide != std::string::npos);
	int svlen = insertion->ins->ins_seq.length();
	if (incomplete_assembly || svlen >= min_disc_pairs_by_size.size()) svlen = min_disc_pairs_by_size.size()-1;
	if (insertion->r_disc_pairs + insertion->l_disc_pairs < min_disc_pairs_by_size[svlen]) filters.push_back("NOT_ENOUGH_DISC_PAIRS");
	if (incomplete_assembly && (divide < stats.read_len && insertion->ins->ins_seq.length()-divide < stats.read_len)) filters.push_back("LOW_QUAL_ASSEMBLY");
	if (filters.empty()) filters.push_back("PASS");
	return filters;
}


int main(int argc, char* argv[]) {
    std::string workdir = argv[1];
    std::string reference_fname = argv[2];
    min_ptn_ratio = std::stod(argv[3]);

    config.parse(workdir + "/config.txt");
	stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contig_map.load(workdir);

	chr_seqs_map_t contigs;
	contigs.read_fasta_into_map(reference_fname);

	std::ifstream mdpbs_fin(workdir + "/min_disc_pairs_by_size.txt");
	int i, min_disc_pairs;
	while (mdpbs_fin >> i >> min_disc_pairs) min_disc_pairs_by_size.push_back(min_disc_pairs);
	mdpbs_fin.close();

	std::vector<bcf1_t*> final_insertions_set;

	std::unordered_map<std::string, std::vector<inss_insertion_t*> > final_insertions_by_contig;

	/* == Insertions predicted by TranSurVeyor == */
	std::string transurveyor_ins_vcf_fname = workdir + "/transurveyor_ins.annotated.vcf.gz";
	htsFile* transurveyor_ins_vcf_file = bcf_open(transurveyor_ins_vcf_fname.c_str(), "r");
	if (!transurveyor_ins_vcf_file) {
		throw std::runtime_error("Unable to open file " + transurveyor_ins_vcf_fname + ".");
	}
	bcf_hdr_t* transurveyor_ins_hdr = bcf_hdr_read(transurveyor_ins_vcf_file);
	bcf_hdr_t* out_vcf_header = bcf_hdr_dup(transurveyor_ins_hdr);
	
	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(transurveyor_ins_vcf_file, transurveyor_ins_hdr, bcf_entry) == 0) {
		std::string contig_name = bcf_seqname_safe(transurveyor_ins_hdr, bcf_entry);

		std::string ins_seq = get_ins_seq(bcf_entry, transurveyor_ins_hdr);

		std::vector<inss_insertion_t*>& dst_contig_insertions = final_insertions_by_contig[contig_name];
		inss_insertion_t* insertion = new inss_insertion_t(contig_name, bcf_entry->pos, get_sv_end(bcf_entry, transurveyor_ins_hdr), 0, 0, 0, 0, 0, 0, 0, ins_seq);

		int* stable_depths = NULL;
		int size = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "STABLE_DEPTHS", &stable_depths, &size);

		int* discordants = NULL;
		size = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "DISCORDANT", &discordants, &size);

		int* split_reads = NULL;
		size = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "SPLIT_READS", &split_reads, &size);

		int* fwd_split_reads = NULL;
		size = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "FWD_SPLIT_READS", &fwd_split_reads, &size);

		int* rev_split_reads = NULL;
		size = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "REV_SPLIT_READS", &rev_split_reads, &size);

		int* spanning_reads = NULL;
		size = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "SPANNING_READS", &spanning_reads, &size);

		insertion->median_lf_cov = stable_depths[0], insertion->median_rf_cov = stable_depths[1];
		insertion->r_disc_pairs = discordants[0], insertion->l_disc_pairs = discordants[1];
		insertion->rc_fwd_reads = fwd_split_reads[0], insertion->lc_fwd_reads = fwd_split_reads[1];
		insertion->rc_rev_reads = rev_split_reads[0], insertion->lc_rev_reads = rev_split_reads[1];
		insertion->r_conc_pairs = spanning_reads[0], insertion->l_conc_pairs = spanning_reads[1];

		std::string mh_seq;
		int mh_len = 0;
		// extract microhomology, if present
		if (insertion->ins->start > insertion->ins->end) {
			mh_len = insertion->ins->start - insertion->ins->end;
			char* mh_seq_cstr = new char[mh_len+1];
			strncpy(mh_seq_cstr, contigs.get_seq(contig_name)+insertion->ins->end, mh_len);
			mh_seq_cstr[mh_len] = '\0';
			for (int i = 0; i < mh_len; i++) {
				mh_seq_cstr[i] = toupper(mh_seq_cstr[i]);
			}
			mh_seq = mh_seq_cstr;
			delete[] mh_seq_cstr;
			insertion->ins->start = insertion->ins->end;
			insertion->mh_len = mh_len;
			bcf_entry->pos = get_sv_end(bcf_entry, transurveyor_ins_hdr);
		}
		std::string ins_seq_w_mh = mh_seq + get_ins_seq(bcf_entry, transurveyor_ins_hdr);
		insertion->ins->ins_seq = ins_seq_w_mh;
		bcf_update_info_string(out_vcf_header, bcf_entry, "SVINSSEQ", ins_seq_w_mh.c_str());
		int int_conv = ins_seq_w_mh.length();
		bcf_update_info_int32(out_vcf_header, bcf_entry, "SVINSLEN", &int_conv, 1);
		int_conv = insertion->ins->ins_seq.length() - (insertion->ins->end-insertion->ins->start);
		bcf_update_info_int32(out_vcf_header, bcf_entry, "SVLEN", &int_conv, 1);
		bcf_update_info_int32(out_vcf_header, bcf_entry, "MHLEN", &insertion->mh_len, 1);

		auto scores = ptn_score(insertion);
		float float2_conv[2];
		float2_conv[0] = scores.first, float2_conv[1] = scores.second;
		bcf_update_info_float(out_vcf_header, bcf_entry, "SCORES", float2_conv, 2);

		std::vector<std::string> filters = get_transurveyor_insertions_filterlist(insertion);
		for (std::string filter : filters) {
			int filter_id = bcf_hdr_id2int(out_vcf_header, BCF_DT_ID, filter.c_str());
			bcf_add_filter(out_vcf_header, bcf_entry, filter_id);
		}

		if (bcf_has_filter(out_vcf_header, bcf_entry, (char*) "PASS")) {
			final_insertions_by_contig[contig_name].push_back(insertion);
		}
		bcf_unpack(bcf_entry, BCF_UN_ALL);
		bcf_translate(out_vcf_header, transurveyor_ins_hdr, bcf_entry);
		final_insertions_set.push_back(bcf_dup(bcf_entry));
	}

	/* == Assembled insertions == */
	std::string assembled_ins_vcf_fname = workdir + "/assembled_ins.annotated.vcf.gz";
	htsFile* assembled_ins_vcf_file = bcf_open(assembled_ins_vcf_fname.c_str(), "r");
	if (!assembled_ins_vcf_file) {
		throw std::runtime_error("Unable to open file " + assembled_ins_vcf_fname + ".");
	}
	bcf_hdr_t* assembled_ins_hdr = bcf_hdr_read(assembled_ins_vcf_file);
	while (bcf_read(assembled_ins_vcf_file, assembled_ins_hdr, bcf_entry) == 0) {
		std::string contig_name = bcf_seqname_safe(assembled_ins_hdr, bcf_entry);

		std::string ins_seq = get_ins_seq(bcf_entry, assembled_ins_hdr);

		std::vector<inss_insertion_t*>& dst_contig_insertions = final_insertions_by_contig[contig_name];
		inss_insertion_t* insertion = new inss_insertion_t(contig_name, bcf_entry->pos, get_sv_end(bcf_entry, assembled_ins_hdr), 0, 0, 0, 0, 0, 0, 0, ins_seq);

		int* stable_depths = NULL;
		int size = 0;
		bcf_get_info_int32(assembled_ins_hdr, bcf_entry, "STABLE_DEPTHS", &stable_depths, &size);

		int* discordants = NULL;
		size = 0;
		bcf_get_info_int32(assembled_ins_hdr, bcf_entry, "DISCORDANT", &discordants, &size);

		insertion->median_lf_cov = stable_depths[0], insertion->median_rf_cov = stable_depths[1];
		insertion->r_disc_pairs = discordants[0], insertion->l_disc_pairs = discordants[1];

		std::vector<std::string> filters = get_assembled_insertions_filterlist(insertion);
		for (std::string filter : filters) {
			int filter_id = bcf_hdr_id2int(out_vcf_header, BCF_DT_ID, filter.c_str());
			bcf_add_filter(out_vcf_header, bcf_entry, filter_id);
		}

		if (bcf_has_filter(out_vcf_header, bcf_entry, (char*) "PASS")) {
			final_insertions_by_contig[contig_name].push_back(insertion);
		}
		bcf_unpack(bcf_entry, BCF_UN_ALL);
		bcf_translate(out_vcf_header, assembled_ins_hdr, bcf_entry);
		final_insertions_set.push_back(bcf_dup(bcf_entry));
	}

	// sort and write final set of insertions to file
	for (bcf1_t* b : final_insertions_set) bcf_unpack(b, BCF_UN_ALL);
	std::sort(final_insertions_set.begin(), final_insertions_set.end(), [](const bcf1_t* b1, const bcf1_t* b2) {
		std::string id1 = std::string(b1->d.id), id2 = std::string(b2->d.id);
		return std::tie(b1->rid, b1->pos, id1) < std::tie(b2->rid, b2->pos, id2);
	});

	std::string out_vcf_fname = workdir + "/out.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + out_vcf_fname + ".");
	}
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	std::string out_pass_vcf_fname = workdir + "/out.pass.vcf.gz";
	htsFile* out_pass_vcf_file = bcf_open(out_pass_vcf_fname.c_str(), "wz");
	if (out_pass_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + out_pass_vcf_fname + ".");
	}
	if (bcf_hdr_write(out_pass_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_pass_vcf_fname + ".");
	}

	for (bcf1_t* bcf_entry : final_insertions_set) {
		if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
			throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
		}
		if (bcf_has_filter(out_vcf_header, bcf_entry, (char*) "PASS")) {
			if (bcf_write(out_pass_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_pass_vcf_fname + ".");
			}
		}
	}
	bcf_close(out_vcf_file);
	bcf_close(out_pass_vcf_file);

	tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
	tbx_index_build(out_pass_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
