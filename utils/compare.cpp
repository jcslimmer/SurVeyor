#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <mutex>

#include "../libs/cxxopts.h"
#include "../libs/IntervalTree.h"
#include "../libs/ssw_cpp.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw.h"
#include "common.h"
#include "htslib/vcf.h"

chr_seqs_map_t chr_seqs;
int max_prec_dist, max_imprec_dist, max_repeat_dist;
double min_prec_frac_overlap, min_imprec_frac_overlap;
int max_prec_len_diff, max_imprec_len_diff;
std::unordered_set<std::string> bdup_ids, cdup_ids;

// StripedSmithWaterman::Aligner aligner(1,4,6,1,false);
std::vector<sv_t*> benchmark_svs;
std::unordered_map<std::string, std::vector<sv_t*> > called_dels_by_chr, called_inss_by_chr, called_invs_by_chr;
std::unordered_map<std::string, IntervalTree<repeat_t>*> reps_i;

bool ignore_seq = false;

bool compatible_gts(sv_t* sv1, sv_t* sv2) {
	int min_sv1_ac = sv1->allele_count(1), min_sv2_ac = sv2->allele_count(1);
	int max_sv1_ac = min_sv1_ac + sv1->missing_alleles(), max_sv2_ac = min_sv2_ac + sv2->missing_alleles();
	bool ac_iv_overlap = !(max_sv1_ac < min_sv2_ac || max_sv2_ac < min_sv1_ac);
	return ac_iv_overlap;
}

std::string get_bsv_type(sv_t* sv) {
	if (bdup_ids.count(sv->id)) return "DUP";
	else return sv->svtype();
}
std::string get_csv_type(sv_t* sv) {
	if (cdup_ids.count(sv->id)) return "DUP";
	else return sv->svtype();
}

int len_diff(sv_t* sv1, sv_t* sv2) {
	if (sv1->svtype() == "DEL" && sv2->svtype() == "DEL") return abs(sv1->svlen()-sv2->svlen());
	else if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") return abs(sv1->svlen()-sv2->svlen());
	else if (sv1->svtype() == "INS" && sv2->svtype() == "INS") return abs((int) (sv1->svlen()-sv2->svlen()));
	else if ((sv1->svtype() == "INS" && sv2->svtype() == "DUP") || (sv1->svtype() == "DUP" && sv2->svtype() == "INS")) return 0;
	else if (sv1->svtype() == "INV" && sv2->svtype() == "INV") return abs(sv1->svlen()-sv2->svlen());
	else return INT32_MAX;
}

bool is_compatible_del_del(sv_t* sv1, sv_t* sv2) {
	if (sv1->imprecise || sv2->imprecise) {
		return  distance(sv1, sv2) <= max_imprec_dist &&
				overlap(sv1, sv2) >= min_imprec_frac_overlap &&
				len_diff(sv1, sv2) <= max_imprec_len_diff;
	} else {
		return  distance(sv1, sv2) <= max_prec_dist &&
				overlap(sv1, sv2) >= min_prec_frac_overlap &&
				len_diff(sv1, sv2) <= max_prec_len_diff;
	}
}
bool is_compatible_dup_dup(sv_t* sv1, sv_t* sv2) {
	return is_compatible_del_del(sv1, sv2);
}

bool check_ins_dup_seq(sv_t* ins_sv, sv_t* dup_sv, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	if (ignore_seq) return true;

	int max_len_diff = (ins_sv->imprecise || dup_sv->imprecise) ? max_imprec_len_diff : max_prec_len_diff;
	if (dup_sv->svlen() > ins_sv->ins_seq.length()+max_len_diff) return false;

	char* dup_seq = new char[dup_sv->end-dup_sv->start+1];
	strncpy(dup_seq, chr_seqs.get_seq(dup_sv->chr)+dup_sv->start, dup_sv->end-dup_sv->start);
	dup_seq[dup_sv->end-dup_sv->start] = '\0';
	std::string ext_dup_seq = dup_seq;
	while (ext_dup_seq.length() <= ins_sv->ins_seq.length() && ext_dup_seq.length() < UINT16_MAX) {
		ext_dup_seq += dup_sv->ins_seq + dup_seq;
	}
	delete[] dup_seq;

	if (ins_sv->start < dup_sv->start) {
		size_t rotation = dup_sv->start - ins_sv->start;
		rotation %= dup_sv->end-dup_sv->start;
		ext_dup_seq = ext_dup_seq.substr(rotation) + ext_dup_seq.substr(0, rotation);
	} else if (ins_sv->start > dup_sv->start) {
		size_t rotation = ins_sv->start - dup_sv->start;
		rotation %= dup_sv->end-dup_sv->start;
		ext_dup_seq = ext_dup_seq.substr(dup_sv->end-dup_sv->start-rotation) + ext_dup_seq.substr(0, dup_sv->end-dup_sv->start-rotation);
	}

	// ssw score is a uint16_t, so we cannot compare strings longer than that
	StripedSmithWaterman::Filter filter;
	if (ins_sv->ins_seq.length() > UINT16_MAX || ext_dup_seq.length() > UINT16_MAX) {
		std::string ins_seq_prefix = ins_sv->ins_seq.substr(0, UINT16_MAX-10);
		std::string ins_seq_suffix = ins_sv->ins_seq.substr(std::max(0, (int) ins_sv->ins_seq.length()-(UINT16_MAX-10)));
		std::string ext_dup_seq_prefix = ext_dup_seq.substr(0, UINT16_MAX-10);
		std::string ext_dup_seq_suffix = ext_dup_seq.substr(std::max(0, (int) ext_dup_seq.length()-(UINT16_MAX-10)));
		aligner.Align(ins_seq_prefix.data(), ext_dup_seq_prefix.data(), ext_dup_seq_prefix.length(), filter, &alignment, 0);
		if (alignment.query_end-alignment.query_begin < ins_seq_prefix.length()*0.8) return false;
		aligner.Align(ins_seq_suffix.data(), ext_dup_seq_suffix.data(), ext_dup_seq_suffix.length(), filter, &alignment, 0);
		if (alignment.query_end-alignment.query_begin < ins_seq_suffix.length()*0.8) return false;
		return true;
	}

	aligner.Align(ins_sv->ins_seq.data(), ext_dup_seq.data(), ext_dup_seq.length(), filter, &alignment, 0);
	return alignment.query_end-alignment.query_begin >= ins_sv->ins_seq.length()*0.8;
}
bool is_compatible_ins_dup(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Alignment alignment;
	if (sv1->imprecise || sv2->imprecise) {
		return distance(sv1, sv2) <= max_imprec_dist && check_ins_dup_seq(sv1, sv2, aligner, alignment);
	} else {
		return distance(sv1, sv2) <= max_prec_dist && check_ins_dup_seq(sv1, sv2, aligner, alignment);
	}
}

bool check_ins_ins_seq(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	if (ignore_seq) return true;
	if (sv1->ins_seq.length() > UINT16_MAX || sv2->ins_seq.length() > UINT16_MAX) return true; // TODO: ssw score is a uint16_t, so we cannot compare strings longer than that

	if (sv1->incomplete_ins_seq() || sv2->incomplete_ins_seq()) return true; // do not compare incomplete assemblies

	int max_len_diff = (sv1->imprecise || sv2->imprecise) ? max_imprec_len_diff : max_prec_len_diff;
	if (abs((int) (sv1->ins_seq.length()-sv2->ins_seq.length())) > max_len_diff) return false;

	sv_t* l_sv = sv1->start < sv2->start ? sv1 : sv2;
	sv_t* r_sv = sv1->start < sv2->start ? sv2 : sv1;
	char* extra_seq = new char[r_sv->start-l_sv->start+1];
	strncpy(extra_seq, chr_seqs.get_seq(sv1->chr)+l_sv->start, r_sv->start-l_sv->start);
	extra_seq[r_sv->start-l_sv->start] = '\0';

	std::string lsv_seq = l_sv->ins_seq + extra_seq;
	std::string rsv_seq = extra_seq + r_sv->ins_seq;
	delete[] extra_seq;

	std::string& query = lsv_seq.length() < rsv_seq.length() ? lsv_seq : rsv_seq;
	std::string& ref = lsv_seq.length() < rsv_seq.length() ? rsv_seq : lsv_seq;
	StripedSmithWaterman::Filter filter;
	aligner.Align(query.data(), ref.data(), ref.length(), filter, &alignment, 0);
	return alignment.query_end-alignment.query_begin >= query.length()*0.8;
}
bool is_compatible_ins_ins(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Alignment alignment;
	if (sv1->imprecise || sv2->imprecise) {
		return distance(sv1, sv2) <= max_imprec_dist && check_ins_ins_seq(sv1, sv2, aligner, alignment);
	} else {
		return distance(sv1, sv2) <= max_prec_dist && check_ins_ins_seq(sv1, sv2, aligner, alignment);
	}
}

bool check_ins_seq(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Alignment alignment;
	if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") {
		return true;
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "INS") {
		return check_ins_dup_seq(sv2, sv1, aligner, alignment);
	} else if (sv1->svtype() == "INS" && sv2->svtype() == "DUP") {
		return check_ins_dup_seq(sv1, sv2, aligner, alignment);
	} else {
		return check_ins_ins_seq(sv1, sv2, aligner, alignment);
	}
}

bool is_compatible_inv_inv(sv_t* sv1, sv_t* sv2) {
	return is_compatible_del_del(sv1, sv2);
}

bool is_compatible(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner) {
	if (sv1->svtype() == "DEL" && sv2->svtype() == "DEL") {
		return is_compatible_del_del(sv1, sv2);
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") {
		return is_compatible_dup_dup(sv1, sv2);
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "INS") {
		return is_compatible_ins_dup(sv2, sv1, aligner);
	} else if (sv1->svtype() == "INS" && sv2->svtype() == "DUP") {
		return is_compatible_ins_dup(sv1, sv2, aligner);
	} else if (sv1->svtype() == "INS" && sv2->svtype() == "INS") {
		return is_compatible_ins_ins(sv1, sv2, aligner);
	} else if (sv1->svtype() == "INV" && sv2->svtype() == "INV") {
		return is_compatible_inv_inv(sv1, sv2);
	} else {
		return false;
	}
}

struct sv_match_t {
	sv_t* b_sv, * c_sv;
	int score = 0;
	bool rep;

	sv_match_t(sv_t* b_sv, sv_t* c_sv, bool rep, StripedSmithWaterman::Aligner& aligner) : b_sv(b_sv), c_sv(c_sv) {
		if (b_sv == NULL || c_sv == NULL) {
			this->score = 0;
			this->rep = rep;
			return;
		}

		StripedSmithWaterman::Alignment alignment;
		alignment.Clear();
		int len_diff = 0;
		if (b_sv->svtype() == "INS" && c_sv->svtype() == "INS") {
			len_diff = abs(b_sv->svlen()-c_sv->svlen());
			check_ins_ins_seq(b_sv, c_sv, aligner, alignment);
		} else if (b_sv->svtype() == "INS" && c_sv->svtype() == "DUP") {
			check_ins_dup_seq(b_sv, c_sv, aligner, alignment);
		} else if (b_sv->svtype() == "DUP" && c_sv->svtype() == "INS") {
			check_ins_dup_seq(c_sv, b_sv, aligner, alignment);
		} else {
			len_diff = abs(b_sv->svlen()-c_sv->svlen());
		}

		int right_clip = 0;
		if (alignment.cigar.size() > 0) {
			uint32_t c = alignment.cigar[alignment.cigar.size()-1];
			if (cigar_int_to_op(c) == 'S') right_clip = cigar_int_to_len(c);
		}

		int dist_log = log(distance(b_sv, c_sv)+1);
		int aln_score = alignment.sw_score/std::max(1.0, double(alignment.query_end+right_clip)) * 100;
		this->score = -len_diff - dist_log + aln_score;
		this->rep = rep;
	}
};
std::vector<sv_match_t> matches;
std::mutex mtx;

void find_match(int id, int start_idx, int end_idx) {

	StripedSmithWaterman::Aligner aligner(1,4,6,1,false);

	for (int i = start_idx; i < end_idx; i++) {
		sv_t* bsv = benchmark_svs[i];

		std::vector<sv_t*>* called_svs_chr_type;
		if (bsv->svtype() == "DEL") {
			called_svs_chr_type = &called_dels_by_chr[bsv->chr];
		} else if (bsv->svtype() == "INS" || bsv->svtype() == "DUP") {
			called_svs_chr_type = &called_inss_by_chr[bsv->chr];
		} else if (bsv->svtype() == "INV") {
			called_svs_chr_type = &called_invs_by_chr[bsv->chr];
		} else {
			continue;
		}

		std::vector<repeat_t> reps_containing_bsv;
		if (reps_i[bsv->chr] != NULL) {
			std::vector<Interval<repeat_t>> intervals_temp = reps_i[bsv->chr]->findOverlapping(bsv->start, bsv->end);
			for (auto &iv : intervals_temp) {
				repeat_t rep = iv.value;
				if (rep.contains(bsv)) {
					reps_containing_bsv.push_back(rep);
				}
			}
		}

		for (sv_t* csv : *called_svs_chr_type) {
			if (is_compatible(bsv, csv, aligner)) {
				sv_match_t match(bsv, csv, false, aligner);
				mtx.lock();
				matches.push_back(match);
				mtx.unlock();
				continue;
			}

			if (distance(bsv, csv) > max_repeat_dist) continue;

			int max_len_diff = (bsv->imprecise || csv->imprecise) ? max_imprec_len_diff : max_prec_len_diff;
			if (len_diff(bsv, csv) <= max_len_diff) {
				bool same_tr = false;
				for (repeat_t& rep : reps_containing_bsv) {
					if (rep.intersects(csv)) {
						same_tr = true;
						break;
					}
				}

				if (same_tr && (bsv->svtype() == "DEL" || 
					((bsv->svtype() == "DUP" || bsv->svtype() == "INS") && check_ins_seq(bsv, csv, aligner)))) {
					sv_match_t match(bsv, csv, true, aligner);
					mtx.lock();
					matches.push_back(match);
					mtx.unlock();
				}
			}
		}
	}
}

int main(int argc, char* argv[]) {

	cxxopts::Options options("compare", "Given a VCF or SV file with benchmark deletions and one with the called ones, reports for "
											"each benchmark deletion whether it was called.");

	options.add_options()
		("benchmark_file", "VCF or SV file with the benchmark calls.", cxxopts::value<std::string>())
		("called_file", "VCF or SV file with the calls to benchmark.", cxxopts::value<std::string>())
		("R,reference", "Reference file in fasta format.", cxxopts::value<std::string>())
		("T,tandem-repeats", "Tandem repeat file in BED format.", cxxopts::value<std::string>())
		("d,max_dist_precise", "Maximum distance allowed between the breakpoints of two precise variants.",
				cxxopts::value<int>()->default_value("100"))
		("D,max_dist_imprecise", "Maximum distance allowed between the breakpoints when at least one of the variants is imprecise.",
				cxxopts::value<int>()->default_value("500"))
		("v,min_overlap_precise", "Minimum overlap required between two precise variants, expressed as the fraction of the shortest variant"
				" that is covered by the longest variant.", cxxopts::value<double>()->default_value("0.8"))
		("V,min_overlap_imprecise", "Minimum overlap required between two variants when at least one is imprecise, "
				"expressed as the fraction of the shortest variant that is covered by the longest variant.",
				cxxopts::value<double>()->default_value("0.5"))
		("s,max_len_diff_precise", "Maximum length difference allowed between two precise variants.",
				cxxopts::value<int>()->default_value("100"))
		("S,max_len_diff_imprecise", "Maximum length difference allowed when at least one variant is imprecise.",
				cxxopts::value<int>()->default_value("500"))
		("max-repeat-dist", "Maximum distance between two variant in the same tandem repeat.", cxxopts::value<int>()->default_value("1000"))
		("r,report", "Print report only", cxxopts::value<bool>()->default_value("false"))
		("f,fps", "Print false positive SVs to file.", cxxopts::value<std::string>())
		("I,ignore-seq", "Only compare insertions by position.", cxxopts::value<bool>()->default_value("false"))
		("i,force-ids", "Generate new IDs for the variants, required when variants do not have IDs or have duplicated IDs.",
				cxxopts::value<bool>()->default_value("false"))
		("a,all-imprecise", "Treat all deletions as imprecise.", cxxopts::value<bool>()->default_value("false"))
		("bdup-ids", "ID of benchmark SVs to be considered duplicatons. Note this only affects the final report, not how SVs are compared.", cxxopts::value<std::string>())
		("cdup-ids", "ID of called SVs to be considered duplicatons. Note this only affects the final report, not how SVs are compared.", cxxopts::value<std::string>())
		("wrong-gts", "Print pairs of matching SVs that have discordant genotypes.", cxxopts::value<std::string>())
		("keep-all-called", "Keep all variants in the called file, even if no alternative allele", cxxopts::value<bool>()->default_value("false"))
		("c,called-to-benchmark-gts", "For each called SV matching a benchmark SV, report their genotype according to the benchmark dataset.", cxxopts::value<std::string>())
		("e,exclusive", "SV cannot be used in multiple matches.", cxxopts::value<bool>()->default_value("false"))
		("t,threads", "Number of threads to use.", cxxopts::value<int>()->default_value("1"))
		("h,help", "Print usage");

	options.parse_positional({"benchmark_file", "called_file"});
	options.positional_help("benchmark_file called_file");
	options.show_positional_help();
	auto parsed_args = options.parse(argc, argv);

	if (parsed_args.count("help") || !parsed_args.count("benchmark_file") || !parsed_args.count("called_file")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}

    std::string benchmark_fname = parsed_args["benchmark_file"].as<std::string>();
    benchmark_svs = read_sv_list(benchmark_fname.c_str());
    std::string called_fname = parsed_args["called_file"].as<std::string>();
	std::vector<sv_t*> called_svs = read_sv_list(called_fname.c_str());
	max_prec_dist = parsed_args["max_dist_precise"].as<int>();
	max_imprec_dist = parsed_args["max_dist_imprecise"].as<int>();
	max_repeat_dist = parsed_args["max-repeat-dist"].as<int>();
	min_prec_frac_overlap = parsed_args["min_overlap_precise"].as<double>();
	min_imprec_frac_overlap = parsed_args["min_overlap_imprecise"].as<double>();
	max_prec_len_diff = parsed_args["max_len_diff_precise"].as<int>();
	max_imprec_len_diff = parsed_args["max_len_diff_imprecise"].as<int>();
    bool report = parsed_args["report"].as<bool>();
	bool exclusive = parsed_args["exclusive"].as<bool>();
	int threads = parsed_args["threads"].as<int>();
	ignore_seq = parsed_args["ignore-seq"].as<bool>();
    if (parsed_args["all-imprecise"].as<bool>()) {
    	max_prec_dist = max_imprec_dist;
    	min_prec_frac_overlap = min_imprec_frac_overlap;
    	max_prec_len_diff = max_imprec_len_diff;
    }

	if (parsed_args.count("bdup-ids")) {
		std::string line;
		std::ifstream dup_f(parsed_args["bdup-ids"].as<std::string>());
		while (getline(dup_f, line)) {
			bdup_ids.insert(line);
		}
	}
	if (parsed_args.count("cdup-ids")) {
		std::string line;
		std::ifstream dup_f(parsed_args["cdup-ids"].as<std::string>());
		while (getline(dup_f, line)) {
			cdup_ids.insert(line);
		}
	}

	if (!parsed_args.count("reference") && !ignore_seq) {
		std::cerr << "If no reference is provided, the --ignore-seq flag must be provided." << std::endl;
		exit(0);
	}

	htsFile* fps_fout = NULL;
	bcf_hdr_t* fps_hdr = NULL;
	if (parsed_args.count("fps")) {
		std::string fps_fname = parsed_args["fps"].as<std::string>();
		fps_fout = hts_open(fps_fname.c_str(), "wz");
		// read header from called file
		htsFile* called_file = bcf_open(called_fname.c_str(), "r");
		fps_hdr = bcf_hdr_read(called_file);
		fps_hdr = bcf_hdr_dup(fps_hdr);
		if (bcf_hdr_write(fps_fout, fps_hdr) != 0) {
			std::cerr << "Error writing header to file " << fps_fname << std::endl;
			exit(1);
		}
		bcf_close(called_file);
	}

	std::ofstream called_to_benchmark_gts_fout;
	if (parsed_args.count("called-to-benchmark-gts")) {
		std::string called_to_benchmark_gts_fname = parsed_args["called-to-benchmark-gts"].as<std::string>();
		called_to_benchmark_gts_fout.open(called_to_benchmark_gts_fname);
	}

	std::ofstream wrong_gts_fout;
	if (parsed_args.count("wrong-gts")) {
		std::string wrong_gts_fname = parsed_args["wrong-gts"].as<std::string>();
		wrong_gts_fout.open(wrong_gts_fname);
	}

	auto is_unsupported_func = [](sv_t* sv) {return sv->svtype() != "DEL" && sv->svtype() != "INS" && sv->svtype() != "DUP" && sv->svtype() != "INV";};
	// erase and count elements from benchmark_svs that are not supported
	int n_unsupported = std::count_if(benchmark_svs.begin(), benchmark_svs.end(), is_unsupported_func);
	if (n_unsupported > 0) {
		std::cerr << "Warning: only SVTYPE=DEL, DUP, INS and INV are supported. " << n_unsupported << " unsupported variants in benchmark file." << std::endl;
	}
	benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), is_unsupported_func), benchmark_svs.end());

	// erase and count elements from called_svs that are not supported
	n_unsupported = std::count_if(called_svs.begin(), called_svs.end(), is_unsupported_func);
	if (n_unsupported > 0) {
		std::cerr << "Warning: only SVTYPE=DEL, DUP, INS and INV are supported. " << n_unsupported << " unsupported variants in called file." << std::endl;
	}
	called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), is_unsupported_func), called_svs.end());

	int n_ac_0 = std::count_if(benchmark_svs.begin(), benchmark_svs.end(), [](sv_t* sv) {return sv->allele_count(1) == 0;});
	if (n_ac_0 > 0) {
		std::cerr << "Warning: excluded " << n_ac_0 << " variants in benchmark file that have no ALT alleles." << std::endl;
	}
	benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), [](sv_t* sv) {return sv->allele_count(1) == 0;}), benchmark_svs.end());

	if (parsed_args["keep-all-called"].as<bool>()) {
		std::cerr << "Warning: keeping all variants in the called file, even if they have no ALT alleles." << std::endl;
	} else {
		n_ac_0 = std::count_if(called_svs.begin(), called_svs.end(), [](sv_t* sv) {return sv->allele_count(1) == 0;});
		if (n_ac_0 > 0) {
			std::cerr << "Warning: excluded " << n_ac_0 << " variants in called file that have no ALT alleles." << std::endl;
		}
		called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), [](sv_t* sv) {return sv->allele_count(1) == 0;}), called_svs.end());
	}

	if (parsed_args["force-ids"].as<bool>()) {
		if (!bdup_ids.empty()) {
			std::cerr << "Warning: --force-ids and --bdup-ids are both used. IDs of benchmark SVs will not be forced." << std::endl;
		} else {
			for (int j = 0; j < benchmark_svs.size(); j++) benchmark_svs[j]->id = "SV_" + std::to_string(j);
		}
		if (!cdup_ids.empty()) {
			std::cerr << "Warning: --force-ids and --cdup-ids are both used. IDs of called SVs will not be forced." << std::endl;
		} else {
			for (int j = 0; j < called_svs.size(); j++) called_svs[j]->id = "SV_" + std::to_string(j);
		}
	}

	// if two variants have the same ID, give a warning that the user should consider --force-ids
	std::unordered_set<std::string> benchmark_id_count, called_id_count;
	for (sv_t* sv : benchmark_svs) {
		if (benchmark_id_count.count(sv->id)) {
			std::cerr << "Warning: benchmark file contains duplicated IDs. Sensitivity may not be computed correctly. Consider using --force-ids." << std::endl;
			break;
		} else {
			benchmark_id_count.insert(sv->id);
		}
	}
	for (sv_t* sv : called_svs) {
		if (called_id_count.count(sv->id)) {
			std::cerr << "Warning: called file contains duplicated IDs. Precision may not be computed correctly. Consider using --force-ids." << std::endl;
			break;
		} else {
			called_id_count.insert(sv->id);
		}
	}

	if (parsed_args.count("reference")) {
		std::string ref_fname = parsed_args["reference"].as<std::string>();
		chr_seqs.read_fasta_into_map(ref_fname);
	}

    if (parsed_args.count("tandem-repeats")) {
		std::string line;
		std::unordered_map<std::string, std::vector<repeat_t>> reps;
		std::unordered_map<std::string, std::vector<Interval<repeat_t>>> reps_iv;
		std::ifstream rep_f(parsed_args["tandem-repeats"].as<std::string>());
		while (getline(rep_f, line)) {
			if (line[0] == '#') continue;
			repeat_t r(line);
			reps[r.chr].push_back(r);
			reps_iv[r.chr].push_back(Interval<repeat_t>(r.start, r.end, r));
		}

		for (auto it = reps.begin(); it != reps.end(); it++) {
			reps_i[it->first] = new IntervalTree<repeat_t>(reps_iv[it->first]);
		}
    }

	for (sv_t* sv : called_svs) {
		if (sv->svtype() == "DEL") called_dels_by_chr[sv->chr].push_back(sv);
		else if (sv->svtype() == "INS" || sv->svtype() == "DUP") called_inss_by_chr[sv->chr].push_back(sv);
		else if (sv->svtype() == "INV") called_invs_by_chr[sv->chr].push_back(sv);
	}

	ctpl::thread_pool thread_pool(threads);
    std::vector<std::future<void> > futures;
	int BLOCK_SIZE = 1;
	for (int i = 0; i < benchmark_svs.size(); i+=BLOCK_SIZE) {
		sv_t* bsv = benchmark_svs[i];
		std::future<void> future = thread_pool.push(find_match, i, std::min(i+BLOCK_SIZE, (int) benchmark_svs.size()));
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

	std::set<std::string> b_tps, c_tps, b_gt_tps, c_gt_tps;

	// sort matches by rep (false first) and then score in descending order
	std::vector<sv_match_t> accepted_matches;
	std::sort(matches.begin(), matches.end(), [](sv_match_t& a, sv_match_t& b) {
		if (a.rep && !b.rep) return false;
		if (!a.rep && b.rep) return true;
		return a.score > b.score;
	});
	for (sv_match_t& match : matches) {
		if (exclusive && (b_tps.count(match.b_sv->id) || c_tps.count(match.c_sv->id))) continue;

		b_tps.insert(match.b_sv->id);
		c_tps.insert(match.c_sv->id);
		accepted_matches.push_back(match);
		if (compatible_gts(match.b_sv, match.c_sv)) {
			b_gt_tps.insert(match.b_sv->id);
			c_gt_tps.insert(match.c_sv->id);
		} else if (wrong_gts_fout.is_open()) {
			wrong_gts_fout << match.b_sv->id << " " << match.c_sv->id << std::endl;
		}
	}

	StripedSmithWaterman::Aligner aligner(1,4,6,1,false);
	for (sv_t* bsv : benchmark_svs) {
		if (!b_tps.count(bsv->id)) {
			accepted_matches.push_back(sv_match_t(bsv, NULL, false, aligner));
		}
	}

	std::sort(accepted_matches.begin(), accepted_matches.end(), [](sv_match_t& a, sv_match_t& b) {
		return a.b_sv->id < b.b_sv->id;
	});
	if (!report) {
		for (sv_match_t& match : accepted_matches) {
			std::cout << match.b_sv->id << " " << (match.c_sv == NULL ? "NONE" : match.c_sv->id);
			std::cout << (match.rep ? " REP" : "") << " " << (match.rep ? std::to_string(distance(match.b_sv, match.c_sv)) : "") << std::endl;
		}
	}
	if (called_to_benchmark_gts_fout.is_open()) {
		for (sv_match_t& match : accepted_matches) {
			if (match.c_sv != NULL) {
				called_to_benchmark_gts_fout << match.c_sv->id << " " << match.b_sv->print_gt() << std::endl;
			}
		}
		std::set<std::string> c_unchoosen_ids;
		for (sv_match_t& match : matches) {
			if (match.c_sv != NULL && !c_tps.count(match.c_sv->id)) {
				c_unchoosen_ids.insert(match.c_sv->id);
			}
		}
		for (std::string id : c_unchoosen_ids) {
			called_to_benchmark_gts_fout << id << " " << "./." << std::endl;
		}
	}

	// count tp and fn benchmark calls, by sv type
	std::unordered_map<std::string, int> n_benchmark_tp, n_benchmark_fn, n_benchmark_gt_tp, n_benchmark_gt_fn;
	for (sv_t* bsv : benchmark_svs) {
		if (b_tps.count(bsv->id)) {
			n_benchmark_tp[get_bsv_type(bsv)]++;
		} else {
			n_benchmark_fn[get_bsv_type(bsv)]++;
		}
		if (b_gt_tps.count(bsv->id)) {
			n_benchmark_gt_tp[get_bsv_type(bsv)]++;
		} else {
			n_benchmark_gt_fn[get_bsv_type(bsv)]++;
		}
	}

	std::unordered_map<std::string, int> n_called_tp, n_called_fp, n_called_gt_tp, n_called_gt_fp;
	for (sv_t* csv : called_svs) {
		if (c_tps.count(csv->id)) {
			n_called_tp[get_csv_type(csv)]++;
		} else {
			n_called_fp[get_csv_type(csv)]++;
		}
		if (c_gt_tps.count(csv->id)) {
			n_called_gt_tp[get_csv_type(csv)]++;
		} else {
			n_called_gt_fp[get_csv_type(csv)]++;
		}
	}

	if (report) {
		std::cout.precision(2);
		for (std::string svtype : {"DEL", "DUP", "INS", "INV"}) {
			int tp = n_benchmark_tp[svtype], fn = n_benchmark_fn[svtype];
			std::cout << svtype << " SENSITIVITY: " << tp << "/" << (tp+fn) << " = " << tp/std::max(1.0, double(tp+fn)) << " ";

			tp = n_benchmark_gt_tp[svtype]; fn = n_benchmark_gt_fn[svtype];
			std::cout << "(considering GT: " << tp << "/" << (tp+fn) << " = " << tp/std::max(1.0, double(tp+fn)) << ")" << std::endl;
			
			tp = n_called_tp[svtype]; int fp = n_called_fp[svtype];
			std::cout << svtype << " PRECISION: " << tp << "/" << (tp+fp) << " = " << tp/std::max(1.0, double(tp+fp)) << " ";

			tp = n_called_gt_tp[svtype]; fp = n_called_gt_fp[svtype];
			std::cout << "(considering GT: " << tp << "/" << (tp+fp) << " = " << tp/std::max(1.0, double(tp+fp)) << ")" << std::endl;
		}
	}

	if (fps_fout != NULL) {
		for (sv_t* csv : called_svs) {
			if (!c_tps.count(csv->id)) {
				if (bcf_write(fps_fout, fps_hdr, csv->vcf_entry) != 0) {
					std::cerr << "Error writing variant to file." << std::endl;
					exit(1);
				}
			}
		}
		bcf_close(fps_fout);
		bcf_hdr_destroy(fps_hdr);
	}
}
