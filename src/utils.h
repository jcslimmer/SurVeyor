#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <queue>
#include <numeric>
#include <cmath>
#include <random>
#include <unistd.h>
#include <algorithm>

#include <htslib/sam.h>
#include "htslib/hts.h"
#include "htslib/kseq.h"
KSEQ_INIT(int, read)

struct config_t {

    int threads, seed;
    int min_clip_len, min_stable_mapq, min_diff_hsr;
    int min_sv_size, max_trans_size;
    int min_size_for_depth_filtering;
    double max_seq_error;
    int max_clipped_pos_dist;
    bool per_contig_stats;
    std::string sampling_regions, version;

    const int min_score_diff = 15;
    const int high_confidence_mapq = 60;

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
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        min_stable_mapq = std::stoi(config_params["min_stable_mapq"]);
        min_diff_hsr = std::stoi(config_params["min_diff_hsr"]);
        min_sv_size = std::stoi(config_params["min_sv_size"]);
        max_trans_size = std::stoi(config_params["max_trans_size"]);
        min_size_for_depth_filtering = std::stoi(config_params["min_size_for_depth_filtering"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        max_clipped_pos_dist = std::stoi(config_params["max_clipped_pos_dist"]);
        per_contig_stats = std::stoi(config_params["per_contig_stats"]);
        sampling_regions = config_params["sampling_regions"];
        version = config_params["version"];
    }
};

struct stats_t {

    int min_is, max_is;
    int read_len;
    std::unordered_map<std::string, int> min_depths, median_depths, max_depths;
    std::unordered_map<int, int> min_disc_pairs_by_insertion_size, max_disc_pairs_by_insertion_size;
    std::unordered_map<int, int> median_disc_pairs_by_del_size;
    int min_avg_base_qual;
    int pop_avg_crossing_is = 0;
    bool per_contig_stats = false;

    void parse(std::string stats_file, bool per_contig_stats) {
        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> stats_params;
        std::ifstream fin(stats_file);
        std::string stat_name, stat_subname, value;
        while (fin >> stat_name >> stat_subname >> value) {
            stats_params[stat_name][stat_subname] = value;
            if (stat_name == "min_is" && stat_subname == ".") min_is = std::stoi(value);
            if (stat_name == "max_is" && stat_subname == ".") max_is = std::stoi(value);
            if (stat_name == "read_len" && stat_subname == ".") read_len = std::stoi(value);
            if (stat_name == "min_depth") min_depths[stat_subname] = std::stoi(value);
            if (stat_name == "median_depth") median_depths[stat_subname] = std::stoi(value);
            if (stat_name == "max_depth") max_depths[stat_subname] = std::stoi(value);
            if (stat_name == "min_disc_pairs_by_insertion_size") min_disc_pairs_by_insertion_size[std::stoi(stat_subname)] = std::stoi(value);
            if (stat_name == "max_disc_pairs_by_insertion_size") max_disc_pairs_by_insertion_size[std::stoi(stat_subname)] = std::stoi(value);
            if (stat_name == "median_disc_pairs_by_del_size") median_disc_pairs_by_del_size[std::stoi(stat_subname)] = std::stoi(value);
            if (stat_name == "min_avg_base_qual") min_avg_base_qual = std::stoi(value);
            if (stat_name == "pop_avg_crossing_is" && stat_subname == ".") pop_avg_crossing_is = std::stoi(value);
        }
        fin.close();
        this->per_contig_stats = per_contig_stats;
    }

    int get_min_depth(std::string contig_name) {
        if (per_contig_stats && min_depths.count(contig_name))
            return min_depths[contig_name];
        else return min_depths["."];
    }
    int get_median_depth(std::string contig_name) {
        if (per_contig_stats && median_depths.count(contig_name))
            return median_depths[contig_name];
        else return median_depths["."];
    }
    int get_max_depth(std::string contig_name) {
        if (per_contig_stats && max_depths.count(contig_name))
            return max_depths[contig_name];
        else return max_depths["."];
    }

    int get_min_disc_pairs_by_insertion_size(int is) {
        if (min_disc_pairs_by_insertion_size.count(is))
            return min_disc_pairs_by_insertion_size[is];
        else return min_disc_pairs_by_insertion_size[max_is];
    }
    int get_max_disc_pairs_by_insertion_size(int is) {
        if (max_disc_pairs_by_insertion_size.count(is))
            return max_disc_pairs_by_insertion_size[is];
        else return max_disc_pairs_by_insertion_size[max_is];
    }

    int get_median_disc_pairs_by_del_size(int ds) {
        if (min_disc_pairs_by_insertion_size.count(ds))
            return min_disc_pairs_by_insertion_size[ds];
        else return min_disc_pairs_by_insertion_size[max_is];
    }
};

struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    contig_map_t() {}
    contig_map_t(std::string workdir) {
        load(workdir);
    }

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
    size_t get_id(std::string& name) {return name_to_id[name];};
};

template<typename T>
T mean(std::vector<T>& v) {
    return std::accumulate(v.begin(), v.end(), (T)0.0)/v.size();
}

template<typename T>
T stddev(std::vector<T>& v) {
    T m = mean(v);
    T sum = 0;
    for (T& e : v) {
        sum += (e - m)*(e - m);
    }
    return std::sqrt(sum/v.size());
}

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

bool file_exists(std::string& fname) {
	return std::ifstream(fname).good();
}

struct chr_seqs_map_t {
    struct chr_seq_t {
        char* seq;
        hts_pos_t len;

        chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
        ~chr_seq_t() {delete[] seq;}
    };

    std::unordered_map<std::string, chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname) {
        ordered_contigs.clear();
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            if (seqs.count(seq_name)) delete seqs[seq_name];
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    void read_lens_into_map(std::string& reference_fname) {
        ordered_contigs.clear();
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            if (seqs.count(seq_name)) delete seqs[seq_name];
            seqs[seq_name] = new chr_seq_t(NULL, seq->seq.l);
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

struct base_frequencies_t {
    int a = 0, c = 0, g = 0, t = 0;

    base_frequencies_t() {}
    base_frequencies_t(int a, int c, int g, int t) : a(a), c(c), g(g), t(t) {}

    int max() {
        return std::max(std::max(a, c), std::max(g, t));
    }
    int tot() {
        return a + c + g + t;
    }
    double max_freq() {
        return max()/double(std::max(1, tot()));
    }

    bool empty() {
        return a + c + g + t == 0;
    }
};

base_frequencies_t get_base_frequencies(const char* seq, int len) {
    int a = 0, c = 0, g = 0, t = 0;
    for (int i = 0; i < len; i++) {
        char b = seq[i];
        if (b == 'A' || b == 'a') a++;
        else if (b == 'C' || b == 'c') c++;
        else if (b == 'G' || b == 'g') g++;
        else if (b == 'T' || b == 't') t++;
    }
    return base_frequencies_t(a, c, g, t);
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
	return get_base_frequencies(seq, len).max_freq() >= 0.8;
}
bool is_homopolymer(std::string seq) {
	return is_homopolymer(seq.data(), seq.length());
}

void to_uppercase(char* s) {
    for (int i = 0; s[i] != '\0'; i++) {
        s[i] = toupper(s[i]);
    }
}

int find_char_in_str(std::string s, char c, int start = 0, int stop = -1) {
    if (start < 0) start = 0;
    if (stop == -1 || stop > s.length()) stop = s.length();
    for (int i = start; i < stop; i++) {
        if (s[i] == c) return i;
    }
    return -1;
}

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

std::vector<int> find_rev_topological_order(int n, std::vector<int> out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

struct random_pos_generator_t {

    const int NO_SAMPLING_PADDING = 1000;

    chr_seqs_map_t& chr_seqs_map;
    hts_pos_t reference_len = 0;
    std::mt19937 rng;
    std::uniform_int_distribution<hts_pos_t> dist;

    struct region_t {
        std::string chr;
        hts_pos_t start, end;
        region_t(std::string chr, hts_pos_t start, hts_pos_t end) : chr(chr), start(start), end(end) {}
    };
    std::vector<region_t> regions;

    random_pos_generator_t(chr_seqs_map_t& chr_seqs_map, int seed, std::string sampling_regions_fname = "") : chr_seqs_map(chr_seqs_map) {
        rng.seed(seed);

        if (sampling_regions_fname.empty()) {
            for (std::string& chr : chr_seqs_map.ordered_contigs) {
                regions.push_back(region_t(chr, NO_SAMPLING_PADDING, chr_seqs_map.get_len(chr)-NO_SAMPLING_PADDING));
            }
        } else {
            std::ifstream fin(sampling_regions_fname);
            std::string line;
            while (std::getline(fin, line)) {
                std::istringstream iss(line);

                // tokenize iss imitating >>, do not assume format
                std::vector<std::string> tokens;
                std::string token;
                while (iss >> token) {
                    tokens.push_back(token);
                }

                if (tokens.size() == 1) {
                    std::string chr = tokens[0];
                    regions.push_back(region_t(chr, NO_SAMPLING_PADDING, chr_seqs_map.get_len(chr)-NO_SAMPLING_PADDING));
                } else if (tokens.size() == 3) {
                    std::string chr = tokens[0];
                    hts_pos_t start = std::stoll(tokens[1]);
                    hts_pos_t end = std::stoll(tokens[2]);
                    regions.push_back(region_t(chr, start, end));
                } else {
                    std::string error_msg = "Error: invalid sampling region format in " + sampling_regions_fname + "\n";
                    throw std::runtime_error(error_msg);
                }
            }
        }
        regions.erase(std::remove_if(regions.begin(), regions.end(), [](region_t& r) {return r.start >= r.end;}), regions.end());

        for (region_t& r : regions) {
            reference_len += r.end - r.start;
        }

        dist = std::uniform_int_distribution<hts_pos_t>(0, reference_len-1);
    }

    std::pair<std::string, hts_pos_t> get_random_pos() {
        hts_pos_t random_pos = dist(rng);

        for (region_t& r : regions) {
            if (random_pos < r.end - r.start) {
                return {r.chr, r.start + random_pos};
            } else {
                random_pos -= r.end - r.start;
            }
        }
        throw std::runtime_error("Error: random_pos_generator_t::get_random_pos() failed\n");
    }
};

#endif
