#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <unistd.h>

#include <htslib/sam.h>
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
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
    int min_avg_base_qual;
    int pop_avg_crossing_is = 0;
    bool per_contig_stats = false;

    void parse(std::string stats_file, bool per_contig_stats) {
        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> stats_params;
        std::ifstream fin(stats_file);
        std::string name, contig_name, value;
        while (fin >> name >> contig_name >> value) {
            stats_params[name][contig_name] = value;
            if (name == "min_is" && contig_name == ".") min_is = std::stoi(value);
            if (name == "max_is" && contig_name == ".") max_is = std::stoi(value);
            if (name == "read_len" && contig_name == ".") read_len = std::stoi(value);
            if (name == "min_depth") min_depths[contig_name] = std::stoi(value);
            if (name == "median_depth") median_depths[contig_name] = std::stoi(value);
            if (name == "max_depth") max_depths[contig_name] = std::stoi(value);
            if (name == "min_avg_base_qual") min_avg_base_qual = std::stoi(value);
            if (name == "pop_avg_crossing_is" && contig_name == ".") pop_avg_crossing_is = std::stoi(value);
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

void to_uppercase(char* s) {
    for (int i = 0; s[i] != '\0'; i++) {
        s[i] = toupper(s[i]);
    }
}

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

#endif
