#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <unistd.h>

#include "htslib/kseq.h"
KSEQ_INIT(int, read)

struct config_t {
    
    int threads;
    int min_clip_len, min_stable_mapq, min_diff_hsr;
    double max_seq_error;
    int max_clipped_pos_dist;

    int min_score_diff = 15;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        threads = std::stoi(config_params["threads"]);
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        min_stable_mapq = std::stoi(config_params["min_stable_mapq"]);
        min_diff_hsr = std::stoi(config_params["min_diff_hsr"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        max_clipped_pos_dist = std::stoi(config_params["max_clipped_pos_dist"]);
    }
};

struct stats_t {

    int read_len;
    int max_is;

    void parse(std::string stats_file) {
        std::unordered_map<std::string, std::string> stats_params;
        std::ifstream fin(stats_file);
        std::string name, value;
        while (fin >> name >> value) {
            stats_params[name] = value;
        }
        fin.close();

        read_len = std::stoi(stats_params["read_len"]);
        max_is = std::stoi(stats_params["max_is"]);
    }

};

struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    void parse(std::string workdir) {
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

#endif
