#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <unistd.h>

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
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

struct consensus_t {
    bool left_clipped;
    int contig_id;
    hts_pos_t start, breakpoint, end;
    std::string consensus;
    int fwd_clipped, rev_clipped;
    uint8_t max_mapq;
    hts_pos_t remap_boundary;
    int clip_len, lowq_clip_portion;
    bool is_hsr = false;

    static const int LOWER_BOUNDARY_NON_CALCULATED = 0, UPPER_BOUNDARY_NON_CALCULATED = INT32_MAX;
    static const int UNKNOWN_CLIP_LEN = INT16_MAX;

    consensus_t(bool left_clipped, int contig_id, hts_pos_t start, hts_pos_t breakpoint, hts_pos_t end,
                const std::string& consensus, int fwd_clipped, int rev_clipped, uint8_t max_mapq, hts_pos_t remap_boundary,
                int lowq_clip_portion)
                : left_clipped(left_clipped), contig_id(contig_id), start(start), breakpoint(breakpoint), end(end),
                consensus(consensus), fwd_clipped(fwd_clipped), rev_clipped(rev_clipped), max_mapq(max_mapq), 
                remap_boundary(remap_boundary), lowq_clip_portion(lowq_clip_portion) {}

    std::string to_string() {
        std::stringstream ss;
        ss << start << " " << end << " " << breakpoint << (left_clipped ? " L " : " R ") << consensus << " ";
        ss << fwd_clipped << " " << rev_clipped << " " << (int)max_mapq << " " << remap_boundary << " " << lowq_clip_portion;
        return ss.str();
    }

    int supp_clipped_reads() { return fwd_clipped + rev_clipped; }
};

int get_left_clip_size(StripedSmithWaterman::Alignment& aln) {
	return cigar_int_to_op(aln.cigar[0]) == 'S' ? cigar_int_to_len(aln.cigar[0]) : 0;
}
int get_right_clip_size(StripedSmithWaterman::Alignment& aln) {
	return cigar_int_to_op(aln.cigar[aln.cigar.size()-1]) == 'S' ? cigar_int_to_len(aln.cigar[aln.cigar.size()-1]) : 0;
}
bool is_left_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return get_left_clip_size(aln) >= min_clip_len;
}
bool is_right_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return get_right_clip_size(aln) >= min_clip_len;
}
bool is_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return is_left_clipped(aln, min_clip_len) || is_right_clipped(aln, min_clip_len);
}

#endif
