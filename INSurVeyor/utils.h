#ifndef SMALLINSFINDER_UTILS_H
#define SMALLINSFINDER_UTILS_H

#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <emmintrin.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "../src/utils.h"

struct inss_contig_map_t {

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

struct inss_chr_seq_t {
    char* seq;
    hts_pos_t len;

    inss_chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
    ~inss_chr_seq_t() {delete[] seq;}
};
struct inss_chr_seqs_map_t {
    std::unordered_map<std::string, inss_chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            seqs[seq_name] = new inss_chr_seq_t(chr_seq, seq->seq.l);
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

    ~inss_chr_seqs_map_t() {
        clear();
    }
};

struct inss_suffix_prefix_aln_t {
    int overlap, score, mismatches;

    inss_suffix_prefix_aln_t(int overlap, int score, int mismatches) : overlap(overlap), score(score), mismatches(mismatches) {}
};


int popcnt(uint32_t x) {
    x = x - ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

int number_of_mismatches_SIMD(const char* s1, const char* s2, int len) {
    // count number of mismatches between the first len characters of s1 and s2 using SIMD instructions
    int n_mismatches = 0;
    __m128i* s1_128 = (__m128i*) s1;
    __m128i* s2_128 = (__m128i*) s2;
    int n_128 = len/16;
    for (int i = 0; i < n_128; i++) {
        __m128i cmp = _mm_cmpeq_epi8(_mm_loadu_si128(s1_128), _mm_loadu_si128(s2_128));
        n_mismatches += 16-popcnt(_mm_movemask_epi8(cmp));
        s1_128++;
        s2_128++;
    }
    for (int i = n_128*16; i < len; i++) {
        if (s1[i] != s2[i]) n_mismatches++;
    }
    return n_mismatches;
}

// Finds the best alignment between a suffix of s1 and a prefix of s2
// Disallows gaps
inss_suffix_prefix_aln_t inss_aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score, double max_seq_error,
                                      int min_overlap = 1, int max_overlap = INT32_MAX, int max_mismatches = INT32_MAX) {
    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < s1.length()-min_overlap+1; i++) {
        if (i+s2.length() < s1.length()) continue;

        int sp_len = s1.length()-i;
        if (best_score >= sp_len*match_score) break; // current best score is unbeatable

        const char* s1_suffix = s1.data()+i;
        const char* s2_prefix = s2.data();
        int mismatches = 0;
//        mismatches = number_of_mismatches_SIMD(s1.data()+i, s2.data(), sp_len);
        while (*s1_suffix) {
            if (*s1_suffix != *s2_prefix) mismatches++;
            s1_suffix++; s2_prefix++;
        }

        int score = (sp_len-mismatches)*match_score + mismatches*mismatch_score;

        int max_acceptable_mm = max_seq_error == 0.0 ? 0 : std::max(1.0, sp_len*max_seq_error);
        if (best_score < score && mismatches <= max_acceptable_mm && mismatches <= max_mismatches) {
            best_score = score;
            best_aln_mismatches = mismatches;
            overlap = sp_len;
        }
    }
    return inss_suffix_prefix_aln_t(overlap, best_score, best_aln_mismatches);
}

template<typename T>
inline T inss_max(T a, T b, T c) { return std::max(std::max(a,b), c); }

template<typename T>
inline T inss_max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

int inss_max_pos(int* a, int n) {
	int pos = 0;
	for (int i = 1; i < n; i++) {
		if (a[i] > a[pos]) pos = i;
	}
	return pos;
}

bool inss_file_exists(std::string& fname) {
	return std::ifstream(fname).good();
}

int64_t inss_overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

struct inss_insertion_t {
    std::string id;
    std::string chr;
    hts_pos_t start, end;
    int overlap, r_disc_pairs, l_disc_pairs, rc_fwd_reads = 0, rc_rev_reads = 0, lc_fwd_reads = 0, lc_rev_reads = 0;
    int r_conc_pairs = 0, l_conc_pairs = 0, median_lf_cov = 0, median_rf_cov = 0;
    int mh_len = 0;
    std::string ins_seq;
    std::string left_anchor, right_anchor, left_anchor_cigar, right_anchor_cigar;
    int left_seq_cov = 0, right_seq_cov = 0;
    bool left_bp_precise = false, right_bp_precise = false;
    double rc_avg_nm = 0.0, lc_avg_nm = 0.0;

    inss_insertion_t(std::string chr, hts_pos_t start, hts_pos_t end, int r_disc_pairs, int l_disc_pairs,
    		int rc_fwd_reads, int rc_rev_reads, int lc_fwd_reads, int lc_rev_reads, int overlap, std::string ins_seq) :
	chr(chr), start(start), end(end), r_disc_pairs(r_disc_pairs), l_disc_pairs(l_disc_pairs),
	rc_fwd_reads(rc_fwd_reads), rc_rev_reads(rc_rev_reads), lc_fwd_reads(lc_fwd_reads), lc_rev_reads(lc_rev_reads),
	overlap(overlap), ins_seq(ins_seq) {}

    int rc_reads() { return rc_fwd_reads + rc_rev_reads; }
    int lc_reads() { return lc_fwd_reads + lc_rev_reads; }

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + ins_seq;
    }
};

int inss_get_left_clip_size(const StripedSmithWaterman::Alignment& aln) {
    uint32_t l = aln.cigar[0];
    return cigar_int_to_op(l) == 'S' ? cigar_int_to_len(l) : 0;
}
int inss_get_right_clip_size(const StripedSmithWaterman::Alignment& aln) {
    uint32_t r = aln.cigar[aln.cigar.size()-1];
    return cigar_int_to_op(r) == 'S' ? cigar_int_to_len(r) : 0;
}
bool inss_is_left_clipped(const StripedSmithWaterman::Alignment& aln) {
	return inss_get_left_clip_size(aln) > 0;
}
bool inss_is_right_clipped(const StripedSmithWaterman::Alignment& aln) {
    return inss_get_right_clip_size(aln) > 0;
}

// Returns a vector scores s.t. scores[N] contains the score of the alignment between reference[aln.ref_begin:aln.ref_begin+N]
// and query[aln.query_begin:aln.query_begin+M]
std::vector<int> inss_ssw_cigar_to_prefix_ref_scores(uint32_t* cigar, int cigar_len,
		int match = 1, int mismatch = -4, int gap_open = -6, int gap_extend = -1) {

	std::vector<int> scores;
	scores.push_back(0);
	int score = 0;

	for (int i = 0; i < cigar_len; i++) {
		uint32_t c = cigar[i];
		char op = cigar_int_to_op(c);
		int len = cigar_int_to_len(c);
		if (op == 'S') {
			continue;
		} else if (op == '=' || op == 'X') {
			for (int i = 0; i < len; i++) {
				score += (op == '=' ? match : mismatch);
				scores.push_back(score);
			}
		} else if (op == 'D') {
			scores.resize(scores.size()+len, score);
			score += gap_open + len*gap_extend;
		} else if (op == 'I') {
			score += gap_open + len*gap_extend;
		}
	}

	return scores;
}

void inss_to_uppercase(char* s) {
    for (int i = 0; s[i] != '\0'; i++) {
        s[i] = toupper(s[i]);
    }
}

std::string inss_cigar_to_string(std::vector<uint32_t>& cigar) {
	std::stringstream ss;
	for (uint32_t c : cigar) {
		ss << cigar_int_to_len(c) << cigar_int_to_op(c);
	}
	return ss.str();
}

#endif //SMALLINSFINDER_UTILS_H
