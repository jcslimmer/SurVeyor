#ifndef SW_UTILS_H
#define SW_UTILS_H

#include <cstdint>
#include <iostream>

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "simd_macros.h"
#include "sam_utils.h"
#include "types.h"
#include "utils.h"

int get_left_clip_size(StripedSmithWaterman::Alignment& aln) {
	return cigar_int_to_op(aln.cigar[0]) == 'S' ? cigar_int_to_len(aln.cigar[0]) : 0;
}
int get_right_clip_size(StripedSmithWaterman::Alignment& aln) {
	return cigar_int_to_op(aln.cigar[aln.cigar.size()-1]) == 'S' ? cigar_int_to_len(aln.cigar[aln.cigar.size()-1]) : 0;
}
int get_unclipped_start(StripedSmithWaterman::Alignment& aln) {
	return aln.ref_begin - get_left_clip_size(aln);
}
int get_unclipped_end(StripedSmithWaterman::Alignment& aln) {
	return aln.ref_end + get_right_clip_size(aln);
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

size_t gcd(size_t a, size_t b) {
    return (b == 0) ? a : gcd(b, a % b);
}
size_t lcm(size_t a, size_t b) {
    return a * b / gcd(a, b);
}

int* smith_waterman_gotoh(const char* ref, int ref_len, const char* read, int read_len,
		int match_score, int mismatch_penalty, int gap_open, int gap_extend) {

	const int INF = 1000000;

	// turn read_len+1 into a multiple of INT_PER_BLOCK
	int read_len_rounded = (read_len+1+INT_PER_BLOCK-1)/INT_PER_BLOCK*INT_PER_BLOCK;

	const char* alphabet = "NACGT";
	const int alphabet_size = strlen(alphabet);

	int* H = NULL;
	int* E = NULL;
	int* F = NULL;
	int* prefix_scores = NULL;

	int** profile = new int*[alphabet_size];
	size_t alignment = lcm(BYTES_PER_BLOCK, sizeof(void*));
	int p1 = posix_memalign(reinterpret_cast<void**>(&H), alignment, 2*read_len_rounded * sizeof(int));
	int p2 = posix_memalign(reinterpret_cast<void**>(&E), alignment, 2*read_len_rounded * sizeof(int));
	int p3 = posix_memalign(reinterpret_cast<void**>(&F), alignment, 2*read_len_rounded * sizeof(int));
	int p4 = posix_memalign(reinterpret_cast<void**>(&prefix_scores), alignment, read_len_rounded * sizeof(int));
	int p5 = 0;
	for (int i = 0; i < alphabet_size; i++) {
		p5 += posix_memalign(reinterpret_cast<void**>(&profile[i]), alignment, read_len_rounded * sizeof(int));
		profile[i][0] = 0;
		for (int j = 1; j <= read_len; j++) {
			profile[i][j] = (read[j-1] == alphabet[i]) ? match_score : mismatch_penalty;
		}
		for (int j = read_len+1; j < read_len_rounded; j++) {
			profile[i][j] = 0;
		}
	}
	if (p1 || p2 || p3 || p4 || p5) {
		std::cerr << "Error allocating aligned memory of size " << (2*read_len_rounded * sizeof(int)) << std::endl;
	}

	int* H_prev = H, *H_curr = H+read_len_rounded;
	int* E_prev = E, *E_curr = E+read_len_rounded;
	int* F_prev = F, *F_curr = F+read_len_rounded;

	std::fill(H_prev, H_prev+read_len_rounded, 0);
	std::fill(E_prev, E_prev+read_len_rounded, 0);
	std::fill(F_prev, F_prev+read_len_rounded, 0);
	H_curr[0] = F_curr[0] = E_curr[0] = 0;

	SIMD_INT gap_open_v = SET1_INT(gap_open);
	SIMD_INT gap_open_v_pos = SET1_INT(-gap_open);
	SIMD_INT gap_extend_v = SET1_INT(gap_extend);
	SIMD_INT zero_v = SET1_INT(0);

	std::fill(prefix_scores, prefix_scores+read_len, 0);
	for (int i = 1; i <= ref_len; i++) {
		for (int j = 0; j < read_len_rounded; j += INT_PER_BLOCK) {
			SIMD_INT H_up_v = LOAD_INT((SIMD_INT*)&H_prev[j]);
			SIMD_INT E_up_v = LOAD_INT((SIMD_INT*)&E_prev[j]);
			SIMD_INT F_up_v = LOAD_INT((SIMD_INT*)&F_prev[j]);
			SIMD_INT m1 = ADD_INT(gap_open_v, MAX_INT(H_up_v, F_up_v));
			SIMD_INT E_curr_v = MAX_INT(m1, ADD_INT(gap_extend_v, E_up_v));
			STORE_INT((SIMD_INT*)&E_curr[j], E_curr_v);
		}

		int* ref_profile = profile[0];
		switch (ref[i-1]) {
			case 'A': ref_profile = profile[1]; break;
			case 'C': ref_profile = profile[2]; break;
			case 'G': ref_profile = profile[3]; break;
			case 'T': ref_profile = profile[4]; break;
		}
		for (int j = 1; j <= read_len_rounded-INT_PER_BLOCK; j += INT_PER_BLOCK) {
			SIMD_INT H_diag_v = LOAD_INT((SIMD_INT*)&H_prev[j-1]);
			SIMD_INT F_diag_v = LOAD_INT((SIMD_INT*)&F_prev[j-1]);
			SIMD_INT E_diag_v = LOAD_INT((SIMD_INT*)&E_prev[j-1]);
			SIMD_INT m1 = MAX_INT(H_diag_v, F_diag_v);
			m1 = MAX_INT(m1, E_diag_v);
			SIMD_INT H_curr_v = MAX_INT(m1, zero_v);

			SIMD_INT profile_curr = LOADU_INT((SIMD_INT*)&ref_profile[j]);
			H_curr_v = ADD_INT(H_curr_v, profile_curr);

			STOREU_INT((SIMD_INT*)&H_curr[j], H_curr_v);

			SIMD_INT prefix_v = LOAD_INT((SIMD_INT*)&prefix_scores[j-1]);
			prefix_v = MAX_INT(prefix_v, H_curr_v);
			STORE_INT((SIMD_INT*)&prefix_scores[j-1], prefix_v);
		}
		for (int j = read_len_rounded-INT_PER_BLOCK; j < read_len_rounded; j++) {
			H_curr[j] = ref_profile[j] + max(H_prev[j-1], F_prev[j-1], E_prev[j-1], 0);
			prefix_scores[j] = std::max(prefix_scores[j], H_curr[j]);
		}

		int j = 0;
		for (; j < read_len_rounded; j += INT_PER_BLOCK) {
			SIMD_INT H_curr_v = LOAD_INT((SIMD_INT*)&H_curr[j]);
			auto cmp = CMP_GT_INT32(H_curr_v, gap_open_v_pos);
			if (cmp) break;
			STORE_INT((SIMD_INT*)&F_curr[j], zero_v);
		}
		if (j == 0) j = 1;
		for (; j <= read_len; j++) {
			F_curr[j] = std::max(gap_open + H_curr[j-1], gap_extend + F_curr[j-1]);
		}

		std::swap(H_prev, H_curr);
		std::swap(E_prev, E_curr);
		std::swap(F_prev, F_curr);
	}

	free(H);
	free(E);
	free(F);
	for (int i = 0; i < alphabet_size; i++) free(profile[i]);
	delete[] profile;

	return prefix_scores;
}

std::vector<int> compute_prefix_scores(std::vector<uint32_t> cigar, int query_len, int match_score, int mismatch_score, 
						int gap_open_score, int gap_extend_score) {
	
	std::vector<int> prefix_scores(query_len + 1, 0); // Stores scores for each prefix up to query_len.
    int score = 0, query_prefix_len = 0;

	for (int i = 0; i < cigar.size() && query_prefix_len < query_len; ++i) {
        int op = cigar_int_to_op(cigar[i]);
        int len = cigar_int_to_len(cigar[i]);

		for (int k = 0; k < len; ++k) {
            if (op == '=') {
                score += match_score;
                prefix_scores[++query_prefix_len] = score;
            } else if (op == 'X') {
                score += mismatch_score;
                prefix_scores[++query_prefix_len] = score;
            } else if (op == 'I') {
                if (k == 0) {
                    score += gap_open_score;
                } else {
                    score += gap_extend_score;
                }
                prefix_scores[++query_prefix_len] = score;
            } else if (op == 'D') {
                if (k == 0) {
                    score += gap_open_score;
                } else {
                    score += gap_extend_score;
                }
                // No change to query_prefix_len since deletion skips query sequence.
            } else if (op == 'S') { // Soft-clipping
                // Soft-clipped bases don't change the score, but they do increment query_prefix_len.
                prefix_scores[++query_prefix_len] = score;
            }
        }
	}

	return prefix_scores;
}

std::vector<int> compute_suffix_scores(std::vector<uint32_t> cigar, int query_len, int match_score, int mismatch_score, 
						int gap_open_score, int gap_extend_score) {

	std::vector<uint32_t> cigar_rev(cigar.rbegin(), cigar.rend());
	std::vector<int> suffix_scores = compute_prefix_scores(cigar_rev, query_len, match_score, mismatch_score, gap_open_score, gap_extend_score);
	std::reverse(suffix_scores.begin(), suffix_scores.end());
	return suffix_scores;
}

std::pair<int, int> find_aln_prefix_score(std::vector<uint32_t>& cigar, int ref_prefix_len, int match_score, int mismatch_score,
						int gap_open_score, int gap_extend_score) {
	
	if (cigar.empty()) return {0, 0};

	int score = 0, query_prefix_len = 0, i = 0;
	if (cigar_int_to_op(cigar[0]) == 'S') {
		query_prefix_len = cigar_int_to_len(cigar[0]);
		i++;
	}

	for (int j = 0; i < cigar.size() && j < ref_prefix_len; i++) {
		int op = cigar_int_to_op(cigar[i]);
		int len = cigar_int_to_len(cigar[i]);
		if (op == '=') {
			len = std::min(len, ref_prefix_len-j);
			score += len*match_score;
			query_prefix_len += len;
			j += len;
		} else if (op == 'X') {
			len = std::min(len, ref_prefix_len-j);
			score += len*mismatch_score;
			query_prefix_len += len;
			j += len;
		} else if (op == 'I') {
			score += gap_open_score + (len-1)*gap_extend_score;
			query_prefix_len += len;
		} else if (op == 'D') {
			len = std::min(len, ref_prefix_len-j);
			score += gap_open_score + (len-1)*gap_extend_score;
			j += len;
		} else if (op == 'S') {
			query_prefix_len += len;
		}
	}
	return {score, query_prefix_len};
}

std::pair<int, int> find_aln_suffix_score(std::vector<uint32_t> cigar, int ref_suffix_len, int match_score, int mismatch_score,
						  int gap_open_score, int gap_extend_score) {
	std::vector<uint32_t> rev_cigar(cigar.rbegin(), cigar.rend());
	return find_aln_prefix_score(rev_cigar, ref_suffix_len, match_score, mismatch_score, gap_open_score, gap_extend_score);
}


// the (min) length of the query prefix that cover ref_len bases on the reference
int query_prefix_len(std::vector<uint32_t>& cigar, int ref_len, int& extra_deletion) {
    int query_len = 0;
    int ref_count = 0;
    for (auto it = cigar.begin(); it != cigar.end(); ++it) {
        char op = cigar_int_to_op(*it);
        int len = cigar_int_to_len(*it);
        if (op == '=' || op == 'X' || op == 'D') {
            ref_count += len;
        }
        if (op == '=' || op == 'X' || op == 'I' || op == 'S') {
            query_len += len;
        }
        if (ref_count >= ref_len) {
            int overshoot = ref_count - ref_len;
            if (op == '=' || op == 'X') {
                query_len -= overshoot;
            } else if (op == 'D') {
                extra_deletion = overshoot;
            }
            break;
        }
    }
    return query_len;
}
int query_suffix_len(std::vector<uint32_t>& cigar, int ref_len, int& extra_deletion) {
    auto cigar_rev = std::vector<uint32_t>(cigar.rbegin(), cigar.rend());
    return query_prefix_len(cigar_rev, ref_len, extra_deletion);
}


struct suffix_prefix_aln_t {
    int overlap, score, mismatches;

    suffix_prefix_aln_t(int overlap, int score, int mismatches) : overlap(overlap), score(score), mismatches(mismatches) {}

	double mismatch_rate() { return overlap > 0 ? double(mismatches)/overlap : 1; }
};

suffix_prefix_aln_t aln_suffix_prefix_perfect(std::string& s1, std::string& s2, int min_overlap = 1) {
    int best_overlap = 0, overlap = 0;

    const char* _s1 = s1.c_str(),* _s2 = s2.c_str();
    int _s1_len = s1.length(), _s2_len = s2.length();

    int max_overlap = std::min(s1.length(), s2.length());
    for (int i = max_overlap; i >= min_overlap; i--) {
    	if (strncmp(_s1+(_s1_len-i), _s2, i) == 0) {
			return suffix_prefix_aln_t(i, i, 0);
    	}
    }
    return suffix_prefix_aln_t(0, 0, 0);
}

// Finds the best alignment between a suffix of s1 and a prefix of s2
// Disallows gaps
int number_of_mismatches_fast(const char* s1, const char* s2, int len) {
    // count number of mismatches between the first len characters of s1 and s2
    int n_matches = 0;
    SIMD_INT* s1_it = (SIMD_INT*) s1;
    SIMD_INT* s2_it = (SIMD_INT*) s2;
    int scaled_len = len/BYTES_PER_BLOCK;
    for (int i = 0; i < scaled_len; i++) {
        SIMD_INT n1 = LOADU_INT(s1_it);
        SIMD_INT n2 = LOADU_INT(s2_it);
        n_matches += COUNT_EQUAL_BYTES(n1, n2);
        s1_it++;
        s2_it++;
    }
    for (int i = scaled_len*BYTES_PER_BLOCK; i < len; i++) {
        if (s1[i] == s2[i]) n_matches++;
    }
    return len-n_matches;
}
suffix_prefix_aln_t aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score, double max_seq_error,
                                      int min_overlap = 1, int max_overlap = INT32_MAX, int max_mismatches = INT32_MAX) {

	if (max_seq_error == 0.0 || max_mismatches == 0) return aln_suffix_prefix_perfect(s1, s2, min_overlap);

    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < s1.length()-min_overlap+1; i++) {
        if (i+s2.length() < s1.length()) continue;	// this means that the suffix of s1 I am considering is larger than the the whole s2,
                                                    // therefore it wouldn't be a suffix-prefix alignment   

        int sp_len = s1.length()-i;
        if (best_score >= sp_len*match_score) break; // current best score is unbeatable

        const char* s1_suffix = s1.data()+i;
        const char* s2_prefix = s2.data();
        int mismatches = mismatches = number_of_mismatches_fast(s1_suffix, s2_prefix, sp_len);

        int score = (sp_len-mismatches)*match_score + mismatches*mismatch_score;

        int max_acceptable_mm = sp_len * max_seq_error;
        if (best_score < score && mismatches <= max_acceptable_mm && mismatches <= max_mismatches) {
            best_score = score;
            best_aln_mismatches = mismatches;
            overlap = sp_len;
        }
    }
    return suffix_prefix_aln_t(overlap, best_score, best_aln_mismatches);
}


std::vector<StripedSmithWaterman::Alignment> get_best_alns(char* contig_seq, hts_pos_t remap_start, hts_pos_t remap_len, char* query, StripedSmithWaterman::Aligner& aligner) {
	
	contig_seq += remap_start;
	
	std::vector<StripedSmithWaterman::Alignment> best_alns;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	aligner.Align(query, contig_seq, remap_len, filter, &aln, 15);
	best_alns.push_back(aln);
	int best_score = aln.sw_score;
	int new_start = aln.ref_begin+5;

	if (new_start >= remap_len) return best_alns;

	aligner.Align(query, contig_seq+new_start, remap_len-new_start, filter, &aln, 15);
	while (aln.sw_score == best_score) {
		aln.ref_begin += new_start;
		aln.ref_end += new_start;
		best_alns.push_back(aln);
		new_start += aln.ref_begin+5;
		if (new_start >= remap_len) break;
		aligner.Align(query, contig_seq+new_start, remap_len-new_start, filter, &aln, 15);
	}
	return best_alns;
}

std::vector<sv_t*> detect_svs_from_aln(StripedSmithWaterman::Alignment& aln, std::string contig_name, hts_pos_t ref_start,
	std::string junction_seq, int match_score = 1, int mismatch_score = -4, int gap_open_score = -6, int gap_extend_score = -1) {
	std::vector<sv_t*> svs;
    hts_pos_t current_pos = ref_start + aln.ref_begin;
	hts_pos_t junction_pos = 0;

	for (int i = 0; i < aln.cigar.size(); i++) {
		uint32_t c = aln.cigar[i];
        int op_length = cigar_int_to_len(c);
        char op = cigar_int_to_op(c);

		hts_pos_t start = current_pos-1, end = (op == 'D') ? current_pos+op_length-1 : current_pos;
		std::string ins_seq = (op == 'I') ? junction_seq.substr(junction_pos, op_length) : "";
		sv_t::anchor_aln_t* left_part_anchor_aln = new sv_t::anchor_aln_t(ref_start+aln.ref_begin, start, junction_pos, aln.sw_score);
		sv_t::anchor_aln_t* right_part_anchor_aln = new sv_t::anchor_aln_t(end, ref_start+aln.ref_end, junction_seq.length()-junction_pos-ins_seq.length(), aln.sw_score);
		if (op == 'D') {
			sv_t* sv = new deletion_t(contig_name, start, end, "", NULL, NULL, left_part_anchor_aln, right_part_anchor_aln);
			svs.push_back(sv);
		} else if (op == 'I') {
			sv_t* sv = new insertion_t(contig_name, current_pos-1, current_pos-1, ins_seq, NULL, NULL, left_part_anchor_aln, right_part_anchor_aln);
			svs.push_back(sv);
		}

        if (op != 'I' && op != 'S') {
            current_pos += op_length;
        }
		if (op != 'D') {
			junction_pos += op_length;
		}
    }

	return svs;
}

std::vector<sv_t*> detect_svs_from_junction(std::string& contig_name, char* contig_seq, std::string junction_seq, hts_pos_t ref_remap_lh_start, hts_pos_t ref_remap_lh_end,
                    hts_pos_t ref_remap_rh_start, hts_pos_t ref_remap_rh_end, StripedSmithWaterman::Aligner& aligner, int min_clip_len) {
    
    hts_pos_t ref_remap_lh_len = ref_remap_lh_end - ref_remap_lh_start;
    hts_pos_t ref_remap_rh_len = ref_remap_rh_end - ref_remap_rh_start;

    char ref_lh_cstr[100000];
	for (int i = 0; i < ref_remap_lh_len; i++) {
		ref_lh_cstr[i] = toupper(contig_seq[ref_remap_lh_start+i]);
	} ref_lh_cstr[ref_remap_lh_len] = '\0';

	int sep = junction_seq.find("-");
	std::string prefix_junction_seq = sep == std::string::npos ? junction_seq : junction_seq.substr(0, sep);
    int* prefix_scores = smith_waterman_gotoh(ref_lh_cstr, ref_remap_lh_len, prefix_junction_seq.c_str(), prefix_junction_seq.length(), 1, -4, -6, -1);

    char ref_rh_cstr[100000];
    for (int i = 0; i < ref_remap_rh_len; i++) {
        ref_rh_cstr[i] = toupper(contig_seq[ref_remap_rh_start+i]);
    } ref_rh_cstr[ref_remap_rh_len] = '\0';

    char ref_rh_cstr_rev[100000];
    for (int i = 0; i < ref_remap_rh_len; i++) {
        ref_rh_cstr_rev[i] = ref_rh_cstr[ref_remap_rh_len-1-i];
    } ref_rh_cstr_rev[ref_remap_rh_len] = '\0';

	std::string suffix_junction_seq = sep == std::string::npos ? junction_seq : junction_seq.substr(sep+1);
    std::string suffix_junction_seq_rev = std::string(suffix_junction_seq.rbegin(), suffix_junction_seq.rend());
    int* suffix_scores = smith_waterman_gotoh(ref_rh_cstr_rev, ref_remap_rh_len, suffix_junction_seq_rev.c_str(), suffix_junction_seq_rev.length(), 1, -4, -6, -1);
	int suffix_begin = sep == std::string::npos ? 0 : sep+1;

    int max_score = 0, best_i = 0, best_j = 0;
	for (int i = min_clip_len; i < prefix_junction_seq.length()-min_clip_len; i++) {
        int prefix_score = prefix_scores[i-1]; // score of the best aln of [0..i-1]
        for (int j = std::max(i, suffix_begin); j <= junction_seq.length()-min_clip_len; j++) {
            int suffix_score = suffix_scores[junction_seq.length()-j-1]; // score of the best aln of [j..junction_seq.length()-1]
            // note that we want the score of the suffix of length junction_seq.length()-j, 
            // so we need to subtract 1 because suffix_scores[n] is the score of the best suffix of length n+1
            if (prefix_score + suffix_score > max_score) {
                max_score = prefix_score + suffix_score;
                best_i = i, best_j = j;
            }
        }
    }

	free(prefix_scores);
	free(suffix_scores);

    if (max_score == 0) return std::vector<sv_t*>();

    std::string left_part = junction_seq.substr(0, best_i);
    std::string middle_part = junction_seq.substr(best_i, best_j-best_i);
    std::string right_part = junction_seq.substr(best_j);

    std::vector<StripedSmithWaterman::Alignment> left_part_alns = get_best_alns(ref_lh_cstr, 0, ref_remap_lh_len, (char*) left_part.c_str(), aligner);
	std::vector<StripedSmithWaterman::Alignment> right_part_alns = get_best_alns(ref_rh_cstr, 0, ref_remap_rh_len, (char*) right_part.c_str(), aligner);

    StripedSmithWaterman::Alignment left_part_aln, right_part_aln;
	int min_size = INT32_MAX;
	for (StripedSmithWaterman::Alignment& _lh_aln : left_part_alns) {
		for (StripedSmithWaterman::Alignment& _rh_aln : right_part_alns) {
			int size = abs((ref_remap_rh_start + _rh_aln.ref_begin - 1) - (ref_remap_lh_start + _lh_aln.ref_end));
			if (size < min_size) {
				left_part_aln = _lh_aln, right_part_aln = _rh_aln;
				min_size = size;
			}
		}
	}

    hts_pos_t left_anchor_start = ref_remap_lh_start + left_part_aln.ref_begin;
    hts_pos_t left_anchor_end = ref_remap_lh_start + left_part_aln.ref_end;
    hts_pos_t right_anchor_start = ref_remap_rh_start + right_part_aln.ref_begin;
    hts_pos_t right_anchor_end = ref_remap_rh_start + right_part_aln.ref_end;
    hts_pos_t left_bp = left_anchor_end, right_bp = right_anchor_start - 1; 
	if (right_bp < 0) right_bp = 0; // this can be -1, if the duplication starts at the beginning of the contig. 
	/* The VCF specification says: "For simple insertions and deletions in which either the REF or one of the ALT alleles 
		would otherwise be null/empty, the REF and ALT Strings must include the base before the event (which must be reflected 
		in the POS field), unless the event occurs at position 1 on the contig in which case it must include the base after 
		the event"
		However, I am not sure what the base "after the event" is in the case of a duplication. I will report 1 in the VCF (meaning 0 here) for now */
	sv_t::anchor_aln_t* left_part_anchor_aln = new sv_t::anchor_aln_t(left_anchor_start, left_anchor_end, left_part.length(), left_part_aln.sw_score);
	sv_t::anchor_aln_t* right_part_anchor_aln = new sv_t::anchor_aln_t(right_anchor_start, right_anchor_end, right_part.length(), right_part_aln.sw_score);

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment full_aln_lh, full_aln_rh;
	aligner.Align(junction_seq.c_str(), contig_seq+ref_remap_lh_start, ref_remap_lh_len, filter, &full_aln_lh, 0);
	aligner.Align(junction_seq.c_str(), contig_seq+ref_remap_rh_start, ref_remap_rh_len, filter, &full_aln_rh, 0);
	StripedSmithWaterman::Alignment& best_full_aln = full_aln_lh.sw_score >= full_aln_rh.sw_score ? full_aln_lh : full_aln_rh;

	int prefix_mh_len = 0;
    if (left_bp > right_bp) { // there is microhomology in the inserted seq or it's a duplication
        int mh_len = left_bp - right_bp;
        std::pair<int, int> lp_suffix_score = find_aln_suffix_score(left_part_aln.cigar, mh_len, 1, -4, -6, -1);
        std::pair<int, int> rp_prefix_score = find_aln_prefix_score(right_part_aln.cigar, mh_len, 1, -4, -6, -1);

		if (middle_part.size() < left_bp - right_bp) {
			if (right_anchor_end - left_anchor_end < min_clip_len || right_anchor_start - left_anchor_start < min_clip_len ||
				(lp_suffix_score.first == mh_len && rp_prefix_score.first == mh_len && middle_part.empty() &&
				!is_right_clipped(left_part_aln) && !is_left_clipped(right_part_aln))) { // it's a duplication
				duplication_t* sv = new duplication_t(contig_name, right_bp, left_bp, middle_part, NULL, NULL, left_part_anchor_aln, right_part_anchor_aln);
				sv->precompute_base_frequencies(contig_seq);
				return std::vector<sv_t*>({sv});
			}
		}

        // otherwise, it's an insertion
        std::string mh;
        if (lp_suffix_score.first > rp_prefix_score.first) { 
            // the suffix of the left part is more similar to the reference, hence we choose the prefix of the right part 
            //to add as part of the inserted sequence
            int right_bp_adjustment = 0;
            int query_mh_bases = query_prefix_len(right_part_aln.cigar, mh_len, right_bp_adjustment);
            mh = right_part.substr(0, query_mh_bases);
            right_bp = left_bp + right_bp_adjustment;
            middle_part = middle_part + mh;
        } else {
            // the prefix of the right part is more similar to the reference, hence we choose the suffix of the left part
            // to add as part of the inserted sequence
            int left_bp_adjustment = 0;
            int query_mh_bases = query_suffix_len(left_part_aln.cigar, mh_len, left_bp_adjustment);
            mh = left_part.substr(left_part.length() - query_mh_bases);
            left_bp = right_bp - left_bp_adjustment;
            middle_part = mh + middle_part;
			prefix_mh_len = mh.length();
        }
    }

    std::vector<sv_t*> svs;
    if (right_bp - left_bp > middle_part.length()) { // length of ALT < REF, deletion
		sv_t* sv = new deletion_t(contig_name, left_bp, right_bp, middle_part, NULL, NULL, left_part_anchor_aln, right_part_anchor_aln);
		sv->precompute_base_frequencies(contig_seq);
        svs.push_back(sv);
    } else { // length of ALT > REF, insertion
		sv_t* sv = new insertion_t(contig_name, left_bp, right_bp, middle_part, NULL, NULL, left_part_anchor_aln, right_part_anchor_aln);
		sv->precompute_base_frequencies(contig_seq);
        svs.push_back(sv);
    }
	svs[0]->mh_len = prefix_mh_len;

	hts_pos_t forbidden_zone_start = std::min(left_anchor_end, right_anchor_start);
	hts_pos_t forbidden_zone_end = std::max(left_anchor_start, right_anchor_end);
	std::vector<sv_t*> extra_svs = detect_svs_from_aln(left_part_aln, contig_name, ref_remap_lh_start, left_part);
	for (sv_t* sv : extra_svs) {
		if (overlap(forbidden_zone_start, forbidden_zone_end, sv->start, sv->end)) {
			delete sv;
		} else {
			svs.push_back(sv);
		}
	}
	extra_svs = detect_svs_from_aln(right_part_aln, contig_name, ref_remap_rh_start, right_part);
	for (sv_t* sv : extra_svs) {
		if (overlap(forbidden_zone_start, forbidden_zone_end, sv->start, sv->end)) {
			delete sv;
		} else {
			svs.push_back(sv);
		}
	}
    return svs;
}

std::vector<sv_t*> detect_svs(std::string& contig_name, char* contig_seq, hts_pos_t contig_len, 
					consensus_t* rc_consensus, consensus_t* lc_consensus,
					StripedSmithWaterman::Aligner& aligner, int min_overlap, int min_clip_len, double max_seq_error) {

	std::string consensus_junction_seq;
	hts_pos_t ref_remap_lh_start = 0, ref_remap_lh_end = 0, ref_remap_rh_start = 0, ref_remap_rh_end = 0;
	int overlap = 0;
	double mismatch_rate = 0.0;

	if (rc_consensus != NULL && lc_consensus != NULL) {
		std::string rc_consensus_seq = rc_consensus->sequence;
		std::string lc_consensus_seq = lc_consensus->sequence;
		suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus_seq, lc_consensus_seq, 1, -4, max_seq_error, min_overlap);
		if (spa.overlap < min_overlap || is_homopolymer(lc_consensus->sequence.c_str(), spa.overlap)) {
			rc_consensus_seq = rc_consensus->sequence.substr(0, rc_consensus->sequence.length()-rc_consensus->lowq_clip_portion);
			lc_consensus_seq = lc_consensus->sequence.substr(lc_consensus->lowq_clip_portion);
			spa = aln_suffix_prefix(rc_consensus_seq, lc_consensus_seq, 1, -4, max_seq_error, min_overlap);
			if (spa.overlap < min_overlap || is_homopolymer(lc_consensus_seq.c_str(), spa.overlap)) {
				return std::vector<sv_t*>();
			}
		}

		overlap = spa.overlap;
		mismatch_rate = spa.mismatch_rate();
		
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment aln_rh, aln_lh;

		consensus_junction_seq = rc_consensus_seq + lc_consensus_seq.substr(spa.overlap);
		ref_remap_lh_start = rc_consensus->breakpoint - consensus_junction_seq.length();
		ref_remap_lh_end = rc_consensus->breakpoint + consensus_junction_seq.length();
		if (ref_remap_lh_start < 0) ref_remap_lh_start = 0;
		if (ref_remap_lh_end > contig_len) ref_remap_lh_end = contig_len;

		ref_remap_rh_start = lc_consensus->breakpoint - consensus_junction_seq.length();
		ref_remap_rh_end = lc_consensus->breakpoint + consensus_junction_seq.length();
		if (ref_remap_rh_start < 0) ref_remap_rh_start = 0;
		if (ref_remap_rh_end > contig_len) ref_remap_rh_end = contig_len;
	} else if (rc_consensus != NULL) {
		consensus_junction_seq = rc_consensus->sequence;

		ref_remap_lh_start = std::max(hts_pos_t(0), rc_consensus->start - (int) rc_consensus->sequence.length());
		ref_remap_lh_end = std::min(rc_consensus->end + (int) rc_consensus->sequence.length(), contig_len);

		ref_remap_rh_start = rc_consensus->remap_boundary - rc_consensus->sequence.length();
		ref_remap_rh_end = rc_consensus->remap_boundary + rc_consensus->sequence.length();
		if (rc_consensus->remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
			ref_remap_rh_start = rc_consensus->breakpoint - 2*rc_consensus->sequence.length();
			ref_remap_rh_end = rc_consensus->breakpoint + 2*rc_consensus->sequence.length();
		}
		ref_remap_rh_start = std::max(ref_remap_rh_start, hts_pos_t(0));
		ref_remap_rh_end = std::min(ref_remap_rh_end, contig_len);
	} else if (lc_consensus != NULL) {
		consensus_junction_seq = lc_consensus->sequence;

		ref_remap_lh_start = lc_consensus->remap_boundary - lc_consensus->sequence.length();
		ref_remap_lh_end = lc_consensus->remap_boundary + lc_consensus->sequence.length();
		if (lc_consensus->remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
			ref_remap_lh_start = lc_consensus->breakpoint - 2*lc_consensus->sequence.length();
			ref_remap_lh_end = lc_consensus->breakpoint + 2*lc_consensus->sequence.length();
		}
		ref_remap_lh_start = std::max(ref_remap_lh_start, hts_pos_t(0));
		ref_remap_lh_end = std::min(ref_remap_lh_end, contig_len);

		ref_remap_rh_start = std::max(hts_pos_t(0), lc_consensus->start - (int) lc_consensus->sequence.length());
		ref_remap_rh_end = std::min(lc_consensus->end + (int) lc_consensus->sequence.length(), contig_len);
	} else {
		return std::vector<sv_t*>();
	}

    std::vector<sv_t*> svs = detect_svs_from_junction(contig_name, contig_seq, consensus_junction_seq, ref_remap_lh_start, ref_remap_lh_end, ref_remap_rh_start, ref_remap_rh_end, aligner, min_clip_len);
	for (sv_t* sv : svs) {
		sv->rc_consensus = rc_consensus;
		sv->lc_consensus = lc_consensus;

		if (rc_consensus != NULL && lc_consensus == NULL) {
			if (rc_consensus->is_hsr) {
				sv->source = "1HSR_RC";
			} else {
				sv->source = "1SR_RC";
			}
		} else if (rc_consensus == NULL && lc_consensus != NULL) {
			if (lc_consensus->is_hsr) {
				sv->source = "1HSR_LC";
			} else {
				sv->source = "1SR_LC";
			}
		} else if (!rc_consensus->is_hsr && !lc_consensus->is_hsr) {
			sv->source = "2SR";
		} else if (rc_consensus->is_hsr && lc_consensus->is_hsr) {
			sv->source = "2HSR";
		} else if (!rc_consensus->is_hsr && lc_consensus->is_hsr) {
			sv->source = "SR-HSR";
		} else if (rc_consensus->is_hsr && !lc_consensus->is_hsr) {
			sv->source = "HSR-SR";
		}
	}

    return svs;
}

breakend_t* detect_bnd(std::string contig_name, char* contig_seq, hts_pos_t contig_len, consensus_t* leftmost_consensus, consensus_t* rightmost_consensus, 
	suffix_prefix_aln_t& spa, StripedSmithWaterman::Aligner& aligner, int min_clip_len) {

	std::string lm_seq = leftmost_consensus->sequence, rm_seq = rightmost_consensus->sequence;
	if (leftmost_consensus->left_clipped) rc(lm_seq);
	else rc(rm_seq);
	std::string full_junction_seq = lm_seq + rm_seq.substr(spa.overlap);

	hts_pos_t ref_remap_lh_start = leftmost_consensus->breakpoint - full_junction_seq.length();
	if (ref_remap_lh_start < 0) ref_remap_lh_start = 0;
	hts_pos_t ref_remap_lh_end = leftmost_consensus->breakpoint + full_junction_seq.length();
	if (ref_remap_lh_end > contig_len) ref_remap_lh_end = contig_len;

	hts_pos_t ref_remap_rh_start = rightmost_consensus->breakpoint - full_junction_seq.length();
	if (ref_remap_rh_start < 0) ref_remap_rh_start = 0;
	hts_pos_t ref_remap_rh_end = rightmost_consensus->breakpoint + full_junction_seq.length();
	if (ref_remap_rh_end > contig_len) ref_remap_rh_end = contig_len;

	if (!leftmost_consensus->left_clipped) {
		int* fwd_prefix_scores = smith_waterman_gotoh(contig_seq+ref_remap_lh_start, ref_remap_lh_end-ref_remap_lh_start, full_junction_seq.c_str(), full_junction_seq.length(), 1, -4, -6, -1);
		rc(full_junction_seq);
		int* revc_prefix_scores = smith_waterman_gotoh(contig_seq+ref_remap_rh_start, ref_remap_rh_end-ref_remap_rh_start, full_junction_seq.c_str(), full_junction_seq.length(), 1, -4, -6, -1);

		int max_score = 0, best_i = 0, best_j = 0;
		for (int i = min_clip_len; i < full_junction_seq.length()-min_clip_len; i++) {
			int fwd_prefix_score = fwd_prefix_scores[i-1]; // score of the best aln of full_junction_seq[0..i-1]
			for (int j = i; j <= full_junction_seq.length()-min_clip_len; j++) {
				int rev_prefix_score = revc_prefix_scores[full_junction_seq.length()-j-1]; // score of the best aln of RC of full_junction_seq[j..junction_seq.length()-1]
				if (fwd_prefix_score + rev_prefix_score >= max_score) {
					max_score = fwd_prefix_score + rev_prefix_score;
					best_i = i, best_j = j;
				}
			}
		}

		free(fwd_prefix_scores);
		free(revc_prefix_scores);

		if (max_score == 0) return NULL;

		rc(full_junction_seq);
		std::string left_part = full_junction_seq.substr(0, best_i);
		std::string middle_part = full_junction_seq.substr(best_i, best_j-best_i);
		std::string right_part = full_junction_seq.substr(best_j);
		rc(right_part);

		std::vector<StripedSmithWaterman::Alignment> left_part_alns = get_best_alns(contig_seq, ref_remap_lh_start, ref_remap_lh_end-ref_remap_lh_start, (char*) left_part.c_str(), aligner);
		std::vector<StripedSmithWaterman::Alignment> right_part_alns = get_best_alns(contig_seq, ref_remap_rh_start, ref_remap_rh_end-ref_remap_rh_start, (char*) right_part.c_str(), aligner);
		StripedSmithWaterman::Alignment left_part_aln = left_part_alns[left_part_alns.size()-1], right_part_aln = right_part_alns[0];

		sv_t::anchor_aln_t* left_anchor_aln = new sv_t::anchor_aln_t(ref_remap_lh_start+left_part_aln.ref_begin, ref_remap_lh_start+left_part_aln.ref_end, left_part.length(), left_part_aln.sw_score);
		sv_t::anchor_aln_t* right_anchor_aln = new sv_t::anchor_aln_t(ref_remap_rh_start+right_part_aln.ref_begin, ref_remap_rh_start+right_part_aln.ref_end, right_part.length(), right_part_aln.sw_score);

		hts_pos_t start = ref_remap_lh_start + left_part_aln.ref_end, end = ref_remap_rh_start + right_part_aln.ref_end;
		return new breakend_t(contig_name, start, end, middle_part, leftmost_consensus, rightmost_consensus, left_anchor_aln, right_anchor_aln, '-');
	} else {

		rc(full_junction_seq);
		std::string full_junction_seq_rev = std::string(full_junction_seq.rbegin(), full_junction_seq.rend());
		char* ref_remap_lh_rev = new char[ref_remap_lh_end-ref_remap_lh_start+1];
		for (int i = 0; i < ref_remap_lh_end-ref_remap_lh_start; i++) {
			ref_remap_lh_rev[i] = std::toupper(contig_seq[ref_remap_lh_end-1-i]);
		} ref_remap_lh_rev[ref_remap_lh_end-ref_remap_lh_start] = '\0';
		int* revc_suffix_scores = smith_waterman_gotoh(ref_remap_lh_rev, ref_remap_lh_end-ref_remap_lh_start, full_junction_seq_rev.c_str(), full_junction_seq_rev.length(), 1, -4, -6, -1);
		
		rc(full_junction_seq);
		full_junction_seq_rev = std::string(full_junction_seq.rbegin(), full_junction_seq.rend());
		char* ref_remap_rh_rev = new char[ref_remap_rh_end-ref_remap_rh_start+1];
		for (int i = 0; i < ref_remap_rh_end-ref_remap_rh_start; i++) {
			ref_remap_rh_rev[i] = std::toupper(contig_seq[ref_remap_rh_end-1-i]);
		} ref_remap_rh_rev[ref_remap_rh_end-ref_remap_rh_start] = '\0';
		int* fwd_suffix_scores = smith_waterman_gotoh(ref_remap_rh_rev, ref_remap_rh_end-ref_remap_rh_start, full_junction_seq_rev.c_str(), full_junction_seq_rev.length(), 1, -4, -6, -1);

		int max_score = 0, best_i = 0, best_j = 0;
		for (int i = min_clip_len; i < full_junction_seq.length()-min_clip_len; i++) {
			int revc_suffix_score = revc_suffix_scores[i-1]; // score of the best aln of the RC of full_junction_seq[i..n]
			for (int j = i; j <= full_junction_seq.length()-min_clip_len; j++) {
				int fwd_suffix_score = fwd_suffix_scores[full_junction_seq.length()-j-1]; // score of the best aln of [j..junction_seq.length()-1]
				if (fwd_suffix_score + revc_suffix_score > max_score) {
					max_score = fwd_suffix_score + revc_suffix_score;
					best_i = i, best_j = j;
				}
			}
		}

		free(fwd_suffix_scores);
		free(revc_suffix_scores);
		delete[] ref_remap_lh_rev;
		delete[] ref_remap_rh_rev;

		if (max_score == 0) return NULL;

		std::string left_part = full_junction_seq.substr(0, best_i);
		rc(left_part);
		std::string middle_part = full_junction_seq.substr(best_i, best_j-best_i);
		std::string right_part = full_junction_seq.substr(best_j);

		std::vector<StripedSmithWaterman::Alignment> left_part_alns = get_best_alns(contig_seq, ref_remap_lh_start, ref_remap_lh_end-ref_remap_lh_start, (char*) left_part.c_str(), aligner);
		std::vector<StripedSmithWaterman::Alignment> right_part_alns = get_best_alns(contig_seq, ref_remap_rh_start, ref_remap_rh_end-ref_remap_rh_start, (char*) right_part.c_str(), aligner);
		StripedSmithWaterman::Alignment left_part_aln = left_part_alns[left_part_alns.size()-1], right_part_aln = right_part_alns[0];

		sv_t::anchor_aln_t* left_anchor_aln = new sv_t::anchor_aln_t(ref_remap_lh_start+left_part_aln.ref_begin, ref_remap_lh_start+left_part_aln.ref_end, left_part.length(), left_part_aln.sw_score);
		sv_t::anchor_aln_t* right_anchor_aln = new sv_t::anchor_aln_t(ref_remap_rh_start+right_part_aln.ref_begin, ref_remap_rh_start+right_part_aln.ref_end, right_part.length(), right_part_aln.sw_score);

		hts_pos_t start = ref_remap_lh_start + left_part_aln.ref_begin-1, end = ref_remap_rh_start + right_part_aln.ref_begin-1;
		return new breakend_t(contig_name, start, end, middle_part, leftmost_consensus, rightmost_consensus, left_anchor_aln, right_anchor_aln, '+');
	}
}

bool accept(StripedSmithWaterman::Alignment& aln, int min_clip_len, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);

	bool lc_is_noise = false, rc_is_noise = false;
	if (!qual_ascii.empty()) {
		if (lc_size > min_clip_len && lc_size < qual_ascii.length()/2) {
			std::string lc_qual = qual_ascii.substr(0, lc_size);
			std::string not_lc_qual = qual_ascii.substr(lc_size);
			lc_is_noise = avg_qual(lc_qual) < min_avg_base_qual && avg_qual(not_lc_qual) >= min_avg_base_qual;
		}
		if (rc_size > min_clip_len && rc_size < qual_ascii.length()/2) {
			std::string not_rc_qual = qual_ascii.substr(0, qual_ascii.length()-rc_size);
			std::string rc_qual = qual_ascii.substr(qual_ascii.length()-rc_size);
			rc_is_noise = avg_qual(rc_qual) < min_avg_base_qual && avg_qual(not_rc_qual) >= min_avg_base_qual;
		}
	}

	double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
	return (lc_size <= min_clip_len || lc_is_noise) && (rc_size <= min_clip_len || rc_is_noise) && mismatch_rate <= max_seq_error;
}

std::pair<int, int> get_dels_ins_in_first_n_chars(StripedSmithWaterman::Alignment& aln, int n) {
    int dels = 0, inss = 0;
    int offset = 0;
    for (uint32_t c : aln.cigar) {
        int len = cigar_int_to_len(c);
        char op = cigar_int_to_op(c);

        // since we are "unrolling" soft-clipped bases, they must be accounted for
        bool consumes_ref = op == '=' || op == 'X' || op == 'D' || op == 'S';
        if (consumes_ref && offset + len > n) {
            len = n-offset;
        }

        if (op == 'D') {
            dels += len;
        } else if (op == 'I') {
            inss += len;
        }

        if (consumes_ref) {
            offset += len;
            if (offset == n) break;
        }
    }
    return {dels, inss};
}

hts_pos_t get_start_offset(StripedSmithWaterman::Alignment& r1, StripedSmithWaterman::Alignment& r2) {
	int r1_uc_start = r1.ref_begin - get_left_clip_size(r1);
	int r2_uc_start = r2.ref_begin - get_left_clip_size(r2);
    hts_pos_t offset = r2_uc_start - r1_uc_start;
    /* suppose the difference in starting position between R1 and R2 is N, we may be tempted to align the start
     * of R2 to position N+1 of R1
     * However, if R1 has insertions or deletions in the first N bps, the difference in starting positions
     * is no longer accurate to decide where the start of R2 aligns on R1.
     * If R1 has I bps inserted in the first N bps, then R2 will align to position N+1+I.
     * Conversely, if D bps are deleted, R2 will align to position N+1-D
     */
    std::pair<int, int> del_ins = get_dels_ins_in_first_n_chars(r1, offset);
    return offset + del_ins.second - del_ins.first;
}

#endif // SW_UTILS_H
