#ifndef SW_UTILS_H
#define SW_UTILS_H

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "simd_macros.h"

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


size_t gcd(size_t a, size_t b) {
    return (b == 0) ? a : gcd(b, a % b);
}
size_t lcm(size_t a, size_t b) {
    return a * b / gcd(a, b);
}

int* smith_waterman_gotoh(const char* ref, int ref_len, const char* read, int read_len,
		int match_score, int mismatch_penalty, int gap_open, int gap_extend) {

	const int INF = 1000000;

	const int BYTES_PER_BLOCK = INT_PER_BLOCK * 4;

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
			auto cmp = CMP_INT(H_curr_v, gap_open_v_pos);
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

std::pair<int, int> find_aln_prefix_score(std::vector<uint32_t> cigar, int ref_prefix_len, int match_score, int mismatch_score,
						  int gap_open_score, int gap_extend_score) {
	int score = 0, query_prefix_len = 0;
	for (int i = 0, j = 0; i < cigar.size() && j < ref_prefix_len; i++) {
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

#endif // SW_UTILS_H