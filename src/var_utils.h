#ifndef VAR_UTILS_H
#define VAR_UTILS_H

#include "types.h"

bool is_homopolymer_indel(sv_t* sv, char* chr_seq) {
	if (sv->svtype() == "DEL") {
		// for deletions, we need
		// 1) empty inserted sequence
		// 2) all deleted bases are the same, and A/C/G/T
		if (!sv->ins_seq.empty() || sv->start >= sv->end) return false;

		char target_base = toupper(chr_seq[sv->start+1]);
		if (target_base != 'A' && target_base != 'C' && target_base != 'G' && target_base != 'T') {
			return false;
		}

		for (int i = sv->start+2; i <= sv->end; i++) {
			if (toupper(chr_seq[i]) != target_base) {
				return false;
			}
		}
		return true;
	} else if (sv->svtype() == "INS") {
		// for insertions, we need
		// 1) start == end
		// 2) all inserted bases are the same, and A/C/G/T
		if (sv->start != sv->end || sv->ins_seq.empty()) return false;
		
		char target_base = toupper(sv->ins_seq[0]);
		if (target_base != 'A' && target_base != 'C' && target_base != 'G' && target_base != 'T') {
			return false;
		}

		for (char c : sv->ins_seq) {
			if (toupper(c) != target_base) {
				return false;
			}
		}
		return true;
	}
	return false;
}

char get_homopolymer_base(sv_t* sv, char* chr_seq) {
	if (!is_homopolymer_indel(sv, chr_seq)) {
		return '\0';
	}
	if (sv->svtype() == "DEL") {
		return toupper(chr_seq[sv->start+1]);
	} else if (sv->svtype() == "INS") {
		return toupper(sv->ins_seq[0]);
	}
	return '\0';
}

// The returned range is 0-based, half-open [beg, end)
hts_pair_pos_t find_ref_hp_range_for_indel(sv_t* sv, char* chr_seq, hts_pos_t chr_len) {
	char hp_base = get_homopolymer_base(sv, chr_seq);
	if (hp_base == '\0') {
		throw std::runtime_error("SV " + sv->id + " is not a homopolymer indel.");
	}

	hts_pair_pos_t hp_range;
	hp_range.beg = sv->start;
	while (hp_range.beg >= 0 && toupper(chr_seq[hp_range.beg]) == hp_base) {
		hp_range.beg--;
	}
	hp_range.beg++; // move back to the first base of the homopolymer run

	hp_range.end = sv->end + 1;
	while (hp_range.end < chr_len && toupper(chr_seq[hp_range.end]) == hp_base) {
		hp_range.end++;
	}

	return hp_range;
}

#endif // VAR_UTILS_H