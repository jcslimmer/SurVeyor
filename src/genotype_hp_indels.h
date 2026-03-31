#ifndef GENOTYPE_HP_INDELS_H
#define GENOTYPE_HP_INDELS_H

#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "utils.h"
#include "stat_tests.h"
#include "../libs/ssw_cpp.h"

#include "genotype.h"
#include "var_utils.h"

struct hp_read_info_t {
    int hp_len;
    int tail_5p_len, tail_3p_len;
    int tail_5p_mismatches, tail_3p_mismatches;
};


// Find the mode of the HP lengths among the good reads
// Good reads = <20% mismatches on both tails
int find_good_reads_hp_len_mode(const std::vector<hp_read_info_t>& hp_read_infos) {
    std::unordered_map<int, int> hp_len_counts;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if ((hp_read_info.tail_5p_len == 0 || double(hp_read_info.tail_5p_mismatches)/hp_read_info.tail_5p_len <= 0.2) &&
            (hp_read_info.tail_3p_len == 0 || double(hp_read_info.tail_3p_mismatches)/hp_read_info.tail_3p_len <= 0.2)) {
            hp_len_counts[hp_read_info.hp_len]++;
        }
    }

    int mode_hp_len = -1, max_count = 0;
    for (const auto& kv : hp_len_counts) {
        if (kv.second > max_count) {
            mode_hp_len = kv.first;
            max_count = kv.second;
        }
    }
    return mode_hp_len;
}

std::vector<std::vector<hp_read_info_t>> cluster_reads_by_two_modes(const std::vector<hp_read_info_t>& hp_read_infos) {
    int n = hp_read_infos.size();
    if (n == 0) {
        return {{}, {}};
    }

    if (n == 1) { // Return a single cluster with the one read
        return {{hp_read_infos[0]}};
    }
    
    std::unordered_map<int, int> counts;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        counts[hp_read_info.hp_len]++;
    }

    if (counts.size() == 1) { // All reads have the same HP length, so we form a single cluster
        return {hp_read_infos};
    }
    
    std::vector<std::vector<hp_read_info_t>> clusters;

    // Pick the two most common observed HP lengths, breaking count ties toward the smaller length.
    std::vector<std::pair<int, int>> count_items(counts.begin(), counts.end());
    std::sort(count_items.begin(), count_items.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });

    int a_mode = count_items[0].first, b_mode = count_items[1].first;
    if (a_mode > b_mode) std::swap(a_mode, b_mode);

    // Absorb every read into the nearest mode; distance ties go to the larger cluster.
    std::vector<int> a_idx, b_idx;
    for (int i = 0; i < n; i++) {
        int x = hp_read_infos[i].hp_len;
        int a_dist = abs(x - a_mode);
        int b_dist = abs(x - b_mode);
        if (a_dist < b_dist) {
            a_idx.push_back(i);
        } else if (b_dist < a_dist) {
            b_idx.push_back(i);
        } else if (counts[a_mode] >= counts[b_mode]) {
            a_idx.push_back(i);
        } else {
            b_idx.push_back(i);
        }
    }

    std::vector<hp_read_info_t> a_reads, b_reads;
    for (int i : a_idx) a_reads.push_back(hp_read_infos[i]);
    for (int i : b_idx) b_reads.push_back(hp_read_infos[i]);
    return {a_reads, b_reads};
}


// For clipped tails, place the tail at its best ungapped offset in the local reference window.
int best_ungapped_mismatch_count(const std::string& query, char* ref, hts_pos_t ref_len) {
    if (query.empty() || ref_len < (hts_pos_t) query.length()) {
        return query.length();
    }

    int best_mismatches = query.length();
    for (hts_pos_t i = 0; i + query.length() <= ref_len; i++) {
        int mismatches = 0;
        for (int j = 0; j < query.length(); j++) {
            if (std::toupper(query[j]) != std::toupper(ref[i + j])) {
                mismatches++;
            }
        }
        if (mismatches < best_mismatches) {
            best_mismatches = mismatches;
            if (best_mismatches == 0) break;
        }
    }
    return best_mismatches;
}

// If the tail is unclipped, calculate mismatches by simply counting mismatches in its alignment
// If the tail is clipped, force the whole tail to align next to the HP ref region (allowing some extra leeway)
int tail_mismatch_count_simple(bam1_t* read, const std::string& read_seq, const std::vector<hts_pos_t>& qpos_to_rpos,
    char* contig_seq, hts_pos_t contig_len, int q_lo, int q_hi, bool left_side, int leeway, hts_pos_t ref_boundary) {

    if (q_hi <= q_lo || read_seq.empty()) return 0;

    int left_clip = get_left_clip_size(read), right_clip = get_right_clip_size(read);
    if (left_side && left_clip > 0 && q_lo == 0) {
        hts_pos_t ref_hi = std::min(ref_boundary, contig_len);
        hts_pos_t ref_lo = std::max((hts_pos_t) 0, ref_hi - (q_hi - q_lo) - leeway);
        return best_ungapped_mismatch_count(read_seq.substr(q_lo, q_hi - q_lo), contig_seq+ref_lo, ref_hi-ref_lo);
    }
    if (!left_side && right_clip > 0 && q_hi == (int) read_seq.length()) {
        hts_pos_t ref_lo = std::max((hts_pos_t) 0, std::min(ref_boundary, contig_len));
        hts_pos_t ref_hi = std::min(contig_len, ref_lo + (q_hi - q_lo) + leeway);
        return best_ungapped_mismatch_count(read_seq.substr(q_lo, q_hi - q_lo), contig_seq+ref_lo, ref_hi-ref_lo);
    }

    int mismatches = 0;
    for (int qpos = q_lo; qpos < q_hi; qpos++) {
        hts_pos_t rpos = qpos_to_rpos[qpos];
        if (rpos == -1 || rpos >= contig_len) continue;
        if (std::toupper(read_seq[qpos]) != std::toupper(contig_seq[rpos])) {
            mismatches++;
        }
    }
    return mismatches;
}


hp_read_info_t calculate_hp_read_info(bam1_t* read, hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len) {
    if (read == NULL || is_unmapped(read) || !is_primary(read) || read->core.l_qseq <= 0) {
        return {0, 0, 0, 0, 0};
    }

    const int tail_align_leeway = 10;
    std::string read_seq = get_sequence(read);
    std::vector<hts_pos_t> qpos_to_rpos(read->core.l_qseq, -1);
    std::vector<int> anchors;

    uint32_t* cigar = bam_get_cigar(read);
    int qpos = 0;
    hts_pos_t rpos = read->core.pos;
    hts_pos_t last_qpos_before_ref_hp = -1, first_qpos_after_ref_hp = -1;
    for (uint32_t i = 0; i < read->core.n_cigar; i++) {
        char opchar = bam_cigar_opchr(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);

        if (opchar == 'M' || opchar == '=' || opchar == 'X') {
            for (int j = 0; j < oplen; j++) {
                qpos_to_rpos[qpos] = rpos;
                if (rpos < ref_hp_range.beg) {
                    last_qpos_before_ref_hp = qpos;
                } else if (rpos >= ref_hp_range.end && first_qpos_after_ref_hp == -1) {
                    first_qpos_after_ref_hp = qpos;
                } else if (ref_hp_range.beg <= rpos && rpos < ref_hp_range.end) {
                    anchors.push_back(qpos);
                }
                qpos++;
                rpos++;
            }
        } else if (opchar == 'I') {
            if (rpos < ref_hp_range.beg) {
                last_qpos_before_ref_hp = qpos + oplen - 1;
            } else if (rpos >= ref_hp_range.end && first_qpos_after_ref_hp == -1) {
                first_qpos_after_ref_hp = qpos;
            }
            qpos += oplen;
        } else if (opchar == 'S') {
            qpos += oplen;
        } else if (opchar == 'D' || opchar == 'N') {
            rpos += oplen;
        }
    }

    // If no query base lands inside the reference HP, treat this as a 0-bp run
    // and derive the tails from the nearest mapped flanks on each side.
    if (anchors.empty()) {
        int left_tail_len, right_tail_len;
        if (last_qpos_before_ref_hp != -1 && first_qpos_after_ref_hp == -1) {
            left_tail_len = last_qpos_before_ref_hp + 1;
            right_tail_len = read->core.l_qseq - left_tail_len;
        } else if (last_qpos_before_ref_hp == -1 && first_qpos_after_ref_hp != -1) {
            left_tail_len = first_qpos_after_ref_hp;
            right_tail_len = read->core.l_qseq - first_qpos_after_ref_hp;
        } else if (last_qpos_before_ref_hp != -1 && first_qpos_after_ref_hp != -1) {
            left_tail_len = last_qpos_before_ref_hp + 1;
            right_tail_len = read->core.l_qseq - first_qpos_after_ref_hp;
        } else {
            return {0, 0, 0, 0, 0};
        }

        int left_mismatches = tail_mismatch_count_simple(read, read_seq, qpos_to_rpos,
            contig_seq, contig_len, 0, left_tail_len, true, tail_align_leeway, ref_hp_range.beg);
        int right_mismatches = tail_mismatch_count_simple(read, read_seq, qpos_to_rpos,
            contig_seq, contig_len, read->core.l_qseq - right_tail_len, read->core.l_qseq, false, tail_align_leeway, ref_hp_range.end);

        if (!bam_is_rev(read)) {
            return {0, left_tail_len, right_tail_len, left_mismatches, right_mismatches};
        } else {
            return {0, right_tail_len, left_tail_len, right_mismatches, left_mismatches};
        }
    }

    // Start from the aligned bases inside the reference HP and expand through
    // adjacent query bases that still look like part of the same HP run.
    int left = anchors.front();
    while (left > 0 && read_seq[left-1] == hp_base) {
        hts_pos_t mapped_rpos = qpos_to_rpos[left-1];
        if (mapped_rpos != -1 && mapped_rpos < ref_hp_range.beg) break;
        left--;
    }

    int right = anchors.back() + 1;
    while (right < read->core.l_qseq && read_seq[right] == hp_base) {
        hts_pos_t mapped_rpos = qpos_to_rpos[right];
        if (mapped_rpos != -1 && mapped_rpos >= ref_hp_range.end) break;
        right++;
    }

    int left_len = left;
    int right_len = read->core.l_qseq - right;

    hp_read_info_t hp_read_info;
    hp_read_info.hp_len = right - left;
    int left_mismatches = tail_mismatch_count_simple(read, read_seq, qpos_to_rpos,
        contig_seq, contig_len, 0, left, true, tail_align_leeway, ref_hp_range.beg);
    int right_mismatches = tail_mismatch_count_simple(read, read_seq, qpos_to_rpos,
        contig_seq, contig_len, right, read->core.l_qseq, false, tail_align_leeway, ref_hp_range.end);
    if (bam_is_rev(read)) {
        hp_read_info.tail_5p_len = right_len;
        hp_read_info.tail_3p_len = left_len;
        hp_read_info.tail_5p_mismatches = right_mismatches;
        hp_read_info.tail_3p_mismatches = left_mismatches;
    } else {
        hp_read_info.tail_5p_len = left_len;
        hp_read_info.tail_3p_len = right_len;
        hp_read_info.tail_5p_mismatches = left_mismatches;
        hp_read_info.tail_3p_mismatches = right_mismatches;
    }

    return hp_read_info;
}


void genotype_hp_indels_group(std::vector<sv_t*>& hp_indels, hts_pair_pos_t ref_hp_range, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len, 
    stats_t& stats, config_t& config) {

    if (hp_indels.empty()) return;

    char hp_base = get_homopolymer_base(hp_indels[0], contig_seq);

    std::vector<hp_read_info_t> hp_read_infos;

    std::stringstream region_ss;
    region_ss << hp_indels[0]->chr << ":" << ref_hp_range.beg-1 << "-" << ref_hp_range.end+1;
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region_ss.str().c_str());

    // For each read overlapping the reference HP, calculate its observed HP length,
    // 5'/3' tail lengths and 5'/3' tail mismatch counts.
    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (!is_proper_pair(read, stats.min_is, stats.max_is)) continue;

        hp_read_info_t hp_read_info = calculate_hp_read_info(read, ref_hp_range, hp_base, contig_seq, contig_len);
        if (hp_read_info.tail_3p_len < config.min_clip_len || hp_read_info.tail_5p_len < config.min_clip_len) continue;
        hp_read_infos.push_back(hp_read_info);
    }
    bam_destroy1(read);
    hts_itr_destroy(iter);

    // Cluster reads into up to two clusters by their observed HP lengths
    std::vector<std::vector<hp_read_info_t>> clusters = cluster_reads_by_two_modes(hp_read_infos);

    std::vector<int> allele_lens;
    for (sv_t* hp_indel : hp_indels) {
        allele_lens.push_back(ref_hp_range.end - ref_hp_range.beg + hp_indel->svlen());
    }
    allele_lens.push_back(ref_hp_range.end - ref_hp_range.beg); // allele_lens[hp_indels.size()] corresponds to the reference allele

    // Associate each cluster to the most likely allele based on its mode HP length
    
}

// hp_indels are guaranteed to be on the same chromosome
void genotype_hp_indels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<sv_t*> hp_indels, 
    stats_t stats, config_t config, bam_pool_t* bam_pool) {

    std::unordered_map<hts_pos_t, std::vector<sv_t*>> hp_indels_by_ref_hp_range;
    for (sv_t* hp_indel : hp_indels) {
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indel, contig_seq, contig_len);
        hp_indels_by_ref_hp_range[ref_hp_range.beg].push_back(hp_indel);
    }

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    for (auto& kv : hp_indels_by_ref_hp_range) {
        std::vector<sv_t*>& hp_indels_in_range = kv.second;
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indels_in_range[0], contig_seq, contig_len);
        genotype_hp_indels_group(hp_indels_in_range, ref_hp_range, bam_file, contig_seq, contig_len, stats, config);
    }
}

#endif // GENOTYPE_HP_INDELS_H
