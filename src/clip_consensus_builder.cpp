#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#include "htslib/sam.h"
#include "sam_utils.h"
#include "hsr_utils.h"
#include "utils.h"
#include "assemble.h"
#include "../libs/cptl_stl.h"

std::mutex mtx;
config_t config;
std::string workdir, workspace;
chr_seqs_map_t contigs;

struct bam_redux_t {
    static const uint8_t IS_REV = 1, IS_MREV = 2, IS_INTER_CHR = 4, IS_MATE_LC = 8, IS_MATE_RC = 16;

    hts_pos_t start, end, mstart, isize;
    int left_clip_size, right_clip_size;
    int nm = 0;
    uint8_t flag = 0, mapq = 0;
    std::vector<uint8_t> seq;
    std::vector<uint8_t> qual;
    std::vector<uint32_t> cigar;
    int as = 0;

    bam_redux_t() {}
    bam_redux_t(bam1_t* read) : start(read->core.pos), end(bam_endpos(read)), mstart(read->core.mpos),
        isize(read->core.isize), mapq(read->core.qual), left_clip_size(get_left_clip_size(read)), right_clip_size(get_right_clip_size(read)),
		nm(bam_aux2i(bam_aux_get(read, "NM"))), as(bam_aux2i(bam_aux_get(read, "AS"))) {

        if (bam_is_rev(read)) flag |= IS_REV;
        if (bam_is_mrev(read)) flag |= IS_MREV;
        if (!is_samechr(read) || is_unmapped(read) != is_mate_unmapped(read)) flag |= IS_INTER_CHR;
        if (is_mate_left_clipped(read)) flag |= IS_MATE_LC;
        if (is_mate_right_clipped(read)) flag |= IS_MATE_RC;

        uint8_t* seq_array = bam_get_seq(read);
        seq = std::vector<uint8_t>(seq_array, seq_array+(read->core.l_qseq+1)/2);

        uint8_t* qual_array = bam_get_qual(read);
        qual = std::vector<uint8_t>(qual_array, qual_array+read->core.l_qseq);

        uint32_t* cigar_array = bam_get_cigar(read);
        cigar = std::vector<uint32_t>(cigar_array, cigar_array+read->core.n_cigar);
    }

    int seq_len() {
        return qual.size();
    }

    bool is_rev() {
        return flag & IS_REV;
    }
    bool is_mrev() {
        return flag & IS_MREV;
    }
    bool is_inter_chr() {
        return flag & IS_INTER_CHR;
    }
    bool mate_left_clipped() {
        return flag & IS_MATE_LC;
    }
    bool mate_right_clipped() {
        return flag & IS_MATE_RC;
    }

    hts_pos_t unclipped_start() {
        return start-left_clip_size;
    }
    hts_pos_t unclipped_end() {
        return end+right_clip_size;
    }

    std::string get_sequence() {
        std::string seq_str;
        for (int i = 0; i < seq_len(); i++) {
            seq_str += get_base(seq.data(), i);
        }
        return seq_str;
    }

    std::string cigar_string() {
        std::stringstream ss;
        for (uint32_t c : cigar) ss << bam_cigar_oplen(c) << bam_cigar_opchr(c);
        return ss.str();
    }
};

std::pair<int, int> get_dels_ins_in_first_n_chars(std::vector<uint32_t>& cigar, int n) {

    if (n < 0) return {0, 0};

    int dels = 0, inss = 0;
    int offset = 0;
    for (uint32_t c : cigar) {
        int len = bam_cigar_oplen(c);
        char op = bam_cigar_opchr(c);

        // since we are "unrolling" soft-clipped bases, they must be accounted for
        bool consumes_ref = bam_cigar_type(c) & 2 || op == 'S';
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

hts_pos_t get_start_offset(bam_redux_t* r1, bam_redux_t* r2) {
    hts_pos_t offset = r2->unclipped_start() - r1->unclipped_start();
    /* suppose the difference in starting position between R1 and R2 is N, we may be tempted to align the start
     * of R2 to position N+1 of R1
     * However, if R1 has insertions or deletions in the first N bps, the difference in starting positions
     * is no longer accurate to decide where the start of R2 aligns on R1.
     * If R1 has I bps inserted in the first N bps, then R2 will align to position N+1+I.
     * Conversely, if D bps are deleted, R2 will align to position N+1-D
     */
    std::pair<int, int> del_ins = get_dels_ins_in_first_n_chars(r1->cigar, offset);
    return offset + del_ins.second - del_ins.first;
}

hts_pos_t get_end_offset(bam_redux_t* r1, bam_redux_t* r2) {
    hts_pos_t offset = r2->unclipped_end() - r1->unclipped_end();
    std::vector<uint32_t> rev_cigar(r2->cigar.rbegin(), r2->cigar.rend());
    std::pair<int, int> del_ins = get_dels_ins_in_first_n_chars(rev_cigar, offset);
    return offset + del_ins.second - del_ins.first;
}

// mismatch_score, gap_open_score and gap_extend_score must be negative
int compute_read_score(bam_redux_t* r, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
	std::vector<uint32_t>& cigar = r->cigar;
	int matches = 0;
	int mismatches = r->nm;
	int score = 0;
	for (int i = 0; i < cigar.size(); i++) {
		char op = bam_cigar_opchr(cigar[i]);
		int len = bam_cigar_oplen(cigar[i]);
		if (op == 'M') matches += bam_cigar_oplen(cigar[i]);
		else if (op == 'I' || op == 'D') {
			score += gap_open_score + (len-1)*gap_extend_score;
			mismatches -= len;
		}
	}
	score += match_score*matches + mismatch_score*mismatches - match_score*mismatches;
	return score;
}

std::vector<hts_pos_t> get_read_start_offsets(std::vector<bam_redux_t*>& reads, bool left_clipped) {
    std::vector<hts_pos_t> read_start_offsets;
    int smallest_unclipped_len_i = 0;
    for (int i = 1; i < reads.size(); i++) {
        if (reads[i]->unclipped_end() < reads[smallest_unclipped_len_i]->unclipped_end()) {
            smallest_unclipped_len_i = i;
        }
    }
    for (bam_redux_t* r : reads) {
        if (left_clipped) {
            hts_pos_t offset = get_end_offset(reads[smallest_unclipped_len_i], r);
            offset -= r->seq_len() - reads[smallest_unclipped_len_i]->seq_len();
            read_start_offsets.push_back(offset);
        } else {
            read_start_offsets.push_back(get_start_offset(reads[0], r));
        }
    }
    return read_start_offsets;
}

struct base_score_t {
    int freq = 0, qual = 0;
    char base;

    base_score_t(char base) : base(base) {}
};
bool operator < (const base_score_t& bs1, const base_score_t& bs2) {
    if (bs1.freq != bs2.freq) return bs1.freq < bs2.freq;
    return bs1.qual < bs2.qual;
}

std::string build_full_consensus_seq(std::vector<std::string>& seqs, std::vector<uint8_t*>& quals, std::vector<hts_pos_t>& read_start_offsets) {

    // if not already sorted, sort by start offset
    if (!std::is_sorted(read_start_offsets.begin(), read_start_offsets.end())) {
        std::vector<int> order(seqs.size());
        for (int i = 0; i < order.size(); i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&read_start_offsets](int i1, int i2) {
            return read_start_offsets[i1] < read_start_offsets[i2];
        });
        std::vector<std::string> sorted_seqs;
        std::vector<uint8_t*> sorted_quals;
        std::vector<hts_pos_t> sorted_read_start_offsets;
        for (int i : order) {
            sorted_seqs.push_back(seqs[i]);
            sorted_quals.push_back(quals[i]);
            sorted_read_start_offsets.push_back(read_start_offsets[i]);
        }
        seqs = sorted_seqs;
        quals = sorted_quals;
        read_start_offsets = sorted_read_start_offsets;
    }

    hts_pos_t consensus_len = 0;
    for (int i = 0; i < read_start_offsets.size(); i++) {
        if (consensus_len < read_start_offsets[i] + seqs[i].length()) {
            consensus_len = read_start_offsets[i] + seqs[i].length();
        }
    }
    if (consensus_len > 10000) {
        std::cerr << "WARNING: consensus sequence too long: " << consensus_len << " ";
        std::cerr << "with " << seqs.size() << " reads" << std::endl;
        return "";
    }
    std::string consensus(consensus_len, 'N');

    std::vector<hts_pos_t> read_end_offsets;
    for (int i = 0; i < seqs.size(); i++) {
        hts_pos_t start_offset = read_start_offsets[i];
        hts_pos_t end_offset = start_offset + seqs[i].length() - 1;
        read_end_offsets.push_back(end_offset);
    }

    int s = 0;
    for (int i = 0; i < consensus_len; i++) {
        while (s < seqs.size() && read_end_offsets[s] < i) s++;

        base_score_t base_scores[4] = { base_score_t('A'), base_score_t('C'), base_score_t('G'), base_score_t('T') };
        for (int j = s; j < seqs.size() && read_start_offsets[j] <= i; j++) {
            if (read_end_offsets[j] < i) continue;

            char nucl = seqs[j][i - read_start_offsets[j]];
            uint8_t qual = quals[j][i - read_start_offsets[j]];
            if (nucl == 'A') {
                base_scores[0].freq++;
                base_scores[0].qual += qual;
            } else if (nucl == 'C') {
                base_scores[1].freq++;
                base_scores[1].qual += qual;
            } else if (nucl == 'G') {
                base_scores[2].freq++;
                base_scores[2].qual += qual;
            } else if (nucl == 'T') {
                base_scores[3].freq++;
                base_scores[3].qual += qual;
            }
        }

        base_score_t best_base_score = max(base_scores[0], base_scores[1], base_scores[2], base_scores[3]);

        consensus[i] = best_base_score.base;
    }
    return consensus;
}

// Use kmers to select reads that are likely to be part of the same haplotype
std::vector<int> select_reads_by_kmer(std::vector<std::string>& seqs, std::vector<hts_pos_t>& read_start_offsets) {

    const int K = sizeof(uint32_t)*8/2; // 16-mers

    uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

    int consensus_len = 0;
    for (int i = 0; i < read_start_offsets.size(); i++) {
        if (consensus_len < read_start_offsets[i] + seqs[i].length()) {
            consensus_len = read_start_offsets[i] + seqs[i].length();
        }
    }
    std::vector<std::vector<std::pair<uint32_t, int>>> kmer_counts_by_pos(consensus_len);

    for (int i = 0; i < seqs.size(); i++) {
        std::string& seq = seqs[i];
        if (seq.length() < K) continue;

        uint32_t kmer = 0;
        for (int j = 0; j < seq.length(); j++) {
            kmer = ((kmer << 2) | nucl_bm[seq[j]]);

            if (j >= K-1) {
                std::vector<std::pair<uint32_t, int>>& kmer_counts = kmer_counts_by_pos[read_start_offsets[i]+j];
                bool found = false;
                for (int i = 0; i < kmer_counts.size(); i++) {
                    if (kmer_counts[i].first == kmer) {
                        kmer_counts[i].second++;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    kmer_counts.push_back({kmer, 1});
                }
            }
        }
    }

    // select pos that maximizes the product of the frequencies of the 1st and 2nd most frequent kmers
    // select most frequent kmer at that pos as mandatory kmer
    uint32_t chosen_kmer = 0, chosen_pos = 0, chosen_freq = 0;
    for (int i = 0; i < kmer_counts_by_pos.size(); i++) {
        // find 1st and 2nd most frequent kmers
        std::sort(kmer_counts_by_pos[i].begin(), kmer_counts_by_pos[i].end(), [](const std::pair<uint32_t, int>& p1, const std::pair<uint32_t, int>& p2) {
            return p1.second > p2.second;
        });

        int kmer1_freq = kmer_counts_by_pos[i].empty() ? 0 : kmer_counts_by_pos[i][0].second;
        int kmer2_freq = kmer_counts_by_pos[i].size() <= 1 ? 0 : kmer_counts_by_pos[i][1].second;
        if (kmer1_freq*kmer2_freq > chosen_freq) {
            chosen_freq = kmer1_freq*kmer2_freq;
            chosen_kmer = kmer_counts_by_pos[i][0].first;
            chosen_pos = i;
        }
    }

    std::vector<int> selected_idxs;
    for (int i = 0; i < seqs.size(); i++) {
        std::string& seq = seqs[i];
        if (seq.length() < K) continue;

        uint32_t kmer = 0;
        int kmer_start = chosen_pos - K + 1 - read_start_offsets[i], kmer_end = chosen_pos - read_start_offsets[i];
        if (kmer_start < 0 || kmer_end >= seq.length()) continue;
        for (int j = kmer_start; j <= kmer_end; j++) {
            kmer = ((kmer << 2) | nucl_bm[seq[j]]);
        }

        if (kmer == chosen_kmer) {
            selected_idxs.push_back(i);
        }
    }
    return selected_idxs;
}

std::vector<bool> find_accepted_reads(std::string& consensus_seq, std::vector<bam_redux_t*>& clipped, std::vector<hts_pos_t>& read_start_offsets,
                         bool left_clipped) {

    std::vector<bool> accepted(clipped.size(), false);
    for (int i = 0; i < clipped.size(); i++) {
        bam_redux_t* r = clipped[i];
        hts_pos_t offset = read_start_offsets[i];
        int mm = 0;
        int mm_clip = 0; // mismatches in the clipped portion only

        // filter reads with too many differences from the consensus_seq
        hts_pos_t clip_start = left_clipped ? 0 : r->seq_len() - r->right_clip_size;
        hts_pos_t clip_end = left_clipped ? r->left_clip_size : r->seq_len();
        for (int i = 0; i < r->seq_len(); i++) {
            if (i+offset >= consensus_seq.length()) {
                std::cerr << "WARNING: consensus_seq out of boundary." << std::endl;
            }
            if (consensus_seq[i + offset] != get_base(r->seq.data(), i)) {
                mm++;
                if (i >= clip_start && i < clip_end) mm_clip++;
            }
        }
        if (mm <= std::ceil(config.max_seq_error * consensus_seq.length()) &&
            mm_clip <= (clip_end-clip_start)/2) {
            /* this is meant to fix a corner case where the clip is very short, and completely different from
            * the consensus_seq. However, the rest of the reads is the same as the consensus_seq, therefore the mismatches
            * accumulated in the clip are not enough to discard the read, despite the read not belonging to the cluster.
            * We are being very permissive (i.e. we allow up to 50% mismatches in the clip) because
            * 1. Illumina is known to accumulate sequencing errors at the 3' tail, which is the clipped portion in about
            * 50% of the cases
            * 2. Especially for short clips, if we used config.max_seq_error as for the rest of the consensus_seq, 2-3 mismatches
            * would already discard them
            */
            // the read should be much better when mapped to the consensus than when mapped to the reference
            int orig_score = compute_read_score(r, 1, -4, -6, -1);
            int new_score = (r->seq_len()-mm)*1 - mm*4;

            if (new_score - orig_score >= config.min_diff_hsr*5) { // each mismatch costs 5 points
                accepted[i] = true;
            }
        }
    }
    return accepted;
}

std::string build_full_consensus_seq(std::vector<bam_redux_t*>& clipped, bool left_clipped, bool use_kmer_selection,
                                     std::vector<bool>& accepted) {
    std::vector<std::string> seqs;
    std::vector<uint8_t*> quals;
    std::vector<hts_pos_t> read_start_offsets = get_read_start_offsets(clipped, left_clipped);
    for (bam_redux_t* r : clipped) {
        std::string seq(r->seq_len(), 'N');
        for (int i = 0; i < r->seq_len(); i++) {
            seq[i] = get_base(r->seq.data(), i);
        }
        seqs.push_back(seq);
        quals.push_back(r->qual.data());
    }

    std::vector<bam_redux_t*> selected_clipped;
    std::vector<int> selected_idxs;
    if (use_kmer_selection) { // let's try partitioning the sequences according to kmer
        selected_idxs = select_reads_by_kmer(seqs, read_start_offsets);
        if (selected_idxs.size() >= 3) {
            std::vector<std::string> selected_seqs;
            std::vector<uint8_t*> selected_quals;
            std::vector<hts_pos_t> selected_read_start_offsets;
            int min_offset = INT32_MAX;
            for (int i : selected_idxs) {
                selected_seqs.push_back(seqs[i]);
                selected_quals.push_back(quals[i]);
                selected_read_start_offsets.push_back(read_start_offsets[i]);
                selected_clipped.push_back(clipped[i]);
                if (min_offset > read_start_offsets[i]) min_offset = read_start_offsets[i];
            }
            seqs = selected_seqs;
            quals = selected_quals;
            read_start_offsets = selected_read_start_offsets;
            for (int i = 0; i < read_start_offsets.size(); i++) {
                read_start_offsets[i] -= min_offset;
            }
        }
    } else {
        selected_clipped = clipped;
        selected_idxs.resize(clipped.size());
        for (int i = 0; i < clipped.size(); i++) selected_idxs[i] = i;
    }

    std::string consensus_seq = build_full_consensus_seq(seqs, quals, read_start_offsets);
    if (consensus_seq == "") return consensus_seq;
    
    std::vector<bool> selected_accepted = find_accepted_reads(consensus_seq, selected_clipped, read_start_offsets, left_clipped);

    accepted = std::vector<bool>(clipped.size(), false);
    for (int i = 0; i < selected_accepted.size(); i++) {
        if (selected_accepted[i]) {
            accepted[selected_idxs[i]] = true;
        }
    }

    return consensus_seq;
}

std::vector<consensus_t*> build_full_consensus(int contig_id, std::vector<bam_redux_t*> clipped, bool left_clipped) {

    std::vector<consensus_t*> consensuses;

    if (clipped.size() <= 2) return consensuses;

    std::sort(clipped.begin(), clipped.end(), [](bam_redux_t* r1, bam_redux_t* r2) {
        return r1->unclipped_start() < r2->unclipped_start();
    });
    if (clipped[0]->unclipped_start() < 0) return consensuses;

    while (clipped.size() >= 3) {
        std::vector<bool> accepted;
        std::string consensus_seq = build_full_consensus_seq(clipped, left_clipped, false, accepted);
        if (consensus_seq == "") return consensuses;

        int accepted_reads_n = std::count(accepted.begin(), accepted.end(), true);
        if (accepted_reads_n < 3) {
            consensus_seq = build_full_consensus_seq(clipped, left_clipped, true, accepted);
        }

        std::vector<bam_redux_t*> accepted_reads, rejected_reads;
        for (size_t i = 0; i < clipped.size(); i++) {
            if (accepted[i]) accepted_reads.push_back(clipped[i]);
            else rejected_reads.push_back(clipped[i]);
        }

        if (accepted_reads.empty()) return consensuses;

        if (accepted_reads.size() >= 3) {
            hts_pos_t start = accepted_reads[0]->unclipped_start(), end = 0;
            hts_pos_t breakpoint = left_clipped ? INT32_MAX : 0; // the current HTS_POS_MAX does not compile on some compilers
            hts_pos_t remap_boundary = left_clipped ? consensus_t::LOWER_BOUNDARY_NON_CALCULATED : consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
            int fwd_clipped = 0, rev_clipped = 0;
            uint8_t max_mapq = 0;
            for (bam_redux_t* r : accepted_reads) {
                breakpoint = left_clipped ? std::min(breakpoint, r->start) : std::max(breakpoint, r->end);
                end = std::max(end, r->unclipped_end());
                if (r->is_rev()) rev_clipped++;
                else fwd_clipped++;

                max_mapq = std::max(max_mapq, r->mapq);

                if (r->is_inter_chr() || r->is_rev() == r->is_mrev()) continue;
                if (left_clipped && r->is_rev() && !r->mate_left_clipped()) {
                    remap_boundary = std::max(remap_boundary, r->mstart);
                } else if (!left_clipped && !r->is_rev() && !r->mate_right_clipped()) {
                    remap_boundary = std::min(remap_boundary, r->start+r->isize);
                }
            }

            // rebuild consensus sequence using only accepted reads
            consensus_seq = build_full_consensus_seq(accepted_reads, left_clipped, false, accepted); 

            // these bps have support from only one or two reads, so they are prone to errors
            hts_pos_t remove_from_start = get_start_offset(accepted_reads[0], accepted_reads[2]);

            sort(accepted_reads.begin(), accepted_reads.end(), [](bam_redux_t* r1, bam_redux_t* r2) {
                return r1->unclipped_end() > r2->unclipped_end();
            });
            hts_pos_t remove_from_end = get_end_offset(accepted_reads[2], accepted_reads[0]);

            int clip_len = left_clipped ? breakpoint-start : end-breakpoint;
            int lowq_clip_portion = left_clipped ? remove_from_start : remove_from_end;

            consensus_t* consensus = new consensus_t(left_clipped, start, breakpoint, end, consensus_seq, 
                fwd_clipped, rev_clipped, clip_len, max_mapq, remap_boundary, lowq_clip_portion);
            consensuses.push_back(consensus);
        }
        clipped.swap(rejected_reads);
    }
    return consensuses;
}

void build_sr_consensuses(int id, int contig_id, std::string contig_name, hts_pos_t contig_len) {

    std::string clip_fname = workspace + "/clipped/" + std::to_string(contig_id) + ".bam";
    if (!file_exists(clip_fname)) return;

    open_samFile_t* bam_file = open_samFile(clip_fname, true);
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::vector<bam_redux_t*> lc_reads, rc_reads;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_left_clipped(read, config.min_clip_len) && !is_right_clipped(read, 1)) {
            bam_redux_t* bam_redux = new bam_redux_t(read);
        	if (bam_redux->unclipped_end() < contig_len) {
				lc_reads.push_back(bam_redux);
        	} else {
        		delete bam_redux;
        	}
        }
        if (is_right_clipped(read, config.min_clip_len) && !is_left_clipped(read, 1)) {
        	bam_redux_t* bam_redux = new bam_redux_t(read);
        	if (bam_redux->unclipped_start() >= 0) {
				rc_reads.push_back(bam_redux);
        	} else {
        		delete bam_redux;
        	}
        }
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);
    close_samFile(bam_file);

    std::vector<consensus_t*> full_consensuses;
    std::vector<bam_redux_t*> curr_candidate_cluster;

    auto lc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->start-r2->start) <= config.max_clipped_pos_dist;};
    if (!lc_reads.empty()) {
        for (bam_redux_t* lc_read : lc_reads) {
            if (!curr_candidate_cluster.empty() &&
            !lc_same_cluster(curr_candidate_cluster[0], lc_read)) { // candidate cluster complete
                std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, curr_candidate_cluster, true);
                curr_candidate_cluster.clear();
                full_consensuses.insert(full_consensuses.end(), consensuses.begin(), consensuses.end());            
            }
            curr_candidate_cluster.push_back(lc_read);
        }
        // process last cluster
        std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, curr_candidate_cluster, true);
        curr_candidate_cluster.clear();
        full_consensuses.insert(full_consensuses.end(), consensuses.begin(), consensuses.end());
    }
    auto rc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->end-r2->end) <= config.max_clipped_pos_dist;};
    sort(rc_reads.begin(), rc_reads.end(), [](bam_redux_t* r1, bam_redux_t* r2) { return r1->end < r2->end; });
    if (!rc_reads.empty()) {
        for (bam_redux_t* rc_read : rc_reads) {
            if (!curr_candidate_cluster.empty() &&
            !rc_same_cluster(curr_candidate_cluster[0], rc_read)) { // candidate cluster complete
                std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, curr_candidate_cluster, false);
                curr_candidate_cluster.clear();
                full_consensuses.insert(full_consensuses.end(), consensuses.begin(), consensuses.end());
            }
            curr_candidate_cluster.push_back(rc_read);
        }
        // process last cluster
        std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, curr_candidate_cluster, false);
        curr_candidate_cluster.clear();
        full_consensuses.insert(full_consensuses.end(), consensuses.begin(), consensuses.end());
    }

    for (bam_redux_t* r : lc_reads) delete r;
    for (bam_redux_t* r : rc_reads) delete r;

    if (full_consensuses.empty()) return;

    std::ofstream clip_fout(workspace + "/sr_consensuses/" + std::to_string(contig_id) + ".txt");
    for (consensus_t* consensus : full_consensuses) {
        clip_fout << consensus->to_string() << std::endl;
        delete consensus;
    }
    clip_fout.close();
}

void build_hsr_consensuses(int id, int contig_id, std::string contig_name, hts_pos_t contig_len) {

    std::string bam_fname = workdir + "/workspace/hsr/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(bam_fname)) return;

    open_samFile_t* bam_file = open_samFile(bam_fname, true);

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::deque<bam1_t*> lc_cluster, rc_cluster;
    std::vector<consensus_t*> rc_consensuses, lc_consensuses;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {

        if (is_mate_unmapped(read)) continue; // TODO: should I keep this?

        std::pair<int, int> left_and_right_diffs = compute_left_and_right_differences(read, true);
        if (left_and_right_diffs.first == left_and_right_diffs.second) {
            left_and_right_diffs = compute_left_and_right_differences(read, false);
        }

        bool lc_clipped = left_and_right_diffs.first > left_and_right_diffs.second;
        std::deque<bam1_t*>& cluster = lc_clipped ? lc_cluster : rc_cluster;

        if (cluster.size() >= 3 && bam_endpos(cluster.front())-read->core.pos < read->core.l_qseq/2) {

            std::vector<bam_redux_t*> cluster_v;
			for (bam1_t* r : cluster) cluster_v.push_back(new bam_redux_t(r));
			std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, cluster_v, lc_clipped);
			for (consensus_t* consensus : consensuses) {
				consensus->is_hsr = true;
				consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
				if (lc_clipped) lc_consensuses.push_back(consensus);
                else rc_consensuses.push_back(consensus);
			}
			for (bam_redux_t* r : cluster_v) delete r;
        }
        while (!cluster.empty() && bam_endpos(cluster.front())-read->core.pos < read->core.l_qseq/2) {
            bam_destroy1(cluster.front());
            cluster.pop_front();
        }
        cluster.push_back(bam_dup1(read));
    }
    if (rc_cluster.size() >= 3) {
        std::vector<bam_redux_t*> cluster_v;
        for (bam1_t* r : rc_cluster) cluster_v.push_back(new bam_redux_t(r));
        std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, cluster_v, false);
        for (auto consensus : consensuses) {
            consensus->is_hsr = true;
            consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
            rc_consensuses.push_back(consensus);
        }
        for (bam_redux_t* r : cluster_v) delete r;
    }
    for (bam1_t* r : rc_cluster) bam_destroy1(r);
    if (lc_cluster.size() >= 3) {
        std::vector<bam_redux_t*> cluster_v;
        for (bam1_t* r : lc_cluster) cluster_v.push_back(new bam_redux_t(r));
        std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, cluster_v, true);
        for (auto consensus : consensuses) {
            consensus->is_hsr = true;
            consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
            lc_consensuses.push_back(consensus);
        }
        for (bam_redux_t* r : cluster_v) delete r;
    }
    for (bam1_t* r : lc_cluster) bam_destroy1(r);
    
    close_samFile(bam_file);
    bam_destroy1(read);
    hts_itr_destroy(iter);

    if (lc_consensuses.empty() && rc_consensuses.empty()) return;

    filter_well_aligned_to_ref(contigs.get_seq(contig_name), contigs.get_len(contig_name), rc_consensuses, config);
    filter_well_aligned_to_ref(contigs.get_seq(contig_name), contigs.get_len(contig_name), lc_consensuses, config);
    select_nonoverlapping_clusters(rc_consensuses);
    select_nonoverlapping_clusters(lc_consensuses);
    enforce_max_ploidy(rc_consensuses, 2);
    enforce_max_ploidy(lc_consensuses, 2);

    std::vector<consensus_t*> all_consensuses;
    all_consensuses.insert(all_consensuses.end(), rc_consensuses.begin(), rc_consensuses.end());
    all_consensuses.insert(all_consensuses.end(), lc_consensuses.begin(), lc_consensuses.end());

    std::ofstream clip_fout(workspace + "/hsr_consensuses/" + std::to_string(contig_id) + ".txt");
    for (consensus_t* consensus : all_consensuses) {
        clip_fout << consensus->to_string() << std::endl;
        delete consensus;
    }
    clip_fout.close();
}

int main(int argc, char* argv[]) {
    workdir = argv[1];
    workspace = workdir + "/workspace";

    std::string reference_fname = argv[2];

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");

    contigs.read_fasta_into_map(reference_fname);

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        hts_pos_t contig_len = contigs.get_len(contig_name);
        std::future<void> future = thread_pool.push(build_sr_consensuses, contig_id, contig_name, contig_len);
        futures.push_back(std::move(future));
        future = thread_pool.push(build_hsr_consensuses, contig_id, contig_name, contig_len);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();
}
