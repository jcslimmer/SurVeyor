#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <algorithm>
#include <queue>
#include <unistd.h>
#include <climits>
#include <random>
#include <cassert>

#include <htslib/sam.h>
#include <htslib/kseq.h>

#include "dc_remapper.h"
#include "../src/sam_utils.h"
#include "../src/clustering_utils.h"
#include "utils.h"
#include "vcf_utils.h"
#include "assemble.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "../libs/IntervalTree.h"


config_t config;
stats_t stats;
std::string workdir;
std::mutex mtx;

contig_map_t contig_map;
chr_seqs_map_t contigs;

std::ofstream assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer, assembly_failed_bad_anchors_writer;
std::ofstream assembly_failed_lt50bp, assembly_failed_no_seq, assembly_succeeded;

std::vector<inss_insertion_t*> assembled_insertions, trans_insertions;

const int SMALL_SAMPLE_SIZE = 15;
const int CLUSTER_CANDIDATES = 3;
const double BASE_ACCEPTANCE_THRESHOLD = 0.95;
const int TOO_MANY_READS = 1000;


struct region_score_t {
    int total_score = 0;
    int l_clip_score = 0, r_clip_score = 0;
    int remap_start = 0, remap_end = 0;
};
bool operator > (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score > s2.total_score;
}
bool operator >= (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score >= s2.total_score;
}
bool operator < (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score < s2.total_score;
}

struct region_t {
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    region_score_t score;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}
};
bool operator > (const region_t& r1, const region_t& r2) {
	if (r1.score.total_score != r2.score.total_score) {
		return r1.score.total_score > r2.score.total_score;
	}
	return r1.end-r1.start < r2.end-r2.start;
}

struct cc_v_distance_t {
    reads_cluster_t* c1, * c2;
    int distance;

    cc_v_distance_t(reads_cluster_t* c1, reads_cluster_t* c2, int distance) : c1(c1), c2(c2), distance(distance) {}
};
bool operator < (const cc_v_distance_t& ccd1, const cc_v_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance < ccd2.distance;
}

struct remap_info_t {
    int start, end;
    int score;
    std::string cigar;
    bool accepted, left_clipped, right_clipped;

    remap_info_t() : start(0), end(0), score(0), cigar(""), accepted(false), left_clipped(false), right_clipped(false) {}

    remap_info_t(StripedSmithWaterman::Alignment& aln, bool accepted) : start(aln.ref_begin), end(aln.ref_end),
    score(aln.sw_score), cigar(aln.cigar_string), accepted(accepted) {
        uint32_t f = aln.cigar[0], l = aln.cigar[aln.cigar.size()-1];
        left_clipped = cigar_int_to_op(f) == 'S' && cigar_int_to_len(f) >= config.min_clip_len;
        right_clipped = cigar_int_to_op(l) == 'S' && cigar_int_to_len(l) >= config.min_clip_len;
    }
};

region_t get_region(std::vector<bam1_t*> subcluster, std::string& m_contig_name) {
    bam1_t* leftmost_reverse_mate = NULL, * rightmost_forward_mate = NULL;
    /* gets two regions
     * 1. from left-most reverse read, max_is to the left and max_insertion_size to the right
     * 2. from right-most forward read, max_is to the right and max_insertion_size to the left
     */
    for (bam1_t* r : subcluster) { // get leftmost reverse read
        if (bam_is_mrev(r) && (leftmost_reverse_mate == NULL || leftmost_reverse_mate->core.mpos > r->core.mpos)) {
            leftmost_reverse_mate = r;
        }
    }
    for (bam1_t* r : subcluster) { // get rightmost forward read
        if (!bam_is_mrev(r) && (rightmost_forward_mate == NULL || rightmost_forward_mate->core.mpos < r->core.mpos)) {
            rightmost_forward_mate = r;
        }
    }

    hts_pos_t start = INT_MAX;
    hts_pos_t end = 0;
    if (leftmost_reverse_mate != NULL) {
        start = std::min(start, leftmost_reverse_mate->core.mpos-stats.max_is);
        end = std::max(end, leftmost_reverse_mate->core.mpos+config.max_trans_size);
    }
    if (rightmost_forward_mate != NULL) {
        start = std::min(start, rightmost_forward_mate->core.mpos-config.max_trans_size);
        end = std::max(end, rightmost_forward_mate->core.mpos+stats.max_is);
    }

    hts_pos_t contig_len = contigs.get_len(m_contig_name);
    return region_t(contig_map.get_id(m_contig_name), subcluster[0]->core.mtid, std::max(hts_pos_t(0),start), std::min(end,contig_len));
}

remap_info_t remap_mate(char* region, int region_len, std::string& mate_seq, StripedSmithWaterman::Aligner& aligner,
                StripedSmithWaterman::Filter& filter) {
    if (region_len <= 0) {
        return remap_info_t();
    }

    StripedSmithWaterman::Alignment alignment;
    aligner.Align(mate_seq.c_str(), region, region_len, filter, &alignment, 0);
    bool accepted = alignment.sw_score >= 30; //s.length(); TODO

    return remap_info_t(alignment, accepted);
}

int get_strict_region_start(std::vector<remap_info_t>& rc_remap_infos, std::vector<remap_info_t>& lc_remap_infos) {
    std::vector<remap_info_t> remap_infos;
    for (remap_info_t& ri : rc_remap_infos) {
        if (ri.left_clipped) remap_infos.push_back(ri);
    }
    for (remap_info_t& ri : lc_remap_infos) {
        if (ri.left_clipped) remap_infos.push_back(ri);
    }
    if (remap_infos.empty()) return 0;

    remap_info_t guard; guard.start = 1000000;
    remap_infos.push_back(guard);

    std::sort(remap_infos.begin(), remap_infos.end(),
              [](remap_info_t& ri1, remap_info_t& ri2) { return ri1.start < ri2.start; });

    int strict_region_start = remap_infos[0].start, max_freq = 1, curr_freq = 1;
    for (int i = 1; i < remap_infos.size(); i++) {
        if (remap_infos[i-1].start != remap_infos[i].start) {
            if (curr_freq > max_freq) {
                strict_region_start = remap_infos[i-1].start;
                max_freq = curr_freq;
            }
            curr_freq = 0;
        }
        curr_freq++;
    }

    return strict_region_start;
}
int get_strict_region_end(std::vector<remap_info_t>& rc_remap_infos, std::vector<remap_info_t>& lc_remap_infos) {
    std::vector<remap_info_t> remap_infos;
    for (remap_info_t& ri : rc_remap_infos) {
        if (ri.right_clipped) remap_infos.push_back(ri);
    }
    for (remap_info_t& ri : lc_remap_infos) {
        if (ri.right_clipped) remap_infos.push_back(ri);
    }
    if (remap_infos.empty()) return 0;

    remap_info_t guard; guard.end = 0;
    remap_infos.push_back(guard);

    std::sort(remap_infos.begin(), remap_infos.end(),
              [](remap_info_t& ri1, remap_info_t& ri2) { return ri1.end > ri2.end; });

    int strict_region_end = remap_infos[0].end, max_freq = 1, curr_freq = 1;
    for (int i = 1; i < remap_infos.size(); i++) {
        if (remap_infos[i-1].end != remap_infos[i].end) {
            if (curr_freq > max_freq) {
                strict_region_end = remap_infos[i-1].end;
                max_freq = curr_freq;
            }
            curr_freq = 0;
        }
        curr_freq++;
    }

    return strict_region_end;
}

// if remapped reads are too far away, find the max_is bp window where the sum of read scores is maximised, and remap all the reads there
void remap_to_maxis_window(std::vector<remap_info_t>& remap_infos, remap_info_t& clip_remap_info, char* region_ptr,
		int remapped_region_start, int remapped_region_end, std::vector<std::string>& read_seqs, std::string& clip_seq, char clip_dir,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
	std::sort(remap_infos.begin(), remap_infos.end(),
			[](const remap_info_t& ri1, const remap_info_t& ri2) { return ri1.start < ri2.start; });
	int s = 0, score = 0;
	int best_start = 0, best_end = 0, best_score = 0;
	for (int e = 0; e < remap_infos.size(); e++) {
		while (remap_infos[e].start-remap_infos[s].start > stats.max_is) {
			score -= remap_infos[s].score;
			s++;
		}
		score += remap_infos[e].score;
		if (score > best_score) {
			best_score = score;
			best_start = s, best_end = e;
		}
	}

	int strictest_region_start = remap_infos[best_start].start, strictest_region_end = remap_infos[best_end].end;
	int strictest_region_len = strictest_region_end - strictest_region_start;
	remap_infos.clear();


	// add some additional space to remap the clip
	int padding = stats.max_is - strictest_region_len;
	if (!clip_seq.empty() && padding > 0) {
		if (clip_dir == 'L') {
			strictest_region_start -= padding;
			// we do not allow the new strictest region to cross the original remapped region - this is the simplest way to make sure
			// it does not go out of the contig boundary and crash
			strictest_region_start = std::max(strictest_region_start, remapped_region_start);
		}
		else if (clip_dir == 'R') {
			strictest_region_end += padding;
			strictest_region_end = std::min(strictest_region_end, remapped_region_end);
		}
		strictest_region_len = strictest_region_end - strictest_region_start;
	}

	for (std::string& s : read_seqs) {
		remap_info_t remap_info = remap_mate(region_ptr+strictest_region_start, strictest_region_len, s, aligner, filter);
		remap_info.start += strictest_region_start;
		remap_info.end += strictest_region_start;
		remap_infos.push_back(remap_info);
	}
	if (!clip_seq.empty()) {
		clip_remap_info = remap_mate(region_ptr+strictest_region_start, strictest_region_len, clip_seq, aligner, filter);
		clip_remap_info.start += strictest_region_start;
		clip_remap_info.end += strictest_region_start;
	}
}

region_score_t compute_score_supp(region_t& region, reads_cluster_t* r_cluster, reads_cluster_t* l_cluster,
                       std::unordered_map<std::string, std::string>& mateseqs,
                       std::vector<remap_info_t>* rc_remap_infos_, std::vector<remap_info_t>* lc_remap_infos_,
                       StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                       StripedSmithWaterman::Filter& filter, bool do_rc) {

    // mateseqs contains seqs are stored as in fasta/q file
    std::vector<std::string> r_mates, l_mates;
    for (bam1_t* read : r_cluster->reads) {
        std::string s = get_mate_seq(read, mateseqs);
        if (!do_rc) rc(s);   // by default, we map mates of fwd reads on the rev strand
        r_mates.push_back(s);
    }
    for (bam1_t* read : l_cluster->reads) {
        std::string s = get_mate_seq(read, mateseqs);
        if (do_rc) rc(s);   // by default, we map mates of rev reads on the fwd strand
        l_mates.push_back(s);
    }

    char* region_ptr = contigs.get_seq(contig_map.get_name(region.contig_id)) + region.start;
    int strict_region_start = 0, strict_region_end = region.end - region.start;

    region_score_t score;
    remap_info_t rc_remap_info, lc_remap_info;

    // map the clips (if available)
    if (r_cluster->clip_consensus) {
        std::string s = r_cluster->clip_consensus->clip_sequence();
        if (do_rc) rc(s);
        remap_info_t remap_info = remap_mate(region_ptr, region.end-region.start, s, permissive_aligner, filter);
        if (remap_info.accepted) {
            score.total_score += remap_info.score;
            score.r_clip_score = remap_info.score;
            if (!do_rc && remap_info.accepted) {
                strict_region_start = remap_info.start;
            } else if (do_rc && remap_info.accepted) {
                strict_region_end = remap_info.end;
            }
        }
        rc_remap_info = remap_info;
    }
    if (l_cluster->clip_consensus) {
        std::string s = l_cluster->clip_consensus->clip_sequence();
        if (do_rc) rc(s);
        remap_info_t remap_info = remap_mate(region_ptr, region.end-region.start, s, permissive_aligner, filter);
        if (remap_info.accepted) {
            score.total_score += remap_info.score;
            score.r_clip_score = remap_info.score;
            if (do_rc && remap_info.accepted) {
                strict_region_start = remap_info.start;
            } else if (!do_rc && remap_info.accepted) {
                strict_region_end = remap_info.end;
            }
        }
        lc_remap_info = remap_info;
    }

    // map the reads
    std::vector<remap_info_t> rc_remap_infos, lc_remap_infos;
    int strict_region_len = strict_region_end - strict_region_start;

    int min_start = INT32_MAX, max_end = 0;
    for (std::string& s : r_mates) {
        remap_info_t remap_info = remap_mate(region_ptr+strict_region_start, strict_region_len, s, aligner, filter);
        remap_info.start += strict_region_start;
        min_start = std::min(min_start, remap_info.start);
        remap_info.end += strict_region_start;
        max_end = std::max(max_end, remap_info.end);
        rc_remap_infos.push_back(remap_info);
    }
    if (max_end-min_start > stats.max_is) { // all the mappings need to be within a config.max_is bp window
    	std::string clip_seq = r_cluster->clip_consensus ? r_cluster->clip_consensus->clip_sequence() : "";
    	char clip_dir = do_rc ? 'R' : 'L'; // R means that the clip determines the right-most boundary of the inserted region/sequence, L the left-most
    	remap_to_maxis_window(rc_remap_infos, rc_remap_info, region_ptr, strict_region_start, strict_region_end,
    			r_mates, clip_seq, clip_dir, aligner, filter);
    }

    min_start = INT32_MAX, max_end = 0;
    for (std::string& s : l_mates) {
        remap_info_t remap_info = remap_mate(region_ptr+strict_region_start, strict_region_len, s, aligner, filter);
        remap_info.start += strict_region_start;
        min_start = std::min(min_start, remap_info.start);
        remap_info.end += strict_region_start;
        max_end = std::max(max_end, remap_info.end);
        lc_remap_infos.push_back(remap_info);
    }
    if (max_end-min_start > stats.max_is) { // all the mappings need to be within a config.max_is bp window
    	std::string clip_seq = l_cluster->clip_consensus ? l_cluster->clip_consensus->clip_sequence() : "";
    	char clip_dir = do_rc ? 'L' : 'R';
    	remap_to_maxis_window(lc_remap_infos, lc_remap_info, region_ptr, strict_region_start, strict_region_end,
    			l_mates, clip_seq, clip_dir, aligner, filter);
    }

    // if strict limits not identified through remapping of the clips, choose most common start and most common ends
    if (strict_region_start == 0) {
        strict_region_start = get_strict_region_start(rc_remap_infos, lc_remap_infos);
    }
    if (strict_region_end == region.end-region.start) {
        int temp = get_strict_region_end(rc_remap_infos, lc_remap_infos);
        if (temp > 0) strict_region_end = temp;
    }

    // un-accept all remappings that fall outside the strict limits or they are not clipped correctly
    for (remap_info_t& remap_info : rc_remap_infos) {
        if (remap_info.left_clipped && strict_region_start < remap_info.start-5) remap_info.accepted = false;
        if (remap_info.end < strict_region_start) remap_info.accepted = false;
        if (remap_info.right_clipped && remap_info.end+5 < strict_region_end) remap_info.accepted = false;
        if (strict_region_end < remap_info.start) remap_info.accepted = false;
    }
    for (remap_info_t& remap_info : lc_remap_infos) {
        if (remap_info.left_clipped && strict_region_start < remap_info.start-5) remap_info.accepted = false;
        if (remap_info.end < strict_region_start) remap_info.accepted = false;
        if (remap_info.right_clipped && remap_info.end+5 < strict_region_end) remap_info.accepted = false;
        if (strict_region_end < remap_info.start) remap_info.accepted = false;
    }

    for (remap_info_t& remap_info : rc_remap_infos) {
        if (remap_info.accepted) score.total_score += remap_info.score;
    }
    for (remap_info_t& remap_info : lc_remap_infos) {
        if (remap_info.accepted) score.total_score += remap_info.score;
    }

    // find inserted sequence coordinates
    int rc_start = INT32_MAX, rc_end = 0;
	if (rc_remap_info.accepted) {
		rc_start = rc_remap_info.start;
		rc_end = rc_remap_info.end;
	} else {
		for (remap_info_t& ri : rc_remap_infos) {
			if (ri.accepted) {
				rc_start = std::min(rc_start, ri.start);
				rc_end = std::max(rc_end, ri.end);
			}
		}
	}
	int lc_start = INT32_MAX, lc_end = 0;
	if (lc_remap_info.accepted) {
		lc_start = lc_remap_info.start;
		lc_end = lc_remap_info.end;
	} else {
		for (remap_info_t& ri : lc_remap_infos) {
			if (ri.accepted) {
				lc_start = std::min(lc_start, ri.start);
				lc_end = std::max(lc_end, ri.end);
			}
		}
	}

	if (!do_rc) {
		score.remap_start = rc_start;
		score.remap_end = lc_end;
	} else {
		score.remap_start = lc_start;
		score.remap_end = rc_end;
	}

    if (rc_remap_infos_ && r_cluster->clip_consensus) {
        rc_remap_infos.push_back(rc_remap_info);
    }
    if (lc_remap_infos_ && l_cluster->clip_consensus) {
        lc_remap_infos.push_back(lc_remap_info);
    }

    if (rc_remap_infos_) rc_remap_infos_->swap(rc_remap_infos);
    if (lc_remap_infos_) lc_remap_infos_->swap(lc_remap_infos);

    return score;
}

void compute_score(region_t& region, reads_cluster_t* r_cluster, reads_cluster_t* l_cluster,
                   std::unordered_map<std::string, std::string>& mateseqs,
                   std::vector<remap_info_t>* rc_remap_infos, std::vector<remap_info_t>* lc_remap_infos,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                   StripedSmithWaterman::Filter& filter, bool& is_rc) {

    region_score_t score = compute_score_supp(region, r_cluster, l_cluster, mateseqs, NULL, NULL,
                                              aligner,permissive_aligner, filter, false);
    region_score_t rc_score = compute_score_supp(region, r_cluster, l_cluster, mateseqs, NULL, NULL,
                                                 aligner, permissive_aligner, filter, true);
    region.score = std::max(score, rc_score);
    if (score >= rc_score) {
        is_rc = false;
        if (rc_remap_infos != NULL) {
            compute_score_supp(region, r_cluster, l_cluster, mateseqs, rc_remap_infos, lc_remap_infos, aligner,
                               permissive_aligner, filter, false);
        }
    } else {
        is_rc = true;
        if (rc_remap_infos != NULL) {
            compute_score_supp(region, r_cluster, l_cluster, mateseqs, rc_remap_infos, lc_remap_infos, aligner,
                               permissive_aligner, filter, true);
        }
    }
}

reads_cluster_t* subsample_cluster(reads_cluster_t* reads_cluster, int size, std::default_random_engine& rng) {
    std::vector<bam1_t*> subset(reads_cluster->reads);
    if (subset.size() > size) {
        std::shuffle(subset.begin(), subset.end(), rng);
        subset.erase(subset.begin() + size, subset.end());
    }
    reads_cluster_t* subsampled_cluster = new reads_cluster_t();
    for (bam1_t* read : subset) subsampled_cluster->add_read(read);
    subsampled_cluster->add_clip_cluster(reads_cluster->clip_consensus);
    return subsampled_cluster;
}

void update_read(bam1_t* read, region_t& chosen_region, remap_info_t& remap_info, bool region_is_rc) {
    read->core.mtid = chosen_region.original_bam_id;
    read->core.mpos = chosen_region.start + remap_info.start;
    if (region_is_rc == bam_is_rev(read)) {
        read->core.flag |= BAM_FMREVERSE; //sets flag to true
    } else {
        read->core.flag &= ~BAM_FMREVERSE; //sets flag to false
    }
    bam_aux_update_str(read, "MC", remap_info.cigar.length() + 1, remap_info.cigar.c_str());
    char ok = remap_info.accepted ? 'T' : 'F';
    bam_aux_append(read, "OK", 'A', 1, (const uint8_t*) &ok);
    bam_aux_append(read, "MS", 'i', 4, (const uint8_t*) &remap_info.score);
}

std::string generate_consensus_sequences(std::string contig_name, reads_cluster_t* r_cluster, reads_cluster_t* l_cluster, region_t& best_region,
		bool is_rc, bool& left_bp_precise, bool& right_bp_precise, std::unordered_map<std::string, std::string>& mateseqs,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner) {

	left_bp_precise = false, right_bp_precise = false;

	if (best_region.score.remap_start >= best_region.score.remap_end) return "";

	if (r_cluster->end()-r_cluster->start() <= 0 || l_cluster->end()-l_cluster->start() <= 0) return "";

	char* left_flanking = new char[r_cluster->end()-r_cluster->start()+2];
	strncpy(left_flanking, contigs.get_seq(contig_name)+r_cluster->start(), r_cluster->end()-r_cluster->start()+1);
	left_flanking[r_cluster->end()-r_cluster->start()+1] = '\0';
	to_uppercase(left_flanking);

	int ins_seq_start = best_region.start + best_region.score.remap_start,
		ins_seq_len = best_region.score.remap_end - best_region.score.remap_start + 1;
	char* pred_ins_seq = new char[ins_seq_len+1];
	strncpy(pred_ins_seq, contigs.get_seq(contig_map.get_name(best_region.contig_id))+ins_seq_start, ins_seq_len);
	pred_ins_seq[ins_seq_len] = '\0';
	if (is_rc) rc(pred_ins_seq);
	to_uppercase(pred_ins_seq);

	char* right_flanking = new char[l_cluster->end()-l_cluster->start()+2];
	strncpy(right_flanking, contigs.get_seq(contig_name)+l_cluster->start(), l_cluster->end()-l_cluster->start()+1);
	right_flanking[l_cluster->end()-l_cluster->start()+1] = '\0';
	to_uppercase(right_flanking);

	// NOTE: prefixes and suffixes may be duplicated - i.e., suffix of left flanking also appearing as prefix of predicted ins seq
	// This is a problem because it predicts an insertion where there is none. One such case is fragment LINE insertions.
	// There are small regions in LINE/L1P* that accumulate a lot of mutations compared to reference, and it may be predicted as a block-swap
	// (which is not an insertion as the length of alt allele is not changed)
	// However, because of the duplicated suffix/prefix, a small insertion is predicted
	std::string left_flanking_str = left_flanking, pred_ins_seq_str = pred_ins_seq, right_flanking_str = right_flanking;

	suffix_prefix_aln_t spa12 = aln_suffix_prefix(left_flanking_str, pred_ins_seq_str, 1, -4, 1.0, config.min_clip_len);
	int overlap12 = spa12.overlap;
	if (!overlap12) {
		StripedSmithWaterman::Alignment aln;
		StripedSmithWaterman::Filter filter;
		aligner.Align(pred_ins_seq_str.c_str(), left_flanking_str.c_str(), left_flanking_str.length(), filter, &aln, 0);
		if (aln.query_begin == 0 && aln.ref_end == left_flanking_str.length()-1) overlap12 = aln.query_end;
	}

	suffix_prefix_aln_t spa23 = aln_suffix_prefix(pred_ins_seq_str, right_flanking_str, 1, -4, 1.0, config.min_clip_len);
	int overlap23 = spa23.overlap;
	if (!overlap23) {
		StripedSmithWaterman::Alignment aln;
		StripedSmithWaterman::Filter filter;
		aligner.Align(right_flanking_str.c_str(), pred_ins_seq_str.c_str(), pred_ins_seq_str.length(), filter, &aln, 0);
		if (aln.query_begin == 0 && aln.ref_end == pred_ins_seq_str.length()-1) overlap23 = aln.query_end;
	}

	std::string full_junction_sequence = left_flanking_str;
	full_junction_sequence += pred_ins_seq_str.substr(overlap12);
	full_junction_sequence += right_flanking_str.substr(overlap23);

	delete[] left_flanking;
	delete[] pred_ins_seq;
	delete[] right_flanking;


	// assembled contigs
	std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
	std::vector<std::string> consensus_contigs = generate_reference_guided_consensus(full_junction_sequence, r_cluster, l_cluster, mateseqs,
			aligner, harsh_aligner, consensus_contigs_alns, config, stats);

	if (consensus_contigs.empty()) return full_junction_sequence;

	int l_bp_in_seq = r_cluster->end()-r_cluster->start(), r_bp_in_seq = full_junction_sequence.length()-(l_cluster->end()-l_cluster->start());

	const int OFFSET = 50;
	for (StripedSmithWaterman::Alignment& aln : consensus_contigs_alns) {
		left_bp_precise |= (aln.ref_begin < l_bp_in_seq-OFFSET && aln.ref_end >= l_bp_in_seq+OFFSET);
		right_bp_precise |= (aln.ref_begin < r_bp_in_seq-OFFSET && aln.ref_end >= r_bp_in_seq+OFFSET);
	}

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
	for (int i = 0; i < consensus_contigs.size(); i++) {
		contigs_sorted_by_pos.push_back({consensus_contigs[i], consensus_contigs_alns[i]});
	}
	std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	// complement holes in contigs with reference
	std::string corrected_junction_seq = full_junction_sequence.substr(0, contigs_sorted_by_pos[0].second.ref_begin);
	corrected_junction_seq += contigs_sorted_by_pos[0].first;
	for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
		corrected_junction_seq += full_junction_sequence.substr(contigs_sorted_by_pos[i-1].second.ref_end+1, contigs_sorted_by_pos[i].second.ref_begin-contigs_sorted_by_pos[i-1].second.ref_end);
		corrected_junction_seq += contigs_sorted_by_pos[i].first;
	}
	corrected_junction_seq += full_junction_sequence.substr(contigs_sorted_by_pos[contigs_sorted_by_pos.size()-1].second.ref_end+1);

	return corrected_junction_seq;
}

void remap_cluster(reads_cluster_t* r_cluster, reads_cluster_t* l_cluster, std::vector<bam1_t*>& kept,
                   int contig_id, bam_hdr_t* header, std::unordered_map<std::string, std::string>& mateseqs,
				   std::unordered_map<std::string, std::string>& matequals,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                   StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
                   stats_t& stats) {

	if (l_cluster->reads.size() + r_cluster->reads.size() > TOO_MANY_READS) return;

    std::vector<region_t> regions;

    /* == Find candidate regions for remapping == */
    std::vector<bam1_t*> full_cluster;
    full_cluster.insert(full_cluster.end(), l_cluster->reads.begin(), l_cluster->reads.end());
    full_cluster.insert(full_cluster.end(), r_cluster->reads.begin(), r_cluster->reads.end());
    std::sort(full_cluster.begin(), full_cluster.end(), [] (bam1_t* r1, bam1_t* r2) {
        if (r1->core.mtid != r2->core.mtid) return r1->core.mtid < r2->core.mtid;
        else return r1->core.mpos < r2->core.mpos;
    });

    std::vector<bam1_t*> subcluster;
    for (bam1_t* r : full_cluster) {
        if (!subcluster.empty() && (subcluster[0]->core.mtid != r->core.mtid ||
                r->core.mpos-subcluster[0]->core.mpos > stats.max_is)) {
            std::string m_contig_name = header->target_name[subcluster[0]->core.mtid];
            regions.push_back(get_region(subcluster, m_contig_name));
            subcluster.clear();
        }
        if (!is_mate_unmapped(r)) subcluster.push_back(r);
    }
    if (!subcluster.empty()) {
        std::string m_contig_name = header->target_name[subcluster[0]->core.mtid];
        regions.push_back(get_region(subcluster, m_contig_name));
    }

    std::string contig_name = contig_map.get_name(contig_id);
    if (regions.empty()) {
		inss_insertion_t* ins = assemble_insertion(contig_name, contigs, r_cluster, l_cluster, mateseqs, matequals,
				aligner_to_base, harsh_aligner, kept, config, stats);

		if (ins != NULL) {
			mtx.lock();
			assembled_insertions.push_back(ins);
			mtx.unlock();
		}
		return;
    }

    /* == Find best candidate region == */
    StripedSmithWaterman::Filter filter, filter_w_cigar;
    filter_w_cigar.report_cigar = true;
    bool is_rc;

    // if too many regions and too many reads, subsample
    if (full_cluster.size() > 2*SMALL_SAMPLE_SIZE && regions.size() > CLUSTER_CANDIDATES) {
    	std::default_random_engine rng(config.seed);
        reads_cluster_t* subsampled_lc = subsample_cluster(l_cluster, SMALL_SAMPLE_SIZE, rng);
        reads_cluster_t* subsampled_rc = subsample_cluster(r_cluster, SMALL_SAMPLE_SIZE, rng);

        // compute best score
        for (int i = 0; i < regions.size(); i++) {
            compute_score(regions[i], subsampled_rc, subsampled_lc, mateseqs, NULL, NULL,
                          aligner, permissive_aligner, filter, is_rc);
        }
        std::sort(regions.begin(), regions.end(), std::greater<region_t>());

        regions.erase(regions.begin()+CLUSTER_CANDIDATES, regions.end());
        delete subsampled_lc;
        delete subsampled_rc;
    }

    // compute best score
    for (int i = 0; i < regions.size(); i++) {
        compute_score(regions[i], r_cluster, l_cluster, mateseqs, NULL, NULL,
                      aligner, permissive_aligner, filter, is_rc);
    }
    std::sort(regions.begin(), regions.end(), std::greater<region_t>());

    region_t best_region = regions[0];

    // get base region
    int start = r_cluster->start() - stats.max_is;
    int end = l_cluster->end() + stats.max_is;
    int contig_len = contigs.get_len(contig_name);
    region_t base_region(contig_id, full_cluster[0]->core.tid, std::max(0, start), std::min(end, contig_len));

    compute_score(base_region, r_cluster, l_cluster, mateseqs, NULL, NULL,
                  aligner_to_base, permissive_aligner, filter, is_rc);

    if (base_region.score.total_score >= best_region.score.total_score*BASE_ACCEPTANCE_THRESHOLD) {
    	std::vector<remap_info_t> rc_remap_infos, lc_remap_infos;
    	compute_score(base_region, r_cluster, l_cluster, mateseqs, &rc_remap_infos, &lc_remap_infos,
					  aligner_to_base, permissive_aligner, filter, is_rc);
    	int accepted_reads = 0;
		for (remap_info_t& remap_info : rc_remap_infos) accepted_reads += remap_info.accepted;
		for (remap_info_t& remap_info : lc_remap_infos) accepted_reads += remap_info.accepted;
		int tot_reads = r_cluster->reads.size() + l_cluster->reads.size();
		if (double(accepted_reads)/tot_reads >= 0.5) {
			return;
		}
    }

    std::vector<remap_info_t> rc_remap_infos, lc_remap_infos;
    compute_score(best_region, r_cluster, l_cluster, mateseqs, &rc_remap_infos, &lc_remap_infos,
                  aligner, permissive_aligner, filter_w_cigar, is_rc);

    // in order to construct the corrected consensus, we rely on read clusters
    // we refine the reads by eliminating reads eexcluded by the remapping, in order to get the most precise possible locations
    reads_cluster_t* refined_r_cluster = new reads_cluster_t();
	reads_cluster_t* refined_l_cluster = new reads_cluster_t();
	for (int i = 0; i < r_cluster->reads.size(); i++) {
		if (rc_remap_infos[i].accepted) {
			refined_r_cluster->add_read(r_cluster->reads[i]);
		}
	}
	for (int i = 0; i < l_cluster->reads.size(); i++) {
		if (lc_remap_infos[i].accepted) {
			refined_l_cluster->add_read(l_cluster->reads[i]);
		}
	}
	if (r_cluster->clip_consensus && rc_remap_infos.rbegin()->accepted) {
		refined_r_cluster->add_clip_cluster(r_cluster->clip_consensus);
	}
	if (l_cluster->clip_consensus && lc_remap_infos.rbegin()->accepted) {
		refined_l_cluster->add_clip_cluster(l_cluster->clip_consensus);
	}

    // build corrected consensus sequence and realign reads to it
    std::string corrected_consensus_sequence;
	std::vector<std::pair<int,int>> covered_segments;
    bool left_bp_precise = false, right_bp_precise = false;
    if (!refined_r_cluster->empty() && !refined_l_cluster->empty()) {
		corrected_consensus_sequence = generate_consensus_sequences(contig_name, refined_r_cluster, refined_l_cluster,
				best_region, is_rc, left_bp_precise, right_bp_precise, mateseqs, aligner_to_base, harsh_aligner);
    }

    if (!corrected_consensus_sequence.empty()) {
		StripedSmithWaterman::Alignment aln;
		for (int i = 0; i < r_cluster->reads.size(); i++) {
			std::string mate_seq = get_mate_seq(r_cluster->reads[i], mateseqs);
			std::string mate_qual = get_mate_qual(r_cluster->reads[i], matequals);
			rc(mate_seq);
			mate_qual = std::string(mate_qual.rbegin(), mate_qual.rend());
			harsh_aligner.Align(mate_seq.c_str(), corrected_consensus_sequence.c_str(), corrected_consensus_sequence.length(), filter, &aln, 0);
			rc_remap_infos[i].accepted = accept(aln, config.min_clip_len, config.max_seq_error, mate_qual, stats.min_avg_base_qual);
			if (rc_remap_infos[i].accepted) {
				covered_segments.push_back({aln.ref_begin, aln.ref_end});
			}
		}
		for (int i = 0; i < l_cluster->reads.size(); i++) {
			std::string mate_seq = get_mate_seq(l_cluster->reads[i], mateseqs);
			std::string mate_qual = get_mate_qual(l_cluster->reads[i], matequals);
			harsh_aligner.Align(mate_seq.c_str(), corrected_consensus_sequence.c_str(), corrected_consensus_sequence.length(), filter, &aln, 0);
			lc_remap_infos[i].accepted = accept(aln, config.min_clip_len, config.max_seq_error, mate_qual, stats.min_avg_base_qual);
			if (lc_remap_infos[i].accepted) {
				covered_segments.push_back({aln.ref_begin, aln.ref_end});
			}
		}
		if (r_cluster->clip_consensus) {
			harsh_aligner.Align(r_cluster->clip_consensus->sequence.c_str(), corrected_consensus_sequence.c_str(), corrected_consensus_sequence.length(), filter, &aln, 0);
			rc_remap_infos.rbegin()->accepted = accept(aln, config.min_clip_len, config.max_seq_error);
			if (rc_remap_infos.rbegin()->accepted) {
				covered_segments.push_back({aln.ref_begin, aln.ref_end});
			}
		}
		if (l_cluster->clip_consensus) {
			harsh_aligner.Align(l_cluster->clip_consensus->sequence.c_str(), corrected_consensus_sequence.c_str(), corrected_consensus_sequence.length(), filter, &aln, 0);
			lc_remap_infos.rbegin()->accepted = accept(aln, config.min_clip_len, config.max_seq_error);
			if (lc_remap_infos.rbegin()->accepted) {
				covered_segments.push_back({aln.ref_begin, aln.ref_end});
			}
		}
    }


    /* == Try assembly == */
	int rc_accepted_reads = 0, lc_accepted_reads = 0;
	for (remap_info_t& remap_info : rc_remap_infos) rc_accepted_reads += remap_info.accepted;
	for (remap_info_t& remap_info : lc_remap_infos) lc_accepted_reads += remap_info.accepted;
	int tot_reads = r_cluster->reads.size() + l_cluster->reads.size();
	if (rc_accepted_reads == 0 || lc_accepted_reads == 0 || double(rc_accepted_reads+lc_accepted_reads)/tot_reads < 0.5) {
		inss_insertion_t* ins = assemble_insertion(contig_name, contigs, r_cluster, l_cluster, mateseqs, matequals,
				aligner_to_base, harsh_aligner, kept, config, stats);

		if (ins != NULL) {
			mtx.lock();
			assembled_insertions.push_back(ins);
			mtx.unlock();
			return;
		}
	}

	if (corrected_consensus_sequence.empty()) return;

	reads_cluster_t* pos_cluster = new reads_cluster_t();
    for (int i = 0; i < r_cluster->reads.size(); i++) {
        bam1_t* r = r_cluster->reads[i];
        update_read(r, best_region, rc_remap_infos[i], is_rc);
        kept.push_back(r);

        if (rc_remap_infos[i].accepted) {
            pos_cluster->add_read(r);
        }
    }
	if (r_cluster->clip_consensus && rc_remap_infos.rbegin()->accepted) {
        pos_cluster->add_clip_cluster(r_cluster->clip_consensus);
    }
    
	reads_cluster_t* neg_cluster = new reads_cluster_t();
    for (int i = 0; i < l_cluster->reads.size(); i++) {
        bam1_t* r = l_cluster->reads[i];
        update_read(r, best_region, lc_remap_infos[i], is_rc);
        kept.push_back(r);

        if (lc_remap_infos[i].accepted) {
            neg_cluster->add_read(r);
        }
    }
    if (l_cluster->clip_consensus && lc_remap_infos.rbegin()->accepted) {
        neg_cluster->add_clip_cluster(l_cluster->clip_consensus);
    }

    if (!pos_cluster->empty() && !neg_cluster->empty()) {
    	// realign corrected consensus to find insertion
		std::string full_assembled_seq, ins_seq;
		bool good_left_anchor, good_right_anchor;
		int remap_start = std::max(hts_pos_t(0), r_cluster->start()-50);
		int remap_end = std::min(l_cluster->end()+50, contigs.get_len(contig_name)-1);
		if (remap_start < remap_end) {
			std::vector<std::string> v = { corrected_consensus_sequence };
			auto p = remap_assembled_sequence(v, contigs.get_seq(contig_name)+remap_start, remap_end-remap_start, aligner_to_base,
					good_left_anchor, good_right_anchor, ins_seq, full_assembled_seq, config);
            StripedSmithWaterman::Alignment lh_aln = p.first, rh_aln = p.second;
			if (good_left_anchor && good_right_anchor && ins_seq.find("N") == std::string::npos) {
				// Get transposed sequence coverage
				int ins_seq_start = lh_aln.query_end;
				int ins_seq_end = corrected_consensus_sequence.length() - (rh_aln.query_end - rh_aln.query_begin + 1);

				int i = 0;
				int left_seq_cov = ins_seq_start;
				std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.first < p2.first;});
				while (i < covered_segments.size() && covered_segments[i].first <= left_seq_cov) {
					left_seq_cov = std::max(left_seq_cov, covered_segments[i].second);
					i++;
				}

				i = covered_segments.size()-1;
				int right_seq_cov = ins_seq_end;
				std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.second < p2.second;});
				while (i >= 0 && covered_segments[i].second >= right_seq_cov) {
					right_seq_cov = std::min(right_seq_cov, covered_segments[i].first);
					i--;
				}

				int rc_fwd_reads = 0, rc_rev_reads = 0, lc_fwd_reads = 0, lc_rev_reads = 0;
				if (pos_cluster->clip_consensus) {
					rc_fwd_reads = pos_cluster->clip_consensus->fwd_clipped;
					rc_rev_reads = pos_cluster->clip_consensus->rev_clipped;
				}
				if (neg_cluster->clip_consensus) {
					lc_fwd_reads = neg_cluster->clip_consensus->fwd_clipped;
					lc_rev_reads = neg_cluster->clip_consensus->rev_clipped;
				}

				inss_insertion_t* insertion = new inss_insertion_t(contig_name, remap_start + lh_aln.ref_end, remap_start + rh_aln.ref_begin-1,
						pos_cluster->reads.size(), neg_cluster->reads.size(), rc_fwd_reads, rc_rev_reads, lc_fwd_reads, lc_rev_reads, 0, ins_seq);
				insertion->left_anchor = std::to_string(remap_start + lh_aln.ref_begin) + "-" + std::to_string(remap_start + lh_aln.ref_end);
				insertion->right_anchor = std::to_string(remap_start + rh_aln.ref_begin) + "-" + std::to_string(remap_start + rh_aln.ref_end);
				insertion->left_bp_precise = left_bp_precise, insertion->right_bp_precise = right_bp_precise;
				for (bam1_t* read : pos_cluster->reads) {
					insertion->rc_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
				}
				insertion->rc_avg_nm /= pos_cluster->reads.size();
				for (bam1_t* read : neg_cluster->reads) {
					insertion->lc_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
				}
				insertion->lc_avg_nm /= neg_cluster->reads.size();

				insertion->left_anchor_cigar = lh_aln.cigar_string;
				insertion->right_anchor_cigar = rh_aln.cigar_string;

				insertion->left_seq_cov = std::min(left_seq_cov - ins_seq_start, (int) ins_seq.length());
				insertion->right_seq_cov = std::min(ins_seq_end - right_seq_cov, (int) ins_seq.length());
				mtx.lock();
				trans_insertions.push_back(insertion);
				mtx.unlock();
			}
		}
    }

    delete pos_cluster;
    delete neg_cluster;

    std::sort(kept.begin(), kept.end(), [] (bam1_t* r1, bam1_t* r2) {return get_endpoint(r1) < get_endpoint(r2);});
}

int find(int* parents, int i) {
    int root = i;
    while (root != parents[root]) {
        root = parents[root];
    }
    while (i != root) {
        int newp = parents[i];
        parents[i] = root;
        i = newp;
    }
    return root;
}
void merge(int* parents, int* sizes, int x, int y) {
    int i = find(parents, x);
    int j = find(parents, y);
    if (i == j) return;

    if (sizes[i] < sizes[j]) {
        parents[i] = j;
        sizes[j] += sizes[i];
    } else {
        parents[j] = i;
        sizes[i] += sizes[j];
    }
}

void remove_cluster_from_mm(std::multimap<int, cluster_t*>& mm, cluster_t* c, int pos) {
    auto bounds = mm.equal_range(pos);
    for (auto it = bounds.first; it != bounds.second; it++) {
        if (it->second == c) {
            mm.erase(it);
            break;
        }
    }
}
void remove_cluster_from_mm(std::multimap<int, cluster_t*>& mm, cluster_t* c) {
    remove_cluster_from_mm(mm, c, c->la_start);
    remove_cluster_from_mm(mm, c, c->la_end);
}

std::vector<reads_cluster_t*> cluster_reads(open_samFile_t* dc_file, int contig_id, std::unordered_map<std::string, std::string>& mateseqs,
		std::unordered_map<std::string, std::string>& matequals, std::vector<consensus_t*>& clip_consensuses) {

	std::string contig_name = contig_map.get_name(contig_id);
    hts_itr_t* iter = sam_itr_querys(dc_file->idx, dc_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    std::vector<bam1_t*> reads;
    std::vector<cluster_t*> clusters;
    while (sam_itr_next(dc_file->file, iter, read) >= 0) {
        std::string mate_seq = get_mate_seq(read, mateseqs);
        std::string mate_qual = get_mate_qual(read, matequals);

        if (mate_seq.empty() || mate_seq.find("N") != std::string::npos) continue;

        cluster_t* cluster = new cluster_t(read, config.high_confidence_mapq);
        cluster->id = clusters.size();
        cluster->ra_start = cluster->la_start;
        cluster->ra_end = cluster->la_end;

        clusters.push_back(cluster);
        reads.push_back(bam_dup1(read));
    }
    sam_itr_destroy(iter);
    bam_destroy1(read);

    if (clusters.empty()) return std::vector<reads_cluster_t*>();

    // union-find data structure
    int n_reads = clusters.size();
    int* parents = new int[n_reads], * sizes = new int[n_reads];
    for (int i = 0; i < n_reads; i++) {
        parents[i] = i;
        sizes[i] = 1;
    }
    
    // merge first equal clusters
    int curr = 0;
    for (int i = 1; i < clusters.size(); i++) {
        if (*clusters[curr] == *clusters[i]) {
            cluster_t* merged = cluster_t::merge(clusters[curr], clusters[i]);
            merged->id = curr;
            merge(parents, sizes, curr, i);
            delete clusters[curr];
            delete clusters[i];
            clusters[curr] = merged;
            clusters[i] = NULL;
        } else {
            curr = i;
        }
    }
    clusters.erase(std::remove(clusters.begin(), clusters.end(), (cluster_t*) NULL), clusters.end());

    std::multimap<int, cluster_t*> clusters_map;
    for (cluster_t* c : clusters) {
        clusters_map.insert(std::make_pair(c->la_start, c));
        clusters_map.insert(std::make_pair(c->la_end, c));
    }

    std::priority_queue<cc_distance_t> pq;
    for (cluster_t* c1 : clusters) {
        auto end = clusters_map.upper_bound(c1->la_end+stats.max_is);
        for (auto map_it = clusters_map.lower_bound(c1->la_start); map_it != end; map_it++) {
            cluster_t* c2 = map_it->second;
            if (c1 != c2 && cluster_t::can_merge(c1, c2, stats.max_is) && c1->la_start <= c2->la_start) {
                pq.push(cc_distance_t(c1, c2));
            }
        }
    }

    while (!pq.empty()) {
        cc_distance_t ccd = pq.top();
        pq.pop();

        if (ccd.c1->used || ccd.c2->used) continue;

        cluster_t* new_cluster = cluster_t::merge(ccd.c1, ccd.c2);
        new_cluster->id = std::min(ccd.c1->id, ccd.c2->id);
        merge(parents, sizes, ccd.c1->id, ccd.c2->id);
        clusters.push_back(new_cluster);

        ccd.c1->used = true;
        remove_cluster_from_mm(clusters_map, ccd.c1);
        ccd.c2->used = true;
        remove_cluster_from_mm(clusters_map, ccd.c2);

        auto end = clusters_map.upper_bound(new_cluster->la_end + stats.max_is);
        for (auto map_it = clusters_map.lower_bound(new_cluster->la_start - stats.max_is);
             map_it != end; map_it++) {
            if (!map_it->second->used && cluster_t::can_merge(new_cluster, map_it->second, stats.max_is)) {
                pq.push(cc_distance_t(new_cluster, map_it->second));
            }
        }
        clusters_map.insert(std::make_pair(new_cluster->la_start, new_cluster));
        clusters_map.insert(std::make_pair(new_cluster->la_end, new_cluster));
    }

    // for each set of reads, make a cluster-vector
    std::vector<reads_cluster_t*> read_clusters(n_reads);
    for (int i = 0; i < n_reads; i++) read_clusters[i] = new reads_cluster_t();
    for (int i = 0; i < n_reads; i++) {
        read_clusters[find(parents, i)]->add_read(reads[i]);
    }

    if (!clip_consensuses.empty()) { // TODO: temporary, enhance logic
        bool left_facing = clip_consensuses[0]->left_clipped;

        std::sort(clip_consensuses.begin(), clip_consensuses.end(), [](consensus_t* c1, consensus_t* c2) {
            return c1->breakpoint < c2->breakpoint;
        });
        
        int curr_j = 0;
        std::vector<consensus_t*> clip_consensus_per_cluster(read_clusters.size(), NULL);
        if (left_facing) {
            std::sort(read_clusters.begin(), read_clusters.end(), [](reads_cluster_t* rc1, reads_cluster_t* rc2) {
                return rc1->end() < rc2->end();
            });

            for (consensus_t* clip_consensus : clip_consensuses) {
                while (curr_j < read_clusters.size() && read_clusters[curr_j]->end() < clip_consensus->breakpoint) curr_j++;

                for (int j = curr_j; j < read_clusters.size(); j++) {
                    if (read_clusters[j]->end() - clip_consensus->breakpoint <= stats.max_is) {
                        if (clip_consensus_per_cluster[j] == NULL || clip_consensus_per_cluster[j]->supp_clipped_reads() < clip_consensus->supp_clipped_reads()) {
                            clip_consensus_per_cluster[j] = clip_consensus;
                        }
                    } else break;
                }
            }
        } else {
            std::sort(read_clusters.begin(), read_clusters.end(), [](reads_cluster_t* rc1, reads_cluster_t* rc2) {
                return rc1->start() < rc2->start();
            });

            for (consensus_t* clip_consensus : clip_consensuses) {
                while (curr_j < read_clusters.size() && clip_consensus->breakpoint-read_clusters[curr_j]->start() > stats.max_is) curr_j++;

                for (int j = curr_j; j < read_clusters.size(); j++) {
                    if (clip_consensus->breakpoint - read_clusters[j]->start() >= 0) {
                        if (clip_consensus_per_cluster[j] == NULL || clip_consensus_per_cluster[j]->supp_clipped_reads() < clip_consensus->supp_clipped_reads()) {
                            clip_consensus_per_cluster[j] = clip_consensus;
                        }
                    } else break;
                }
            }
        }

        for (int j = 0; j < read_clusters.size(); j++) {
            read_clusters[j]->add_clip_cluster(clip_consensus_per_cluster[j]);
        }
    }

    for (int i = 0; i < read_clusters.size(); i++) {
        if (read_clusters[i]->reads.size() <= 1) {
            read_clusters[i]->deallocate_reads();
            delete read_clusters[i];
            read_clusters[i] = nullptr;
        }
    }
    read_clusters.erase(std::remove(read_clusters.begin(), read_clusters.end(), nullptr), read_clusters.end());

    std::sort(read_clusters.begin(), read_clusters.end(), [](const reads_cluster_t* rc1, const reads_cluster_t* rc2) {
    	return rc1->reads.size() > rc2->reads.size();
    });

    for (cluster_t* c : clusters) delete c;
    delete[] parents;
    delete[] sizes;

    return read_clusters;
}

bool is_semi_mapped(bam1_t* read) {
	return is_proper_pair(read, stats.max_is) && !is_mate_clipped(read);
}

void add_semi_mapped(std::string clipped_fname, int contig_id, std::vector<reads_cluster_t*>& r_clusters, std::vector<reads_cluster_t*>& l_clusters) {

	if (!file_exists(clipped_fname)) return;

	std::vector<bam1_t*> l_semi_mapped_pairs, r_semi_mapped_pairs;

	open_samFile_t* clipped_file = open_samFile(clipped_fname.c_str(), true);
	std::string contig_name = contig_map.get_name(contig_id);
	hts_itr_t* iter = sam_itr_querys(clipped_file->idx, clipped_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(clipped_file->file, iter, read) >= 0) {
		if (is_semi_mapped(read) && get_sequence(read).find("N") == std::string::npos) {
			if (read->core.flag & BAM_FMREVERSE) {
				l_semi_mapped_pairs.push_back(bam_dup1(read));
			} else {
				r_semi_mapped_pairs.push_back(bam_dup1(read));
			}
		}
	}
	close_samFile(clipped_file);

	// add semimapped pairs
	std::vector<Interval<bam1_t*>> r_iv, l_iv;
	for (bam1_t* r : r_semi_mapped_pairs) {
		r_iv.push_back(Interval<bam1_t*>(r->core.mpos, get_mate_endpos(r), r));
	}
	for (bam1_t* r : l_semi_mapped_pairs) {
		l_iv.push_back(Interval<bam1_t*>(r->core.mpos, get_mate_endpos(r), r));
	}
	IntervalTree<bam1_t*> r_it(r_iv), l_it(l_iv);
	std::unordered_set<bam1_t*> used_sm_reads;
	for (reads_cluster_t* rc : r_clusters) {
		int end = rc->end(), start = end - stats.max_is;
		auto sm_reads = r_it.findContained(start, end);
		for (auto i : sm_reads) {
			bam1_t* r = i.value;
			if (!used_sm_reads.count(r)) {
				rc->add_semi_mapped_reads(bam_dup1(r));
				used_sm_reads.insert(r);
			}
		}
	}
	for (reads_cluster_t* rc : l_clusters) {
		int start, end;
		start = rc->start(); end = start + stats.max_is;
		auto sm_reads = l_it.findContained(start, end);
		for (auto i : sm_reads) {
			bam1_t* r = i.value;
			if (!used_sm_reads.count(r)) {
				rc->add_semi_mapped_reads(bam_dup1(r));
				used_sm_reads.insert(r);
			}
		}
	}
	for (bam1_t* r : r_semi_mapped_pairs) {
		bam_destroy1(r);
	}
	for (bam1_t* r : l_semi_mapped_pairs) {
		bam_destroy1(r);
	}
}

void remap(int id, int contig_id) {

    std::string r_dc_fname = workdir + "/workspace/R/" + std::to_string(contig_id) + ".noremap.bam";
    std::string l_dc_fname = workdir + "/workspace/L/" + std::to_string(contig_id) + ".noremap.bam";
    if (!file_exists(r_dc_fname) || !file_exists(l_dc_fname)) return;

    std::unordered_map<std::string, std::string> mateseqs;
    std::unordered_map<std::string, std::string> matequals;
    std::ifstream mateseqs_fin(workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt");
    std::string qname, seq, qual; int mapq;
    while (mateseqs_fin >> qname >> seq >> qual >> mapq) {
        mateseqs[qname] = seq;
        matequals[qname] = qual;
    }
    mateseqs_fin.close();

    open_samFile_t* r_dc_file = open_samFile(r_dc_fname.c_str(), true);
    open_samFile_t* l_dc_file = open_samFile(l_dc_fname.c_str(), true);

    std::string clip_consensus_fname = workdir + "/workspace/sr_consensuses/" + std::to_string(contig_id) + ".txt";
    std::ifstream clipped_fin(clip_consensus_fname);
    std::vector<consensus_t*> rc_consensuses, lc_consensuses;
    if (clipped_fin.good()) {
    	std::string contig_name, seq;
    	char dir;
    	hts_pos_t start, end, breakpoint;
    	int fwd_clipped, rev_clipped;
        int max_mapq, lowq_clip_portion;
        hts_pos_t remap_boundary;
    	while (clipped_fin >> contig_name >> start >> end >> breakpoint >> dir >> seq >> 
            fwd_clipped >> rev_clipped >> max_mapq >> remap_boundary >> lowq_clip_portion) {
			if (dir == 'L') {
                start += lowq_clip_portion;
                seq = seq.substr(lowq_clip_portion);
                consensus_t* consensus = new consensus_t(true, start, breakpoint, end, seq, fwd_clipped, rev_clipped, breakpoint-start, max_mapq, remap_boundary, lowq_clip_portion);
                lc_consensuses.push_back(consensus);
			} else {
                end -= lowq_clip_portion;
                seq = seq.substr(0, seq.length()-lowq_clip_portion);
                consensus_t* consensus = new consensus_t(false, start, breakpoint, end, seq, fwd_clipped, rev_clipped, end-breakpoint, max_mapq, remap_boundary, lowq_clip_portion);
                rc_consensuses.push_back(consensus);
			}
    	}
    }

	std::vector<reads_cluster_t*> r_clusters = cluster_reads(r_dc_file, contig_id, mateseqs, matequals, rc_consensuses);
	std::vector<reads_cluster_t*> l_clusters = cluster_reads(l_dc_file, contig_id, mateseqs, matequals, lc_consensuses);

	std::string clipped_fname = workdir + "/workspace/clipped/" + std::to_string(contig_id) + ".bam";
	add_semi_mapped(clipped_fname, contig_id, r_clusters, l_clusters);

    std::vector<bam1_t*> r_reads_to_write, l_reads_to_write;

    std::priority_queue<cc_v_distance_t> pq;
    auto score_f = [](const reads_cluster_t* v1, const reads_cluster_t* v2) {return v1->reads.size()*v2->reads.size();};

    std::multimap<int, reads_cluster_t*> l_clusters_map;
    for (reads_cluster_t* l_cluster : l_clusters) {
        hts_pos_t cl_start = l_cluster->start();
        l_clusters_map.insert({cl_start, l_cluster});
    }
    for (reads_cluster_t* r_cluster : r_clusters) {
        hts_pos_t cl_end = r_cluster->end();
        auto begin = l_clusters_map.lower_bound(cl_end - stats.max_is);
        auto end = l_clusters_map.upper_bound(cl_end + stats.max_is);

        reads_cluster_t* l_cluster = NULL;
        for (auto it = begin; it != end; it++) {
            if (l_cluster == NULL || it->second->reads.size() > l_cluster->reads.size()) {
                l_cluster = it->second;
            }
        }
        if (l_cluster == NULL) continue;

        pq.push(cc_v_distance_t(r_cluster, l_cluster, score_f(r_cluster, l_cluster)));
    }

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Aligner permissive_aligner(2, 2, 4, 1, false);
    StripedSmithWaterman::Aligner aligner_to_base(1, 4, 6, 1, true);
    StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, true);
    while (!pq.empty()) {
        cc_v_distance_t cc_v_distance = pq.top();
        pq.pop();

        reads_cluster_t* c1 = cc_v_distance.c1;
        reads_cluster_t* c2 = cc_v_distance.c2;
        if (c1->used || c2->used) continue;

        // remap clusters
        std::vector<bam1_t*> to_write;
		remap_cluster(c1, c2, to_write, contig_id, r_dc_file->header, mateseqs, matequals, aligner, permissive_aligner,
				aligner_to_base, harsh_aligner, stats);

        for (bam1_t* r : to_write) {
            if (bam_is_rev(r)) {
                l_reads_to_write.push_back(bam_dup1(r));
            } else {
                r_reads_to_write.push_back(bam_dup1(r));
            }
        };
        c1->used = true; c2->used = true;
        c1->deallocate_reads();
        c2->deallocate_reads();
    }

    std::string r_dc_remapped_fname = workdir + "/workspace/R/" + std::to_string(contig_id) + ".remap.bam";
    write_and_index_file(r_reads_to_write, r_dc_remapped_fname, r_dc_file->header);
    for (bam1_t* r : r_reads_to_write) bam_destroy1(r);

    std::string l_dc_remapped_fname = workdir + "/workspace/L/" + std::to_string(contig_id) + ".remap.bam";
    write_and_index_file(l_reads_to_write, l_dc_remapped_fname, l_dc_file->header);
    for (bam1_t* r : l_reads_to_write) bam_destroy1(r);

    for (reads_cluster_t* c : r_clusters) c->deallocate_reads();
    for (reads_cluster_t* c : l_clusters) c->deallocate_reads();

    close_samFile(l_dc_file);
    close_samFile(r_dc_file);
}

int main(int argc, char* argv[]) {

    workdir = std::string(argv[1]);
    std::string workspace = workdir + "/workspace";
    std::string reference_fname  = argv[2];
    std::string sample_name  = argv[3];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
    std::ifstream full_cmd_fin(full_cmd_fname);
    std::string full_cmd_str;
    std::getline(full_cmd_fin, full_cmd_str);

    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contigs.read_fasta_into_map(reference_fname);
    contig_map.load(workdir);

    assembly_failed_no_seq.open(workdir + "/assembly_failed.no_seq.sv");
    assembly_failed_cycle_writer.open(workdir + "/assembly_failed.w_cycle.sv");
    assembly_failed_too_many_reads_writer.open(workdir + "/assembly_failed.too_many_reads.sv");
    assembly_failed_bad_anchors_writer.open(workdir + "/assembly_failed.bad_anchors.sv");
    assembly_failed_lt50bp.open(workdir + "/assembly_failed.lt50bp.sv");
    assembly_succeeded.open(workdir + "/assembly_succeeded.sv");

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::future<void> future = thread_pool.push(remap, contig_id);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    // out assembled file
	std::string assembled_out_vcf_fname = workdir + "/assembled_ins.vcf.gz";
	htsFile* assembled_out_vcf_file = bcf_open(assembled_out_vcf_fname.c_str(), "wz");
	if (assembled_out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + assembled_out_vcf_fname + ".");
	}
	bcf_hdr_t* out_vcf_header = generate_vcf_header(contigs, sample_name, config, full_cmd_str);
	if (bcf_hdr_write(assembled_out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + assembled_out_vcf_fname + ".");
	}

	std::sort(assembled_insertions.begin(), assembled_insertions.end(), [&out_vcf_header](const inss_insertion_t* i1, const inss_insertion_t* i2) {
		int contig_id1 = bcf_hdr_name2id(out_vcf_header, i1->ins->chr.c_str());
		int contig_id2 = bcf_hdr_name2id(out_vcf_header, i2->ins->chr.c_str());
		int disc_score1 = -(i1->r_disc_pairs*i1->l_disc_pairs), disc_score2 = -(i2->r_disc_pairs*i2->l_disc_pairs);
		// negative because we want descending order
		int disc_score_mul1 = -(i1->r_disc_pairs*i1->l_disc_pairs), disc_score_mul2 = -(i2->r_disc_pairs*i2->l_disc_pairs);
		int disc_score_sum1 = -(i1->r_disc_pairs+i1->l_disc_pairs), disc_score_sum2 = -(i2->r_disc_pairs+i2->l_disc_pairs);
		return std::tie(contig_id1, i1->ins->start, i1->ins->end, i1->ins->ins_seq, disc_score_mul1, disc_score_sum1) <
			   std::tie(contig_id2, i2->ins->start, i2->ins->end, i2->ins->ins_seq, disc_score_mul2, disc_score_sum2);
	});

	bcf1_t* bcf_entry = bcf_init();
	int a_id = 0;
	for (inss_insertion_t* insertion : assembled_insertions) {
		std::string id = "A_INS_" + std::to_string(a_id++);
		insertion_to_bcf_entry(insertion, out_vcf_header, bcf_entry, id, contigs);

		int int2_conv[2];
		int2_conv[0] = insertion->r_disc_pairs, int2_conv[1] = insertion->l_disc_pairs;
		bcf_update_info_int32(out_vcf_header, bcf_entry, "DISCORDANT", int2_conv, 2);
		bcf_update_info_string(out_vcf_header, bcf_entry, "ALGORITHM", "assembly");
		bcf_update_info_string(out_vcf_header, bcf_entry, "LEFT_ANCHOR", insertion->left_anchor.c_str());
		bcf_update_info_string(out_vcf_header, bcf_entry, "RIGHT_ANCHOR", insertion->right_anchor.c_str());
		float f2_conv[2];
		f2_conv[0] = insertion->rc_avg_nm, f2_conv[1] = insertion->lc_avg_nm;
		bcf_update_info_float(out_vcf_header, bcf_entry, "AVG_STABLE_NM", f2_conv, 2);

		if (bcf_write(assembled_out_vcf_file, out_vcf_header, bcf_entry) != 0) {
			throw std::runtime_error("Failed to write to " + assembled_out_vcf_fname + ".");
		}
	}

	// out transurveyor insertions
	std::string transurveyor_out_vcf_fname = workdir + "/transurveyor_ins.vcf.gz";
	htsFile* transurveyor_out_vcf_file = bcf_open(transurveyor_out_vcf_fname.c_str(), "wz");
	if (transurveyor_out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + transurveyor_out_vcf_fname + ".");
	}
	if (bcf_hdr_write(transurveyor_out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + transurveyor_out_vcf_fname + ".");
	}

	std::sort(trans_insertions.begin(), trans_insertions.end(), [&out_vcf_header](const inss_insertion_t* i1, const inss_insertion_t* i2) {
		int contig_id1 = bcf_hdr_name2id(out_vcf_header, i1->ins->chr.c_str());
		int contig_id2 = bcf_hdr_name2id(out_vcf_header, i2->ins->chr.c_str());
		// negative because we want descending order
		int disc_score_mul1 = -(i1->r_disc_pairs*i1->l_disc_pairs), disc_score_mul2 = -(i2->r_disc_pairs*i2->l_disc_pairs);
		int disc_score_sum1 = -(i1->r_disc_pairs+i1->l_disc_pairs), disc_score_sum2 = -(i2->r_disc_pairs+i2->l_disc_pairs);
		return std::tie(contig_id1, i1->ins->start, i1->ins->end, i1->ins->ins_seq, disc_score_mul1, disc_score_sum1) <
			   std::tie(contig_id2, i2->ins->start, i2->ins->end, i2->ins->ins_seq, disc_score_mul2, disc_score_sum2);
	});

	int t_id = 0;
	for (inss_insertion_t* insertion : trans_insertions) {
		std::string id = "T_INS_" + std::to_string(t_id++);
		insertion_to_bcf_entry(insertion, out_vcf_header, bcf_entry, id, contigs);

		int int2_conv[2];
		int2_conv[0] = insertion->r_disc_pairs, int2_conv[1] = insertion->l_disc_pairs;
		bcf_update_info_int32(out_vcf_header, bcf_entry, "DISCORDANT", int2_conv, 2);
		bcf_update_info_string(out_vcf_header, bcf_entry, "ALGORITHM", "transurveyor");
		bcf_update_info_flag(out_vcf_header, bcf_entry, "IMPRECISE", "", !insertion->left_bp_precise || !insertion->right_bp_precise);
		bcf_update_info_string(out_vcf_header, bcf_entry, "LEFT_ANCHOR", insertion->left_anchor.c_str());
		bcf_update_info_string(out_vcf_header, bcf_entry, "RIGHT_ANCHOR", insertion->right_anchor.c_str());
		bcf_update_info_string(out_vcf_header, bcf_entry, "LEFT_ANCHOR_CIGAR", insertion->left_anchor_cigar.c_str());
		bcf_update_info_string(out_vcf_header, bcf_entry, "RIGHT_ANCHOR_CIGAR", insertion->right_anchor_cigar.c_str());
		float f2_conv[2];
		f2_conv[0] = insertion->rc_avg_nm, f2_conv[1] = insertion->lc_avg_nm;
		bcf_update_info_float(out_vcf_header, bcf_entry, "AVG_STABLE_NM", f2_conv, 2);
		int2_conv[0] = insertion->left_seq_cov, int2_conv[1] = insertion->right_seq_cov;
		bcf_update_info_int32(out_vcf_header, bcf_entry, "TRANS_QUERY_COV", int2_conv, 2);

		if (bcf_write(transurveyor_out_vcf_file, out_vcf_header, bcf_entry) != 0) {
			throw std::runtime_error("Failed to write to " + transurveyor_out_vcf_fname + ".");
		}
	}

    bcf_close(assembled_out_vcf_file);
    bcf_close(transurveyor_out_vcf_file);
}
