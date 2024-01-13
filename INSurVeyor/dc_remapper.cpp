#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <queue>
#include <unistd.h>
#include <climits>
#include <random>
#include <cassert>

#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <vector>

#include "../libs/cptl_stl.h"
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "../libs/IntervalTree.h"
#include "../src/dc_remapper.h"
#include "../src/sam_utils.h"
#include "../src/clustering_utils.h"
#include "remapping.h"
#include "utils.h"
#include "vcf_utils.h"
#include "../src/assemble.h"
#include "reference_guided_assembly.h"
#include "../src/vcf_utils.h"


config_t config;
stats_t stats;
std::string workdir;
std::mutex mtx;

contig_map_t contig_map;
chr_seqs_map_t contigs;

std::ofstream assembly_failed_no_seq, assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer;

std::vector<sv_t*> insertions;

const int SMALL_SAMPLE_SIZE = 15;
const int CLUSTER_CANDIDATES = 3;
const double BASE_ACCEPTANCE_THRESHOLD = 0.95;
const int TOO_MANY_READS = 1000;

struct cc_v_distance_t {
    insertion_cluster_t* c1, * c2;
    int distance;

    cc_v_distance_t(insertion_cluster_t* c1, insertion_cluster_t* c2, int distance) : c1(c1), c2(c2), distance(distance) {}
};
bool operator < (const cc_v_distance_t& ccd1, const cc_v_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance < ccd2.distance;
}

remap_info_t remap_seq(char* region, int region_len, std::string& seq, StripedSmithWaterman::Aligner& aligner,
                StripedSmithWaterman::Filter& filter, int min_clip_len) {
    if (region_len <= 0) {
        return remap_info_t();
    }

    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), region, region_len, filter, &alignment, 0);
    bool accepted = alignment.sw_score >= 30; //s.length(); TODO

    return remap_info_t(alignment, accepted, min_clip_len);
}
remap_info_t remap_seq2(char* region, int region_len, std::string& seq, StripedSmithWaterman::Aligner& aligner,
                StripedSmithWaterman::Filter& filter, int min_clip_len) {
    if (region_len <= 0) {
        return remap_info_t();
    }

    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), region, region_len, filter, &alignment, 0);
    bool accepted = 2*alignment.sw_score >= seq.length() && alignment.query_end-alignment.query_begin >= min_clip_len; 

    return remap_info_t(alignment, accepted, min_clip_len);
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
		remap_info_t remap_info = remap_seq(region_ptr+strictest_region_start, strictest_region_len, s, aligner, filter, config.min_clip_len);
		remap_info.start += strictest_region_start;
		remap_info.end += strictest_region_start;
		remap_infos.push_back(remap_info);
	}
	if (!clip_seq.empty()) {
		clip_remap_info = remap_seq(region_ptr+strictest_region_start, strictest_region_len, clip_seq, aligner, filter, config.min_clip_len);
		clip_remap_info.start += strictest_region_start;
		clip_remap_info.end += strictest_region_start;
	}
}

region_score_t compute_score_supp(region_t& region, insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster,
                       std::unordered_map<std::string, std::string>& mateseqs,
                       std::vector<remap_info_t>* rc_remap_infos_, std::vector<remap_info_t>* lc_remap_infos_,
                       StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                       StripedSmithWaterman::Filter& filter, bool do_rc) {

    // mateseqs contains seqs are stored as in fasta/q file
    std::vector<std::string> r_mates, l_mates;
    for (bam1_t* read : r_cluster->cluster->reads) {
        std::string s = get_mate_seq(read, mateseqs);
        if (!do_rc) rc(s);   // by default, we map mates of fwd reads on the rev strand
        r_mates.push_back(s);
    }
    for (bam1_t* read : l_cluster->cluster->reads) {
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
        remap_info_t remap_info = remap_seq2(region_ptr, region.end-region.start, s, aligner, filter, config.min_clip_len);
        if (remap_info.accepted) {
            score.total_score += remap_info.score;
            if (!do_rc) {
                strict_region_start = remap_info.start;
            } else if (do_rc) {
                strict_region_end = remap_info.end;
            }
        }
        rc_remap_info = remap_info;
    }
    if (l_cluster->clip_consensus) {
        std::string s = l_cluster->clip_consensus->clip_sequence();
        if (do_rc) rc(s);
        remap_info_t remap_info = remap_seq2(region_ptr, region.end-region.start, s, aligner, filter, config.min_clip_len);
        if (remap_info.accepted) {
            score.total_score += remap_info.score;
            if (do_rc) {
                strict_region_start = remap_info.start;
            } else if (!do_rc) {
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
        remap_info_t remap_info = remap_seq(region_ptr+strict_region_start, strict_region_len, s, aligner, filter, config.min_clip_len);
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
        remap_info_t remap_info = remap_seq(region_ptr+strict_region_start, strict_region_len, s, aligner, filter, config.min_clip_len);
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
        if (remap_info.left_clipped && strict_region_start < remap_info.start-config.max_clipped_pos_dist) remap_info.accepted = false;
        if (remap_info.start < strict_region_start-config.max_clipped_pos_dist) remap_info.accepted = false;
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

void compute_score(region_t& region, insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster,
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

insertion_cluster_t* subsample_cluster(insertion_cluster_t* reads_cluster, int size, std::default_random_engine& rng) {
    std::vector<bam1_t*> subset(reads_cluster->cluster->reads);
    if (subset.size() > size) {
        std::shuffle(subset.begin(), subset.end(), rng);
        subset.erase(subset.begin() + size, subset.end());
    }
    insertion_cluster_t* subsampled_cluster = new insertion_cluster_t(new cluster_t());
    for (bam1_t* read : subset) subsampled_cluster->add_stable_read(read);
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

bool remap_cluster(insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster, std::vector<bam1_t*>& kept,
                   int contig_id, bam_hdr_t* header, std::unordered_map<std::string, std::string>& mateseqs,
				   std::unordered_map<std::string, std::string>& matequals,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                   StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
                   stats_t& stats) {

	if (l_cluster->cluster->reads.size() + r_cluster->cluster->reads.size() > TOO_MANY_READS) return false;

    std::vector<region_t> regions;

    /* == Find candidate regions for remapping == */
    std::vector<bam1_t*> full_cluster;
    full_cluster.insert(full_cluster.end(), l_cluster->cluster->reads.begin(), l_cluster->cluster->reads.end());
    full_cluster.insert(full_cluster.end(), r_cluster->cluster->reads.begin(), r_cluster->cluster->reads.end());
    std::sort(full_cluster.begin(), full_cluster.end(), [] (bam1_t* r1, bam1_t* r2) {
        if (r1->core.mtid != r2->core.mtid) return r1->core.mtid < r2->core.mtid;
        else return r1->core.mpos < r2->core.mpos;
    });

    std::vector<bam1_t*> subcluster;
    for (bam1_t* r : full_cluster) {
        if (!subcluster.empty() && (subcluster[0]->core.mtid != r->core.mtid ||
                r->core.mpos-subcluster[0]->core.mpos > stats.max_is)) {
            std::string m_contig_name = header->target_name[subcluster[0]->core.mtid];
            regions.push_back(get_candidate_region(subcluster, m_contig_name, contig_map.get_id(m_contig_name), contigs.get_len(m_contig_name), stats.max_is, config.max_trans_size));
            subcluster.clear();
        }
        if (!is_mate_unmapped(r)) subcluster.push_back(r);
    }
    if (!subcluster.empty()) {
        std::string m_contig_name = header->target_name[subcluster[0]->core.mtid];
        regions.push_back(get_candidate_region(subcluster, m_contig_name, contig_map.get_id(m_contig_name), contigs.get_len(m_contig_name), stats.max_is, config.max_trans_size));
    }

    std::string contig_name = contig_map.get_name(contig_id);
    if (regions.empty()) {
		sv_t* ins = detect_de_novo_insertion(contig_name, contigs, r_cluster, l_cluster, mateseqs, matequals,
				aligner_to_base, harsh_aligner, kept, config, stats);

		if (ins != NULL) {
			mtx.lock();
            insertions.push_back(ins);
			mtx.unlock();
            return true;
		}
		return false;
    }

    /* == Find best candidate region == */
    StripedSmithWaterman::Filter filter, filter_w_cigar;
    filter_w_cigar.report_cigar = true;
    bool is_rc;

    // if too many regions and too many reads, subsample
    if (full_cluster.size() > 2*SMALL_SAMPLE_SIZE && regions.size() > CLUSTER_CANDIDATES) {
    	std::default_random_engine rng(config.seed);
        insertion_cluster_t* subsampled_lc = subsample_cluster(l_cluster, SMALL_SAMPLE_SIZE, rng);
        insertion_cluster_t* subsampled_rc = subsample_cluster(r_cluster, SMALL_SAMPLE_SIZE, rng);

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
    int start = r_cluster->start - stats.max_is;
    int end = l_cluster->end + stats.max_is;
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
		int tot_reads = r_cluster->cluster->reads.size() + l_cluster->cluster->reads.size();
		if (double(accepted_reads)/tot_reads >= 0.5) {
			return false;
		}
    }

    std::vector<remap_info_t> rc_remap_infos, lc_remap_infos;
    compute_score(best_region, r_cluster, l_cluster, mateseqs, &rc_remap_infos, &lc_remap_infos,
                  aligner, permissive_aligner, filter_w_cigar, is_rc);

    // build corrected consensus sequence and realign reads to it
    bool left_bp_precise = false, right_bp_precise = false;
    std::string corrected_consensus_sequence = generate_consensus_sequences(contig_name, contigs, contig_map, r_cluster, l_cluster, rc_remap_infos, lc_remap_infos,
				best_region, is_rc, left_bp_precise, right_bp_precise, mateseqs, matequals, aligner, harsh_aligner, config, stats);

    /* == Try assembly == */
	int rc_accepted_reads = 0, lc_accepted_reads = 0;
	for (remap_info_t& remap_info : rc_remap_infos) rc_accepted_reads += remap_info.accepted;
	for (remap_info_t& remap_info : lc_remap_infos) lc_accepted_reads += remap_info.accepted;
	int tot_reads = r_cluster->cluster->reads.size() + l_cluster->cluster->reads.size();
	if (rc_accepted_reads == 0 || lc_accepted_reads == 0 || double(rc_accepted_reads+lc_accepted_reads)/tot_reads < 0.5) {
		sv_t* ins = detect_de_novo_insertion(contig_name, contigs, r_cluster, l_cluster, mateseqs, matequals,
				aligner_to_base, harsh_aligner, kept, config, stats);

		if (ins != NULL) {
			mtx.lock();
            insertions.push_back(ins);
			mtx.unlock();
			return true;
		}
	}

	if (corrected_consensus_sequence.empty()) return false;

	insertion_cluster_t* pos_cluster = new insertion_cluster_t(new cluster_t());
    for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
        bam1_t* r = r_cluster->cluster->reads[i];
        if (rc_remap_infos[i].accepted) {
            pos_cluster->add_stable_read(r);
        }
    }
	if (r_cluster->clip_consensus && rc_remap_infos.rbegin()->accepted) {
        pos_cluster->add_clip_cluster(r_cluster->clip_consensus);
    }

	insertion_cluster_t* neg_cluster = new insertion_cluster_t(new cluster_t());
    for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
        bam1_t* r = l_cluster->cluster->reads[i];
        if (lc_remap_infos[i].accepted) {
            neg_cluster->add_stable_read(r);
        }
    }
    if (l_cluster->clip_consensus && lc_remap_infos.rbegin()->accepted) {
        neg_cluster->add_clip_cluster(l_cluster->clip_consensus);
    }

    bool success = false;
    if (!pos_cluster->empty() && !neg_cluster->empty()) {
		int remap_start = std::max(hts_pos_t(0), r_cluster->start-50);
		int remap_end = std::min(l_cluster->end+50, contigs.get_len(contig_name)-1);
		if (remap_start < remap_end) {
            std::vector<sv_t*> inss = detect_svs_from_junction(contig_name, contigs.get_seq(contig_name), corrected_consensus_sequence, 
                remap_start, remap_end, remap_start, remap_end, aligner, config.min_clip_len);
			if (!inss.empty() && inss[0]->svtype() == "INS") {
                sv_t* main_sv = inss[0];
                std::vector<std::pair<int,int>> covered_segments;
                StripedSmithWaterman::Alignment aln;
                for (remap_info_t ri : rc_remap_infos) {
                    if (ri.accepted) {
                        covered_segments.push_back({ri.start, ri.end});
                    }
                }
                for (remap_info_t ri : lc_remap_infos) {
                    if (ri.accepted) {
                        covered_segments.push_back({ri.start, ri.end});
                    }
                }

                // Get transposed sequence coverage
                int ins_seq_start = main_sv->left_anchor_aln->seq_len - main_sv->prefix_mh_len;
                int ins_seq_end = corrected_consensus_sequence.length() - (main_sv->right_anchor_aln->seq_len - main_sv->suffix_mh_len);

                int i = 0;
                std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.first < p2.first;});
                while (i < covered_segments.size() && covered_segments[i].second <= ins_seq_start) i++; // find left-most segment that either overlaps or is fully contained in the inserted sequence
                int prefix_cov_start = 0, prefix_cov_end = 0;
                if (i < covered_segments.size()) {
                    prefix_cov_start = std::max(ins_seq_start, covered_segments[i].first);
                    prefix_cov_end = covered_segments[i].second;
                    while (i < covered_segments.size() && covered_segments[i].first <= prefix_cov_end) {
                        prefix_cov_end = std::max(prefix_cov_end, covered_segments[i].second);
                        i++;
                    }
                    prefix_cov_start -= ins_seq_start;
                    prefix_cov_end -= ins_seq_start;
                    if (prefix_cov_end >= main_sv->ins_seq.length()) prefix_cov_end = main_sv->ins_seq.length()-1;
                }

                i = covered_segments.size()-1;
                std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.second < p2.second;});
                while (i >= 0 && covered_segments[i].first >= ins_seq_end) i--; // find right-most segment that either overlaps or is fully contained in the inserted sequence
                int suffix_cov_start = 0, suffix_cov_end = 0;
                if (i >= 0) {
                    suffix_cov_end = std::min(ins_seq_end-1, covered_segments[i].second);
                    suffix_cov_start = covered_segments[i].first;
                    while (i >= 0 && covered_segments[i].second >= suffix_cov_start) {
                        suffix_cov_start = std::min(suffix_cov_start, covered_segments[i].first);
                        i--;
                    }
                    suffix_cov_start -= ins_seq_start;
                    if (suffix_cov_start < 0) suffix_cov_start = 0;
                    suffix_cov_end -= ins_seq_start;
                }

                if (pos_cluster->clip_consensus) main_sv->rc_consensus = pos_cluster->clip_consensus;
                if (neg_cluster->clip_consensus) main_sv->lc_consensus = neg_cluster->clip_consensus;
                main_sv->disc_pairs_lf = pos_cluster->cluster->reads.size();
                main_sv->disc_pairs_rf = neg_cluster->cluster->reads.size();
                for (bam1_t* read : pos_cluster->cluster->reads) {
                    if (read->core.qual >= config.high_confidence_mapq) main_sv->disc_pairs_lf_high_mapq++;
                }
                for (bam1_t* read : neg_cluster->cluster->reads) {
                    if (read->core.qual >= config.high_confidence_mapq) main_sv->disc_pairs_rf_high_mapq++;
                }
                main_sv->disc_pairs_lf_maxmapq = pos_cluster->cluster->max_mapq;
                main_sv->disc_pairs_rf_maxmapq = neg_cluster->cluster->max_mapq;
                for (bam1_t* read : pos_cluster->cluster->reads) {
                    main_sv->disc_pairs_lf_avg_nm += get_nm(read);
                }
                main_sv->disc_pairs_lf_avg_nm /= std::max(1, (int) pos_cluster->cluster->reads.size());
                for (bam1_t* read : neg_cluster->cluster->reads) {
                    main_sv->disc_pairs_rf_avg_nm += get_nm(read);
                }
                main_sv->disc_pairs_rf_avg_nm /= std::max(1, (int) neg_cluster->cluster->reads.size());

                ((insertion_t*) main_sv)->imprecise_bp = !left_bp_precise || !right_bp_precise;
                ((insertion_t*) main_sv)->prefix_cov_start = prefix_cov_start;
                ((insertion_t*) main_sv)->prefix_cov_end = prefix_cov_end;
                ((insertion_t*) main_sv)->suffix_cov_start = suffix_cov_start;
                ((insertion_t*) main_sv)->suffix_cov_end = suffix_cov_end;
                main_sv->source = "REFERENCE_GUIDED_ASSEMBLY";

                mtx.lock();
                insertions.push_back(main_sv);
                mtx.unlock();
                success = true;

                for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
                    bam1_t* r = r_cluster->cluster->reads[i];
                    update_read(r, best_region, rc_remap_infos[i], is_rc);
                    kept.push_back(r);
                }
                for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
                    bam1_t* r = l_cluster->cluster->reads[i];
                    update_read(r, best_region, lc_remap_infos[i], is_rc);
                    kept.push_back(r);
                }
            }
		}
    }

    delete pos_cluster;
    delete neg_cluster;

    std::sort(kept.begin(), kept.end(), [] (bam1_t* r1, bam1_t* r2) {return get_endpoint(r1) < get_endpoint(r2);});

    return success;
}

std::vector<insertion_cluster_t*> cluster_reads(open_samFile_t* dc_file, int contig_id, std::unordered_map<std::string, std::string>& mateseqs,
		std::unordered_map<std::string, std::string>& matequals, std::vector<consensus_t*>& clip_consensuses) {

	std::string contig_name = contig_map.get_name(contig_id);
    hts_itr_t* iter = sam_itr_querys(dc_file->idx, dc_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    std::vector<bam1_t*> reads;
    std::vector<cluster_t*> clusters;
    while (sam_itr_next(dc_file->file, iter, read) >= 0) {
        std::string mate_seq = get_mate_seq(read, mateseqs);

        if (mate_seq.empty() || mate_seq.find("N") != std::string::npos) continue;

        cluster_t* cluster = new cluster_t(read, 0, config.high_confidence_mapq);
        cluster->id = clusters.size();
        cluster->ra_start = cluster->la_start;
        cluster->ra_end = cluster->la_end;

        clusters.push_back(cluster);
        reads.push_back(bam_dup1(read));
    }
    sam_itr_destroy(iter);
    bam_destroy1(read);

    if (clusters.empty()) return std::vector<insertion_cluster_t*>();

    // union-find data structure
    int n_reads = clusters.size();    
    int max_cluster_size = (stats.get_max_depth(contig_name) * stats.max_is)/stats.read_len;
    cluster_clusters(clusters, reads, stats.max_is, max_cluster_size, true);

    std::vector<insertion_cluster_t*> read_clusters;
    for (int i = 0; i < n_reads; i++) {
        if (clusters[i] != NULL) read_clusters.push_back(new insertion_cluster_t(clusters[i]));
    }

    if (!clip_consensuses.empty()) { // TODO: temporary, enhance logic
        std::sort(clip_consensuses.begin(), clip_consensuses.end(), [](consensus_t* c1, consensus_t* c2) {
            return c1->breakpoint < c2->breakpoint;
        });

        int curr_j = 0;
        std::vector<consensus_t*> clip_consensus_per_cluster(read_clusters.size(), NULL);

        bool left_facing = clip_consensuses[0]->left_clipped;
        if (left_facing) {
            std::sort(read_clusters.begin(), read_clusters.end(), [](insertion_cluster_t* rc1, insertion_cluster_t* rc2) {
                return rc1->end < rc2->end;
            });

            for (consensus_t* clip_consensus : clip_consensuses) {
                while (curr_j < read_clusters.size() && read_clusters[curr_j]->end < clip_consensus->breakpoint) curr_j++;

                for (int j = curr_j; j < read_clusters.size(); j++) {
                    if (read_clusters[j]->end - clip_consensus->breakpoint <= stats.max_is) {
                        if (clip_consensus_per_cluster[j] == NULL || clip_consensus_per_cluster[j]->supp_clipped_reads() < clip_consensus->supp_clipped_reads()) {
                            clip_consensus_per_cluster[j] = clip_consensus;
                        }
                    } else break;
                }
            }
        } else {
            std::sort(read_clusters.begin(), read_clusters.end(), [](insertion_cluster_t* rc1, insertion_cluster_t* rc2) {
                return rc1->start < rc2->start;
            });

            for (consensus_t* clip_consensus : clip_consensuses) {
                while (curr_j < read_clusters.size() && clip_consensus->breakpoint-read_clusters[curr_j]->start > stats.max_is) curr_j++;

                for (int j = curr_j; j < read_clusters.size(); j++) {
                    if (clip_consensus->breakpoint - read_clusters[j]->start >= 0) {
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
        if (read_clusters[i]->cluster->reads.size() <= 1) {
            read_clusters[i]->deallocate_reads();
            delete read_clusters[i];
            read_clusters[i] = nullptr;
        }
    }
    read_clusters.erase(std::remove(read_clusters.begin(), read_clusters.end(), nullptr), read_clusters.end());

    std::sort(read_clusters.begin(), read_clusters.end(), [](const insertion_cluster_t* rc1, const insertion_cluster_t* rc2) {
    	return rc1->cluster->reads.size() > rc2->cluster->reads.size();
    });

    return read_clusters;
}

bool is_semi_mapped(bam1_t* read) {
	return is_proper_pair(read, stats.max_is) && !is_mate_clipped(read);
}

void add_semi_mapped_pairs(std::string clipped_fname, int contig_id, std::vector<insertion_cluster_t*>& r_clusters, std::vector<insertion_cluster_t*>& l_clusters) {

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
	for (insertion_cluster_t* rc : r_clusters) {
		int end = rc->end, start = end - stats.max_is;
		auto sm_reads = r_it.findContained(start, end);
		for (auto i : sm_reads) {
			bam1_t* r = i.value;
			if (!used_sm_reads.count(r)) {
				rc->add_semi_mapped_reads(bam_dup1(r));
				used_sm_reads.insert(r);
			}
		}
	}
	for (insertion_cluster_t* rc : l_clusters) {
		int start, end;
		start = rc->start; end = start + stats.max_is;
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
    std::vector<consensus_t*> rc_consensuses, lc_consensuses;
    if (file_exists(clip_consensus_fname)) {
        std::ifstream clipped_fin(clip_consensus_fname);
        std::string line;
        while (std::getline(clipped_fin, line)) {
            consensus_t* consensus = new consensus_t(line, false);
            if (consensus->left_clipped) {
                consensus->start += consensus->lowq_clip_portion;
                consensus->clip_len -= consensus->lowq_clip_portion;
                consensus->sequence = consensus->sequence.substr(consensus->lowq_clip_portion);
                lc_consensuses.push_back(consensus);
            } else {
                consensus->end -= consensus->lowq_clip_portion;
                consensus->clip_len -= consensus->lowq_clip_portion;
                consensus->sequence = consensus->sequence.substr(0, consensus->sequence.length()-consensus->lowq_clip_portion);
                rc_consensuses.push_back(consensus);
            }
        }
    }

	std::vector<insertion_cluster_t*> r_clusters = cluster_reads(r_dc_file, contig_id, mateseqs, matequals, rc_consensuses);
	std::vector<insertion_cluster_t*> l_clusters = cluster_reads(l_dc_file, contig_id, mateseqs, matequals, lc_consensuses);

	std::string clipped_fname = workdir + "/workspace/clipped/" + std::to_string(contig_id) + ".bam";
	add_semi_mapped_pairs(clipped_fname, contig_id, r_clusters, l_clusters);

    std::priority_queue<cc_v_distance_t> pq;
    std::multimap<int, insertion_cluster_t*> l_clusters_map;
    for (insertion_cluster_t* l_cluster : l_clusters) {
        hts_pos_t cl_start = l_cluster->start;
        l_clusters_map.insert({cl_start, l_cluster});
    }
    for (insertion_cluster_t* r_cluster : r_clusters) {
        hts_pos_t cl_end = r_cluster->end;
        auto begin = l_clusters_map.lower_bound(cl_end - stats.max_is);
        auto end = l_clusters_map.upper_bound(cl_end + stats.max_is);

        for (auto it = begin; it != end; it++) {
            insertion_cluster_t* l_cluster = it->second;
            pq.push(cc_v_distance_t(r_cluster, l_cluster, r_cluster->size()*l_cluster->size()));
        }
    }

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Aligner permissive_aligner(2, 2, 4, 1, false);
    StripedSmithWaterman::Aligner aligner_to_base(1, 4, 6, 1, true);
    StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, true);
    std::vector<bam1_t*> r_reads_to_write, l_reads_to_write;
    while (!pq.empty()) {
        cc_v_distance_t cc_v_distance = pq.top();
        pq.pop();

        insertion_cluster_t* c1 = cc_v_distance.c1;
        insertion_cluster_t* c2 = cc_v_distance.c2;
        if (c1->cluster->used || c2->cluster->used) continue;

        // remap clusters
        std::vector<bam1_t*> to_write;
		bool success = remap_cluster(c1, c2, to_write, contig_id, r_dc_file->header, mateseqs, matequals, 
                                aligner, permissive_aligner, aligner_to_base, harsh_aligner, stats);
        if (!success) continue;

        for (bam1_t* r : to_write) {
            if (bam_is_rev(r)) {
                l_reads_to_write.push_back(bam_dup1(r));
            } else {
                r_reads_to_write.push_back(bam_dup1(r));
            }
        };
        c1->cluster->used = true; c2->cluster->used = true;
        c1->deallocate_reads();
        c2->deallocate_reads();
    }

    std::string r_dc_remapped_fname = workdir + "/workspace/R/" + std::to_string(contig_id) + ".remap.bam";
    write_and_index_file(r_reads_to_write, r_dc_remapped_fname, r_dc_file->header);
    for (bam1_t* r : r_reads_to_write) bam_destroy1(r);

    std::string l_dc_remapped_fname = workdir + "/workspace/L/" + std::to_string(contig_id) + ".remap.bam";
    write_and_index_file(l_reads_to_write, l_dc_remapped_fname, l_dc_file->header);
    for (bam1_t* r : l_reads_to_write) bam_destroy1(r);

    for (insertion_cluster_t* c : r_clusters) c->deallocate_reads();
    for (insertion_cluster_t* c : l_clusters) c->deallocate_reads();

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

	bcf1_t* bcf_entry = bcf_init();

	// out transurveyor insertions
	std::string out_vcf_fname = workdir + "/assembled_ins.vcf.gz";
    bcf_hdr_t* out_vcf_header = generate_vcf_header(contigs, sample_name, config, full_cmd_str);
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + out_vcf_fname + ".");
	}
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	std::sort(insertions.begin(), insertions.end(), [&out_vcf_header](const sv_t* i1, const sv_t* i2) {
		int contig_id1 = bcf_hdr_name2id(out_vcf_header, i1->chr.c_str());
		int contig_id2 = bcf_hdr_name2id(out_vcf_header, i2->chr.c_str());
		// negative because we want descending order
		int disc_score_mul1 = -(i1->disc_pairs_lf*i1->disc_pairs_rf), disc_score_mul2 = -(i2->disc_pairs_lf*i2->disc_pairs_rf);
		int disc_score_sum1 = -(i1->disc_pairs_lf+i1->disc_pairs_rf), disc_score_sum2 = -(i2->disc_pairs_lf+i2->disc_pairs_rf);
		return std::tie(contig_id1, i1->start, i1->end, i1->ins_seq, disc_score_mul1, disc_score_sum1) <
			   std::tie(contig_id2, i2->start, i2->end, i2->ins_seq, disc_score_mul2, disc_score_sum2);
	});

	int a_id = 0, t_id = 0;
	for (sv_t* insertion : insertions) {
        if (insertion->source == "REFERENCE_GUIDED_ASSEMBLY") insertion->id = "T_INS_" + std::to_string(t_id++);
        else insertion->id = "A_INS_" + std::to_string(a_id++);
        sv2bcf(out_vcf_header, bcf_entry, insertion, contigs.get_seq(insertion->chr));
		if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
			throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
		}
	}

    // bcf_close(assembled_out_vcf_file);
    bcf_close(out_vcf_file);
}
