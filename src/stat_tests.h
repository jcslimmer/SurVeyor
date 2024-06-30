#ifndef STAT_TESTS_H
#define STAT_TESTS_H

#include <cstring>
#include <set>
#include <cmath>
#include <unordered_set>
#include <random>
#include <htslib/sam.h>

#include "../libs/ks-test.h"
#include "../libs/IntervalTree.h"
#include "htslib/hts.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"

struct read_w_cached_info_t {
    bam1_t* read;
    hts_pos_t start, end, isize;
    int64_t as;
    int aln_len;
    bool left_clipped, right_clipped, is_rev, is_mrev;
    int references = 0;

    read_w_cached_info_t(bam1_t* read) : read(bam_dup1(read)), as(get_AS_tag(read)), start(read->core.pos), end(bam_endpos(read)),
                                         aln_len(get_aligned_portion_len(read)), left_clipped(is_left_clipped(read, 0)),
                                         right_clipped(is_right_clipped(read, 0)), is_rev(bam_is_rev(read)), is_mrev(bam_is_mrev(read)),
                                         isize(read->core.isize) {}

    ~read_w_cached_info_t() { bam_destroy1(read); }
};

const int FLANKING_SIZE = 5000, INDEL_TESTED_REGION_SIZE = 10000;

struct region_depth_t {
	hts_pos_t start, end;
	std::vector<uint32_t> depths, highmq_depths;

	region_depth_t(hts_pos_t start, hts_pos_t end) : start(start), end(end), depths(std::vector<uint32_t>(end-start+1)),
			highmq_depths(std::vector<uint32_t>(end-start+1)) {}
};
struct region_w_median_t {
	hts_pos_t start, end;
	int* median_target;
	int* median_target_highmq_only;

	region_w_median_t(hts_pos_t start, hts_pos_t end, int* median_target, int* median_target_highmq_only) :
		start(start), end(end), median_target(median_target), median_target_highmq_only(median_target_highmq_only) {}
};
void depth_filter_indel(std::string contig_name, std::vector<sv_t*>& svs, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

	if (svs.empty()) return;
	std::sort(svs.begin(), svs.end(), [](const sv_t* s1, const sv_t* s2) {
		return s1->start < s2->start;
	});

    std::vector<char*> regions;
    for (sv_t* sv : svs) {
        std::stringstream ss;
        ss << contig_name << ":" << std::max(hts_pos_t(1), sv->start - FLANKING_SIZE) << "-" << std::min(sv->start + INDEL_TESTED_REGION_SIZE, sv->end);
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);

        ss.str(std::string());
        ss << contig_name << ":" << std::max(sv->end - INDEL_TESTED_REGION_SIZE, sv->start) << "-" << sv->end + FLANKING_SIZE;
        region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);
    }

    std::vector<region_w_median_t> regions_of_interest;
    for (sv_t* sv : svs) {
    	regions_of_interest.emplace_back(std::max(hts_pos_t(0), sv->start - FLANKING_SIZE), sv->start, &(sv->median_left_flanking_cov), &(sv->median_left_flanking_cov_highmq));
    	regions_of_interest.emplace_back(sv->start, std::min(sv->start + INDEL_TESTED_REGION_SIZE, sv->end), &(sv->median_indel_left_cov), &(sv->median_indel_left_cov_highmq));
    	regions_of_interest.emplace_back(std::max(sv->end - INDEL_TESTED_REGION_SIZE, sv->start), sv->end, &(sv->median_indel_right_cov), &(sv->median_indel_right_cov_highmq));
    	regions_of_interest.emplace_back(sv->end, sv->end + FLANKING_SIZE, &(sv->median_right_flanking_cov), &(sv->median_right_flanking_cov_highmq));
    	if (sv->svtype() == "DEL") {
			hts_pos_t l_cluster_start = std::min(sv->left_anchor_aln->start, sv->start-stats.read_len);
			l_cluster_start = std::max(hts_pos_t(0), l_cluster_start);
			regions_of_interest.emplace_back(l_cluster_start, sv->start, &(sv->median_left_cluster_cov), &(sv->median_left_cluster_cov_highmq));
			hts_pos_t r_cluster_end = std::max(sv->right_anchor_aln->end, sv->end+stats.read_len);
	    	regions_of_interest.emplace_back(sv->end, r_cluster_end, &(sv->median_right_cluster_cov), &(sv->median_right_cluster_cov_highmq));
    	}
    }
    std::sort(regions_of_interest.begin(), regions_of_interest.end(), [](region_w_median_t& r1, region_w_median_t& r2) {return r1.start < r2.start;});

    std::vector<region_depth_t> merged_regions_of_interest;
    auto curr_region = regions_of_interest[0];
    for (const auto& region : regions_of_interest) {
    	if (region.start - curr_region.end <= stats.read_len) {
    		curr_region.end = std::max(region.end, curr_region.end);
    	} else {
    		merged_regions_of_interest.emplace_back(curr_region.start, curr_region.end);
    		curr_region = region;
        }
    }
    merged_regions_of_interest.emplace_back(curr_region.start, curr_region.end);

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	int curr_pos = 0;
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        hts_pos_t rs = read->core.pos, re = bam_endpos(read);
        while (curr_pos < merged_regions_of_interest.size() && rs > merged_regions_of_interest[curr_pos].end) curr_pos++;
        if (curr_pos == merged_regions_of_interest.size()) break;

        auto& curr_region = merged_regions_of_interest[curr_pos];

        if (overlap(curr_region.start, curr_region.end, rs, re)) {
			int ov_start = std::max(0, int(rs-curr_region.start));
			int ov_end = std::min(int(curr_region.depths.size()), int(re-curr_region.start));
			for (int i = ov_start; i < ov_end; i++) {
				curr_region.depths[i]++;
			}
		}

        if (read->core.qual < config.high_confidence_mapq) continue;

        if (overlap(curr_region.start, curr_region.end, rs, re)) {
        	int ov_start = std::max(0, int(rs-curr_region.start));
        	int ov_end = std::min(int(curr_region.depths.size()), int(re-curr_region.start));
        	for (int i = ov_start; i < ov_end; i++) {
        		curr_region.highmq_depths[i]++;
        	}
        }
	}

	curr_pos = 0;
	for (auto& region : regions_of_interest) {
		while (curr_pos < merged_regions_of_interest.size() && region.start > merged_regions_of_interest[curr_pos].end) curr_pos++;

		auto& merged_region = merged_regions_of_interest[curr_pos];
		std::vector<uint32_t> depths(merged_region.depths.begin() + (region.start - merged_region.start),
		                             merged_region.depths.begin() + (region.end - merged_region.start));
		if (depths.empty()) continue;
		std::sort(depths.begin(), depths.end());
		*region.median_target = depths[depths.size()/2];

		std::vector<uint32_t> highmq_depths(merged_region.highmq_depths.begin() + (region.start - merged_region.start),
		                                    merged_region.highmq_depths.begin() + (region.end - merged_region.start));
		std::sort(highmq_depths.begin(), highmq_depths.end());
		*region.median_target_highmq_only = highmq_depths[highmq_depths.size()/2];
	}

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void depth_filter_del(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	if (deletions.empty()) return;
	std::vector<sv_t*> testable_dels(deletions.begin(), deletions.end());
    depth_filter_indel(contig_name, testable_dels, bam_file, config, stats);
}
void depth_filter_dup(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	if (duplications.empty()) return;
	std::vector<sv_t*> testable_dups(duplications.begin(), duplications.end());
    depth_filter_indel(contig_name, testable_dups, bam_file, config, stats);
}
void depth_filter_ins(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	if (insertions.empty()) return;
	std::vector<sv_t*> testable_ins(insertions.begin(), insertions.end());
	depth_filter_indel(contig_name, testable_ins, bam_file, config, stats);
}

int find_smallest_range_start(std::vector<int>& v, int range_size, int& min_cum) {
	int cum = 0;
	for (int j = 0; j < range_size && j < v.size(); j++) {
		cum += v[j];
	}
	int best_start = 0, best_cum = cum;
	for (int j = range_size; j < v.size(); j++) {
		cum = cum - v[j-range_size] + v[j];
		if (cum < best_cum) {
			best_cum = cum;
			best_start = j - range_size + 1;
		}
	}
	min_cum = best_cum;
	return best_start;
}

void calculate_cluster_region_disc(std::string contig_name, std::vector<deletion_t*> deletions, open_samFile_t* bam_file) {

	if (deletions.empty()) return;

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->start < d2->start;
	});

	std::vector<char*> l_cluster_regions, r_cluster_regions;
	for (deletion_t* deletion : deletions) {
		std::stringstream ss;
		ss << contig_name << ":" << deletion->left_anchor_aln->start << "-" << deletion->left_anchor_aln->end;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		l_cluster_regions.push_back(region);

		ss.str("");
		ss << contig_name << ":" << deletion->right_anchor_aln->start << "-" << deletion->right_anchor_aln->end;
		region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		r_cluster_regions.push_back(region);
	}

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->left_anchor_aln->start < d2->left_anchor_aln->start;
	});

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, l_cluster_regions.data(), l_cluster_regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		while (curr_pos < deletions.size() && deletions[curr_pos]->left_anchor_aln->end < read->core.pos) curr_pos++;

		// if the pair is discordant and it overlaps left_anchor_aln, increase l_cluster_region_disc_pair
		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_outward(read)) {
			for (int i = curr_pos; i < deletions.size() && deletions[i]->left_anchor_aln->start <= bam_endpos(read); i++) {
				if (read->core.pos <= deletions[i]->left_anchor_aln->end) deletions[i]->l_cluster_region_disc_pairs++;
				// TODO: try to remap the mate instead?
			}
		}
	}

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->end < d2->end;
	});

	curr_pos = 0;
	iter = sam_itr_regarray(bam_file->idx, bam_file->header, r_cluster_regions.data(), r_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < deletions.size() && deletions[curr_pos]->right_anchor_aln->end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_outward(read)) {
			for (int i = curr_pos; i < deletions.size() && deletions[i]->right_anchor_aln->start <= bam_endpos(read); i++) {
				if (read->core.pos <= deletions[i]->right_anchor_aln->end) deletions[i]->r_cluster_region_disc_pairs++;
			}
		}
	}

	for (char* region : l_cluster_regions) {
		delete[] region;
	}
	for (char* region : r_cluster_regions) {
		delete[] region;
	}
}

void calculate_cluster_region_disc(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, stats_t& stats) {

	if (duplications.empty()) return;

	std::vector<char*> lc_cluster_regions, rc_cluster_regions;
	for (duplication_t* dup : duplications) {
		std::stringstream ss;
		ss << contig_name << ":" << dup->left_anchor_aln->start << "-" << dup->left_anchor_aln->end;
		rc_cluster_regions.push_back(strdup(ss.str().c_str()));

		ss.str("");
		ss << contig_name << ":" << dup->right_anchor_aln->start << "-" << dup->right_anchor_aln->end;
		lc_cluster_regions.push_back(strdup(ss.str().c_str()));
	}

	std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->left_anchor_aln->start < d2->left_anchor_aln->start;
	});

	bam1_t* read = bam_init1();
	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, rc_cluster_regions.data(), rc_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		while (curr_pos < duplications.size() && duplications[curr_pos]->left_anchor_aln->end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_long(read, stats.max_is)) {
			for (int i = curr_pos; i < duplications.size() && duplications[i]->left_anchor_aln->start <= bam_endpos(read); i++) {
				if (read->core.pos <= duplications[i]->left_anchor_aln->end) {
					duplications[i]->l_cluster_region_disc_pairs++;
				}
			}
		}
	}
	hts_itr_destroy(iter);
	
	for (char* region : rc_cluster_regions) {
		free(region);
	}

	std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->start < d2->start;
	});

	curr_pos = 0;
	iter = sam_itr_regarray(bam_file->idx, bam_file->header, lc_cluster_regions.data(), lc_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		while (curr_pos < duplications.size() && duplications[curr_pos]->right_anchor_aln->end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_long(read, stats.max_is)) {
			for (int i = curr_pos; i < duplications.size() && duplications[i]->right_anchor_aln->start <= bam_endpos(read); i++) {
				if (read->core.pos <= duplications[i]->right_anchor_aln->end) duplications[i]->r_cluster_region_disc_pairs++;
			}
		}
	}
	hts_itr_destroy(iter);

	for (char* region : lc_cluster_regions) {
		free(region);
	}
}

void calculate_confidence_interval_size(std::string contig_name, std::vector<double>& global_crossing_isize_dist,
										std::vector<sv_t*>& svs, open_samFile_t* bam_file, 
										config_t& config, stats_t& stats, int min_sv_size,
										bool disallow_changes = false) {

	if (svs.empty()) return;

	std::sort(svs.begin(), svs.end(), [](const sv_t* d1, const sv_t* d2) {
		return (d1->start+d1->end)/2 < (d2->start+d2->end)/2;
	});

    std::vector<hts_pos_t> midpoints, sizes;
    std::vector<uint64_t> sums(svs.size()), sums_highmq(svs.size());
	std::vector<uint64_t> sq_sums(svs.size()), sq_sums_highmq(svs.size());
    std::vector<uint32_t> ns(svs.size()), ns_highmq(svs.size());
    std::vector<char*> regions;
    std::vector<std::pair<hts_pos_t, hts_pos_t> > regions_coos;
    for (sv_t* sv : svs) {
        hts_pos_t midpoint = (sv->start+sv->end)/2, size = sv->end-sv->start;
        midpoints.push_back(midpoint);
        sizes.push_back(size);

        regions_coos.push_back({std::max(hts_pos_t(1), sv->start-stats.max_is), sv->start});
        regions_coos.push_back({std::max(hts_pos_t(1), midpoint-stats.max_is), midpoint});
    }
	std::sort(regions_coos.begin(), regions_coos.end());
    for (auto coos : regions_coos) {
    	std::stringstream ss;
    	ss << contig_name << ":" << coos.first << "-" << coos.second;
    	char* region = new char[1000];
    	strcpy(region, ss.str().c_str());
		regions.push_back(region);
    }

	int curr_pos = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    std::vector<std::vector<double> > local_dists(svs.size()), local_dists_highmq(svs.size());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {

    	if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
        if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0) continue;

        while (curr_pos < svs.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

        hts_pos_t start = read->core.pos + read->core.l_qseq/2;
        hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;
        for (int i = curr_pos; i < midpoints.size() && midpoints[i] < read->core.pos+read->core.isize; i++) {
            if (start <= midpoints[i] && midpoints[i] <= end && read->core.isize <= stats.max_is+sizes[i]) {
                sums[i] += read->core.isize;
                sq_sums[i] += read->core.isize*read->core.isize;
                ns[i]++;
                local_dists[i].push_back(read->core.isize);
				if (read->core.qual >= config.high_confidence_mapq) {
					sums_highmq[i] += read->core.isize;
                	sq_sums_highmq[i] += read->core.isize*read->core.isize;
					ns_highmq[i]++;
					local_dists_highmq[i].push_back(read->core.isize);
				}
            }
        }
    }

    for (int i = 0; i < svs.size(); i++) {
		sv_t* sv = svs[i];
        uint32_t n = ns[i], n_highmq = ns_highmq[i];
        uint64_t sum = sums[i], sum_highmq = sums_highmq[i];
		uint64_t sq_sum = sq_sums[i], sq_sum_highmq = sq_sums_highmq[i];
        if (n >= 4) {
            int avg_is = sum/n;
            int var_is = (sq_sum - sum*sum/n)/(n-1);
            int confidence_ival = 2.576 * sqrt(var_is/n);
			sv->min_conf_size = abs(avg_is - stats.pop_avg_crossing_is) - confidence_ival;
            sv->max_conf_size = abs(avg_is - stats.pop_avg_crossing_is) + confidence_ival;

			if (n_highmq >= 4) {
				int avg_is_highmq = sum_highmq/n_highmq;
				int var_is_highmq = (sq_sum_highmq - sum_highmq*sum_highmq/n_highmq)/(n_highmq-1);
				int confidence_ival_highmq = 2.576 * sqrt(var_is_highmq/n_highmq);
				sv->min_conf_size_highmq = abs(avg_is_highmq - stats.pop_avg_crossing_is) - confidence_ival_highmq;
				sv->max_conf_size_highmq = abs(avg_is_highmq - stats.pop_avg_crossing_is) + confidence_ival_highmq;
			}

            if (!global_crossing_isize_dist.empty()) {
				sv->ks_pval = ks_test(global_crossing_isize_dist, local_dists[i]);
				if (n_highmq >= 4) {
					sv->ks_pval_highmq = ks_test(global_crossing_isize_dist, local_dists_highmq[i]);
				}
            
				if (sv->svtype() != "DEL" || !sv->imprecise || disallow_changes) continue;

				int est_size = avg_is - stats.pop_avg_crossing_is;

				// compute depth base by base in the imprecise deleted regions
				// TODO: there are some faster data structures out there for this - e.g., Fenwick tree
				hts_pos_t range_start = std::max(hts_pos_t(1), sv->start-stats.read_len), range_end = sv->end + stats.read_len;
				if (est_size > range_end-range_start || est_size < min_sv_size) continue; // estimated size is grossly off - ignore

				std::vector<int> depth_by_base(range_end-range_start+1);
				std::stringstream ss;
				ss << contig_name << ":" << range_start << "-" << range_end;
				iter = sam_itr_querys(bam_file->idx, bam_file->header, ss.str().c_str());
				while (sam_itr_next(bam_file->file, iter, read) >= 0) {
					if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read) || read->core.qual < 0) continue;
					if (!is_samechr(read) || is_samestr(read) || read->core.isize < -stats.max_is || read->core.isize > stats.max_is) continue;

					int start = std::max(hts_pos_t(0), read->core.pos-range_start);
					int end = std::min(bam_endpos(read)-range_start, range_end-range_start);
					for (int i = start; i <= end; i++) depth_by_base[i]++;
				}

				// find window of est_size bp with minimum depth
				int min_cum_depth;
				int best_start = find_smallest_range_start(depth_by_base, est_size, min_cum_depth);
				int del_cov = min_cum_depth/est_size;

				std::vector<double> homalt_local_dist, het_local_dist;
				for (int i = 0; i < global_crossing_isize_dist.size(); i++) {
					double d = global_crossing_isize_dist[i];
					homalt_local_dist.push_back(d+est_size);
					if (i%2) het_local_dist.push_back(d+2*est_size);
					else het_local_dist.push_back(d);
				}

				int homalt_evidence = 0, het_evidence = 0;

				// whether the distribution supports a het or a homalt deletion
				double homalt_pval = ks_test(local_dists[i], homalt_local_dist);
				double het_pval = ks_test(local_dists[i], het_local_dist);
				if (homalt_pval < het_pval) { // deletion is het
					het_evidence++;
				} else { // deletion is homalt
					homalt_evidence++;
				}

				// whether the discordant pairs are closer to what we expect with a het or a hom alt deletion
				// TODO: this is better than random but not very accurate. If we could have a test on standard deviations (het will have larger
				// than the global, hom alt should be the same) as like Bartlett as the third vote it might be better
				int homalt_idx = std::max(0, stats.max_is - est_size);
				int het_idx = std::max(0, stats.max_is - 2*est_size);
				int dev_if_homalt = abs((int) (sv->disc_pairs_lf-stats.get_median_disc_pairs_by_del_size(homalt_idx)));
				int dev_if_het = abs((int) (sv->disc_pairs_lf-stats.get_median_disc_pairs_by_del_size(het_idx)/2));
				if (dev_if_het <= dev_if_homalt) {
					het_evidence++;
				} else { // deletion is homalt
					homalt_evidence++;
				}

				// depth supports hom alt or het
				if (sv->median_left_flanking_cov*0.25 < del_cov || sv->median_right_flanking_cov*0.25 < del_cov) {
					het_evidence++;
				} else {
					homalt_evidence++;
				}

				std::string genotype = het_evidence > homalt_evidence ? "0/1" : "1/1";

				// if deletion is estimated het, double the size (unless it is bigger than the original range - then assume GT is wrong)
				if (genotype == "0/1" && est_size*2 <= range_end-range_start) {
					est_size *= 2;
					best_start = find_smallest_range_start(depth_by_base, est_size, min_cum_depth);
				}

				int new_start = range_start + best_start, new_end = new_start + est_size;
				if (new_start >= sv->start-stats.read_len/2 && new_end <= sv->end+stats.read_len/2) {
					((deletion_t*) sv)->original_range = std::to_string(sv->start) + "-" + std::to_string(sv->end);
					sv->start = new_start; sv->end = new_end;
					sv->start = new_start; sv->end = new_end;
				}
            }
        }
    }

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void calculate_ptn_ratio(std::string contig_name, std::vector<sv_t*>& svs, open_samFile_t* bam_file, stats_t& stats) {

	if (svs.empty()) return;

	std::sort(svs.begin(), svs.end(), [](const sv_t* sv1, const sv_t* sv2) {
		return (sv1->start+sv1->end)/2 < (sv2->start+sv2->end)/2;
	});

	std::vector<hts_pos_t> midpoints;
	std::vector<char*> mid_regions;
	for (sv_t* sv : svs) {
		hts_pos_t midpoint = (sv->start+sv->end)/2;
		midpoints.push_back(midpoint);

		std::stringstream ss;
		ss << contig_name << ":" << std::max(hts_pos_t(1), midpoint-stats.max_is) << "-" << midpoint;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		mid_regions.push_back(region);
	}

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, mid_regions.data(), mid_regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < svs.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
		if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0 || read->core.isize > stats.max_is) continue;

		hts_pos_t start = read->core.pos + read->core.l_qseq/2;
		hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;
		for (int i = curr_pos; i < svs.size() && midpoints[i] <= end; i++) {
			if (start <= midpoints[i] && midpoints[i] <= end) svs[i]->conc_pairs++;
		}
	}

	for (char* region : mid_regions) {
		delete[] region;
	}
}

void calculate_ptn_ratio(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t& config, stats_t& stats, bool find_disc_pairs = false) {
	if (deletions.empty()) return;
	std::vector<sv_t*> svs(deletions.begin(), deletions.end());
	calculate_ptn_ratio(contig_name, svs, bam_file, stats);

	if (find_disc_pairs) {
		std::vector<char*> regions;
		for (deletion_t* del : deletions) {
			std::stringstream ss;
			ss << contig_name << ":" << std::max(hts_pos_t(1), del->start-stats.max_is) << "-" << del->start;
			char* region = new char[ss.str().length()+1];
			strcpy(region, ss.str().c_str());
			regions.push_back(region);
		}

		std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
			return d1->start < d2->start;
		});

		int curr_pos = 0;
		hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
		bam1_t* read = bam_init1();
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			while (curr_pos < deletions.size() && deletions[curr_pos]->start < read->core.pos) curr_pos++;

			if (read->core.isize > stats.max_is) {
				hts_pos_t pair_start = read->core.pos + read->core.l_qseq/2, pair_end = read->core.pos + read->core.isize - read->core.l_qseq/2;
				for (int i = curr_pos; i < deletions.size() && deletions[i]->start <= pair_end; i++) {
					if (pair_start <= deletions[i]->start && deletions[i]->end <= pair_end) {
						deletions[i]->disc_pairs_lf++;
						deletions[i]->disc_pairs_rf++;
						if (read->core.qual >= config.high_confidence_mapq) {
							deletions[i]->disc_pairs_lf_high_mapq++;
							deletions[i]->disc_pairs_rf_high_mapq++;
						}
						if (read->core.qual > deletions[i]->disc_pairs_lf_maxmapq) {
							deletions[i]->disc_pairs_lf_maxmapq = read->core.qual;
							deletions[i]->disc_pairs_rf_maxmapq = read->core.qual;
						}
						deletions[i]->disc_pairs_lf_avg_nm += get_nm(read);
						// TODO: is it worth it to spend the extra time to calculate the rf avg NM?
					}
				}
			}

		}
		for (deletion_t* del : deletions) {
			if (del->disc_pairs_lf > 0) del->disc_pairs_lf_avg_nm /= del->disc_pairs_lf;
		}
	}
}
void calculate_ptn_ratio(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, stats_t& stats) {
	if (insertions.empty()) return;
	std::vector<sv_t*> svs(insertions.begin(), insertions.end());
	calculate_ptn_ratio(contig_name, svs, bam_file, stats);
}


#endif //STAT_TESTS_H
