#ifndef STAT_TESTS_H
#define STAT_TESTS_H

#include <cstdint>
#include <cstring>
#include <cmath>
#include <htslib/sam.h>

#include "../libs/ks-test.h"
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
		regions_of_interest.emplace_back(sv->left_anchor_aln->start, sv->left_anchor_aln->end, &(sv->median_left_cluster_cov), &(sv->median_left_cluster_cov_highmq));
		regions_of_interest.emplace_back(sv->right_anchor_aln->start, sv->right_anchor_aln->end, &(sv->median_right_cluster_cov), &(sv->median_right_cluster_cov_highmq));
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
	std::vector<sv_t*> testable_dels(deletions.begin(), deletions.end());
    depth_filter_indel(contig_name, testable_dels, bam_file, config, stats);
}
void depth_filter_dup(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_dups(duplications.begin(), duplications.end());
    depth_filter_indel(contig_name, testable_dups, bam_file, config, stats);
}
void depth_filter_ins(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_inss(insertions.begin(), insertions.end());
	depth_filter_indel(contig_name, testable_inss, bam_file, config, stats);
}
void depth_filter_inv(std::string contig_name, std::vector<inversion_t*>& invertions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_invs(invertions.begin(), invertions.end());
	depth_filter_indel(contig_name, testable_invs, bam_file, config, stats);
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

void calculate_cluster_region_disc(std::string contig_name, std::vector<deletion_t*> deletions, open_samFile_t* bam_file, config_t& config) {

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
				if (read->core.pos <= deletions[i]->left_anchor_aln->end) {
					deletions[i]->l_cluster_region_disc_pairs++;
					if (read->core.qual >= config.high_confidence_mapq) deletions[i]->l_cluster_region_disc_pairs_high_mapq++;
				}
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
				if (read->core.pos <= deletions[i]->right_anchor_aln->end) {
					deletions[i]->r_cluster_region_disc_pairs++;
					if (read->core.qual >= config.high_confidence_mapq) deletions[i]->r_cluster_region_disc_pairs_high_mapq++;
				}
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

void calculate_cluster_region_disc(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

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
					if (read->core.qual >= config.high_confidence_mapq) duplications[i]->l_cluster_region_disc_pairs_high_mapq++;
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
				if (read->core.pos <= duplications[i]->right_anchor_aln->end) {
					duplications[i]->r_cluster_region_disc_pairs++;
					if (read->core.qual >= config.high_confidence_mapq) duplications[i]->r_cluster_region_disc_pairs_high_mapq++;
				}
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
			}

            if (!global_crossing_isize_dist.empty()) {
				double ks_pval = ks_test(global_crossing_isize_dist, local_dists[i]);
				if (std::isfinite(ks_pval)) sv->ks_pval = ks_pval;
				else sv->ks_pval = sv_t::KS_PVAL_NOT_COMPUTED;

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
				int dev_if_homalt = abs((int) (sv->sample_info.alt_bp1.supp_pairs-stats.get_median_disc_pairs_by_del_size(homalt_idx)));
				int dev_if_het = abs((int) (sv->sample_info.alt_bp1.supp_pairs-stats.get_median_disc_pairs_by_del_size(het_idx)/2));
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

void calculate_ptn_ratio(std::string contig_name, std::vector<sv_t*>& svs, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

	if (svs.empty()) return;

	struct conc_pairs_count_t {
		hts_pos_t pos;
		int* conc_pairs_ptr, *conc_pairs_highmq_ptr;

		conc_pairs_count_t(hts_pos_t pos, int* conc_pairs_ptr, int* conc_pairs_highmq_ptr) :
			pos(pos), conc_pairs_ptr(conc_pairs_ptr), conc_pairs_highmq_ptr(conc_pairs_highmq_ptr) {}
	};

	std::vector<conc_pairs_count_t> bkp_with_conc_pairs_count;
	for (sv_t* sv : svs) {
		bkp_with_conc_pairs_count.push_back({sv->start, &(sv->conc_pairs_lbp), &(sv->conc_pairs_lbp_high_mapq)});
		bkp_with_conc_pairs_count.push_back({(sv->start+sv->end)/2, &(sv->conc_pairs_midp), &(sv->conc_pairs_midp_high_mapq)});
		bkp_with_conc_pairs_count.push_back({sv->end, &(sv->conc_pairs_rbp), &(sv->conc_pairs_rbp_high_mapq)});
	}

	std::sort(bkp_with_conc_pairs_count.begin(), bkp_with_conc_pairs_count.end(), [](const conc_pairs_count_t& p1, const conc_pairs_count_t& p2) {
		return p1.pos < p2.pos;
	});

	std::vector<char*> regions;
	for (auto& b : bkp_with_conc_pairs_count) {
		std::stringstream ss;
		ss << contig_name << ":" << std::max(hts_pos_t(1), b.pos-stats.max_is) << "-" << b.pos;
		char* region = strdup(ss.str().c_str());
		regions.push_back(region);
	}

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < bkp_with_conc_pairs_count.size() && bkp_with_conc_pairs_count[curr_pos].pos < read->core.pos) curr_pos++;

		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
		if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0 || read->core.isize > stats.max_is) continue;

		hts_pos_t start = read->core.pos + read->core.l_qseq/2;
		hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;
		for (int i = curr_pos; i < bkp_with_conc_pairs_count.size() && bkp_with_conc_pairs_count[i].pos <= end; i++) {
			if (start <= bkp_with_conc_pairs_count[i].pos && bkp_with_conc_pairs_count[i].pos <= end) {
				(*bkp_with_conc_pairs_count[i].conc_pairs_ptr)++;
				if (read->core.qual >= config.high_confidence_mapq) (*bkp_with_conc_pairs_count[i].conc_pairs_highmq_ptr)++;
			}
		}
	}

	for (char* region : regions) {
		delete[] region;
	}
}

void calculate_ptn_ratio(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t& config, stats_t& stats, bool find_disc_pairs = false,
	std::string nm_file = "") {
	
	if (deletions.empty()) return;
	std::vector<sv_t*> svs(deletions.begin(), deletions.end());
	calculate_ptn_ratio(contig_name, svs, bam_file, config, stats);

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

		std::unordered_map<std::string, int64_t> qname_to_mate_nm; 
		std::ifstream mateseqs_fin(nm_file);
		std::string qname, seq;
		int64_t nm;
		while (mateseqs_fin >> qname >> seq >> nm) {
			qname_to_mate_nm[qname] = nm;
		}

		std::vector<hts_pos_t> dp1_start(deletions.size(), INT32_MAX), dp1_end(deletions.size(), 0);
		std::vector<hts_pos_t> dp2_start(deletions.size(), INT32_MAX), dp2_end(deletions.size(), 0);

		int curr_pos = 0;
		hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
		bam1_t* read = bam_init1();
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			if (is_unmapped(read) || !is_primary(read) || bam_is_rev(read)) continue;

			while (curr_pos < deletions.size() && deletions[curr_pos]->start < read->core.pos) curr_pos++;

			if (read->core.isize > stats.max_is) {
				hts_pos_t pair_start = std::min(read->core.pos + read->core.l_qseq/2, bam_endpos(read));
				hts_pos_t pair_end = std::max(read->core.pos + read->core.isize - read->core.l_qseq/2, read->core.mpos);
				for (int i = curr_pos; i < deletions.size() && deletions[i]->start <= pair_end; i++) {
					if (pair_start <= deletions[i]->start+5 && deletions[i]->end-5 <= pair_end && 
						deletions[i]->start-pair_start + pair_end-deletions[i]->end <= stats.max_is) {
						deletions[i]->sample_info.alt_bp1.supp_pairs++;
						deletions[i]->sample_info.alt_bp2.supp_pairs++;

						if (read->core.pos < dp1_start[i]) dp1_start[i] = read->core.pos;
						if (bam_endpos(read) > dp1_end[i]) dp1_end[i] = bam_endpos(read);
						if (read->core.mpos < dp2_start[i]) dp2_start[i] = read->core.mpos;
						if (get_mate_endpos(read) > dp2_end[i]) dp2_end[i] = get_mate_endpos(read);

						int64_t mq = get_mq(read);
						if (read->core.qual >= config.high_confidence_mapq) {
							deletions[i]->sample_info.alt_bp1.supp_pairs_pos_high_mapq++;
						}
						if (mq >= config.high_confidence_mapq) {
							deletions[i]->sample_info.alt_bp2.supp_pairs_neg_high_mapq++;
						}
						if (read->core.qual > deletions[i]->sample_info.alt_bp1.supp_pairs_max_mq) {
							deletions[i]->sample_info.alt_bp1.supp_pairs_max_mq = read->core.qual;
						}
						if (mq > deletions[i]->sample_info.alt_bp2.supp_pairs_max_mq) {
							deletions[i]->sample_info.alt_bp2.supp_pairs_max_mq = mq;
						}
						deletions[i]->disc_pairs_lf_avg_nm += get_nm(read);
						deletions[i]->disc_pairs_rf_avg_nm += qname_to_mate_nm[std::string(bam_get_qname(read))];
					}
				}
			}

		}
		for (int i = 0; i < deletions.size(); i++) {
			if (deletions[i]->sample_info.alt_bp1.supp_pairs > 0) deletions[i]->disc_pairs_lf_avg_nm /= deletions[i]->sample_info.alt_bp1.supp_pairs;
			if (deletions[i]->sample_info.alt_bp2.supp_pairs > 0) deletions[i]->disc_pairs_rf_avg_nm /= deletions[i]->sample_info.alt_bp2.supp_pairs;
			deletions[i]->disc_pairs_lf_span = std::max(hts_pos_t(0), dp1_end[i] - dp1_start[i]);
			deletions[i]->disc_pairs_rf_span = std::max(hts_pos_t(0), dp2_end[i] - dp2_start[i]);
		}
	}
}
void calculate_ptn_ratio(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, config_t& config, stats_t& stats, bool find_disc_pairs = false,
	std::string nm_file = "") {
	
	if (duplications.empty()) return;
	std::vector<sv_t*> svs(duplications.begin(), duplications.end());
	calculate_ptn_ratio(contig_name, svs, bam_file, config, stats);

	if (find_disc_pairs) {
		std::vector<char*> regions;
		for (duplication_t* dup : duplications) {
			std::stringstream ss;
			ss << contig_name << ":" << std::max(hts_pos_t(1), dup->start) << "-" << dup->start+stats.max_is;
			char* region = new char[ss.str().length()+1];
			strcpy(region, ss.str().c_str());
			regions.push_back(region);
		}

		std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
			return d1->start < d2->start;
		});

		std::unordered_map<std::string, int64_t> qname_to_mate_nm; 
		std::ifstream mateseqs_fin(nm_file);
		std::string qname, seq;
		int64_t nm;
		while (mateseqs_fin >> qname >> seq >> nm) {
			qname_to_mate_nm[qname] = nm;
		}

		std::vector<hts_pos_t> dp1_start(duplications.size(), INT32_MAX), dp1_end(duplications.size(), 0);
		std::vector<hts_pos_t> dp2_start(duplications.size(), INT32_MAX), dp2_end(duplications.size(), 0);

		int curr_pos = 0;
		hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
		bam1_t* read = bam_init1();
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			if (!bam_is_rev(read) || is_unmapped(read) || !is_primary(read) || !is_outward(read)) continue;

			while (curr_pos < duplications.size() && duplications[curr_pos]->start+stats.max_is < read->core.pos) curr_pos++;

			for (int i = curr_pos; i < duplications.size() && bam_endpos(read) > duplications[i]->start; i++) {
				duplication_t* dup = duplications[i];
				hts_pos_t pair_start = read->core.pos + read->core.l_qseq/2, pair_end = read->core.mpos + read->core.l_qseq/2;
				if (dup->start < pair_start && pair_start < dup->start+stats.max_is && dup->end-stats.max_is < pair_end && pair_end < dup->end) {
					dup->sample_info.alt_bp1.supp_pairs++;
					dup->sample_info.alt_bp2.supp_pairs++;

					if (read->core.pos < dp1_start[i]) dp1_start[i] = read->core.pos;
					if (bam_endpos(read) > dp1_end[i]) dp1_end[i] = bam_endpos(read);
					if (read->core.mpos < dp2_start[i]) dp2_start[i] = read->core.mpos;
					if (get_mate_endpos(read) > dp2_end[i]) dp2_end[i] = get_mate_endpos(read);

					int64_t mq = get_mq(read);
					if (read->core.qual >= config.high_confidence_mapq) {
						dup->sample_info.alt_bp1.supp_pairs_neg_high_mapq++;
					}
					if (mq >= config.high_confidence_mapq) {
						dup->sample_info.alt_bp1.supp_pairs_pos_high_mapq++;
					}
					if (read->core.qual > duplications[i]->sample_info.alt_bp1.supp_pairs_max_mq) {
						duplications[i]->sample_info.alt_bp1.supp_pairs_max_mq = read->core.qual;
					}
					if (mq > duplications[i]->sample_info.alt_bp2.supp_pairs_max_mq) {
						duplications[i]->sample_info.alt_bp2.supp_pairs_max_mq = mq;
					}
					dup->disc_pairs_lf_avg_nm += get_nm(read);
					dup->disc_pairs_rf_avg_nm += qname_to_mate_nm[std::string(bam_get_qname(read))];
				}
			}
		}
		for (int i = 0; i < duplications.size(); i++) {
			if (duplications[i]->sample_info.alt_bp1.supp_pairs > 0) duplications[i]->disc_pairs_lf_avg_nm /= duplications[i]->sample_info.alt_bp1.supp_pairs;
			if (duplications[i]->sample_info.alt_bp2.supp_pairs > 0) duplications[i]->disc_pairs_rf_avg_nm /= duplications[i]->sample_info.alt_bp2.supp_pairs;
			duplications[i]->disc_pairs_lf_span = std::max(hts_pos_t(0), dp1_end[i] - dp1_start[i]);
			duplications[i]->disc_pairs_rf_span = std::max(hts_pos_t(0), dp2_end[i] - dp2_start[i]);
		}
	}
}
void calculate_ptn_ratio(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, config_t& config, stats_t& stats, bool find_disc_pairs = false) {
	if (insertions.empty()) return;
	std::vector<sv_t*> svs(insertions.begin(), insertions.end());
	calculate_ptn_ratio(contig_name, svs, bam_file, config, stats);
}
void calculate_ptn_ratio(std::string contig_name, std::vector<inversion_t*>& inversions, open_samFile_t* bam_file, config_t& config, stats_t& stats, bool find_disc_pairs = false) {
	if (inversions.empty()) return;

	std::vector<std::pair<hts_pos_t, inversion_t*> > bkp_with_conc_pairs_count;
	for (inversion_t* inv : inversions) {
		bkp_with_conc_pairs_count.push_back({inv->start, inv});
		bkp_with_conc_pairs_count.push_back({inv->end, inv});
	}

	std::sort(bkp_with_conc_pairs_count.begin(), bkp_with_conc_pairs_count.end(), [](const std::pair<hts_pos_t, inversion_t*>& p1, const std::pair<hts_pos_t, inversion_t*>& p2) {
		return p1.first < p2.first;
	});

	std::vector<char*> regions;
	for (auto& b : bkp_with_conc_pairs_count) {
		std::stringstream ss;
		ss << contig_name << ":" << std::max(hts_pos_t(1), b.first-stats.max_is) << "-" << b.first;
		char* region = strdup(ss.str().c_str());
		regions.push_back(region);
	}

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < bkp_with_conc_pairs_count.size() && bkp_with_conc_pairs_count[curr_pos].first < read->core.pos) curr_pos++;

		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
		if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0 || read->core.isize > stats.max_is) continue;

		hts_pos_t start = read->core.pos + read->core.l_qseq/2;
		hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;
		for (int i = curr_pos; i < bkp_with_conc_pairs_count.size() && bkp_with_conc_pairs_count[i].first <= end; i++) {
			if (start <= bkp_with_conc_pairs_count[i].first) { // pair is concordant and crosses the breakpoint
				inversion_t* inv = bkp_with_conc_pairs_count[i].second;
				if (start <= inv->start && end <= inv->end) {
					inv->conc_pairs_lbp++;
					if (read->core.qual >= config.high_confidence_mapq) inv->conc_pairs_lbp_high_mapq++;
				} else if (start >= inv->start && end >= inv->end) {
					inv->conc_pairs_rbp++;
					if (read->core.qual >= config.high_confidence_mapq) inv->conc_pairs_rbp_high_mapq++;
				}
			}
		}
	}

	for (char* region : regions) {
		delete[] region;
	}

	if (find_disc_pairs) {
		std::vector<char*> regions;
		for (inversion_t* inv : inversions) {
			std::stringstream ss;
			ss << contig_name << ":" << std::max(hts_pos_t(1), inv->start-stats.max_is) << "-" << inv->start+stats.max_is;
			char* region = new char[ss.str().length()+1];
			strcpy(region, ss.str().c_str());
			regions.push_back(region);
		}

		std::sort(inversions.begin(), inversions.end(), [](const inversion_t* i1, const inversion_t* i2) {
			return i1->start < i2->start;
		});
		
		int curr_pos = 0;
		hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
		bam1_t* read = bam_init1();
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			if (is_unmapped(read) || !is_primary(read) || !is_samechr(read) || !is_samestr(read)) continue;

			while (curr_pos < inversions.size() && inversions[curr_pos]->start+stats.max_is < bam_endpos(read)) curr_pos++;

			for (int i = curr_pos; i < inversions.size() && read->core.pos >= inversions[i]->start-stats.max_is; i++) {
				inversion_t* inv = inversions[i];
				
				if (!bam_is_rev(read) && read->core.pos+stats.read_len/2 > inv->start) continue;
				if (bam_is_rev(read) && bam_endpos(read)-stats.read_len/2 < inv->start) break;

				hts_pos_t mate_startpos = read->core.mpos, mate_endpos = get_mate_endpos(read);
				if (!bam_is_mrev(read) && mate_startpos >= inv->end-stats.max_is && mate_startpos+stats.read_len/2 <= inv->end) {
					inv->sample_info.alt_bp1.supp_pairs++;
					
					int64_t mq = get_mq(read);
					if (read->core.qual >= config.high_confidence_mapq) {
						inv->sample_info.alt_bp1.supp_pairs_pos_high_mapq++;
					}
					if (mq >= config.high_confidence_mapq) {
						inv->sample_info.alt_bp1.supp_pairs_neg_high_mapq++;
					}
					if (inv->sample_info.alt_bp1.supp_pairs_max_mq < read->core.qual) {
						inv->sample_info.alt_bp1.supp_pairs_max_mq = read->core.qual;
					}
				}
				if (bam_is_mrev(read) && mate_endpos-stats.read_len/2 >= inv->end && mate_endpos <= inv->end+stats.max_is) {
					inv->sample_info.alt_bp2.supp_pairs++;

					int64_t mq = get_mq(read);
					if (read->core.qual >= config.high_confidence_mapq) {
						inv->sample_info.alt_bp2.supp_pairs_pos_high_mapq++;
					}
					if (mq >= config.high_confidence_mapq) {
						inv->sample_info.alt_bp2.supp_pairs_neg_high_mapq++;
					}
					if (inv->sample_info.alt_bp2.supp_pairs_max_mq < mq) {
						inv->sample_info.alt_bp2.supp_pairs_max_mq = mq;
					}
				}
			}
		}
	}
}


#endif //STAT_TESTS_H
