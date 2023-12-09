#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <sstream>
#include <unordered_map>

#include "utils.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "../libs/ssw.h"
#include "../libs/IntervalTree.h"
#include "../src/sam_utils.h"
#include "stat_tests.h"
#include "../src/sw_utils.h"
#include "../src/sam_utils.h"

config_t config;
stats_t stats;
std::string workdir, complete_bam_fname, reference_fname;
std::mutex mtx;

chr_seqs_map_t chr_seqs;
contig_map_t contig_map;

std::vector<double> global_crossing_isize_dist;
std::vector<uint32_t> median_crossing_count_geqi_by_isize;

std::mutex bam_pool_mtx;
std::queue<open_samFile_t*> bam_pool;
open_samFile_t* get_bam_reader(std::string bam_fname) {
    open_samFile_t* o;
    bam_pool_mtx.lock();
    if (!bam_pool.empty()) {
        o = bam_pool.front();
        bam_pool.pop();
    } else {
        o = open_samFile(bam_fname.c_str());
        hts_set_fai_filename(o->file, fai_path(reference_fname.c_str()));
    }
    bam_pool_mtx.unlock();
    return o;
}
void release_bam_reader(open_samFile_t* reader) {
    bam_pool_mtx.lock();
    bam_pool.push(reader);
    bam_pool_mtx.unlock();
}

std::unordered_map<std::string, std::vector<deletion_t*> > deletions_by_chr;
std::unordered_map<std::string, std::vector<duplication_t*> > duplications_by_chr;
std::unordered_map<std::string, std::vector<insertion_t*> > insertions_by_chr;
std::unordered_map<std::string, std::vector<sv2_deletion_t*> > sv2_deletions_by_chr;
std::unordered_map<std::string, std::vector<sv2_duplication_t*> > sv2_duplications_by_chr;

void size_and_depth_filtering_del(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    mtx.lock();
    std::vector<sv2_deletion_t*>& deletions = sv2_deletions_by_chr[contig_name];
    mtx.unlock();
    // TODO: we are only computing the ks-test p-value for large deletions and DP deletions for efficienty reasons. Investigate why it takes so long calculating it for all deletions.
    std::vector<sv2_deletion_t*> large_dp_deletions, small_dp_deletions, sr_deletions;
    for (sv2_deletion_t* del : deletions) {
        if (-del->sv->svlen() >= stats.max_is && del->sv->source == "DP") {
            large_dp_deletions.push_back(del);
        } else if (del->sv->source == "DP") {
            small_dp_deletions.push_back(del);
        } else {
            sr_deletions.push_back(del);
        }
    }
    std::vector<double> v1;
    std::vector<uint32_t> v2;
    calculate_confidence_interval_size(contig_name, v1, v2, sr_deletions, bam_file, stats, config.min_sv_size);
    calculate_confidence_interval_size(contig_name, global_crossing_isize_dist, median_crossing_count_geqi_by_isize, small_dp_deletions, bam_file, stats, config.min_sv_size);
    calculate_ptn_ratio(contig_name, large_dp_deletions, bam_file, stats);
    depth_filter_del(contig_name, deletions, bam_file, config.min_size_for_depth_filtering, stats);

    release_bam_reader(bam_file);
}

void size_and_depth_filtering_dup(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    mtx.lock();
    std::vector<sv2_duplication_t*>& duplications = sv2_duplications_by_chr[contig_name];
    mtx.unlock();
    depth_filter_dup(contig_name, duplications, bam_file, config.min_size_for_depth_filtering, stats);

    release_bam_reader(bam_file);
}

int main(int argc, char* argv[]) {

    complete_bam_fname = argv[1];
	std::string in_vcf_fname = argv[2];
    std::string out_vcf_fname = argv[3];
    workdir = argv[4];
    reference_fname = argv[5];
    std::string sample_name = argv[6];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt");

    chr_seqs.read_fasta_into_map(reference_fname);

    std::ifstream crossing_isizes_dist_fin(workdir + "/crossing_isizes.txt");
	int isize, count;
	while (crossing_isizes_dist_fin >> isize >> count) {
		for (int i = 0; i < count; i++) global_crossing_isize_dist.push_back(isize);
	}
	std::random_shuffle(global_crossing_isize_dist.begin(), global_crossing_isize_dist.end());
	global_crossing_isize_dist.resize(100000);
	crossing_isizes_dist_fin.close();

    std::ifstream crossing_isizes_count_geq_i_fin(workdir + "/crossing_isizes_count_geq_i.txt");
	int median;
	while (crossing_isizes_count_geq_i_fin >> isize >> median) {
		median_crossing_count_geqi_by_isize.push_back(median);
	}
	crossing_isizes_count_geq_i_fin.close();

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	if (in_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
	}
	bcf_hdr_t* sr_vcf_header = bcf_hdr_read(in_vcf_file);
	if (sr_vcf_header == NULL) {
		throw std::runtime_error("Unable to read the header of file " + in_vcf_fname + ".");
	}

	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(in_vcf_file, sr_vcf_header, bcf_entry) == 0) {
		sv_t* sv = bcf_to_sv(sr_vcf_header, bcf_entry);
		if (sv->svtype() == "DEL") {
			deletions_by_chr[sv->chr].push_back((deletion_t*) sv);
		} else if (sv->svtype() == "DUP") {
			duplications_by_chr[sv->chr].push_back((duplication_t*) sv);
		} else {
			insertions_by_chr[sv->chr].push_back((insertion_t*) sv);
		}
	}	

	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
		for (deletion_t* del : deletions) {
			sv2_deletions_by_chr[contig_name].push_back(new sv2_deletion_t(del));
		}

		std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
		for (duplication_t* dup : duplications) {
			sv2_duplications_by_chr[contig_name].push_back(new sv2_duplication_t(dup));
		}
		for (insertion_t* ins : insertions_by_chr[contig_name]) {
			sv2_duplications_by_chr[contig_name].push_back(new sv2_duplication_t(ins));
		}
	}

	std::cout << "Computing statistics for indels." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();

    ctpl::thread_pool thread_pool3(config.threads);
	std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool3.push(size_and_depth_filtering_del, contig_name);
        futures.push_back(std::move(future));
        future = thread_pool3.push(size_and_depth_filtering_dup, contig_name);
        futures.push_back(std::move(future));
    }
    thread_pool3.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

	auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cout << "Statistics computed in " << elapsed_time << " seconds" << std::endl;

	// create VCF out files
	bcf_hdr_t* out_vcf_header = sv2_generate_vcf_header(chr_seqs, sample_name, config, full_cmd_str);
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}


    std::unordered_map<std::string, std::vector<bcf1_t*>> bcf_entries;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!sv2_deletions_by_chr.count(contig_name)) continue;
		std::vector<sv2_deletion_t*>& dels = sv2_deletions_by_chr[contig_name];
		std::sort(dels.begin(), dels.end(), [](sv2_deletion_t* del1, sv2_deletion_t* del2) {
			return std::tie(del1->sv->start, del1->sv->end) < std::tie(del2->sv->start, del2->sv->end);
		});
        for (sv2_deletion_t* del : dels) {
            std::vector<std::string> filters;

            if (-del->sv->svlen() < stats.max_is) {
            	if (-del->sv->svlen()/2 > del->max_conf_size) filters.push_back("SIZE_FILTER");
            }
            if (del->remap_boundary_lower > del->sv->start) {
                filters.push_back("REMAP_BOUNDARY_FILTER");
            } else if (del->remap_boundary_upper < del->sv->end) {
                filters.push_back("REMAP_BOUNDARY_FILTER");
            }
            if (-del->sv->svlen() >= config.min_size_for_depth_filtering &&
                    (del->med_left_flanking_cov*0.74<=del->med_indel_left_cov || del->med_right_flanking_cov*0.74<=del->med_indel_right_cov)) {
            	filters.push_back("DEPTH_FILTER");
            }
            if (del->med_left_flanking_cov > stats.max_depth || del->med_right_flanking_cov > stats.max_depth ||
            	del->med_left_flanking_cov < stats.min_depth || del->med_right_flanking_cov < stats.min_depth ||
				del->med_left_cluster_cov > stats.max_depth || del->med_right_cluster_cov > stats.max_depth) {
                filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (del->med_indel_left_cov > stats.max_depth || del->med_indel_right_cov > stats.max_depth) {
                filters.push_back("ANOMALOUS_DEL_DEPTH");
            }

            if (del->sv->full_junction_aln != NULL && del->sv->left_anchor_aln->best_score+del->sv->right_anchor_aln->best_score-del->sv->full_junction_aln->best_score < config.min_score_diff) {
            	filters.push_back("WEAK_SPLIT_ALIGNMENT");
            }
            if ((del->sv->lc_consensus == NULL || del->sv->lc_consensus->max_mapq < config.high_confidence_mapq) &&
            	(del->sv->rc_consensus == NULL || del->sv->rc_consensus->max_mapq < config.high_confidence_mapq) &&
                 del->sv->source != "DP") {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
            }
            if (del->sv->source == "1SR_RC" || del->sv->source == "1HSR_RC") {
            	if (del->sv->rc_consensus->right_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            } else if (del->sv->source == "1SR_LC" || del->sv->source == "1HSR_LC") {
            	if (del->sv->lc_consensus->left_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            }
            
            if (del->sv->source == "DP") {
                if (del->ks_pval > 0.01) {
                    filters.push_back("KS_FILTER");
                }
                if (-del->sv->svlen() >= stats.max_is && double(del->sv->disc_pairs)/(del->sv->disc_pairs+del->conc_pairs) < 0.25) {
                    filters.push_back("LOW_PTN_RATIO");
                }
                if (-del->sv->svlen() > 10000 && (del->l_cluster_region_disc_pairs >= del->sv->disc_pairs || del->r_cluster_region_disc_pairs >= del->sv->disc_pairs)) {
                    filters.push_back("AMBIGUOUS_REGION");
                }
            }

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            sv2_del2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, del, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
            delete del;
        }
    }


    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!sv2_duplications_by_chr.count(contig_name)) continue;
    	std::vector<sv2_duplication_t*>& dups = sv2_duplications_by_chr[contig_name];
    	std::sort(dups.begin(), dups.end(), [](sv2_duplication_t* dup1, sv2_duplication_t* dup2) {
    		return std::tie(dup1->sv->start, dup1->sv->end) < std::tie(dup2->sv->start, dup2->sv->end);
    	});
        for (sv2_duplication_t* dup : dups) {
            std::vector<std::string> filters;

            if (dup->sv->svlen() >= config.min_size_for_depth_filtering &&
                    (dup->med_left_flanking_cov*1.26>=dup->med_indel_left_cov || dup->med_indel_right_cov<=dup->med_right_flanking_cov*1.26)) {
                // note: using >= so that a 0 0 0 0 depth will not be accepted
                filters.push_back("DEPTH_FILTER");
            }

            if (dup->med_left_flanking_cov > stats.max_depth || dup->med_right_flanking_cov > stats.max_depth ||
                dup->med_left_flanking_cov < stats.min_depth || dup->med_right_flanking_cov < stats.min_depth) {
                filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (dup->sv->full_junction_aln != NULL && dup->sv->left_anchor_aln->best_score+dup->sv->right_anchor_aln->best_score-dup->sv->full_junction_aln->best_score < config.min_score_diff) {
				filters.push_back("WEAK_SPLIT_ALIGNMENT");
			}
            if ((dup->sv->lc_consensus == NULL || dup->sv->lc_consensus->max_mapq < config.high_confidence_mapq) &&
				(dup->sv->rc_consensus == NULL || dup->sv->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
			}
            if (dup->sv->svlen() >= config.min_size_for_depth_filtering && dup->sv->disc_pairs < 3) {
                filters.push_back("NOT_ENOUGH_OW_PAIRS");
            }
            if (dup->sv->source == "1SR_RC" || dup->sv->source == "1HSR_RC") {
				if (dup->sv->rc_consensus->right_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            } else if (dup->sv->source == "1SR_LC" || dup->sv->source == "1HSR_LC") {
            	if (dup->sv->lc_consensus->left_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
			}

            if (filters.empty()) {
                filters.push_back("PASS");
            }

			sv2_dup2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, dup, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
            delete dup;
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	auto& bcf_entries_contig = bcf_entries[contig_name];
    	std::sort(bcf_entries_contig.begin(), bcf_entries_contig.end(), [](const bcf1_t* b1, const bcf1_t* b2) {return b1->pos < b2->pos;});
		for (bcf1_t* bcf_entry : bcf_entries_contig) {
 			if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
			}
		}
    }

    chr_seqs.clear();

    bcf_hdr_destroy(out_vcf_header);
    bcf_close(out_vcf_file);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
