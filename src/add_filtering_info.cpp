#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <sstream>
#include <unordered_map>

#include "../libs/cptl_stl.h"
#include "stat_tests.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "vcf_utils.h"

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

void size_and_depth_filtering_del(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    mtx.lock();
    std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
    mtx.unlock();
    // TODO: we are only computing the ks-test p-value for large deletions and DP deletions for efficienty reasons. Investigate why it takes so long calculating it for all deletions.
    std::vector<deletion_t*> large_dp_deletions, small_dp_deletions, sr_deletions;
    for (deletion_t* del : deletions) {
        if (-del->svlen() >= stats.max_is && del->source == "DP") {
            large_dp_deletions.push_back(del);
        } else if (del->source == "DP") {
            small_dp_deletions.push_back(del);
        } else {
            sr_deletions.push_back(del);
        }
    }
    std::vector<double> v1;
    std::vector<uint32_t> v2;
    depth_filter_del(contig_name, deletions, bam_file, config.min_size_for_depth_filtering, stats);
    calculate_confidence_interval_size(contig_name, v1, v2, sr_deletions, bam_file, stats, config.min_sv_size);
    calculate_confidence_interval_size(contig_name, global_crossing_isize_dist, median_crossing_count_geqi_by_isize, small_dp_deletions, bam_file, stats, config.min_sv_size);
    calculate_ptn_ratio(contig_name, large_dp_deletions, bam_file, stats);

    release_bam_reader(bam_file);
}

void size_and_depth_filtering_dup(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    mtx.lock();
    std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
    mtx.unlock();
    depth_filter_dup(contig_name, duplications, bam_file, config.min_size_for_depth_filtering, stats);
    calculate_cluster_region_disc(contig_name, duplications, bam_file, stats);

    release_bam_reader(bam_file);
}

void size_and_depth_filtering_ins(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    mtx.lock();
    std::vector<insertion_t*>& insertions = insertions_by_chr[contig_name];
    mtx.unlock();
    depth_filter_ins(contig_name, insertions, bam_file, config.min_size_for_depth_filtering, stats);

    release_bam_reader(bam_file);
}

int main(int argc, char* argv[]) {

    complete_bam_fname = argv[1];
	std::string in_vcf_fname = argv[2];
    workdir = argv[3];
    reference_fname = argv[4];
    std::string sample_name = argv[5];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

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
	bcf_hdr_t* in_vcf_header = bcf_hdr_read(in_vcf_file);
	if (in_vcf_header == NULL) {
		throw std::runtime_error("Unable to read the header of file " + in_vcf_fname + ".");
	}

	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(in_vcf_file, in_vcf_header, bcf_entry) == 0) {
		sv_t* sv = bcf_to_sv(in_vcf_header, bcf_entry);
		if (sv->svtype() == "DEL") {
			deletions_by_chr[sv->chr].push_back((deletion_t*) sv);
		} else if (sv->svtype() == "DUP") {
			duplications_by_chr[sv->chr].push_back((duplication_t*) sv);
		} else {
			insertions_by_chr[sv->chr].push_back((insertion_t*) sv);
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
        future = thread_pool3.push(size_and_depth_filtering_ins, contig_name);
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
    std::string in_prefix = in_vcf_fname.substr(0, in_vcf_fname.size()-7);
	bcf_hdr_t* out_vcf_header = in_vcf_header;
    
    std::string out_vcf_fname = in_prefix + ".annotated.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

    std::string out_pass_vcf_fname = in_prefix + ".annotated.pass.vcf.gz";
    htsFile* out_pass_vcf_file = bcf_open(out_pass_vcf_fname.c_str(), "wz");
    if (bcf_hdr_write(out_pass_vcf_file, out_vcf_header) != 0) {
        throw std::runtime_error("Failed to write the VCF header to " + out_pass_vcf_fname + ".");
    }

    std::unordered_map<std::string, std::vector<sv_t*>> sv_entries;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!deletions_by_chr.count(contig_name)) continue;
		std::vector<deletion_t*>& dels = deletions_by_chr[contig_name];
        for (deletion_t* del : dels) {
            if (-del->svlen() < stats.max_is) {
            	if (-del->svlen()/2 > del->max_conf_size) del->filters.push_back("SIZE_FILTER");
            }
            if (del->remap_boundary_lower() > del->start) {
                del->filters.push_back("REMAP_BOUNDARY_FILTER");
            } else if (del->remap_boundary_upper() < del->end) {
                del->filters.push_back("REMAP_BOUNDARY_FILTER");
            }
            if (-del->svlen() >= config.min_size_for_depth_filtering &&
                    (del->median_left_flanking_cov*0.74<=del->median_indel_left_cov || del->median_right_flanking_cov*0.74<=del->median_indel_right_cov)) {
            	del->filters.push_back("DEPTH_FILTER");
            }
            if (del->median_left_flanking_cov > stats.get_max_depth(contig_name) || del->median_right_flanking_cov > stats.get_max_depth(contig_name) ||
            	del->median_left_flanking_cov < stats.get_min_depth(contig_name) || del->median_right_flanking_cov < stats.get_min_depth(contig_name) ||
				del->median_left_cluster_cov > stats.get_max_depth(contig_name) || del->median_right_cluster_cov > stats.get_max_depth(contig_name)) {
                del->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (del->median_indel_left_cov > stats.get_max_depth(contig_name) || del->median_indel_right_cov > stats.get_max_depth(contig_name)) {
                del->filters.push_back("ANOMALOUS_DEL_DEPTH");
            }

            if (del->full_junction_aln != NULL && del->left_anchor_aln->best_score+del->right_anchor_aln->best_score-del->full_junction_aln->best_score < config.min_score_diff) {
            	del->filters.push_back("WEAK_SPLIT_ALIGNMENT");
            }
            if ((del->lc_consensus == NULL || del->lc_consensus->max_mapq < config.high_confidence_mapq) &&
            	(del->rc_consensus == NULL || del->rc_consensus->max_mapq < config.high_confidence_mapq) &&
                 del->source != "DP") {
				del->filters.push_back("LOW_MAPQ_CONSENSUSES");
            }
            if (del->source == "1SR_RC" || del->source == "1HSR_RC") {
            	if (del->rc_consensus->right_ext_reads < 3) del->filters.push_back("FAILED_TO_EXTEND");
            } else if (del->source == "1SR_LC" || del->source == "1HSR_LC") {
            	if (del->lc_consensus->left_ext_reads < 3) del->filters.push_back("FAILED_TO_EXTEND");
            }
            
            if (del->source == "DP") {
                if (del->ks_pval > 0.01) {
                    del->filters.push_back("KS_FILTER");
                }
                if (-del->svlen() >= stats.max_is && double(del->disc_pairs)/(del->disc_pairs+del->conc_pairs) < 0.25) {
                    del->filters.push_back("LOW_PTN_RATIO");
                }
                if (-del->svlen() > 10000 && (del->l_cluster_region_disc_pairs >= del->disc_pairs || del->r_cluster_region_disc_pairs >= del->disc_pairs)) {
                    del->filters.push_back("AMBIGUOUS_REGION");
                }
            }

            if (del->filters.empty()) {
                del->filters.push_back("PASS");
            }

            sv_entries[contig_name].push_back(del);
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!duplications_by_chr.count(contig_name)) continue;
    	std::vector<duplication_t*>& dups = duplications_by_chr[contig_name];
        for (duplication_t* dup : dups) {
            if (dup->svlen() >= config.min_size_for_depth_filtering &&
                    (dup->median_left_flanking_cov*1.26>=dup->median_indel_left_cov || dup->median_indel_right_cov<=dup->median_right_flanking_cov*1.26)) {
                // note: using >= so that a 0 0 0 0 depth will not be accepted
                dup->filters.push_back("DEPTH_FILTER");
            }

            if (dup->median_left_flanking_cov > stats.get_max_depth(contig_name) || dup->median_right_flanking_cov > stats.get_max_depth(contig_name) ||
                dup->median_left_flanking_cov < stats.get_min_depth(contig_name) || dup->median_right_flanking_cov < stats.get_min_depth(contig_name)) {
                dup->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (dup->full_junction_aln != NULL && dup->left_anchor_aln->best_score+dup->right_anchor_aln->best_score-dup->full_junction_aln->best_score < config.min_score_diff) {
				dup->filters.push_back("WEAK_SPLIT_ALIGNMENT");
			}
            if ((dup->lc_consensus == NULL || dup->lc_consensus->max_mapq < config.high_confidence_mapq) &&
				(dup->rc_consensus == NULL || dup->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				dup->filters.push_back("LOW_MAPQ_CONSENSUSES");
			}
            if (dup->svlen() >= config.min_size_for_depth_filtering && dup->disc_pairs < 3) {
                dup->filters.push_back("NOT_ENOUGH_OW_PAIRS");
            }
            if (dup->source == "1SR_RC" || dup->source == "1HSR_RC") {
				if (dup->rc_consensus->right_ext_reads < 3) dup->filters.push_back("FAILED_TO_EXTEND");
            } else if (dup->source == "1SR_LC" || dup->source == "1HSR_LC") {
            	if (dup->lc_consensus->left_ext_reads < 3) dup->filters.push_back("FAILED_TO_EXTEND");
			}

            if (dup->filters.empty()) {
                dup->filters.push_back("PASS");
            }

            sv_entries[contig_name].push_back(dup);
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
        if (!insertions_by_chr.count(contig_name)) continue;
        std::vector<insertion_t*>& insertions = insertions_by_chr[contig_name];
        for (insertion_t* ins : insertions) {
            if (ins->median_left_flanking_cov > stats.get_max_depth(contig_name) || ins->median_right_flanking_cov > stats.get_max_depth(contig_name) ||
                ins->median_left_flanking_cov < stats.get_min_depth(contig_name) || ins->median_right_flanking_cov < stats.get_min_depth(contig_name)) {
                ins->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (ins->full_junction_aln != NULL && ins->left_anchor_aln->best_score+ins->right_anchor_aln->best_score-ins->full_junction_aln->best_score < config.min_score_diff) {
                ins->filters.push_back("WEAK_SPLIT_ALIGNMENT");
            }
            if ((ins->lc_consensus == NULL || ins->lc_consensus->max_mapq < config.high_confidence_mapq) &&
                (ins->rc_consensus == NULL || ins->rc_consensus->max_mapq < config.high_confidence_mapq)) {
                ins->filters.push_back("LOW_MAPQ_CONSENSUSES");
            }
            if (ins->source == "1SR_RC" || ins->source == "1HSR_RC") {
                if (ins->rc_consensus->right_ext_reads < 3) ins->filters.push_back("FAILED_TO_EXTEND");
            } else if (ins->source == "1SR_LC" || ins->source == "1HSR_LC") {
                if (ins->lc_consensus->left_ext_reads < 3) ins->filters.push_back("FAILED_TO_EXTEND");
            }

            if (ins->filters.empty()) {
                ins->filters.push_back("PASS");
            }

            sv_entries[contig_name].push_back(ins);
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	auto& sv_entries_contig = sv_entries[contig_name];
    	std::sort(sv_entries_contig.begin(), sv_entries_contig.end(), [](const sv_t* sv1, const sv_t* sv2) {return sv1->start < sv2->start;});
		for (sv_t* sv : sv_entries_contig) {
            sv2bcf(out_vcf_header, bcf_entry, sv, chr_seqs.get_seq(contig_name));
 			if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
			}
            if (sv->is_pass()) {
                if (bcf_write(out_pass_vcf_file, out_vcf_header, bcf_entry) != 0) {
                    throw std::runtime_error("Failed to write to " + out_pass_vcf_fname + ".");
                }
            }
		}
    }

    chr_seqs.clear();

    bcf_close(out_vcf_file);
    bcf_close(out_pass_vcf_file);
    bcf_hdr_destroy(out_vcf_header);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
    tbx_index_build(out_pass_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
