#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <ostream>
#include <queue>
#include <sstream>
#include <unordered_map>

#include "../libs/cptl_stl.h"
#include "htslib/hts.h"
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

bam_pool_t* del_bam_pool, *dup_bam_pool, *ins_bam_pool;

std::unordered_map<std::string, std::vector<deletion_t*> > deletions_by_chr;
std::unordered_map<std::string, std::vector<duplication_t*> > duplications_by_chr;
std::unordered_map<std::string, std::vector<insertion_t*> > insertions_by_chr;
std::unordered_map<std::string, std::vector<inversion_t*> > inversions_by_chr;

void size_and_depth_filtering_del(int id, std::string contig_name) {
    open_samFile_t* bam_file = del_bam_pool->get_bam_reader(id);

    mtx.lock();
    std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
    mtx.unlock();
    // TODO: we are only computing the ks-test p-value for large deletions and DP deletions for efficienty reasons. Investigate why it takes so long calculating it for all deletions.
    std::vector<deletion_t*> large_deletions, small_deletions;
    std::vector<sv_t*> small_svs;
    for (deletion_t* del : deletions) {
        if (-del->svlen() >= stats.max_is) {
            large_deletions.push_back(del);
        } else {
            small_deletions.push_back(del);
            small_svs.push_back(del);
        }
    }

    depth_filter_del(contig_name, deletions, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, global_crossing_isize_dist, small_svs, bam_file, config, stats, config.min_sv_size);
    calculate_ptn_ratio(contig_name, large_deletions, bam_file, config, stats);
    calculate_cluster_region_disc(contig_name, deletions, bam_file);
}

void size_and_depth_filtering_dup(int id, std::string contig_name) {
    open_samFile_t* bam_file = dup_bam_pool->get_bam_reader(id);

    mtx.lock();
    std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
    mtx.unlock();
    depth_filter_dup(contig_name, duplications, bam_file, config, stats);
    calculate_cluster_region_disc(contig_name, duplications, bam_file, stats);

}

void size_and_depth_filtering_ins(int id, std::string contig_name) {
    open_samFile_t* bam_file = ins_bam_pool->get_bam_reader(id);

    mtx.lock();
    std::vector<insertion_t*>& insertions = insertions_by_chr[contig_name];
    mtx.unlock();
    depth_filter_ins(contig_name, insertions, bam_file, config, stats);
    std::vector<insertion_t*> assembled_insertions;
    for (insertion_t* ins : insertions) {
        if (ins->source == "REFERENCE_GUIDED_ASSEMBLY" || ins->source == "DE_NOVO_ASSEMBLY") {
            assembled_insertions.push_back(ins);
        }
    }
    calculate_ptn_ratio(contig_name, assembled_insertions, bam_file, stats);
}

void size_and_depth_filtering_inv(int id, std::string contig_name) {
    open_samFile_t* bam_file = del_bam_pool->get_bam_reader(id);

    mtx.lock();
    std::vector<inversion_t*>& inversions = inversions_by_chr[contig_name];
    mtx.unlock();
    depth_filter_inv(contig_name, inversions, bam_file, config, stats);
    calculate_ptn_ratio(contig_name, inversions, bam_file, stats);
}

void apply_ALL_filters(sv_t* sv) {
    if (sv->median_left_flanking_cov > stats.get_max_depth(sv->chr) || sv->median_right_flanking_cov > stats.get_max_depth(sv->chr) ||
        sv->median_left_flanking_cov < stats.get_min_depth(sv->chr) || sv->median_right_flanking_cov < stats.get_min_depth(sv->chr)) {
        sv->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
    }
    if (sv->full_junction_aln != NULL && sv->left_anchor_aln->best_score+sv->right_anchor_aln->best_score-sv->full_junction_aln->best_score < config.min_score_diff) {
        sv->filters.push_back("WEAK_SPLIT_ALIGNMENT");
    }
    if ((sv->lc_consensus == NULL || sv->lc_consensus->max_mapq < config.high_confidence_mapq) &&
        (sv->rc_consensus == NULL || sv->rc_consensus->max_mapq < config.high_confidence_mapq) &&
        sv->source != "DP" && sv->source != "REFERENCE_GUIDED_ASSEMBLY" && sv->source != "DE_NOVO_ASSEMBLY") {
        sv->filters.push_back("LOW_MAPQ_CONSENSUSES");
    }
    if (sv->source == "1SR_RC" || sv->source == "1HSR_RC") {
        if (sv->rc_consensus->right_ext_reads < 3) sv->filters.push_back("FAILED_TO_EXTEND");
    } else if (sv->source == "1SR_LC" || sv->source == "1HSR_LC") {
        if (sv->lc_consensus->left_ext_reads < 3) sv->filters.push_back("FAILED_TO_EXTEND");
    }
}

int main(int argc, char* argv[]) {

    complete_bam_fname = argv[1];
	std::string in_vcf_fname = argv[2];
    workdir = argv[3];
    reference_fname = argv[4];
    std::string sample_name = argv[5];

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    chr_seqs.read_fasta_into_map(reference_fname);
    del_bam_pool = new bam_pool_t(config.threads, complete_bam_fname, reference_fname);
    dup_bam_pool = new bam_pool_t(config.threads, complete_bam_fname, reference_fname);
    ins_bam_pool = new bam_pool_t(config.threads, complete_bam_fname, reference_fname);

    std::ifstream crossing_isizes_dist_fin(workdir + "/crossing_isizes.txt");
	int isize, count;
	while (crossing_isizes_dist_fin >> isize >> count) {
		for (int i = 0; i < count; i++) global_crossing_isize_dist.push_back(isize);
	}
	std::random_shuffle(global_crossing_isize_dist.begin(), global_crossing_isize_dist.end());
	global_crossing_isize_dist.resize(100000);
	crossing_isizes_dist_fin.close();

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
		} else if (sv->svtype() == "INS") {
			insertions_by_chr[sv->chr].push_back((insertion_t*) sv);
		} else if (sv->svtype() == "INV") {
            inversions_by_chr[sv->chr].push_back((inversion_t*) sv);
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
        future = thread_pool3.push(size_and_depth_filtering_inv, contig_name);
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
            apply_ALL_filters(del);
            if (del->median_left_cluster_cov > stats.get_max_depth(contig_name) || del->median_right_cluster_cov > stats.get_max_depth(contig_name)) {
                del->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (del->median_indel_left_cov > stats.get_max_depth(contig_name) || del->median_indel_right_cov > stats.get_max_depth(contig_name)) {
                del->filters.push_back("ANOMALOUS_DEL_DEPTH");
            }
            if (del->rc_reads() > stats.get_max_depth(del->chr) || del->lc_reads() > stats.get_max_depth(del->chr)) {
                del->filters.push_back("ANOMALOUS_SC_NUMBER");
            }

            if (del->source == "DP") {
                if (del->ks_pval > 0.01) {
                    del->filters.push_back("KS_FILTER");
                }
                if (-del->svlen() >= stats.max_is && double(del->disc_pairs_lf)/(del->disc_pairs_lf+del->conc_pairs_midp) < 0.25) {
                    del->filters.push_back("LOW_PTN_RATIO");
                }
                if (-del->svlen() > 10000 && (del->l_cluster_region_disc_pairs >= del->disc_pairs_lf || del->r_cluster_region_disc_pairs >= del->disc_pairs_lf)) {
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

            apply_ALL_filters(dup);
            if (dup->svlen() >= config.min_size_for_depth_filtering && dup->disc_pairs_lf < 3) {
                dup->filters.push_back("NOT_ENOUGH_DISC_PAIRS");
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
            apply_ALL_filters(ins);
            if (ins->median_left_cluster_cov > stats.get_max_depth(contig_name) || ins->median_right_cluster_cov > stats.get_max_depth(contig_name)) {
                ins->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (ins->rc_reads() > stats.get_max_depth(ins->chr) || ins->lc_reads() > stats.get_max_depth(ins->chr)) {
                ins->filters.push_back("ANOMALOUS_SC_NUMBER");
            }

            if (ins->left_anchor_aln && double(ins->left_anchor_aln->best_score)/ins->left_anchor_aln->seq_len < 0.5) {
                ins->filters.push_back("WEAK_ANCHOR");
            }
            if (ins->right_anchor_aln && double(ins->right_anchor_aln->best_score)/ins->right_anchor_aln->seq_len < 0.5) {
                ins->filters.push_back("WEAK_ANCHOR");
            }

            if (ins->source == "REFERENCE_GUIDED_ASSEMBLY" || ins->source == "DE_NOVO_ASSEMBLY") {
                auto support = {ins->disc_pairs_rf+ins->rc_reads(), ins->disc_pairs_lf+ins->lc_reads()};
                
                int r_positive = ins->disc_pairs_rf + ins->rc_reads();
                int r_negative = ins->conc_pairs_lbp;
                int l_positive = ins->disc_pairs_lf + ins->lc_reads();
                int l_negative = ins->conc_pairs_rbp;
                if (r_positive < stats.get_median_depth(ins->chr)/5 && l_positive < stats.get_median_depth(ins->chr)/5) {
		            ins->filters.push_back("LOW_SUPPORT");
                }
                if (double(r_positive)/(r_positive+r_negative) < 0.25 || double(l_positive)/(l_positive+l_negative) < 0.25) {
                    ins->filters.push_back("LOW_PTN_RATIO");
                }

                int svinslen = ins->ins_seq.length();
	            if (ins->disc_pairs_rf + ins->disc_pairs_lf < stats.get_min_disc_pairs_by_insertion_size(svinslen)) {
                    ins->filters.push_back("NOT_ENOUGH_DISC_PAIRS");
                }
            }

            int d = ins->ins_seq.find("-");
            std::string ins_seq_fh = ins->ins_seq.substr(0, d);
            std::string ins_seq_sh = ins->ins_seq.substr(d+1);
            if (is_homopolymer(ins_seq_fh) || is_homopolymer(ins_seq_sh)) {
                ins->filters.push_back("HOMOPOLYMER_INSSEQ");
            }

            if (ins->filters.empty()) {
                ins->filters.push_back("PASS");
            }

            sv_entries[contig_name].push_back(ins);
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
        if (!inversions_by_chr.count(contig_name)) continue;
        std::vector<inversion_t*>& inversions = inversions_by_chr[contig_name];
        for (inversion_t* inv : inversions) {
            // apply_ALL_filters(inv);
            if (inv->median_left_flanking_cov > stats.get_max_depth(inv->chr) || inv->median_right_flanking_cov > stats.get_max_depth(inv->chr) ||
                inv->median_left_flanking_cov < stats.get_min_depth(inv->chr) || inv->median_right_flanking_cov < stats.get_min_depth(inv->chr) ||
                inv->median_left_cluster_cov > stats.get_max_depth(inv->chr) || inv->median_right_cluster_cov > stats.get_max_depth(inv->chr)) {
                inv->filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (inv->rc_reads() > stats.get_max_depth(inv->chr) || inv->lc_reads() > stats.get_max_depth(inv->chr)) {
                inv->filters.push_back("ANOMALOUS_SC_NUMBER");
            }
            if (inv->disc_pairs_lf_maxmapq < config.high_confidence_mapq || inv->disc_pairs_rf_maxmapq < config.high_confidence_mapq) {
                inv->filters.push_back("LOW_MAPQ_DISC_PAIRS");
            }

            if (inv->disc_pairs_lf < stats.get_min_disc_pairs_by_insertion_size(inv->svlen())/2) {
                inv->filters.push_back("NOT_ENOUGH_DISC_PAIRS");
            }

            double ptn_ratio_lbp = double(inv->disc_pairs_rf)/(inv->disc_pairs_rf+inv->conc_pairs_lbp);
            double ptn_ratio_rbp = double(inv->disc_pairs_lf)/(inv->disc_pairs_lf+inv->conc_pairs_rbp);
            if (ptn_ratio_lbp < 0.25 || ptn_ratio_rbp < 0.25) {
                inv->filters.push_back("LOW_PTN_RATIO");
            }

            if (inv->right_anchor_aln->end-inv->left_anchor_aln->start < stats.max_is/2 || inv->rbp_right_anchor_aln->end-inv->rbp_left_anchor_aln->start < stats.max_is/2) {
                inv->filters.push_back("SHORT_ANCHOR");
            }

            if (inv->filters.empty()) {
                inv->filters.push_back("PASS");
            }
            sv_entries[contig_name].push_back(inv);
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
		}
    }

    chr_seqs.clear();

    bcf_close(out_vcf_file);
    bcf_hdr_destroy(out_vcf_header);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
