#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <vector>

#include "extend_1sr_consensus.h"
#include "types.h"

std::vector<std::shared_ptr<bam1_t>> gen_consensus_and_find_consistent_seqs_subset(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score);
std::vector<std::shared_ptr<bam1_t>> find_seqs_consistent_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, double& avg_score, double& stddev_score);

void write_reads_associations(std::string id, std::vector<std::shared_ptr<bam1_t>>& reads);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<std::shared_ptr<bam1_t>>& consistent_reads, 
    double consistent_avg_score, double consistent_stddev_score);

void read_mates(int contig_id);
void release_mates(int contig_id);

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension);

#endif // GENOTYPE_H
