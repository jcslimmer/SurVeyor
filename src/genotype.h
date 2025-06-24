#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <fstream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "extend_1sr_consensus.h"
#include "htslib/sam.h"
#include "sam_utils.h"
#include "types.h"

struct evidence_logger_t {

    std::ofstream alt_reads_to_sv_associations;
    std::mutex mtx;

    evidence_logger_t(const std::string& workdir) {
        alt_reads_to_sv_associations.open(workdir + "/alt_reads_to_sv_associations.txt");
    }

    void log_read_association(const std::string& sv_id, bam1_t* pair) {
        std::lock_guard<std::mutex> lock(mtx);
        alt_reads_to_sv_associations << sv_id << " " << bam_get_qname(pair) << "\n";
    }

    void log_reads_associations(std::string sv_id, std::vector<std::shared_ptr<bam1_t>>& reads) {
        std::lock_guard<std::mutex> lock(mtx);
        for (std::shared_ptr<bam1_t> read : reads) {
            alt_reads_to_sv_associations << sv_id << " " << bam_get_qname(read.get()) << std::endl;
        }
    }

    evidence_logger_t(const evidence_logger_t&) = delete;
    evidence_logger_t& operator=(const evidence_logger_t&) = delete;
};

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs);
std::vector<std::shared_ptr<bam1_t>> gen_consensus_and_find_consistent_seqs_subset(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score);
std::vector<std::shared_ptr<bam1_t>> find_seqs_consistent_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, double& avg_score, double& stddev_score);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<std::shared_ptr<bam1_t>>& consistent_reads, 
    double consistent_avg_score, double consistent_stddev_score);

void read_mates(int contig_id);
void release_mates(int contig_id);

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension);

#endif // GENOTYPE_H
