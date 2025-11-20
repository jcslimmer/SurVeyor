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
#include "vcf_utils.h"

struct evidence_logger_t {

    std::ofstream alt_reads_to_sv_associations, alt_pairs_to_sv_associations;
    std::mutex mtx;

    evidence_logger_t(const std::string& workdir) {
        alt_reads_to_sv_associations.open(workdir + "/alt_reads_to_sv_associations.txt");
        alt_pairs_to_sv_associations.open(workdir + "/alt_pairs_to_sv_associations.txt");
    }

    void log_pair_association(const std::string& sv_id, bam1_t* pair) {
        std::lock_guard<std::mutex> lock(mtx);
        alt_pairs_to_sv_associations << sv_id << " " << bam_get_qname(pair) << std::endl;
    }

    void log_reads_associations(std::string sv_id, int bp_n, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<int>& scores) {
        std::lock_guard<std::mutex> lock(mtx);
        for (size_t i = 0; i < reads.size(); i++) {
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << bam_get_qname(reads[i].get()) << " " << scores[i] << std::endl;
        }
    }

    evidence_logger_t(const evidence_logger_t&) = delete;
    evidence_logger_t& operator=(const evidence_logger_t&) = delete;
};

struct evidence_map_t {
    std::unordered_map<std::string, std::string> read_to_sv_map;
    std::unordered_map<std::string, int> read_to_bp_map;
    std::unordered_map<std::string, std::vector<std::string>> read_to_non_chosen_svs_map;

    evidence_map_t() {}

    void load(std::string workdir, std::string vcf_fname) {
        
        htsFile* vcf_file = bcf_open(vcf_fname.c_str(), "r");
        if (vcf_file == NULL) {
            throw std::runtime_error("Unable to open file " + vcf_fname + ".");
        }

        bcf_hdr_t* vcf_header = bcf_hdr_read(vcf_file);
        if (vcf_header == NULL) {
            throw std::runtime_error("Failed to read the VCF header.");
        }

        std::unordered_map<std::string, float> sv_epr_map;

        bcf1_t* vcf_record = bcf_init();
        while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {
            bcf_unpack(vcf_record, BCF_UN_ALL);

            std::string id = vcf_record->d.id;
            float epr = get_sv_epr(vcf_header, vcf_record);
            sv_epr_map[id] = epr;
        }
        hts_close(vcf_file);
        bcf_hdr_destroy(vcf_header);
        bcf_destroy(vcf_record);

        std::string alt_reads_association_fname = workdir + "/alt_reads_to_sv_associations.txt";
        std::ifstream alt_reads_association_fin(alt_reads_association_fname);
        std::string sv_id, read_name;
        int bp, score;
        std::unordered_map<std::string, std::pair<int, float>> read_to_score_epr_map;
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score) {
            float epr = sv_epr_map[sv_id];
            std::pair<int, float> p = {score, epr};
            // remove _DUP suffix if present
            if (sv_id.size() > 4 && sv_id.substr(sv_id.size()-4) == "_DUP") {
                sv_id = sv_id.substr(0, sv_id.size()-4);
            }
            if (p > read_to_score_epr_map[read_name]) {
                read_to_score_epr_map[read_name] = p; // Store the highest score and EPR for the read
                if (read_to_sv_map.count(read_name)) {
                    read_to_non_chosen_svs_map[read_name].push_back(read_to_sv_map[read_name]);
                }
                read_to_sv_map[read_name] = sv_id;
                read_to_bp_map[read_name] = bp;
            } else {
                read_to_non_chosen_svs_map[read_name].push_back(sv_id);
            }
        }
    }

    bool is_read_assigned_to_different_sv(std::string read_name, std::string sv_id) {
        if (!read_to_sv_map.count(read_name)) return false;
        // remove _DUP suffix if present
        if (sv_id.size() > 4 && sv_id.substr(sv_id.size()-4) == "_DUP") {
            sv_id = sv_id.substr(0, sv_id.size()-4);
        }
        return read_to_sv_map[read_name] != sv_id;
    }

    std::vector<std::string> get_non_chosen_svs_for_read(std::string read_name) {
        if (!read_to_non_chosen_svs_map.count(read_name)) return {};
        return read_to_non_chosen_svs_map[read_name];
    }
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
