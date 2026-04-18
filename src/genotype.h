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

constexpr double MIN_EPR = 0.10;

struct bp_support_read_t {
    std::string read_name;
    int64_t mapq, mate_mapq;
    std::string seq;
    bool mate_is_reverse;
    hts_pos_t mate_pos;
    hts_pos_t mate_endpos;
    bool is_first_in_pair;

    bp_support_read_t() : mapq(0), mate_mapq(0), mate_is_reverse(false), mate_pos(-1), mate_endpos(-1), is_first_in_pair(false) {}

    explicit bp_support_read_t(bam1_t* read) :
        read_name(bam_get_qname(read)),
        mapq(read->core.qual),
        mate_mapq(get_mq(read)),
        seq(get_sequence(read)),
        mate_is_reverse(bam_is_mrev(read)),
        mate_pos(read->core.mpos),
        mate_endpos(get_mate_endpos(read)),
        is_first_in_pair(is_first_read(read))
     {}
};


std::string read_name_with_suffix(bp_support_read_t& read) {
    return read.read_name + (read.is_first_in_pair ? "/1" : "/2");
}
std::string read_name_with_suffix(bam1_t* read) {
    return std::string(bam_get_qname(read)) + ((is_first_read(read)) ? "/1" : "/2");
}

std::string remove_svid_dup_suffix(const std::string& sv_id) {
    if (sv_id.size() > 4 && sv_id.substr(sv_id.size()-4) == "_DUP") {
        return sv_id.substr(0, sv_id.size()-4);
    }
    return sv_id;
}

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
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(reads[i].get()) << " " << scores[i] << std::endl;
        }
    }
    void log_reads_associations(std::string sv_id, int bp_n, std::vector<bp_support_read_t>& reads, std::vector<int>& scores) {
        std::lock_guard<std::mutex> lock(mtx);
        for (size_t i = 0; i < reads.size(); i++) {
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(reads[i]) << " " << scores[i] << std::endl;
        }
    }

    evidence_logger_t(const evidence_logger_t&) = delete;
    evidence_logger_t& operator=(const evidence_logger_t&) = delete;
};

struct evidence_map_t {
    std::unordered_map<std::string, std::string> read_to_sv_map;
    std::unordered_map<std::string, int> read_to_bp_map;
    std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> read_to_non_chosen_svs_map;

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
        // For each read, keep only the best association
        // Criteria for best: predicted as existing, highest score, highest EPR
        std::unordered_map<std::string, std::pair<int, float>> read_to_score_epr_map;
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score) {
            float epr = sv_epr_map[sv_id];
            if (epr < MIN_EPR && epr != -1.0) continue; // only consider SVs predicted as existing
            std::pair<int, float> p = {score, epr};
            sv_id = remove_svid_dup_suffix(sv_id); // we avoid INS and INS_TO_DUP from stealing each other's reads
            if (p > read_to_score_epr_map[read_name]) {
                read_to_score_epr_map[read_name] = p; // Store the highest score and EPR for the read
                if (read_to_sv_map.count(read_name)) {
                    read_to_non_chosen_svs_map[read_name].push_back({read_to_sv_map[read_name], read_to_bp_map[read_name]});
                }
                read_to_sv_map[read_name] = sv_id;
                read_to_bp_map[read_name] = bp;
            } else {
                read_to_non_chosen_svs_map[read_name].push_back({sv_id, bp});
            }
        }

        for (auto& kv : read_to_sv_map) {
            std::string read_name = kv.first;
            std::string sv_id = kv.second;
            // remove sv_id from non-chosen svs if it is there
            auto& vec = read_to_non_chosen_svs_map[read_name];
            vec.erase(std::remove_if(vec.begin(), vec.end(), [&](const std::pair<std::string, int>& p) {
                return p.first == sv_id;
            }), vec.end());
        }
    }

    bool is_read_assigned_to_different_sv(bam1_t* read, std::string sv_id) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_sv_map.count(read_name)) return false;
        sv_id = remove_svid_dup_suffix(sv_id);
        return read_to_sv_map[read_name] != sv_id;
    }

    std::vector<std::pair<std::string, int>> get_non_chosen_svs_for_read(bam1_t* read) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_non_chosen_svs_map.count(read_name)) return {};
        return read_to_non_chosen_svs_map[read_name];
    }

    std::vector<std::pair<std::string, int>> get_non_chosen_svs_for_read(const bp_support_read_t& read) {
        std::string read_name = read.read_name + (read.is_first_in_pair ? "/1" : "/2");
        if (!read_to_non_chosen_svs_map.count(read_name)) return {};
        return read_to_non_chosen_svs_map[read_name];
    }
};

std::mutex orc_mtx;

void increase_orc_supp(sv_t::sample_info_t& sample_info, int bp_n, bool hq) {
    std::lock_guard<std::mutex> lock(orc_mtx);
    if (bp_n == 1) {
        sample_info.assigned_to_other_sv_bp1_consistent++;
        if (hq) {
            sample_info.assigned_to_other_sv_bp1_consistent_highmq++;
        }
    } else if (bp_n == 2) {
        sample_info.assigned_to_other_sv_bp2_consistent++;
        if (hq) {
            sample_info.assigned_to_other_sv_bp2_consistent_highmq++;
        }
    }
}

void increase_orc(std::unordered_map<std::string, std::shared_ptr<sv_t>>& sv_map, std::string sv_id, int bp_n, bool hq) {
    if (!sv_map.count(sv_id)) return;
    increase_orc_supp(sv_map[sv_id]->sample_info, bp_n, hq);

    sv_id += "_DUP";
    if (!sv_map.count(sv_id)) return;
    increase_orc_supp(sv_map[sv_id]->sample_info, bp_n, hq);
}

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs);
std::vector<std::shared_ptr<bam1_t>> gen_consensus_and_find_consistent_seqs_subset(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, 
    std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);
std::vector<std::shared_ptr<bam1_t>> find_seqs_consistent_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, 
    double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<bp_support_read_t>& consistent_reads,
    std::vector<bool>& is_exact_read, double consistent_avg_score, double consistent_stddev_score);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<std::shared_ptr<bam1_t>>& consistent_reads, 
    std::vector<bool>& is_exact_read, double consistent_avg_score, double consistent_stddev_score);

void read_mates(int contig_id);
void release_mates(int contig_id);

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension);


char* generate_haplotype_left(char* chrom_seq, hts_pos_t hap_end, hts_pos_t hap_len, 
    std::vector<std::shared_ptr<sv_t>>& aux_indels, std::vector<snp_t>& aux_snps);
char* generate_haplotype_right(char* chrom_seq, hts_pos_t chrom_len, hts_pos_t hap_start, hts_pos_t hap_len,
    std::vector<std::shared_ptr<sv_t>>& aux_indels, std::vector<snp_t>& aux_snps);

// Given a sequence alt_seq, a series of sequences ref_seqs and a read length read_len,
// return all the positions in alt_seq where a read of length read_len that is not present 
// in any of the ref_seqs starts.
std::vector<hts_pos_t> get_diff_reads_expected_positions(std::vector<char*>& ref_seqs, std::vector<hts_pos_t>& ref_lens, char* alt_seq, hts_pos_t alt_len, int read_len) {
    std::vector<hts_pos_t> positions;
    if (read_len > alt_len || read_len <= 0) {
        return positions;
    }
    for (hts_pos_t i = 0; i <= alt_len - read_len; i++) {
        char* read_begin = alt_seq + i;
        char* read_end = read_begin + read_len;
        bool found_in_ref = false;
        for (size_t j = 0; j < ref_seqs.size(); j++) {
            char* ref_begin = ref_seqs[j];
            char* ref_end = ref_begin + ref_lens[j];
            auto it = std::search(ref_begin, ref_end, read_begin, read_end);
            if (it != ref_end) {
                found_in_ref = true;
                break;
            }
        }
        if (!found_in_ref) {
            positions.push_back(i);
        }
    }
    return positions;
}

std::vector<int> get_consistent_reads_start_positions(std::vector<std::shared_ptr<bam1_t>>& consistent_reads,
    std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<int>& start_positions) {
    std::vector<int> consistent_positions;
    std::unordered_set<bam1_t*> consistent_set;
    for (auto& r : consistent_reads) {
        consistent_set.insert(r.get());
    }

    for (int i = 0; i < reads.size(); i++) {
        if (consistent_set.count(reads[i].get())) {
            consistent_positions.push_back(start_positions[i]);
        }
    }
    return consistent_positions;
}

// Given a vector of observed positions (with duplicates allowed) and the number of distinct valid positions
// observable, compute the occupancy ratio as the number of distinct observed positions divided 
// by expected number of observed positions under uniform distribution of observations across valid positions.
// We can assume that all positions in 'positions' are valid.
// If U is the number of unique observed positions, E[U] =  k * (1 - (1 - 1/k)^n), where k is valid_positions and n is positions.size()
double occ_ratio(std::vector<int>& positions, int valid_positions) {
    if (valid_positions <= 0) return sv_t::sample_info_t::NOT_COMPUTED;

    int n = positions.size();
    if (n == 0) return 1.0;

    double k = valid_positions;

    std::unordered_set<int> unique_positions(positions.begin(), positions.end());
    double U = unique_positions.size();
    double EU = k * (1.0 - std::exp((double)n * std::log1p(-1.0 / k)));

    if (EU <= 0.0) return sv_t::sample_info_t::NOT_COMPUTED;
    return U / EU;
}

#endif // GENOTYPE_H
