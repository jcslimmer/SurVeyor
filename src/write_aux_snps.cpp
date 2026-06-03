#include <algorithm>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>
#include <sstream>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "types.h"
#include "vcf_utils.h"

const hts_pos_t MIN_SV_SIZE = 50;
const float MIN_PRIMARY_PPR = 0.5f;

bool is_literal_allele_record(bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele < 2) return false;
    const char* alt = record->d.allele[1];
    return alt[0] != '<' && strchr(alt, '[') == NULL && strchr(alt, ']') == NULL;
}

bool same_record_identity(bcf_hdr_t* hdr, bcf1_t* b1, bcf1_t* b2) {
    bcf_unpack(b1, BCF_UN_STR);
    bcf_unpack(b2, BCF_UN_STR);

    if (b1->rid != b2->rid || b1->pos != b2->pos) return false;

    bool b1_is_literal = is_literal_allele_record(b1);
    bool b2_is_literal = is_literal_allele_record(b2);
    if (b1_is_literal != b2_is_literal) return false;

    if (b1_is_literal) {
        return strcmp(b1->d.allele[0], b2->d.allele[0]) == 0 &&
               strcmp(b1->d.allele[1], b2->d.allele[1]) == 0;
    }

    return get_sv_end(hdr, b1) == get_sv_end(hdr, b2) && get_ins_seq(hdr, b1) == get_ins_seq(hdr, b2) 
        && get_sv_type(hdr, b1) == get_sv_type(hdr, b2);
}

bool is_small_variant_by_svsize(bcf_hdr_t* hdr, bcf1_t* record) {
    std::string svtype = get_sv_type(hdr, record);
    if (svtype == "BND") return false;

    if (svtype == "DEL" || svtype == "DUP" || svtype == "INS" || svtype == "INV") {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, record);
        if (sv == nullptr) return false;
        if (svtype == "INS" && sv->incomplete_ins_seq()) return false;
        return sv->svsize() < MIN_SV_SIZE;
    }

    int len = get_sv_len(hdr, record);
    return len != bcf_int32_missing && len > -MIN_SV_SIZE && len < MIN_SV_SIZE;
}

void divide_variants_by_svsize(bcf_hdr_t* hdr, const std::vector<bcf1_t*>& vcf_records,
                              std::vector<bcf1_t*>& small_variants, std::vector<bcf1_t*>& svs) {
    small_variants.clear();
    svs.clear();

    for (bcf1_t* record : vcf_records) {
        if (is_small_variant_by_svsize(hdr, record)) {
            small_variants.push_back(record);
        } else {
            svs.push_back(record);
        }
    }
}

void write_and_remove_non_primary(bcf_hdr_t* hdr, std::vector<bcf1_t*>& variants, const std::string& out_vcf_fname) {
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    if (out_vcf_file == NULL) {
        throw std::runtime_error("Failed to open " + out_vcf_fname + " for writing.");
    }

    if (bcf_hdr_write(out_vcf_file, hdr) != 0) {
        throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
    }

    for (bcf1_t*& variant : variants) {
        float ppr = get_sv_ppr(hdr, variant);
        if (ppr >= 0.0f && ppr < MIN_PRIMARY_PPR) {
            if (bcf_write(out_vcf_file, hdr, variant) != 0) {
                throw std::runtime_error("Failed to write VCF record to " + out_vcf_fname);
            }
            variant = nullptr;
        }
    }
    variants.erase(std::remove(variants.begin(), variants.end(), nullptr), variants.end());

    hts_close(out_vcf_file);
}

std::string bcf_record_unique_key(bcf_hdr_t* hdr, bcf1_t* record) {
    return std::string(bcf_seqname_safe(hdr, record)) + ":" +
        std::to_string(record->pos) + ":" +
        std::to_string(get_sv_end(hdr, record)) + ":" +
        get_sv_type(hdr, record) + ":" +
        get_ins_seq(hdr, record);
}

void remove_aux_indels_matching_primary(bcf_hdr_t* hdr, std::vector<bcf1_t*>& vcf_records,
                                        const std::unordered_set<bcf1_t*>& aux_indel_records,
                                        const std::unordered_set<std::string>& non_aux_variant_keys) {
    for (bcf1_t*& record : vcf_records) {
        if (aux_indel_records.count(record) == 0) continue;

        if (non_aux_variant_keys.count(bcf_record_unique_key(hdr, record)) > 0) {
            bcf_destroy(record);
            record = nullptr;
        }
    }

    vcf_records.erase(std::remove(vcf_records.begin(), vcf_records.end(), nullptr), vcf_records.end());
}

std::string get_record_id(bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    return record->d.id == NULL ? "." : std::string(record->d.id);
}

void set_gt(bcf_hdr_t* hdr, bcf1_t* record, int allele1, int allele2) {
    int gt[2] = { bcf_gt_unphased(allele1), bcf_gt_unphased(allele2) };
    bcf_update_genotypes(hdr, record, gt, 2);
}

void concat_record_ids(bcf_hdr_t* hdr, bcf1_t* record, bcf1_t* other_record) {
    std::string id = get_record_id(record);
    std::string other_id = get_record_id(other_record);
    if (id == ".") {
        id = other_id;
    } else if (other_id != ".") {
        id += ";" + other_id;
    }
    bcf_update_id(hdr, record, id.c_str());
}

void make_record_multiallelic(bcf_hdr_t* hdr, bcf1_t* record, bcf1_t* second_allele_record) {
    bcf_unpack(record, BCF_UN_STR);
    bcf_unpack(second_allele_record, BCF_UN_STR);
    std::string ref = record->d.allele[0];
    std::string alt1 = record->d.allele[1];
    std::string alt2 = second_allele_record->d.allele[1];

    std::string ref2 = second_allele_record->d.allele[0];
    if (ref2.length() > ref.length() && ref2.compare(0, ref.length(), ref) == 0) {
        alt1 += ref2.substr(ref.length());
        ref = ref2;
    } else if (ref.length() > ref2.length() && ref.compare(0, ref2.length(), ref2) == 0) {
        alt2 += ref.substr(ref2.length());
    }

    std::string alleles = ref + "," + alt1 + "," + alt2;
    bcf_update_alleles_str(hdr, record, alleles.c_str());
    concat_record_ids(hdr, record, second_allele_record);
    set_gt(hdr, record, 1, 2);
}

void add_star_allele(bcf_hdr_t* hdr, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele != 2 || count_alt_alleles(hdr, record) != 1) {
        return;
    }

    std::string alleles = record->d.allele[0];
    for (int i = 1; i < record->n_allele; i++) {
        if (strcmp(record->d.allele[i], "*") == 0) {
            return;
        }
        alleles += ",";
        alleles += record->d.allele[i];
    }

    int star_allele = record->n_allele;
    alleles += ",*";
    bcf_update_alleles_str(hdr, record, alleles.c_str());
    set_gt(hdr, record, 1, star_allele);
}

void add_star_alleles_for_overlapping_del_ins(bcf_hdr_t* hdr, std::vector<bcf1_t*>& small_variants) {
    int curr_rid = -1;
    hts_pos_t max_del_end = -1;

    for (bcf1_t* record : small_variants) {
        if (record->rid != curr_rid) {
            curr_rid = record->rid;
            max_del_end = -1;
        }

        std::string svtype = get_sv_type(hdr, record);
        int alt_alleles = count_alt_alleles(hdr, record);
        if (alt_alleles == 1 && record->pos <= max_del_end) {
            add_star_allele(hdr, record);
        }
        if (svtype == "DEL" && alt_alleles > 0) {
            max_del_end = std::max(max_del_end, (hts_pos_t) get_sv_end(hdr, record));
        }
    }
}

// The group is already sorted by decreasing EPR. Keep at most two ALT alleles at
// each exact position for INS/DEL groups. The best call is always kept. If the first two
// calls are identical hets, merge them into a single 1/1 call. Otherwise, if the
// best call is already 1/1, let the second call replace one allele only when its
// EPR is stronger than the best call's homozygous probability. All later ALT
// calls are suppressed. Accepted second alleles are folded into the first record
// as ALT2, and a 0/0 record ends the group because no later record can add an
// ALT allele.
void apply_multiallelic_logic_to_group(bcf_hdr_t* hdr, std::vector<bcf1_t*>& indel_group,
                                       std::vector<bcf1_t*>& records_to_remove) {
    
    if (indel_group.empty()) return;

    bcf1_t* first = indel_group[0];
    int first_alt_alleles = count_alt_alleles(hdr, first);
    if (first_alt_alleles == 0) return;

    if (indel_group.size() == 1) return;

    bcf1_t* second = indel_group[1];
    int second_alt_alleles = count_alt_alleles(hdr, second);
    if (second_alt_alleles == 0) return;

    size_t suppress_from = 2;

    if (first_alt_alleles > 0 && second_alt_alleles > 0 && same_record_identity(hdr, first, second)) {
        if (first_alt_alleles == 1) concat_record_ids(hdr, first, second);
        set_gt(hdr, first, 1, 1);
        set_gt(hdr, second, 0, 0);
        records_to_remove.push_back(second);
    } else if (first_alt_alleles >= 2) {
        float second_epr = get_sv_epr(hdr, second);
        float first_hopr = get_sv_hopr(hdr, first);
        if (second_epr > first_hopr) {
            set_gt(hdr, first, 0, 1);
            if (second_alt_alleles > 1) {
                set_gt(hdr, second, 0, 1);
            }
            make_record_multiallelic(hdr, first, second);
            records_to_remove.push_back(second);
        } else {
            set_gt(hdr, second, 0, 0);
            records_to_remove.push_back(second);
            suppress_from = 2;
        }
    } else {
        if (second_alt_alleles > 1) {
            set_gt(hdr, second, 0, 1);
        }
        make_record_multiallelic(hdr, first, second);
        records_to_remove.push_back(second);
    }

    for (size_t i = suppress_from; i < indel_group.size(); i++) {
        bcf1_t* record = indel_group[i];
        int alt_alleles = count_alt_alleles(hdr, record);
        if (alt_alleles == 0) {
            return;
        }

        set_gt(hdr, record, 0, 0);
        records_to_remove.push_back(record);
    }
}

void make_multiallelic(bcf_hdr_t* hdr, std::vector<bcf1_t*>& small_variants) {
    int curr_rid = -1;
    hts_pos_t curr_pos = -1;
    std::vector<bcf1_t*> indel_group;
    std::vector<bcf1_t*> records_to_remove;

    for (bcf1_t* record : small_variants) {
        if (!indel_group.empty() && (record->rid != curr_rid || record->pos != curr_pos)) {
            apply_multiallelic_logic_to_group(hdr, indel_group, records_to_remove);
            indel_group.clear();
        }

        if (indel_group.empty()) {
            curr_rid = record->rid;
            curr_pos = record->pos;
        }

        std::string svtype = get_sv_type(hdr, record);
        if (svtype == "INS" || svtype == "DEL") {
            indel_group.push_back(record);
        }
    }

    apply_multiallelic_logic_to_group(hdr, indel_group, records_to_remove);

    std::unordered_set<bcf1_t*> remove_set(records_to_remove.begin(), records_to_remove.end());
    for (bcf1_t*& record : small_variants) {
        if (remove_set.count(record) > 0) {
            record = nullptr;
        }
    }
    small_variants.erase(std::remove(small_variants.begin(), small_variants.end(), nullptr), small_variants.end());
}

int main(int argc, char* argv[]) {
    std::string in_vcf_fname = argv[1];
    std::string out_smvars_prefix = argv[2];
    std::string out_stvars_prefix = argv[3];
    std::string reference_fname = argv[4];
    std::string out_smvars_vcf_fname = out_smvars_prefix + ".vcf.gz";
    std::string out_stvars_vcf_fname = out_stvars_prefix + ".vcf.gz";
    std::string out_smvars_non_primary_vcf_fname = out_smvars_prefix + ".non-primary.vcf.gz";
    std::string out_stvars_non_primary_vcf_fname = out_stvars_prefix + ".non-primary.vcf.gz";

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	bcf_hdr_t* in_vcf_hdr = bcf_hdr_read(in_vcf_file);

    std::vector<bcf1_t*> vcf_records;
    std::unordered_set<bcf1_t*> aux_indel_records;
	bcf1_t* r = bcf_init();
	while (bcf_read(in_vcf_file, in_vcf_hdr, r) == 0) {
        bcf1_t* b = bcf_dup(r);
		vcf_records.push_back(b);

        char* s_data = NULL;
        int len = 0;
        bcf_get_info_string(in_vcf_hdr, b, "AUX_SNPS", (void**) &s_data, &len);
        if (len > 0) {
            bcf_unpack(b, BCF_UN_ALL);
            std::istringstream ss(s_data);
            std::string snp_str;
            int i = 0;
            while (std::getline(ss, snp_str, ',')) {
                snp_t snp(snp_str);
                std::string contig_name = bcf_seqname(in_vcf_hdr, b);
                char* chr_seq = chr_seqs.get_seq(contig_name);
                std::string id = std::string(b->d.id) + ".SNP." + std::to_string(i++);
                std::vector<int> gt = get_bcf_gt(in_vcf_hdr, b);
                bcf1_t* snp_record = generate_snp(in_vcf_hdr, bcf_seqname(in_vcf_hdr, b), snp.pos,
                    chr_seq[snp.pos], snp.alt_base, id, gt);
                copy_all_fmt(in_vcf_hdr, b, snp_record);
                vcf_records.push_back(snp_record);
            }
            // remove INFO/AUX_SNPS from the original record
            bcf_update_info_string(in_vcf_hdr, b, "AUX_SNPS", NULL);
        }
        free(s_data);

        s_data = NULL;
        len = 0;
        bcf_get_info_string(in_vcf_hdr, b, "AUX_INDELS", (void**) &s_data, &len);
        if (len > 0) {
            bcf_unpack(b, BCF_UN_ALL);
            std::istringstream ss(s_data);;
            std::string indel_str;
            int i = 0;
            while (std::getline(ss, indel_str, ',')) {
                std::stringstream indel_ss(indel_str);
                std::string start_str, end_str, ins_seq;
                std::getline(indel_ss, start_str, ':');
                std::getline(indel_ss, end_str, ':');
                std::getline(indel_ss, ins_seq, ':');
                std::string contig_name = bcf_seqname(in_vcf_hdr, b);
                hts_pos_t start = std::stoll(start_str)-1;
                hts_pos_t end = std::stoll(end_str)-1;
                std::shared_ptr<sv_t> indel;
                if (start == end) {
                    // insertion
                    indel = std::make_shared<insertion_t>(contig_name, start, end, ins_seq, nullptr, nullptr, nullptr, nullptr);
                } else {
                    // deletion
                    indel = std::make_shared<deletion_t>(contig_name, start, end, ins_seq, nullptr, nullptr, nullptr, nullptr);
                }
                indel->id = std::string(b->d.id) + ".INDEL." + std::to_string(i++);
                std::vector<int> gt = get_bcf_gt(in_vcf_hdr, b);
                indel->sample_info.gt = gt;
                bcf1_t* indel_record = bcf_init();
                sv2bcf(in_vcf_hdr, indel_record, indel.get(), chr_seqs.get_seq(contig_name));
                copy_all_fmt(in_vcf_hdr, b, indel_record);
                vcf_records.push_back(indel_record);
                aux_indel_records.insert(indel_record);
            }
            // remove INFO/AUX_INDELS from the original record
            bcf_update_info_string(in_vcf_hdr, b, "AUX_INDELS", NULL);
        }
        free(s_data);
    }
    
    // sort by pos and then by non-ascending EPR
    std::stable_sort(vcf_records.begin(), vcf_records.end(),
        [in_vcf_hdr](bcf1_t* b1, bcf1_t* b2) {
            if (std::tie(b1->rid, b1->pos) != std::tie(b2->rid, b2->pos)) {
                return std::tie(b1->rid, b1->pos) < std::tie(b2->rid, b2->pos);
            }
            return get_sv_epr(in_vcf_hdr, b1) > get_sv_epr(in_vcf_hdr, b2);
        });

    std::unordered_set<std::string> non_aux_variant_keys;
    for (int i = 0; i < vcf_records.size(); i++) {
        if (aux_indel_records.count(vcf_records[i]) || count_alt_alleles(in_vcf_hdr, vcf_records[i]) == 0) {
            continue;
        }
        non_aux_variant_keys.insert(bcf_record_unique_key(in_vcf_hdr, vcf_records[i]));
    }
    remove_aux_indels_matching_primary(in_vcf_hdr, vcf_records, aux_indel_records, non_aux_variant_keys);
    
    std::vector<bcf1_t*> small_variants, svs;
    divide_variants_by_svsize(in_vcf_hdr, vcf_records, small_variants, svs);

    write_and_remove_non_primary(in_vcf_hdr, small_variants, out_smvars_non_primary_vcf_fname);
    write_and_remove_non_primary(in_vcf_hdr, svs, out_stvars_non_primary_vcf_fname);

    make_multiallelic(in_vcf_hdr, small_variants);
    add_star_alleles_for_overlapping_del_ins(in_vcf_hdr, small_variants);

    htsFile* out_smvars_vcf_file = bcf_open(out_smvars_vcf_fname.c_str(), "wz");
    if (bcf_hdr_write(out_smvars_vcf_file, in_vcf_hdr) != 0) {
        throw std::runtime_error("Failed to write the VCF header to " + out_smvars_vcf_fname + ".");
    }
    for (bcf1_t* small_variant : small_variants) {
        if (bcf_write(out_smvars_vcf_file, in_vcf_hdr, small_variant) != 0) {
            throw std::runtime_error("Failed to write VCF record to " + out_smvars_vcf_fname);
        }
    }
    hts_close(out_smvars_vcf_file);

    htsFile* out_stvars_vcf_file = bcf_open(out_stvars_vcf_fname.c_str(), "wz");
    if (bcf_hdr_write(out_stvars_vcf_file, in_vcf_hdr) != 0) {
        throw std::runtime_error("Failed to write the VCF header to " + out_stvars_vcf_fname + ".");
    }
    for (bcf1_t* sv : svs) {
        if (bcf_write(out_stvars_vcf_file, in_vcf_hdr, sv) != 0) {
            throw std::runtime_error("Failed to write VCF record to " + out_stvars_vcf_fname);
        }
    }
    hts_close(out_stvars_vcf_file);

    for (bcf1_t* vcf_record : vcf_records) {
        bcf_destroy(vcf_record);
    }
    bcf_destroy(r);

    bcf_hdr_destroy(in_vcf_hdr);
    hts_close(in_vcf_file);
}
