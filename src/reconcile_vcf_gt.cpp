#include "htslib/hts.h"
#include "vcf_utils.h"
#include <htslib/vcf.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>

struct gt_record_t {
    bcf1_t* record;
    float epr;
};

// Function to get FORMAT/EPR from a VCF record
float get_epr(bcf_hdr_t *hdr, bcf1_t *b) {
    int ngt = 0;
    float *epr = NULL;
    if (bcf_get_format_float(hdr, b, "EPR", &epr, &ngt) > 0) {
        float value = epr[0];
        free(epr);
        return value;
    }
    free(epr);
    return -1.0;
}

// Load gt_vcf records, retaining only the one with the highest EPR for each unique ID
std::unordered_map<std::string, gt_record_t> load_gt_vcf(std::string gt_vcf_fname) {
    htsFile* fp = bcf_open(gt_vcf_fname.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening " << gt_vcf_fname << std::endl;
        exit(1);
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);

    bcf1_t* b = bcf_init();
    std::unordered_map<std::string, gt_record_t> gt_records;

    while (bcf_read(fp, hdr, b) == 0) {
        bcf_unpack(b, BCF_UN_STR);
        std::string id = b->d.id;

        float epr = get_epr(hdr, b);
        auto it = gt_records.find(id);

        if (it == gt_records.end() || epr > it->second.epr) {
            if (it != gt_records.end()) {
                bcf_destroy(it->second.record);
            }
            gt_records[id] = { bcf_dup(b), epr };
        }
    }

    bcf_destroy(b);
    return gt_records;
}

// Copy FORMAT/GT from gt_vcf to base_vcf for matching IDs
void copy_all_format_to_base(bcf_hdr_t* base_hdr, htsFile* base_fp, const std::unordered_map<std::string, gt_record_t>& gt_records, bcf_hdr_t* gt_hdr, htsFile* out_fp, bcf_hdr_t* out_hdr) {
    bcf1_t* b = bcf_init();

    while (bcf_read(base_fp, base_hdr, b) == 0) {
        bcf_unpack(b, BCF_UN_STR);
        std::string id = b->d.id;

        auto it = gt_records.find(id);
        if (it != gt_records.end()) {
            bcf1_t* gt_record = it->second.record;
            bcf_unpack(gt_record, BCF_UN_ALL);

            // if SVTYPE of gt_record is DUP continue
            if (get_sv_type(gt_hdr, gt_record) == "DUP" && get_sv_type(base_hdr, b) == "INS") {
                std::stringstream dup_str; 
                dup_str << bcf_seqname(gt_hdr, gt_record) << ":" << std::to_string(gt_record->pos) << "-" << get_sv_end(gt_hdr, gt_record);
                bcf_update_info_string(out_hdr, b, "GT_AS_DUP", dup_str.str().c_str());
            }

            int32_t* gt = NULL, ngt = 0;
            int gt_ret = bcf_get_genotypes(gt_hdr, it->second.record, &gt, &ngt);

            if (gt_ret > 0) {
                bcf_update_genotypes(out_hdr, b, gt, ngt);
            }
            free(gt);

            bcf_fmt_t* fmt = gt_record->d.fmt;
            for (int i = 0; i < gt_record->n_fmt; i++) {
                const char* key = bcf_hdr_int2id(gt_hdr, BCF_DT_ID, fmt[i].id);
                int type = bcf_hdr_id2type(gt_hdr, BCF_HL_FMT, fmt[i].id);

                if (std::string(key) == "GT") {
                    continue;
                }

                if (type == BCF_HT_INT) {
                    int32_t* values = NULL;
                    int nvalues = 0;
                    if (bcf_get_format_int32(gt_hdr, gt_record, key, &values, &nvalues) > 0) {
                        bcf_update_format_int32(out_hdr, b, key, values, nvalues);
                    }
                    free(values);
                } else if (type == BCF_HT_REAL) {
                    float* values = NULL;
                    int nvalues = 0;
                    if (bcf_get_format_float(gt_hdr, gt_record, key, &values, &nvalues) > 0) {
                        bcf_update_format_float(out_hdr, b, key, values, nvalues);
                    }
                    free(values);
                } else if (type == BCF_HT_STR) {
                    char** values = NULL;
                    int nvalues = 0;
                    if (bcf_get_format_string(gt_hdr, gt_record, key, &values, &nvalues) > 0) {
                        bcf_update_format_string(out_hdr, b, key, (const char**)values, nvalues);
                    }
                    for (int j = 0; j < nvalues; ++j) {
                        free(values[j]);
                    }
                    free(values);
                }
            }

        }

        if (bcf_write(out_fp, out_hdr, b) < 0) {
            std::cerr << "Error writing record to base VCF" << std::endl;
            exit(1);
        }
    }

    bcf_destroy(b);
    bcf_close(out_fp);
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <base_vcf> <gt_vcf> <out_vcf> <sample>" << std::endl;
        return 1;
    }

    std::string base_vcf_fname = argv[1];
    std::string gt_vcf_fname = argv[2];
    std::string out_vcf_fname = argv[3];
    std::string sample = argv[4];

    // Load gt_vcf records, keeping the one with the highest EPR per ID
    auto gt_records = load_gt_vcf(gt_vcf_fname);

    htsFile* base_fp = bcf_open(base_vcf_fname.c_str(), "r");
    if (!base_fp) {
        std::cerr << "Error opening base VCF file: " << base_vcf_fname << std::endl;
        return 1;
    }
    bcf_hdr_t* base_hdr = bcf_hdr_read(base_fp);

    // reset samples in out header
    bcf_hdr_t* out_hdr = bcf_hdr_dup(base_hdr);
    add_genotyping_tags(out_hdr);
    bcf_hdr_remove(out_hdr, BCF_HL_INFO, "GT_AS_DUP");
    int len = 0;
    const char* gt_as_dup_tag = "##INFO=<ID=GT_AS_DUP,Number=1,Type=String,Description=\"This insertions was genotyped as the duplication provided in this field.\">";
    bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, gt_as_dup_tag, &len));
    if (bcf_hdr_set_samples(out_hdr, NULL, 0) != 0) {
        throw std::runtime_error("Failed to unset samples in VCF header");
    }
    bcf_hdr_add_sample(out_hdr, sample.c_str());
    if (bcf_hdr_sync(out_hdr) != 0) {
        throw std::runtime_error("Failed to sync header");
    }

    htsFile* out_fp = bcf_open(out_vcf_fname.c_str(), "w");
    if (!out_fp) {
        std::cerr << "Error opening " << out_vcf_fname << std::endl;
        exit(1);
    }
    if (bcf_hdr_write(out_fp, out_hdr) < 0) {
        std::cerr << "Error writing header to output VCF" << std::endl;
        exit(1);
    }

    // Write updated records to base_vcf
    htsFile* gt_fp = bcf_open(gt_vcf_fname.c_str(), "r");
    if (!gt_fp) {
        std::cerr << "Error opening " << gt_vcf_fname << std::endl;
        exit(1);
    }
    bcf_hdr_t* gt_hdr = bcf_hdr_read(gt_fp);
    copy_all_format_to_base(base_hdr, base_fp, gt_records, gt_hdr, out_fp, out_hdr);

    // Cleanup
    for (auto &pair : gt_records) {
        bcf_destroy(pair.second.record);
    }
    bcf_hdr_destroy(base_hdr);
    bcf_close(base_fp);
}
