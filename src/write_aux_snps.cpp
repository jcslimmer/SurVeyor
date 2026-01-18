#include <string>
#include <vector>
#include <sstream>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "types.h"
#include "vcf_utils.h"

int main(int argc, char* argv[]) {
    std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
    std::string reference_fname = argv[3];

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	bcf_hdr_t* in_vcf_hdr = bcf_hdr_read(in_vcf_file);

    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, in_vcf_hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

    std::vector<bcf1_t*> vcf_records;
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
                vcf_records.push_back(indel_record);
            }
            // remove INFO/AUX_INDELS from the original record
            bcf_update_info_string(in_vcf_hdr, b, "AUX_INDELS", NULL);
        }
        free(s_data);
    }

    std::sort(vcf_records.begin(), vcf_records.end(),
        [](const bcf1_t* b1, const bcf1_t* b2) { return std::tie(b1->rid, b1->pos) < std::tie(b2->rid, b2->pos); });

    for (int i = 0; i < vcf_records.size(); i++) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(in_vcf_hdr, vcf_records[i]);
        for (int j = i+1; j < vcf_records.size(); j++) {
            if (sv->start != vcf_records[j]->pos) break;

            std::shared_ptr<sv_t> sv2 = bcf_to_sv(in_vcf_hdr, vcf_records[j]);
            if (sv->chr == sv2->chr && sv->end == sv2->end && sv->ins_seq == sv2->ins_seq) {
                // if both 1/1
                std::vector<int> gt1 = get_bcf_gt(in_vcf_hdr, vcf_records[i]);
                std::vector<int> gt2 = get_bcf_gt(in_vcf_hdr, vcf_records[j]);
                if (gt1.size() == 2 && gt2.size() == 2 && (bcf_gt_allele(gt1[0])+bcf_gt_allele(gt1[1]) >= 1 && bcf_gt_allele(gt2[0])+bcf_gt_allele(gt2[1]) >= 1)) {
                    // set vcf_records[j] to ./.
                    int missing_gt[2] = { bcf_gt_missing, bcf_gt_missing };
                    bcf_update_genotypes(in_vcf_hdr, vcf_records[j], missing_gt, 2);

                    // set vcf_records[i] to 1/1
                    int hom_gt[2] = { bcf_gt_unphased(1), bcf_gt_unphased(1) };
                    bcf_update_genotypes(in_vcf_hdr, vcf_records[i], hom_gt, 2);
                }
            }
        }
    }

    for (bcf1_t* vcf_record : vcf_records) {
        if (bcf_write(out_vcf_file, in_vcf_hdr, vcf_record) != 0) {
            throw std::runtime_error("Failed to write VCF record to " +  out_vcf_fname);
        }
        bcf_destroy(vcf_record);
    }
    bcf_destroy(r);

    bcf_hdr_destroy(in_vcf_hdr);
    hts_close(in_vcf_file);
    hts_close(out_vcf_file);
}
