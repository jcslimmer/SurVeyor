#include <iostream>
#include <unordered_set>

#include "htslib/vcf.h"

int main(int argc, char* argv[]) {

    std::string in_vcf_fname = argv[1];

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (in_vcf_file == NULL) {
        std::cerr << "Error: Unable to open VCF file " << in_vcf_fname << std::endl;
        return 1;
    }
	
    bcf_hdr_t* hdr = bcf_hdr_read(in_vcf_file);
    if (hdr == NULL) {
        std::cerr << "Error: Unable to read header of VCF file " << in_vcf_fname << std::endl;
        return 1;
    }

    std::unordered_set<std::string> ids;

	bcf1_t* b = bcf_init();
    while (bcf_read(in_vcf_file, hdr, b) == 0) {
        bcf_unpack(b, BCF_UN_STR);
        std::string id = b->d.id;
        if (ids.count(id) > 0) {
            std::cerr << "Error: Duplicate ID " << id << std::endl;
            return 1;
        }
        ids.insert(id);
    }

    bcf_destroy(b);
    bcf_hdr_destroy(hdr);
    bcf_close(in_vcf_file);

    return 0;
}