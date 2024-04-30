#ifndef SVCOMPARE_COMMON_H
#define SVCOMPARE_COMMON_H

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

#include "../src/vcf_utils.h"
#include "../src/utils.h"

bool ends_with(const char* str, const char* suffix) {
	int l_string = strlen(str);
	int l_suffix = strlen(suffix);
	return strcmp(str+(l_string-l_suffix), suffix) == 0;
}

std::vector<sv_t*> read_sv_list(const char* filename) {
	std::vector<sv_t*> svs;
	if (ends_with(filename, ".vcf.gz") || ends_with(filename, ".vcf") || ends_with(filename, ".bcf")) { // vcf format
		htsFile* file = bcf_open(filename, "r");
		bcf_hdr_t* hdr = bcf_hdr_read(file);
		bcf1_t* line = bcf_init();
		while (bcf_read(file, hdr, line) == 0) {
            try {
                sv_t* sv = bcf_to_sv(hdr, line);
                sv->vcf_entry = bcf_dup(line);
                if (sv != NULL) {
                    svs.push_back(sv);
                }
            } catch (const std::exception& e) {
                std::cerr << e.what() << std::endl;
            }
		}
		bcf_destroy(line);
		bcf_hdr_destroy(hdr);
		bcf_close(file);
	}
    return svs;
}

struct repeat_t {
    std::string chr;
    int start, end;

    repeat_t() : start(0), end(0) {}
    repeat_t(std::string& line) {
        std::stringstream ss(line);
        ss >> chr >> start >> end;
    }

    bool contains(sv_t* sv) {
        return sv->chr == chr && start-10 <= sv->start && sv->end <= end+10;
    }
    bool intersects(sv_t* sv) {
        return sv->chr == chr && sv->start <= end && sv->end >= start;
    }
};
bool operator == (const repeat_t& r1, const repeat_t& r2) {
    return r1.chr == r2.chr && r1.start == r2.start && r1.end == r2.end;
}

#endif //SVCOMPARE_COMMON_H
