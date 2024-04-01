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


struct general_sv_t {
    std::string id, chr;
    hts_pos_t start, end;
    std::string type, ins_seq;
    int* gt = NULL, ngt;
    bool precise;

    general_sv_t(bcf_hdr_t* hdr, bcf1_t* vcf_sv) {
    	bcf_unpack(vcf_sv, BCF_UN_INFO);
        set_vars(vcf_sv->d.id, bcf_seqname(hdr, vcf_sv), vcf_sv->pos, get_sv_end(hdr, vcf_sv), get_sv_type(hdr, vcf_sv), get_ins_seq(hdr, vcf_sv));
        int n = 0;
        ngt = bcf_get_genotypes(hdr, vcf_sv, &gt, &n);
        std::sort(gt, gt + ngt);
        precise = !(bcf_get_info_flag(hdr, vcf_sv, "IMPRECISE", NULL, 0) == 1);
    }
    general_sv_t(std::string id, std::string chr, hts_pos_t start, hts_pos_t end, std::string type, std::string ins_seq) {
        set_vars(id, chr, start, end, type, ins_seq);
    }

    void set_vars(std::string id, std::string chr, hts_pos_t start, hts_pos_t end, std::string type, std::string ins_seq) {
        this->id = id;
        this->chr = chr;
        this->start = start;
        this->end = end;
        this->type = type;
        this->ins_seq = ins_seq;
    }

    bool incomplete_assembly() {
        return type == "INS" && ins_seq.find('-') != std::string::npos;
    }

    hts_pos_t len() {
        if (type == "DEL") return start - end + ins_seq.length();
        else if (type == "DUP") return end - start + ins_seq.length();
        else if (type == "INS") return ins_seq.length() - (end - start);
        return 0;
    }

    std::string print_gt() {
        std::stringstream ss;
        for (int i = 0; i < ngt; i++) {
            if (i > 0) ss << "/";
            ss << (bcf_gt_is_missing(gt[i]) ? "." : std::to_string(bcf_gt_allele(gt[i])));
        }
        return ss.str();
    
    }
};
bool operator < (const general_sv_t& sv1, const general_sv_t& sv2) {
    if (sv1.chr != sv2.chr) return sv1.chr < sv2.chr;
    return std::tie(sv1.start, sv1.end) < std::tie(sv2.start, sv2.end);
}

std::vector<general_sv_t> read_sv_list(const char* filename) {
	std::vector<general_sv_t> svs;
	if (ends_with(filename, ".vcf.gz") || ends_with(filename, ".vcf") || ends_with(filename, ".bcf")) { // vcf format
		htsFile* file = bcf_open(filename, "r");
		bcf_hdr_t* hdr = bcf_hdr_read(file);
		bcf1_t* line = bcf_init();
		while (bcf_read(file, hdr, line) == 0) {
			svs.push_back(general_sv_t(hdr, line));
		}
		bcf_destroy(line);
		bcf_hdr_destroy(hdr);
		bcf_close(file);
	} else {
		std::ifstream file(filename);
		std::string line;
		while (getline(file, line)) {
			std::stringstream ss(line);
			char strand;
            std::string id, chr, type, ins_seq;
            hts_pos_t start, end;
			ss >> id >> chr >> start >> strand >> chr >> end >> strand >> type >> ins_seq;
			general_sv_t sv(id, chr, start, end, type, ins_seq);
			svs.push_back(sv);
		}
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

    bool contains(general_sv_t& sv) {
        return sv.chr == chr && start <= sv.start && end >= sv.end;
    }
    bool intersects(general_sv_t& sv) {
        return sv.chr == chr && sv.start <= end && sv.end >= start;
    }
};
bool operator == (const repeat_t& r1, const repeat_t& r2) {
    return r1.chr == r2.chr && r1.start == r2.start && r1.end == r2.end;
}

#endif //SVCOMPARE_COMMON_H
