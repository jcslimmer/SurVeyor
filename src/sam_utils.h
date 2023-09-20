#ifndef SAM_UTILS_H
#define SAM_UTILS_H

#include "htslib/sam.h"
#include "utils.h"

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}
bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}
bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid && !is_unmapped(r) && !is_mate_unmapped(r);
}
bool is_samestr(bam1_t* r) {
    return is_samechr(r) && (bam_is_rev(r) == bam_is_mrev(r));
}
bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000 || is_unmapped(r) != is_mate_unmapped(r);
}

bool is_left_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' && bam_cigar_oplen(cigar[0]) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' && bam_cigar_oplen(cigar[r->core.n_cigar-1]) >= min_clip_len;
}

bool is_hidden_split_read(bam1_t* r, config_t config) {
	if (is_left_clipped(r, config.min_clip_len) || is_right_clipped(r, config.min_clip_len)) return false;

    int mismatches = bam_aux2i(bam_aux_get(r, "NM"));

    int indels = 0;
    uint32_t* cigar = bam_get_cigar(r);
    for (uint32_t i = 0; i < r->core.n_cigar; i++) {
        char op_chr = bam_cigar_opchr(cigar[i]);
        if (op_chr == 'D' || op_chr == 'I') {
            mismatches -= bam_cigar_oplen(cigar[i]);
            indels++;
        }
    }

    return (mismatches + indels >= config.min_diff_hsr);
}


int get_mate_endpos(const bam1_t* r) {
    uint8_t* mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0, pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}

int64_t get_mq(bam1_t* r) {
    uint8_t* mq = bam_aux_get(r, "MQ");
    if (mq == NULL) {
    	if ((r->core.flag & BAM_FUNMAP) == 0 && (r->core.flag & BAM_FMUNMAP) == 0) {
			std::cerr << "Warning: read pair " << bam_get_qname(r) << " does not have an MQ tag. Please include it." << std::endl;
			std::cerr << (r->core.flag & BAM_FMUNMAP) << std::endl;
    	}
    	return 0;
    }
    return bam_aux2i(mq);
}

void rc(char* read) {
    int len = strlen(read);
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
    	char c = std::toupper(read[i]);
        if (c == 'A') read[i] = 'T';
        else if (c == 'C') read[i] = 'G';
        else if (c == 'G') read[i] = 'C';
        else if (c == 'T') read[i] = 'A';
        else c = 'N';
    }
}
char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}
std::string get_sequence(bam1_t* r, bool fastq_seq = false) {
    char seq[100000];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (fastq_seq && bam_is_rev(r)) rc(seq);
    return std::string(seq);
}
std::string get_qual_ascii(bam1_t* r, bool fastq_seq = false) {
	uint8_t* qual = bam_get_qual(r);
	std::string qual_ascii(r->core.l_qseq, ' ');
	for (int i = 0; i < r->core.l_qseq; i++) {
		qual_ascii[i] = char(33 + qual[i]);
	}
	if (fastq_seq && bam_is_rev(r)) qual_ascii = std::string(qual_ascii.rbegin(), qual_ascii.rend());
	return qual_ascii;
}

double avg_qual(bam1_t* read) {
	double avg_qual = 0;
	uint8_t* qual = bam_get_qual(read);
	for (int i = 0; i < read->core.l_qseq; i++) {
		avg_qual += qual[i];
	}
	return avg_qual/read->core.l_qseq;
}

samFile* open_writer(std::string filename, bam_hdr_t* header) {
    samFile* writer = sam_open(filename.c_str(), "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + filename;
    }
    return writer;
}

struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() {}

    open_samFile_t(samFile* file, bam_hdr_t* header, hts_idx_t* idx) : file(file), header(header), idx(idx) {}
};

open_samFile_t* open_samFile(std::string fname_str, bool index_file = false) {
    const char* fname = fname_str.c_str();
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
    if (sam_file->idx == NULL) {
        throw "Unable to open index for " + std::string(fname);
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
    hts_idx_destroy(f->idx);
    bam_hdr_destroy(f->header);
    sam_close(f->file);
    delete f;
}

#endif
