#ifndef SAM_UTILS_H
#define SAM_UTILS_H

#include <iostream>
#include <algorithm>
#include <mutex>

#include "htslib/sam.h"
#include "htslib/faidx.h"
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
    return bam_is_rev(r) == bam_is_mrev(r);
}
bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000 || is_unmapped(r) != is_mate_unmapped(r);
}
bool is_outward(bam1_t* r) {
	return is_samechr(r) && !is_samestr(r) && ((!bam_is_rev(r) && r->core.isize < 0) || (bam_is_rev(r) && r->core.isize > 0));
}
bool is_short(bam1_t* r, int min_is) {
    return is_samechr(r) && !is_samestr(r) && ((!bam_is_rev(r) && r->core.isize < min_is) || (bam_is_rev(r) && r->core.isize > -min_is));
}
bool is_long(bam1_t* r, int max_is) {
	return is_samechr(r) && !is_samestr(r) && ((!bam_is_rev(r) && r->core.isize > max_is) || (bam_is_rev(r) && r->core.isize < -max_is));
}
bool is_proper_pair(bam1_t* r, int min_is, int max_is) {
	return is_primary(r) && is_samechr(r) && !is_samestr(r) && !is_dc_pair(r) && !is_outward(r) && !is_short(r, min_is) && !is_long(r, max_is);
}

int get_endpoint(bam1_t* r) {
    return bam_is_rev(r) ? r->core.pos : bam_endpos(r);
}

int get_left_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]): 0;
}
int get_right_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' ? bam_cigar_oplen(cigar[r->core.n_cigar-1]): 0;
}
bool is_left_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    return get_left_clip_size(r) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    return get_right_clip_size(r) >= min_clip_len;
}

hts_pos_t get_unclipped_start(bam1_t* r) {
    return r->core.pos - get_left_clip_size(r);
}
hts_pos_t get_unclipped_end(bam1_t* r) {
    return bam_endpos(r) + get_right_clip_size(r);
}

int get_mate_left_clip_size(bam1_t* r) {
    if (is_mate_unmapped(r)) return 0;
    const uint8_t* mc_tag = bam_aux_get(r, "MC");
    if (!mc_tag) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
    }
    char* mc = bam_aux2Z(mc_tag);
    int i = 0, left_clip = 0;
    while (mc[i] >= '0' && mc[i] <= '9') {
        left_clip = left_clip * 10 + (mc[i] - '0');
        i++;
    }
    return (mc[i] == 'S') ? left_clip : 0;
}

int get_mate_right_clip_size(bam1_t* r) {
    if (is_mate_unmapped(r)) return 0;
    const uint8_t* mc_tag = bam_aux_get(r, "MC");
    if (!mc_tag) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
    }
    char* mc = bam_aux2Z(mc_tag);
    int len = strlen(mc), i = len - 1;
    if (mc[i] != 'S') return 0;
    i--;
    int right_clip = 0, place = 1;
    while (i >= 0 && mc[i] >= '0' && mc[i] <= '9') {
        right_clip += (mc[i] - '0') * place;
        place *= 10;
        i--;
    }
    return right_clip;
}

bool is_mate_left_clipped(bam1_t* r) {
    if (is_mate_unmapped(r)) return false;
    const uint8_t* mc_tag = bam_aux_get(r, "MC");
    if (mc_tag == NULL) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
    }
    char* mc_tag_str = bam_aux2Z(mc_tag);
    int i = 0;
    while (mc_tag_str[i] >= '0' && mc_tag_str[i] <= '9') i++;
    return mc_tag_str[i] == 'S';
}
bool is_mate_right_clipped(bam1_t* r) {
    if (is_mate_unmapped(r)) return false;
    const uint8_t* mc_tag = bam_aux_get(r, "MC");
    if (mc_tag == NULL) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
    }
    char* mc_tag_str = bam_aux2Z(mc_tag);
    int i = strlen(mc_tag_str)-1;
    return mc_tag_str[i] == 'S';
}
bool is_mate_clipped(bam1_t* r) {
	return is_mate_left_clipped(r) || is_mate_right_clipped(r);
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


hts_pos_t get_mate_endpos(const bam1_t* r) {
    uint8_t* mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0;
    hts_pos_t pos = r->core.mpos;
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

hts_pos_t get_mate_unclipped_start(bam1_t* r) {
    return r->core.mpos - get_mate_left_clip_size(r);
}
hts_pos_t get_mate_unclipped_end(bam1_t* r) {
    return get_mate_endpos(r) + get_mate_right_clip_size(r);
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

int64_t get_nm(bam1_t* r) {
    uint8_t* nm = bam_aux_get(r, "NM");
    if (nm == NULL) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the NM tag.";
    }
    return bam_aux2i(nm);
}

std::string get_md(bam1_t* r) {
    uint8_t* md = bam_aux_get(r, "MD");
    if (md == NULL) {
        return "";
    }
    return bam_aux2Z(md);
}

bool has_no_indels(bam1_t* read) {
    uint32_t *cigar = bam_get_cigar(read);
    return read->core.n_cigar == 1 && bam_cigar_op(cigar[0]) == BAM_CMATCH && bam_cigar_oplen(cigar[0]) == read->core.l_qseq;
}

bool is_perfectly_aligned(bam1_t* read) {
    // All M
    if (has_no_indels(read)) {
        uint8_t* nm = bam_aux_get(read, "NM");
        if (nm) {
            return bam_aux2i(nm) == 0; // Whether it is perfectly aligned based on NM
        } 

        std::string md = get_md(read);
        return md == std::to_string(read->core.l_qseq); // Perfectly aligned based on MD
    }
    return false; // Not perfectly aligned or cannot confirm alignment without NM or MD
}

void rc(std::string& read) {
    int len = read.length();
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
double avg_qual(std::string& qual_ascii, int offset = 33) {
	double tot = 0;
	for (char c : qual_ascii) tot += int(c)-offset;
	return tot/qual_ascii.length();
}

int get_aligned_portion_len(bam1_t* read) {
    return read->core.l_qseq - get_left_clip_size(read) - get_right_clip_size(read);
}

int64_t get_AS_tag(bam1_t* read) {
    uint8_t* aux_get = bam_aux_get(read, "AS");
    return bam_aux2i(aux_get);
}

void copy_sequence(bam1_t* r, char* seq, bool fastq_seq = false) { // assumes seq is long enough
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (fastq_seq && bam_is_rev(r)) rc(seq);
}

samFile* open_writer(std::string filename, bam_hdr_t* header) {
    samFile* writer = sam_open(filename.c_str(), "wb");
    if (writer == NULL) {
        throw "Unable to open " + filename;
    }
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + filename;
    }
    return writer;
}

void write_and_index_file(std::vector<bam1_t*>& reads, std::string path, bam_hdr_t* header) {
    samFile* file = open_writer(path, header);
    if (file == NULL) {
        throw "Unable to open " + path;
    }

    // write reads
    std::sort(reads.begin(), reads.end(), [](bam1_t *r1, bam1_t *r2) { return r1->core.pos < r2->core.pos; });
    for (bam1_t* r : reads) {
        int ok = sam_write1(file, header, r);
        if (ok < 0) throw "Unable to write to " + path;
    }

    sam_close(file);

    file = sam_open(path.c_str(), "r");

    int code = sam_index_build(path.c_str(), 0);
    if (code != 0) {
        throw "Cannot index " + path;
    }

    sam_close(file);
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
    if (f) {
        hts_idx_destroy(f->idx);
        bam_hdr_destroy(f->header);
        sam_close(f->file);
        delete f;
    }
}

struct bam_pool_t {
    std::mutex mtx;
    std::vector<open_samFile_t*> pool;
    std::string bam_fname, reference_fname;

    bam_pool_t(int size, std::string bam_fname, std::string reference_fname) : bam_fname(bam_fname), reference_fname(reference_fname) {
        for (int i = 0; i < size; i++) {
            open_samFile_t* o = open_samFile(bam_fname.c_str());
            hts_set_fai_filename(o->file, fai_path(reference_fname.c_str()));
            pool.push_back(o);
        }
    }

    open_samFile_t* get_bam_reader(int i) {
        return pool[i];
    }

    ~bam_pool_t() {
        for (open_samFile_t* o : pool) {
            close_samFile(o);
        }
    }
};


#endif
