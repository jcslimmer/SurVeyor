#ifndef DC_REMAPPER_H_
#define DC_REMAPPER_H_

#include "cluster.h"

struct reads_cluster_t {
    std::vector<bam1_t*> reads;
    std::vector<bam1_t*> semi_mapped_reads;
    consensus_t* clip_consensus = NULL;
    bool used = false;

    void add_read(bam1_t* read) {
        reads.push_back(read);
    }
    void add_semi_mapped_reads(bam1_t* read) {
    	semi_mapped_reads.push_back(read);
    }
    void add_clip_cluster(consensus_t* clip_consensus) {
        if (!clip_consensus) return;
        if (!this->clip_consensus || this->clip_consensus->supp_clipped_reads() < clip_consensus->supp_clipped_reads()) {
            this->clip_consensus = clip_consensus;
        }
    }

    hts_pos_t start() {
        if (clip_consensus && clip_consensus->left_clipped) {
            return clip_consensus->anchor_start();
        } else {
            hts_pos_t _start = INT32_MAX;
            for (bam1_t* read : reads) {
                _start = std::min(_start, read->core.pos);
            }
            for (bam1_t* read : semi_mapped_reads) {
				_start = std::min(_start, read->core.pos);
			}
            if (clip_consensus && clip_consensus->anchor_start() < _start) {
                _start = clip_consensus->anchor_start();
            }
            return _start;
        }
    }
    hts_pos_t end() {
        if (clip_consensus && !clip_consensus->left_clipped) {
            return clip_consensus->anchor_end();
        } else {
            hts_pos_t _end = 0;
            for (bam1_t* read : reads) {
                _end = std::max(_end, bam_endpos(read));
            }
            for (bam1_t* read : semi_mapped_reads) {
				_end = std::max(_end, bam_endpos(read));
			}
            if (clip_consensus && clip_consensus->anchor_end() > _end) {
                _end = clip_consensus->anchor_end();
            }
            return _end;
        }
    }

    void deallocate_reads() {
        for (bam1_t* read : reads) {
            bam_destroy1(read);
        }
        for (bam1_t* read : semi_mapped_reads) {
        	bam_destroy1(read);
        }
        reads.clear();
        semi_mapped_reads.clear();
    }

    bool empty() {
    	return reads.empty() && !clip_consensus;
    }
};

std::string get_mate_seq(bam1_t* read, std::unordered_map<std::string, std::string>& mateseqs) {
    std::string qname = bam_get_qname(read);
    if (is_samechr(read)) {
        if (read->core.flag & BAM_FREAD1) qname += "_2";
        else qname += "_1";
    }
    if (!mateseqs.count(qname)) {
    	std::cerr << "Warning: mate not found for " << qname << std::endl;
    	return "";
    }
    return mateseqs[qname];
}
std::string get_mate_qual(bam1_t* read, std::unordered_map<std::string, std::string>& matequals) {
    std::string qname = bam_get_qname(read);
    if (is_samechr(read)) {
        if (read->core.flag & BAM_FREAD1) qname += "_2";
        else qname += "_1";
    }
    if (!matequals.count(qname)) {
    	std::cerr << "Warning: mate not found for " << qname << std::endl;
    	return "";
    }
    return matequals[qname];
}


#endif /* DC_REMAPPER_H_ */
