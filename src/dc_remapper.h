#ifndef DC_REMAPPER_H_
#define DC_REMAPPER_H_

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "clustering_utils.h"

struct insertion_cluster_t {

    hts_pos_t start, end;
    cluster_t* cluster = NULL;
    std::vector<bam1_t*> semi_mapped_reads;
    consensus_t* clip_consensus = NULL;

    insertion_cluster_t(cluster_t* cluster) : cluster(cluster), start(cluster->la_start), end(cluster->la_end) {
        std::vector<bam1_t*> reads;
        reads.swap(cluster->reads);
        for (bam1_t* read : reads) {
            add_stable_read(read);
        }
    }

    void add_stable_read(bam1_t* read) {
        cluster->reads.push_back(read);
        start = std::min(start, read->core.pos);
        end = std::max(end, bam_endpos(read));
        cluster->la_max_mapq = std::max(cluster->la_max_mapq, read->core.qual);
        cluster->ra_max_mapq = std::max(cluster->ra_max_mapq, read->core.qual);
    }
    void add_semi_mapped_reads(bam1_t* read) {
        if (clip_consensus) {
            bool sm_left_clipped = is_left_clipped(read, 0);
            bool sm_right_clipped = is_right_clipped(read, 0);
            if (sm_left_clipped && clip_consensus->left_clipped && read->core.pos == clip_consensus->breakpoint) {
                return;
            } else if (sm_right_clipped && !clip_consensus->left_clipped && bam_endpos(read) == clip_consensus->breakpoint) {
                return;
            }
        }
    	semi_mapped_reads.push_back(read);
        start = std::min(start, read->core.pos);
        end = std::max(end, bam_endpos(read));
    }
    void add_clip_cluster(consensus_t* clip_consensus) {
        if (!clip_consensus) return;
        if (!this->clip_consensus || this->clip_consensus->reads() < clip_consensus->reads()) {
            this->clip_consensus = clip_consensus;
            if (clip_consensus->left_clipped) start = clip_consensus->anchor_start();
            else start = std::min(start, clip_consensus->anchor_start());
            if (!clip_consensus->left_clipped) end = clip_consensus->anchor_end();
            else end = std::max(end, clip_consensus->anchor_end());
        }
    }

    void deallocate_reads() {
        for (bam1_t* read : cluster->reads) {
            bam_destroy1(read);
        }
        cluster->reads.clear();
        for (bam1_t* read : semi_mapped_reads) {
        	bam_destroy1(read);
        }
        semi_mapped_reads.clear();
    }

    int size() {
        return cluster->reads.size() + semi_mapped_reads.size() + (clip_consensus ? 1 : 0);
    }

    bool empty() {
    	return this->size() == 0;
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
