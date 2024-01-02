#ifndef CLUSTERING_UTILS_H
#define CLUSTERING_UTILS_H

#include "htslib/sam.h"
#include "sam_utils.h"
#include <cstdint>

// IMPORTANT: All of these functions do not consider contig and direction, because they assume that we are trying to cluster pairs 
// from the same contig and compatible direction

struct cluster_t {
	int id = -1;
	hts_pos_t la_start, la_end, ra_start, ra_end;
	int count = 1, confident_count = 0;
	uint8_t max_mapq;
	bool used = false;
	std::string rightmost_lseq, leftmost_rseq; // right-most right-facing sequence and left-most left-facing seq
	
	std::vector<bam1_t*> reads; // reads that are part of this cluster

	cluster_t() : la_start(INT32_MAX), la_end(0), ra_start(INT32_MAX), ra_end(0), max_mapq(0) {}

	cluster_t(bam1_t* read, int high_confidence_mapq, bool store_read = false) : la_start(read->core.pos), la_end(bam_endpos(read)), ra_start(read->core.mpos), ra_end(get_mate_endpos(read)),
			max_mapq(std::min(read->core.qual, (uint8_t) get_mq(read))) {
		confident_count = max_mapq == high_confidence_mapq;
		rightmost_lseq = get_sequence(read);
		leftmost_rseq = bam_get_qname(read);
		if (store_read) reads.push_back(bam_dup1(read));
	}

	cluster_t(hts_pos_t la_start, hts_pos_t la_end, hts_pos_t ra_start, hts_pos_t ra_end, int count) : 
		la_start(la_start), la_end(la_end), ra_start(ra_start), ra_end(ra_end), count(count), max_mapq(0) {}

	static hts_pos_t distance(cluster_t* c1, cluster_t* c2) {
		hts_pos_t la_dist = std::max(c1->la_end, c2->la_end) - std::min(c1->la_start, c2->la_start);
		hts_pos_t ra_dist = std::max(c1->ra_end, c2->ra_end) - std::min(c1->ra_start, c2->ra_start);
		return std::max(la_dist, ra_dist);
	}

	static bool can_merge(cluster_t* c1, cluster_t* c2, int max_distance) {
		return distance(c1, c2) <= max_distance;
	}

	static cluster_t* merge(cluster_t* c1, cluster_t* c2) {
		cluster_t* merged = new cluster_t;
		merged->la_start = std::min(c1->la_start, c2->la_start);
		merged->la_end = std::max(c1->la_end, c2->la_end);
		merged->ra_start = std::min(c1->ra_start, c2->ra_start);
		merged->ra_end = std::max(c1->ra_end, c2->ra_end);
		merged->count = c1->count + c2->count;
		merged->confident_count = c1->confident_count + c2->confident_count;
		merged->max_mapq = std::max(c1->max_mapq, c2->max_mapq);

		// set rightmost_lseq
		if (c1->la_end > c2->la_end) merged->rightmost_lseq = c1->rightmost_lseq;
		else if (c1->la_end < c2->la_end) merged->rightmost_lseq = c2->rightmost_lseq;
		else if (c1->la_start >= c2->la_start) merged->rightmost_lseq = c1->rightmost_lseq; // both reads may be clipped at the same position, then we choose the one that starts later
		else merged->rightmost_lseq = c2->rightmost_lseq;

		// set leftmost_lseq
		if (c1->ra_start < c2->ra_start) merged->leftmost_rseq = c1->leftmost_rseq;
		else if (c1->ra_start > c2->ra_start) merged->leftmost_rseq = c2->leftmost_rseq;
		else if (c1->la_end < c2->la_end) merged->leftmost_rseq = c1->leftmost_rseq;
		else merged->leftmost_rseq = c2->leftmost_rseq;

		merged->reads.insert(merged->reads.end(), c1->reads.begin(), c1->reads.end());
		merged->reads.insert(merged->reads.end(), c2->reads.begin(), c2->reads.end());
		return merged;
	}

	~cluster_t() {
		for (bam1_t* read : reads) {
			bam_destroy1(read);
		}
	}
};

bool operator < (const cluster_t& c1, const cluster_t& c2) {
	return std::make_tuple(c1.la_start, c1.la_end, c1.ra_start, c1.ra_end) < std::make_tuple(c2.la_start, c2.la_end, c2.ra_start, c2.ra_end);
}
bool operator == (const cluster_t& c1, const cluster_t& c2) {
	return c1.la_start == c2.la_start && c1.la_end == c2.la_end && c1.ra_start == c2.ra_start && c1.ra_end == c2.ra_end;
}
bool operator != (const cluster_t& c1, const cluster_t& c2) {
	return !(c1 == c2);
}

struct cc_distance_t {
    int distance;
    cluster_t* c1,* c2;

    cc_distance_t(cluster_t *c1, cluster_t *c2) : distance(cluster_t::distance(c1, c2)), c1(c1), c2(c2) {}
};
bool operator < (const cc_distance_t& ccd1, const cc_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance > ccd2.distance;
}

#endif // CLUSTERING_UTILS_H
