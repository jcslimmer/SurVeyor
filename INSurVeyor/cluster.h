#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include <iostream>
#include <atomic>

#include "../src/sam_utils.h"
#include "utils.h"
#include "htslib/sam.h"
#include "../src/clustering_utils.h"

struct anchor_t {

    char dir;
    int contig_id;
    int start, end;
    int fwd_sc_reads, rev_sc_reads;

    anchor_t() {}
    anchor_t(char dir, int contig_id, int start, int end, int fwd_sc_reads, int rev_sc_reads) :
    	dir(dir), contig_id(contig_id), start(start), end(end),
		fwd_sc_reads(fwd_sc_reads), rev_sc_reads(rev_sc_reads) {}

    int pos() { return dir == 'L' ? start : end; }

    int sc_reads() { return fwd_sc_reads + rev_sc_reads; }

    static bool check_clipped(anchor_t& clipped, anchor_t& other, config_t& config) {
        if (clipped.dir == 'L') {
            return clipped.pos() <= other.pos()+config.max_clipped_pos_dist;
        } else {
            return clipped.pos() >= other.pos()-config.max_clipped_pos_dist;
        }
    }

    static bool can_merge(anchor_t& a1, anchor_t& a2, config_t& config, stats_t& stats) {
        if (a1.sc_reads() > 0 && !check_clipped(a1, a2, config)) return false;
        if (a2.sc_reads() > 0 && !check_clipped(a2, a1, config)) return false;
        return inss_distance(a1, a2) <= stats.max_is;
    }

    static anchor_t merge(anchor_t& a1, anchor_t& a2) {
        return anchor_t(a1.dir, a1.contig_id, std::min(a1.start, a2.start), std::max(a1.end,  a2.end),
        		a1.fwd_sc_reads+a2.fwd_sc_reads, a1.rev_sc_reads+a2.rev_sc_reads);
    }

    static int inss_distance(anchor_t& a1, anchor_t& a2) {
        if (a1.contig_id != a2.contig_id || a1.dir != a2.dir) return INT_MAX;
        return std::max(a1.end, a2.end) - std::min(a1.start, a2.start);
    }

    int size() {
        return end-start+1;
    }
};
bool operator < (const anchor_t& a1, const anchor_t& a2) {
    if (a1.start != a2.start) return a1.start < a2.start;
    return a1.end < a2.end;
}

struct inss_cluster_t {

	int id = 0; // not necessarily set, use if needed
    anchor_t a1, a2;
    int disc_pairs;
    bool dead = false;

    cluster_t* cluster;

    inss_cluster_t(const anchor_t &a1, const anchor_t &a2, int disc_pairs) : a1(a1), a2(a2), disc_pairs(disc_pairs) {
        to_cluster();
    }

    inss_cluster_t(inss_cluster_t* c) {
        this->id = c->id;
        this->a1 = c->a1;
        this->a2 = c->a2;
        this->disc_pairs = c->disc_pairs;
        this->dead = c->dead;
        this->cluster = c->cluster;
    }

    void to_cluster() {
        cluster = new cluster_t();
        cluster->la_start = a1.start;
        cluster->la_end = a1.end;
        cluster->ra_start = a2.start;
        cluster->ra_end = a2.end;
        cluster->count = disc_pairs;
        cluster->confident_count = 0;
        cluster->max_mapq = 0;
        cluster->used = false;
        cluster->rightmost_lseq = "";
        cluster->leftmost_rseq = "";
    }

    static bool can_merge(inss_cluster_t* c1, inss_cluster_t* c2, config_t& config, stats_t& stats) {
        return anchor_t::can_merge(c1->a1, c2->a1, config, stats) && anchor_t::can_merge(c1->a2, c2->a2, config, stats);
    }
    static int inss_distance(inss_cluster_t* c1, inss_cluster_t* c2) {
        if (c1->a1.dir != c2->a1.dir || c1->a2.dir != c2->a2.dir) {
            return INT32_MAX;
        }
        return distance(c1->cluster, c2->cluster);
    }

    static inss_cluster_t* inss_merge(inss_cluster_t* c1, inss_cluster_t* c2) {
        cluster_t* merged = merge(c1->cluster, c2->cluster);
        // return new inss_cluster_t(merged);
        return new inss_cluster_t(anchor_t::merge(c1->a1, c2->a1), anchor_t::merge(c1->a2, c2->a2),
                           c1->disc_pairs+c2->disc_pairs);
    }
};


struct cc_distance_t {
    int distance;
    inss_cluster_t* c1,* c2;

    cc_distance_t(int distance, inss_cluster_t *c1, inss_cluster_t *c2) : distance(distance), c1(c1), c2(c2) {}
};
bool operator < (const cc_distance_t& ccd1, const cc_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance > ccd2.distance;
}


#endif //SURVEYOR_CLUSTER_H
