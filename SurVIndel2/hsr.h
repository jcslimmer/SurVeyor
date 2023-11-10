#ifndef HSR_H_
#define HSR_H_

#include "utils.h"

const int EXTRA_SEQ = 10;

// remove HSR clusters that overlap with a clipped position, i.e., the clipped position of a clipped cluster is contained in the HSR cluster
// direction of clip must be the same
// clipped_consensuses must be sorted by breakpoint
void remove_hsr_overlapping_clipped(std::vector<consensus_t*>& hsr_consensuses, std::vector<consensus_t*>& clipped_consensuses) {
	std::sort(hsr_consensuses.begin(), hsr_consensuses.end(),
			[](const consensus_t* c1, const consensus_t* c2) { return c1->start < c2->start; });
	std::vector<consensus_t*> kept_consensus;
	int i = 0;
	for (consensus_t* c : hsr_consensuses) {
		while (i < clipped_consensuses.size() && clipped_consensuses[i]->breakpoint < c->start) i++;

		// clipped_consensuses[i] must be same clip direction as cluster and breakpoint must be the smallest s.t. >= cluster.start
		if (i >= clipped_consensuses.size() || clipped_consensuses[i]->breakpoint > c->end) {
			kept_consensus.push_back(c);
		} else {
			delete c;
		}
	}
	kept_consensus.swap(hsr_consensuses);
}

#endif /* HSR_H_ */
