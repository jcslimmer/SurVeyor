#ifndef EXTEND_1SR_CONSENSUS_H_
#define EXTEND_1SR_CONSENSUS_H_

#include <set>
#include <unordered_map>
#include <queue>
#include <stack>

#include "../libs/IntervalTree.h"
#include "utils.h"
#include "sw_utils.h"
#include "sam_utils.h"
#include "types.h"

struct ext_read_t {
	std::string qname;
	hts_pos_t start, end;
	uint8_t mapq;
	bool rev, same_chr;
	int sequence_len = 0;
	char* sequence = NULL;

	ext_read_t(bam1_t* read) {
		init(read);
	}

	void init(bam1_t* read) {
		qname = bam_get_qname(read);
		start = read->core.pos;
		end = bam_endpos(read);
		mapq = read->core.qual;
		rev = bam_is_rev(read);
		same_chr = is_samechr(read);
		if (sequence_len < read->core.l_qseq) {
			if (sequence != NULL) delete[] sequence;
			sequence = new char[read->core.l_qseq+1];
		}
		copy_sequence(read, sequence);
		sequence_len = read->core.l_qseq;
	}

	~ext_read_t() {
		if (sequence != NULL) delete[] sequence;
	}
};

struct ext_read_allocator_t {
	std::queue<ext_read_t*> reads;

	ext_read_allocator_t() {}

	ext_read_t* get(bam1_t* read) {
		if (reads.empty()) {
			return new ext_read_t(read);
		} else {
			ext_read_t* r = reads.front();
			reads.pop();
			r->init(read);
			return r;
		}
	}

	void release(ext_read_t* read) {
		reads.push(read);
	}
};

bool is_vertex_in_cycle(std::vector<std::vector<edge_t> >& l_adj, int i) {
	std::stack<int> s;
	for (edge_t& e : l_adj[i]) {
		s.push(e.next);
	}

	// check if node i can reach itself - i.e., it is part of a cycle
	std::vector<bool> visited(l_adj.size(), false);
	bool i_in_cycle = false;
	while (!s.empty()) {
		int curr = s.top();
		s.pop();
		if (visited[curr]) continue;
		visited[curr] = true;

		if (curr == i) return true;
		for (edge_t& e : l_adj[curr]) {
			s.push(e.next);
		}
	}
	return false;
}

std::vector<bool> find_vertices_in_cycles(std::vector<std::vector<edge_t> >& l_adj) {

	int n = l_adj.size();
	std::vector<bool> in_cycle(n);
	for (int i = 0; i < n; i++) {
		in_cycle[i] = is_vertex_in_cycle(l_adj, i);
	}
	return in_cycle;
}

void build_graph_fwd(std::vector<std::string>& read_seqs, std::vector<int> starting_idxs, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, int min_overlap) {

	uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

	int n = read_seqs.size();
	std::vector<bool> visited(n);

	std::unordered_map<uint64_t, std::set<int> > kmer_to_idx;
	for (int i = 0; i < n; i++) {
		uint64_t kmer = 0;

		std::string& seq = read_seqs[i];
		for (int j = 0; j < seq.length()-1; j++) {
			uint64_t nv = nucl_bm[seq[j]];
			kmer = ((kmer << 2) | nv);

			if (j >= 32) {
				kmer_to_idx[kmer].insert(i);
			}
		}
	}

	std::queue<int> bfs;
	for (int i : starting_idxs) bfs.push(i);
	while (!bfs.empty()) {
		int i = bfs.front();
		bfs.pop();

		if (visited[i]) continue;
		visited[i] = true;

		std::string& s1 = read_seqs[i];
		uint64_t kmer = 0;
		for (int j = s1.length()-32; j < s1.length(); j++) {
			uint64_t nv = nucl_bm[s1[j]];
			kmer = ((kmer << 2) | nv);
		}

		for (int j : kmer_to_idx[kmer]) {
			std::string& s2 = read_seqs[j];

			if (s1 == s2) continue;

			suffix_prefix_aln_t spa = aln_suffix_prefix_perfect(s1, s2, min_overlap);
			if (spa.overlap) {
				if (!is_homopolymer(s2.c_str(), spa.overlap)) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa.score, spa.overlap});
					l_adj_rev[j].push_back({i, spa.score, spa.overlap});
					bfs.push(j);
				}
			}
		}
	}
}

void build_graph_rev(std::vector<std::string>& read_seqs, std::vector<int> starting_idxs, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, int min_overlap) {

	uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

	int n = read_seqs.size();
	std::vector<bool> visited(n);

	std::unordered_map<uint64_t, std::vector<int> > kmer_to_idx;

	for (int i = 0; i < n; i++) {
		uint64_t kmer = 0;

		std::string& seq = read_seqs[i];
		for (int j = 0; j < seq.length(); j++) {
			uint64_t nv = nucl_bm[seq[j]];
			kmer = ((kmer << 2) | nv);

			if (j >= 32+1) {
				kmer_to_idx[kmer].push_back(i);
			}
		}
	}

	std::queue<int> bfs;
	for (int i : starting_idxs) bfs.push(i);
	while (!bfs.empty()) {
		int i = bfs.front();
		bfs.pop();

		if (visited[i]) continue;
		visited[i] = true;

		std::string& s1 = read_seqs[i];
		uint64_t kmer = 0;
		for (int j = 0; j < 32; j++) {
			uint64_t nv = nucl_bm[s1[j]];
			kmer = ((kmer << 2) | nv);
		}

		for (int j : kmer_to_idx[kmer]) {
			std::string& s2 = read_seqs[j];

			if (s1 == s2) continue;

			suffix_prefix_aln_t spa2 = aln_suffix_prefix_perfect(s2, s1, min_overlap);
			if (spa2.overlap) {
				bool spa2_homopolymer = is_homopolymer(s1.c_str(), spa2.overlap);
				if (!spa2_homopolymer) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa2.score, spa2.overlap});
					l_adj_rev[j].push_back({i, spa2.score, spa2.overlap});
					bfs.push(j);
				}
			}
		}
	}
}

void get_extension_read_seqs(IntervalTree<ext_read_t*>& candidate_reads_itree, std::vector<std::string>& read_seqs, std::vector<int>& read_mapqs,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len, int high_confidence_mapq, stats_t& stats, int max_reads = INT32_MAX) {

	hts_pos_t fwd_mates_start = std::max(hts_pos_t(1), target_start-stats.max_is+stats.read_len);
	hts_pos_t fwd_mates_end = std::max(hts_pos_t(1), target_end-stats.min_is);

	hts_pos_t rev_mates_start = std::min(target_start+stats.min_is, contig_len);
	hts_pos_t rev_mates_end = std::min(target_end+stats.max_is-stats.read_len, contig_len);

	std::vector<Interval<ext_read_t*> > target_reads = candidate_reads_itree.findOverlapping(
			std::min(fwd_mates_start, rev_mates_start)-10, std::max(fwd_mates_end, rev_mates_end)+10);
	for (Interval<ext_read_t*> i_read : target_reads) {
		ext_read_t* ext_read = i_read.value;
		if (!ext_read->rev && ext_read->mapq >= high_confidence_mapq &&
			!ext_read->same_chr && ext_read->start >= fwd_mates_start && ext_read->start <= fwd_mates_end) {
			if (mateseqs_w_mapq.count(ext_read->qname) == 0) {
				std::cerr << "ERROR: mateseqs_w_mapq does not contain " << ext_read->qname << std::endl;
				exit(1);
			};
			std::pair<std::string, int> mateseq_w_mapq = mateseqs_w_mapq[ext_read->qname];
			rc(mateseq_w_mapq.first);
			read_seqs.push_back(mateseq_w_mapq.first);
			read_mapqs.push_back(mateseq_w_mapq.second);
		} else if (ext_read->rev && ext_read->mapq >= high_confidence_mapq &&
				!ext_read->same_chr && ext_read->end >= rev_mates_start && ext_read->end <= rev_mates_end) {
			if (mateseqs_w_mapq.count(ext_read->qname) == 0) {
				std::cerr << "ERROR: mateseqs_w_mapq does not contain " << ext_read->qname << std::endl;
				exit(1);
			};
			std::pair<std::string, int> mateseq_w_mapq = mateseqs_w_mapq[ext_read->qname];
			read_seqs.push_back(mateseq_w_mapq.first);
			read_mapqs.push_back(mateseq_w_mapq.second);
		}

		if (ext_read->end > target_start && ext_read->start < target_end) {
			read_seqs.push_back(ext_read->sequence);
			read_mapqs.push_back(ext_read->mapq);
		}

		if (read_seqs.size() > max_reads) return;
	}

	// TODO : results very slightly change when sorting, investigate
//	std::vector<std::pair<std::string, int>> temp;
//	for (int i = 1; i < read_seqs.size(); i++) {
//		temp.push_back({read_seqs[i], read_mapqs[i]});
//	}
//	std::sort(temp.begin(), temp.end());
//
//	read_seqs = std::vector<std::string>(1, read_seqs[0]);
//	read_mapqs = std::vector<int>(1, read_mapqs[0]);
//	for (int i = 0; i < temp.size(); i++) {
//		read_seqs.push_back(temp[i].first);
//		read_mapqs.push_back(temp[i].second);
//	}
}
std::vector<ext_read_t*> get_extension_reads(std::string contig_name, std::vector<hts_pair_pos_t>& target_ivals, hts_pos_t contig_len,
		stats_t& stats, open_samFile_t* bam_file) {

	ext_read_allocator_t ext_read_allocator;

	std::sort(target_ivals.begin(), target_ivals.end(), [](hts_pair_pos_t& a, hts_pair_pos_t& b) {return a.beg < b.beg;});

	std::vector<hts_pair_pos_t> merged_target_ivals;
	merged_target_ivals.push_back(target_ivals[0]);
	for (int i = 1; i < target_ivals.size(); i++) {
		if (target_ivals[i].beg-merged_target_ivals.back().end < stats.read_len) {
			merged_target_ivals.back().end = std::max(merged_target_ivals.back().end, target_ivals[i].end);
		} else {
			merged_target_ivals.push_back(target_ivals[i]);
		}
	}

	std::vector<ext_read_t*> reads;
	bam1_t* read = bam_init1();
	for (hts_pair_pos_t target_ival : merged_target_ivals) {
		hts_pos_t fwd_mates_start = std::max(hts_pos_t(1), target_ival.beg-stats.max_is+stats.read_len);
		hts_pos_t rev_mates_end = std::min(target_ival.end+stats.max_is-stats.read_len, contig_len);

		std::stringstream ss;
		ss << contig_name << ":" << fwd_mates_start << "-" << rev_mates_end;

		std::vector<ext_read_t*> region_reads;
		region_reads.reserve(target_ival.end-target_ival.beg+1);
		hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, ss.str().c_str());
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			if (is_unmapped(read) || !is_primary(read)) continue;

			region_reads.push_back(ext_read_allocator.get(read));
			if (region_reads.size() > target_ival.end-target_ival.beg) { // too many reads
				for (ext_read_t* r : region_reads) ext_read_allocator.release(r);
				region_reads.clear();
				break;
			}
		}
		reads.insert(reads.end(), region_reads.begin(), region_reads.end());
		hts_itr_destroy(iter);
	}

	bam_destroy1(read);

	return reads;
}

std::vector<ext_read_t*> get_extension_reads_from_consensuses(std::vector<consensus_t*>& consensuses, std::string contig_name, hts_pos_t contig_len,
		stats_t& stats, open_samFile_t* bam_file) {

	if (consensuses.empty()) return std::vector<ext_read_t*>();

	std::vector<hts_pair_pos_t> target_ivals;
	for (consensus_t* consensus : consensuses) {
		hts_pos_t left_ext_target_start = consensus->left_ext_target_start(stats.max_is, stats.read_len);
		hts_pos_t left_ext_target_end = consensus->left_ext_target_end(stats.max_is, stats.read_len);
		hts_pos_t right_ext_target_start = consensus->right_ext_target_start(stats.max_is, stats.read_len);
		hts_pos_t right_ext_target_end = consensus->right_ext_target_end(stats.max_is, stats.read_len);
		hts_pair_pos_t target_ival;
		target_ival.beg = left_ext_target_start, target_ival.end = left_ext_target_end;
		target_ivals.push_back(target_ival);
		target_ival.beg = right_ext_target_start, target_ival.end = right_ext_target_end;
		target_ivals.push_back(target_ival);
	}
	return get_extension_reads(contig_name, target_ivals, contig_len, stats, bam_file);
}

void break_cycles(std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {
	int n = l_adj.size();
	std::vector<bool> in_cycle = find_vertices_in_cycles(l_adj);
	for (int i = 0; i < n; i++) {
		if (in_cycle[i]) {
			l_adj[i].clear();
			out_edges[i] = 0;
		}
		l_adj_rev[i].erase(
			std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(), [&in_cycle](edge_t& e) {return in_cycle[e.next];}),
			l_adj_rev[i].end()
		);
	}
}

void extend_consensus_to_right(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len,
		const int high_confidence_mapq, stats_t& stats, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	if (consensus->extended_to_right) return;

	std::vector<std::string> read_seqs;
	std::vector<int> read_mapqs;
	read_seqs.push_back(consensus->sequence);
	read_mapqs.push_back(high_confidence_mapq);

	get_extension_read_seqs(candidate_reads_itree, read_seqs, read_mapqs, mateseqs_w_mapq, target_start, target_end, contig_len, high_confidence_mapq, stats, 5000);

	if (read_seqs.size() > 5000) return;

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_graph_fwd(read_seqs, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, stats.read_len/2);

	bool cycle_contains_0 = is_vertex_in_cycle(l_adj, 0);
	if (cycle_contains_0) {
		// remove all edges pointing to 0
		for (std::vector<edge_t>& l_adj_i : l_adj) {
			l_adj_i.erase(std::remove_if(l_adj_i.begin(), l_adj_i.end(), [](edge_t& e) {return e.next == 0;}), l_adj_i.end());
		}
		l_adj_rev[0].clear();
	}

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	bool cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) {
		break_cycles(out_edges, l_adj, l_adj_rev);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		for (edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}


	// 0 is the consensus, we start from there
	edge_t e = best_edges[0];
	std::string ext_consensus = consensus->sequence;
	while (e.overlap) {
		ext_consensus += read_seqs[e.next].substr(e.overlap);
		e = best_edges[e.next];
		consensus->right_ext_reads++;
		if (read_mapqs[e.next] >= high_confidence_mapq) consensus->hq_right_ext_reads++;
	}

	if (!consensus->left_clipped) {
		if (consensus->clip_len != consensus_t::UNKNOWN_CLIP_LEN) {
			consensus->clip_len += ext_consensus.length() - consensus->sequence.length();
		}
	} else {
		consensus->end += ext_consensus.length() - consensus->sequence.length();
	}
	consensus->sequence = ext_consensus;
	consensus->extended_to_right = true;
}

void extend_consensus_to_left(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len,
		const int high_confidence_mapq, stats_t& stats, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	if (consensus->extended_to_left) return;

	std::vector<std::string> read_seqs;
	std::vector<int> read_mapqs;
	read_seqs.push_back(consensus->sequence);
	read_mapqs.push_back(high_confidence_mapq);

	get_extension_read_seqs(candidate_reads_itree, read_seqs, read_mapqs, mateseqs_w_mapq, target_start, target_end, contig_len, high_confidence_mapq, stats, 5000);

	if (read_seqs.size() > 5000) return;

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_graph_rev(read_seqs, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, stats.read_len/2);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	bool cycle_contains_0 = is_vertex_in_cycle(l_adj, 0);
	if (cycle_contains_0) {
		// remove all edges pointing to 0
		for (std::vector<edge_t>& l_adj_i : l_adj) {
			l_adj_i.erase(std::remove_if(l_adj_i.begin(), l_adj_i.end(), [](edge_t& e) {return e.next == 0;}), l_adj_i.end());
		}
		l_adj_rev[0].clear();
	}

	bool cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) {
		break_cycles(out_edges, l_adj, l_adj_rev);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		for (edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}

	// 0 is the consensus, we start from there
	edge_t e = best_edges[0];
	std::string ext_consensus = consensus->sequence;
	while (e.overlap) {
		ext_consensus = read_seqs[e.next].substr(0, read_seqs[e.next].length()-e.overlap) + ext_consensus;
		e = best_edges[e.next];
		consensus->left_ext_reads++;
		if (read_mapqs[e.next] >= high_confidence_mapq) consensus->hq_left_ext_reads++;
	}

	if (consensus->left_clipped) {
		if (consensus->clip_len != consensus_t::UNKNOWN_CLIP_LEN) {
			consensus->clip_len += ext_consensus.length() - consensus->sequence.length();
		}
	} else {
		consensus->start -= ext_consensus.length() - consensus->sequence.length();
	}
	consensus->sequence = ext_consensus;
	consensus->extended_to_left = true;
}


#endif /* EXTEND_1SR_CONSENSUS_H_ */
