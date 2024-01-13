#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_

#include <queue>
#include <unordered_set>
#include <mutex>
#include <htslib/sam.h>

#include "../libs/ssw.h"
#include "sw_utils.h"
#include "types.h"
#include "dc_remapper.h"

extern std::ofstream assembly_failed_no_seq, assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer;
std::mutex failed_assembly_mtx;

extern const int TOO_MANY_READS;

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

struct path_permission_t {
	bool can_start_path, can_end_path;
};
struct seq_w_pp_t {
	std::string seq;
	path_permission_t clip_pair;

	seq_w_pp_t() : seq(), clip_pair() {}
	seq_w_pp_t(std::string& seq, bool can_start_path, bool can_end_path) : seq(seq) {
		clip_pair.can_start_path = can_start_path;
		clip_pair.can_end_path = can_end_path;
	}
};

std::vector<int> find_rev_topological_order(int n, std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

void build_aln_guided_graph(std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& alns, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, config_t& config) {
	std::sort(alns.begin(), alns.end(),
			[](const std::pair<std::string, StripedSmithWaterman::Alignment>& aln1, const std::pair<std::string, StripedSmithWaterman::Alignment>& aln2) {
		return aln1.second.ref_begin < aln2.second.ref_begin;
	});

	for (int i = 0; i < alns.size(); i++) {
		for (int j = i+1; j < alns.size() && alns[i].second.ref_end-alns[j].second.ref_begin >= config.min_clip_len; j++) {
			suffix_prefix_aln_t spa = aln_suffix_prefix(alns[i].first, alns[j].first, 1, -4, config.max_seq_error, config.min_clip_len);
			if (spa.overlap) {
				out_edges[i]++;
				l_adj[i].push_back({j, spa.score, spa.overlap});
				l_adj_rev[j].push_back({i, spa.score, spa.overlap});
			}
		}
	}

	// two major differences with the regular assembly:
	// 1 - by how the graph is defined, no cycle is possible here
	// 2 - in regular assembly, we only report contigs made of at last 2 reads. Here we report even single reads
}

bool accept(StripedSmithWaterman::Alignment& aln, int min_clip_len, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);

	bool lc_is_noise = false, rc_is_noise = false;
	if (!qual_ascii.empty()) {
		if (lc_size > min_clip_len && lc_size < qual_ascii.length()/2) {
			std::string lc_qual = qual_ascii.substr(0, lc_size);
			std::string not_lc_qual = qual_ascii.substr(lc_size);
			lc_is_noise = avg_qual(lc_qual) < min_avg_base_qual && avg_qual(not_lc_qual) >= min_avg_base_qual;
		}
		if (rc_size > min_clip_len && rc_size < qual_ascii.length()/2) {
			std::string not_rc_qual = qual_ascii.substr(0, qual_ascii.length()-rc_size);
			std::string rc_qual = qual_ascii.substr(qual_ascii.length()-rc_size);
			rc_is_noise = avg_qual(rc_qual) < min_avg_base_qual && avg_qual(not_rc_qual) >= min_avg_base_qual;
		}
	}

	double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
	return (lc_size <= min_clip_len || lc_is_noise) && (rc_size <= min_clip_len || rc_is_noise) && mismatch_rate <= max_seq_error;
}
void add_alignment(std::string& reference, std::string& query, std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& accepted_alns,
		std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& rejected_alns, StripedSmithWaterman::Aligner& aligner, config_t& config) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	aligner.Align(query.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
	if (accept(aln, config.min_clip_len)) {
		accepted_alns.push_back({query, aln});
	} else {
		rejected_alns.push_back({query, aln});
	}
}

void correct_contig(std::string& contig, std::vector<std::string>& reads, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config) {
	std::vector<int> As(contig.length()), Cs(contig.length()), Gs(contig.length()), Ts(contig.length());
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	for (std::string& read : reads) {
		harsh_aligner.Align(read.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len)) {
			for (int i = aln.query_begin; i < aln.query_end; i++) {
				char c = read[i];
				if (c == 'A') As[i-aln.query_begin+aln.ref_begin]++;
				else if (c == 'C') Cs[i-aln.query_begin+aln.ref_begin]++;
				else if (c == 'G') Gs[i-aln.query_begin+aln.ref_begin]++;
				else if (c == 'T') Ts[i-aln.query_begin+aln.ref_begin]++;
			}
		}
	}

	for (int i = 0; i < contig.length(); i++) {
		int max_freq = max(As[i], Cs[i], Gs[i], Ts[i]);
		if (max_freq == 0) continue;
		if (max_freq == As[i]) contig[i] = 'A';
		else if (max_freq == Cs[i]) contig[i] = 'C';
		else if (max_freq == Gs[i]) contig[i] = 'G';
		else if (max_freq == Ts[i]) contig[i] = 'T';
	}
}

void build_graph(std::vector<std::string>& read_seqs, std::vector<int>& order, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev,
		int max_mismatches, int min_overlap) {

	int n = read_seqs.size();

	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			std::string& s1 = read_seqs[i];
			std::string& s2 = read_seqs[j];

			int max_overlap = std::min(s1.length(), s2.length())-1;
			suffix_prefix_aln_t spa1 = aln_suffix_prefix(s1, s2, 1, -4, 1.0, min_overlap, max_overlap, max_mismatches);
			bool spa1_homopolymer = is_homopolymer(s2.c_str(), spa1.overlap);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(s2, s1, 1, -4, 1.0, min_overlap, max_overlap, max_mismatches);
			bool spa2_homopolymer = is_homopolymer(s1.c_str(), spa2.overlap);
			if (spa1.overlap && spa2.overlap) {
				if (spa1.score >= spa2.score && order[i] <= order[j]) {
					if (!spa1_homopolymer) {
						out_edges[i]++;
						l_adj[i].push_back({j, spa1.score, spa1.overlap});
						l_adj_rev[j].push_back({i, spa1.score, spa1.overlap});
					}
				} else if (spa1.score < spa2.score && order[j] <= order[i]) {
					if (!spa2_homopolymer) {
						out_edges[j]++;
						l_adj[j].push_back({i, spa2.score, spa2.overlap});
						l_adj_rev[i].push_back({j, spa2.score, spa2.overlap});
					}
				}
			} else if (spa1.overlap && order[i] <= order[j]) {
				if (!spa1_homopolymer) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa1.score, spa1.overlap});
					l_adj_rev[j].push_back({i, spa1.score, spa1.overlap});
				}
			} else if (spa2.overlap && order[j] <= order[i]) {
				if (!spa2_homopolymer) {
					out_edges[j]++;
					l_adj[j].push_back({i, spa2.score, spa2.overlap});
					l_adj_rev[i].push_back({j, spa2.score, spa2.overlap});
				}
			}
		}
	}
}

std::vector<std::string> assemble_reads(std::vector<seq_w_pp_t>& left_stable_read_seqs, std::vector<seq_w_pp_t>& unstable_read_seqs,
		std::vector<seq_w_pp_t>& right_stable_read_seqs, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, stats_t& stats) {

	std::vector<std::string> read_seqs;
	std::vector<path_permission_t> path_permissions;
	std::vector<int> order;
	for (seq_w_pp_t& s : left_stable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(1);
	}
	for (seq_w_pp_t& s : unstable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(2);
	}
	for (seq_w_pp_t& s : right_stable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(3);
	}

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);

	build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 1, config.min_clip_len);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	if (rev_topological_order.size() < n) {
		build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 0.0, config.min_clip_len);

		int min_overlap = config.min_clip_len;
		for (; min_overlap <= stats.read_len/2; min_overlap += 10) {
			for (int i = 0; i < n; i++) {
				l_adj[i].erase(std::remove_if(l_adj[i].begin(), l_adj[i].end(),
						[&min_overlap](edge_t& e) { return e.overlap < min_overlap; }), l_adj[i].end());
				l_adj_rev[i].erase(std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(),
						[&min_overlap](edge_t& e) { return e.overlap < min_overlap; }), l_adj_rev[i].end());
				out_edges[i] = l_adj[i].size();
			}

			rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

			if (rev_topological_order.size() == n) {
				break;
			}
		}
		if (rev_topological_order.size() < n) {
			return {"HAS_CYCLE"};
		}
	}

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			if (best_scores[i] == 0 && !path_permissions[i].can_end_path) continue; // sink and cannot end path => discard
			for (edge_t& e : l_adj_rev[i]) {
				if (best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (path_permissions[i].can_start_path && best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = read_seqs[curr_vertex];
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		used_reads.push_back(read_seqs[curr_vertex]);
		while (best_edges[curr_vertex].overlap) {
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += read_seqs[curr_vertex].substr(overlap);
			used_reads.push_back(read_seqs[curr_vertex]);
		}
		used[curr_vertex] = true;
		correct_contig(assembled_sequence, used_reads, harsh_aligner, config);
		assembled_sequences.push_back(assembled_sequence);
	}

	return assembled_sequences;
}

std::vector<std::string> assemble_sequences(std::string contig_name, insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, stats_t& stats) {
	std::vector<seq_w_pp_t> left_stable_read_seqs, unstable_read_seqs, right_stable_read_seqs;
	std::unordered_set<std::string> used_ls, used_us, used_rs;
	if (r_cluster->clip_consensus) left_stable_read_seqs.push_back({r_cluster->clip_consensus->sequence, true, false});
	for (bam1_t* read : r_cluster->cluster->reads) {
		std::string seq = get_sequence(read);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		rc(mate_seq);
		if (!used_ls.count(seq)) {
			left_stable_read_seqs.push_back({seq, !is_left_clipped(read, 1), true});
			used_ls.insert(seq);
		}
		if (!used_us.count(mate_seq)) {
			unstable_read_seqs.push_back({mate_seq, true, true});
			used_us.insert(mate_seq);
		}
	}
	for (bam1_t* read : l_cluster->cluster->reads) {
		std::string seq = get_sequence(read);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		if (!used_rs.count(seq)) {
			right_stable_read_seqs.push_back({seq, true, !is_right_clipped(read, 1)});
			used_rs.insert(seq);
		}
		if (!used_us.count(mate_seq)) {
			unstable_read_seqs.push_back({mate_seq, true, true});
			used_us.insert(mate_seq);
		}
	}
	if (l_cluster->clip_consensus) right_stable_read_seqs.push_back({l_cluster->clip_consensus->sequence, false, true});

	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		if (!used_us.count(read_seq)) {
			unstable_read_seqs.push_back({read_seq, !is_left_clipped(read, 1), !is_right_clipped(read, 1)});
			used_us.insert(read_seq);
		}
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		if (!used_us.count(read_seq)) {
			unstable_read_seqs.push_back({read_seq, !is_left_clipped(read, 1), !is_right_clipped(read, 1)});
			used_us.insert(read_seq);
		}
	}

	if (unstable_read_seqs.size() + left_stable_read_seqs.size() + right_stable_read_seqs.size() >= TOO_MANY_READS) return {"TOO_MANY_READS"};

	std::vector<std::string> assembled_sequences = assemble_reads(left_stable_read_seqs, unstable_read_seqs, right_stable_read_seqs,
			harsh_aligner, config, stats);

	return assembled_sequences;
}

sv_t* detect_de_novo_insertion(std::string& contig_name, chr_seqs_map_t& contigs,
		insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<bam1_t*>& assembled_reads, config_t& config, stats_t& stats) {

	std::vector<std::string> assembled_sequences = assemble_sequences(contig_name, r_cluster, l_cluster, mateseqs, harsh_aligner, config, stats);

	std::string ins_full_id = "NO_ID";
	if (assembled_sequences.empty()) {
		failed_assembly_mtx.lock();
		assembly_failed_no_seq << ins_full_id << " " << contig_name << " " << r_cluster->end << " + ";
		assembly_failed_no_seq << contig_name << " " << l_cluster->start << " - INS" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	} else if (assembled_sequences[0] == "HAS_CYCLE") {
		failed_assembly_mtx.lock();
		assembly_failed_cycle_writer << ins_full_id << " " << contig_name << " " << r_cluster->end << " + ";
		assembly_failed_cycle_writer << contig_name << " " << l_cluster->start << " - INS HAS_CYCLE" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	} else if (assembled_sequences[0] == "TOO_MANY_READS") {
		failed_assembly_mtx.lock();
		assembly_failed_too_many_reads_writer << ins_full_id << " " << contig_name << " " << r_cluster->end << " + ";
		assembly_failed_too_many_reads_writer << contig_name << " " << l_cluster->start << " - INS TOO_MANY_READS" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	}

	// remap assembled sequences to find breakpoints
	char* contig_seq = contigs.get_seq(contig_name);
	int contig_len = contigs.get_len(contig_name);

	std::string assembled_sequence = assembled_sequences[0];

	// divide the assembled sequence into two and remap them
	int extend = 2*assembled_sequence.length();
	int remap_region_start = std::min(r_cluster->end, l_cluster->start)-extend;
	remap_region_start = std::max(0, remap_region_start);
	int remap_region_end = std::max(r_cluster->end, l_cluster->start)+extend;
	remap_region_end = std::min(remap_region_end, contig_len-1);

	StripedSmithWaterman::Alignment aln;
	std::string full_assembled_seq;
	std::vector<sv_t*> svs = detect_svs_from_junction(contig_name, contig_seq, assembled_sequence, remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, config.min_clip_len);
	sv_t* chosen_ins = NULL;
	if (svs.empty() || svs[0]->svtype() != "INS" || svs[0]->svlen() < config.min_sv_size) { // cannot identify an insertion, check if it is an incomplete assembly
		if (assembled_sequences.size() == 1) return NULL;

		extend = 2*(assembled_sequence.length() + assembled_sequences[1].length());
		remap_region_start = std::min(r_cluster->end, l_cluster->start)-extend;
		if (remap_region_start < 0) remap_region_start = 0;
		remap_region_end = std::max(r_cluster->end, l_cluster->start)+extend;
		if (remap_region_end >= contig_len) remap_region_end = contig_len-1;

		std::vector<sv_t*> svs1 = detect_svs_from_junction(contig_name, contig_seq, assembled_sequence + "-" + assembled_sequences[1], remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, config.min_clip_len);
		std::vector<sv_t*> svs2 = detect_svs_from_junction(contig_name, contig_seq, assembled_sequences[1] + "-" + assembled_sequence, remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, config.min_clip_len);
		bool sv1_is_ins = !svs1.empty() && svs1[0]->svtype() == "INS";
		bool sv2_is_ins = !svs2.empty() && svs2[0]->svtype() == "INS";
		if (sv1_is_ins && !sv2_is_ins) {
			chosen_ins = svs1[0];
		} else if (!sv1_is_ins && sv2_is_ins) {
			chosen_ins = svs2[0];
		} else if (sv1_is_ins && sv2_is_ins) {
			int sv1_score = svs1[0]->left_anchor_aln->best_score + svs1[0]->right_anchor_aln->best_score;
			int sv2_score = svs2[0]->left_anchor_aln->best_score + svs2[0]->right_anchor_aln->best_score;
			if (sv1_score >= sv2_score) chosen_ins = svs1[0];
			else chosen_ins = svs2[0];
		}

		if (chosen_ins == svs1[0]) full_assembled_seq = assembled_sequence + "-" + assembled_sequences[1];
		else if (chosen_ins == svs2[0]) full_assembled_seq = assembled_sequences[1] + "-" + assembled_sequence;
	} else if (assembled_sequences.size() > 1) { // can identify an insertion, but check if we can get a better one as incomplete assembly
		chosen_ins = svs[0];
		full_assembled_seq = assembled_sequence;

		double laa_score_ratio = double(chosen_ins->left_anchor_aln->best_score)/chosen_ins->left_anchor_aln->seq_len;
		double raa_score_ratio = double(chosen_ins->right_anchor_aln->best_score)/chosen_ins->right_anchor_aln->seq_len;
		if (laa_score_ratio >= raa_score_ratio) {
			std::vector<sv_t*> new_svs = detect_svs_from_junction(contig_name, contig_seq, assembled_sequence + "-" + assembled_sequences[1], remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, config.min_clip_len);
			if (!new_svs.empty() && new_svs[0]->svtype() == "INS") {
				sv_t* new_ins = new_svs[0];
				int sep = new_ins->ins_seq.find("-");
				if (sep != std::string::npos && sep >= config.min_clip_len && new_ins->ins_seq.length()-sep >= config.min_clip_len &&
					raa_score_ratio < double(new_ins->right_anchor_aln->best_score)/new_ins->right_anchor_aln->seq_len && new_ins->svlen() >= config.min_sv_size) {
						chosen_ins = new_ins;
						full_assembled_seq = assembled_sequence + "-" + assembled_sequences[1];
					}
			}
		} else {
			std::vector<sv_t*> new_svs = detect_svs_from_junction(contig_name, contig_seq, assembled_sequences[1] + "-" + assembled_sequence, remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, config.min_clip_len);
			if (!new_svs.empty() && new_svs[0]->svtype() == "INS") {
				sv_t* new_ins = new_svs[0];
				int sep = new_ins->ins_seq.find("-");
				if (sep != std::string::npos && sep >= config.min_clip_len && new_ins->ins_seq.length()-sep >= config.min_clip_len &&
					laa_score_ratio < double(new_ins->left_anchor_aln->best_score)/new_ins->left_anchor_aln->seq_len && new_ins->svlen() >= config.min_sv_size) {
						chosen_ins = new_ins;
						full_assembled_seq = assembled_sequences[1] + "-" + assembled_sequence;
					}
			}
		}
	} else {
		chosen_ins = svs[0];
		full_assembled_seq = assembled_sequence;
	}

	if (chosen_ins == NULL) return NULL;

	int sep = chosen_ins->ins_seq.find("-");
	if (sep != std::string::npos && (sep < config.min_clip_len || chosen_ins->ins_seq.length()-sep < config.min_clip_len)) return NULL; // "-" is too close to the edge

	StripedSmithWaterman::Filter filter;
	if (r_cluster->clip_consensus) {
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(r_cluster->clip_consensus->sequence.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len, config.max_seq_error)) {
			chosen_ins->rc_consensus = r_cluster->clip_consensus;
		}
	}
	for (bam1_t* read : r_cluster->cluster->reads) {
		std::string mate_seq = get_mate_seq(read, mateseqs);
		std::string mate_qual = get_mate_qual(read, matequals);
		rc(mate_seq);
		mate_qual = std::string(mate_qual.rbegin(), mate_qual.rend());
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(mate_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len, config.max_seq_error, mate_qual, stats.min_avg_base_qual)) {
			chosen_ins->disc_pairs_lf++;
			chosen_ins->disc_pairs_lf_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
			bam1_t* d = bam_dup1(read);
			d->core.mpos = 0;
			d->core.mtid = d->core.tid;
			d->core.flag |= BAM_FMUNMAP;
			assembled_reads.push_back(d);
		}
	}
	chosen_ins->disc_pairs_lf_avg_nm /= chosen_ins->disc_pairs_lf;
	for (bam1_t* read : l_cluster->cluster->reads) {
		std::string mate_seq = get_mate_seq(read, mateseqs);
		std::string mate_qual = get_mate_qual(read, matequals);
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(mate_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len, config.max_seq_error, mate_qual, stats.min_avg_base_qual)) {
			chosen_ins->disc_pairs_rf++;
			chosen_ins->disc_pairs_rf_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
			bam1_t* d = bam_dup1(read);
			d->core.mpos = 0;
			d->core.mtid = d->core.tid;
			d->core.flag |= BAM_FMUNMAP;
			assembled_reads.push_back(d);
		}
	}
	chosen_ins->disc_pairs_rf_avg_nm /= chosen_ins->disc_pairs_rf++;
	if (l_cluster->clip_consensus) {
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(l_cluster->clip_consensus->sequence.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len, config.max_seq_error)) {
			chosen_ins->lc_consensus = l_cluster->clip_consensus;
		}
	}

	// start and end of inserted sequence within the full assembled sequence
	int ins_seq_start = chosen_ins->left_anchor_aln->seq_len - chosen_ins->prefix_mh_len;
	int ins_seq_end = ins_seq_start + chosen_ins->ins_seq.length(); // corrected_consensus_sequence.length() - (chosen_ins->right_anchor_aln->seq_len - chosen_ins->suffix_mh_len);
	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		StripedSmithWaterman::Alignment aln;
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		harsh_aligner.Align(read_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		std::string qual_ascii = get_qual_ascii(read, true);
		qual_ascii = std::string(qual_ascii.rbegin(), qual_ascii.rend());
		if (overlap(ins_seq_start, ins_seq_end, aln.ref_begin, aln.ref_end) >= config.min_clip_len
				&& accept(aln, config.min_clip_len, config.max_seq_error, qual_ascii, stats.min_avg_base_qual)) {
			chosen_ins->disc_pairs_lf++;
			assembled_reads.push_back(bam_dup1(read));
		}
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		StripedSmithWaterman::Alignment aln;
		std::string read_seq = get_sequence(read, true);
		harsh_aligner.Align(read_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		std::string qual_ascii = get_qual_ascii(read, true);
		if (overlap(ins_seq_start, ins_seq_end, aln.ref_begin, aln.ref_end) >= config.min_clip_len
				&& accept(aln, config.min_clip_len, config.max_seq_error, qual_ascii, stats.min_avg_base_qual)) {
			chosen_ins->disc_pairs_rf++;
			assembled_reads.push_back(bam_dup1(read));
		}
	}

	for (bam1_t* read : assembled_reads) {
		bam_aux_update_str(read, "ID", chosen_ins->id.length(), chosen_ins->id.c_str());
	}

	chosen_ins->source = "DE_NOVO_ASSEMBLY";
	return chosen_ins;
}

#endif /* ASSEMBLE_H_ */
