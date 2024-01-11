#ifndef GUIDED_REFERENCE_ASSEMBLE_H
#define GUIDED_REFERENCE_ASSEMBLE_H

#include <iostream>

#include "../libs/ssw_cpp.h"
#include "../src/utils.h"
#include "../src/dc_remapper.h"
#include "assemble.h"

std::vector<std::string> generate_reference_guided_consensus(std::string reference, insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, config_t& config, stats_t& stats) {

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > accepted_alns, _rejected_alns_lf, _rejected_alns_is, _rejected_alns_rf;
	for (bam1_t* read : r_cluster->cluster->reads) {
		std::string read_seq = get_sequence(read);
		add_alignment(reference, read_seq, accepted_alns, _rejected_alns_lf, aligner, config);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		rc(mate_seq);
		add_alignment(reference, mate_seq, accepted_alns, _rejected_alns_is, aligner, config);
	}
	for (bam1_t* read : l_cluster->cluster->reads) {
		std::string read_seq = get_sequence(read);
		add_alignment(reference, read_seq, accepted_alns, _rejected_alns_rf, aligner, config);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		add_alignment(reference, mate_seq, accepted_alns, _rejected_alns_is, aligner, config);
	}
	if (r_cluster->clip_consensus) {
		add_alignment(reference, r_cluster->clip_consensus->sequence, accepted_alns, _rejected_alns_lf, aligner, config);
	}
	if (l_cluster->clip_consensus) {
		add_alignment(reference, l_cluster->clip_consensus->sequence, accepted_alns, _rejected_alns_rf, aligner, config);
	}

	int n = accepted_alns.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_aln_guided_graph(accepted_alns, out_edges, l_adj, l_adj_rev, config);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			for (edge_t& e : l_adj_rev[i]) {
				if (best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = accepted_alns[curr_vertex].first;
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		used_reads.push_back(accepted_alns[curr_vertex].first);
		while (best_edges[curr_vertex].overlap) {
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += accepted_alns[curr_vertex].first.substr(overlap);
			used_reads.push_back(accepted_alns[curr_vertex].first);
		}
		used[curr_vertex] = true;

		std::string corrected_assembled_sequence = assembled_sequence;
		correct_contig(corrected_assembled_sequence, used_reads, harsh_aligner, config);
		assembled_sequences.push_back(corrected_assembled_sequence);
	}
	for (int i = 0; i < n; i++) {
		if (!used[i]) assembled_sequences.push_back(accepted_alns[i].first);
	}

	// retain assembled sequences that align without clipping and do not overlap a higher rated sequence
	std::vector<std::string> retained_assembled_sequences;
	for (std::string& assembled_sequence : assembled_sequences) {
		aligner.Align(assembled_sequence.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		if (!accept(aln, config.min_clip_len)) continue;

		bool overlaps = false;
		for (StripedSmithWaterman::Alignment& existing_aln : consensus_contigs_alns) {
			if (std::max(existing_aln.ref_begin, aln.ref_begin) <= std::min(existing_aln.ref_end, aln.ref_end)) {
				overlaps = true;
				break;
			}
		}
		if (!overlaps) {
			consensus_contigs_alns.push_back(aln);
			retained_assembled_sequences.push_back(assembled_sequence);
		}
	}

	if (retained_assembled_sequences.empty()) return {};

	int l_bp_in_seq = r_cluster->end-r_cluster->start, r_bp_in_seq = reference.length()-(l_cluster->end-l_cluster->start);

	// try scaffolding using rejected reads
	std::vector<seq_w_pp_t> rejected_alns_lf, rejected_alns_is, rejected_alns_rf;
	for (auto& e : _rejected_alns_lf) rejected_alns_lf.push_back({e.first, true, true});
	for (auto& e : _rejected_alns_is) rejected_alns_is.push_back({e.first, true, true});
	for (auto& e : _rejected_alns_rf) rejected_alns_rf.push_back({e.first, true, true});
	std::vector<char> rejected_alns_lf_clipped(rejected_alns_lf.size(), 'N'), rejected_alns_is_clipped(rejected_alns_is.size(), 'N'),
			rejected_alns_rf_clipped(rejected_alns_rf.size(), 'N');
	std::vector<std::string> scaffolds = assemble_reads(rejected_alns_lf, rejected_alns_is, rejected_alns_rf,
			harsh_aligner, config, stats);

	for (std::string a : assembled_sequences) {
		aligner.Align(a.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
	}
	for (int i = 0; i < retained_assembled_sequences.size(); i++) {
		StripedSmithWaterman::Alignment& aln = consensus_contigs_alns[i];
	}

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
	for (int i = 0; i < consensus_contigs_alns.size(); i++) {
		contigs_sorted_by_pos.push_back({retained_assembled_sequences[i], consensus_contigs_alns[i]});
	}
	std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	std::vector<int> linked(contigs_sorted_by_pos.size()-1, -1);
	std::vector<std::pair<int, int> > link_overlaps(contigs_sorted_by_pos.size()-1);
	for (int i = 0; i < scaffolds.size(); i++) {
		std::string& scaffold = scaffolds[i];
		int best_link = -1, best_link_w = 0;
		std::pair<int, int> best_link_overlap;
		for (int j = 0; j < contigs_sorted_by_pos.size()-1; j++) {
			if (linked[j] >= 0) continue;

			suffix_prefix_aln_t spa1 = aln_suffix_prefix(contigs_sorted_by_pos[j].first, scaffold, 1, -4, config.max_seq_error, config.min_clip_len);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(scaffold, contigs_sorted_by_pos[j+1].first, 1, -4, config.max_seq_error, config.min_clip_len);
			if (spa1.overlap && spa2.overlap && best_link_w < spa1.score+spa2.score) {
				best_link = j, best_link_w = spa1.score+spa2.score;
				best_link_overlap = {spa1.overlap, spa2.overlap};
			}
		}

		if (best_link >= 0) {
			linked[best_link] = i;
			link_overlaps[best_link] = best_link_overlap;
		}
	}

	std::vector<std::string> scaffolded_sequences;
	std::string curr_seq = contigs_sorted_by_pos[0].first;
	for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
		if (linked[i-1] == -1) {
			scaffolded_sequences.push_back(curr_seq);
			curr_seq = contigs_sorted_by_pos[i].first;
		} else {
			std::string link = scaffolds[linked[i-1]];
			auto& lo = link_overlaps[i-1];
			link = link.substr(lo.first, link.length()-lo.first-lo.second);
			curr_seq += link + contigs_sorted_by_pos[i].first;
		}
	}
	scaffolded_sequences.push_back(curr_seq);

	std::vector<StripedSmithWaterman::Alignment> scaffolded_seqs_alns;
	bool scaffolding_failed = false;
	for (std::string s : scaffolded_sequences) {
		aligner.Align(s.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		if (!accept(aln, 0)) {
			scaffolding_failed = true;
			break;
		}
		scaffolded_seqs_alns.push_back(aln);
	}

	if (!scaffolding_failed)
	for (int i = 0; i < scaffolded_seqs_alns.size(); i++) {
		for (int j = i+1; j < scaffolded_seqs_alns.size(); j++) {
			if (std::max(scaffolded_seqs_alns[i].ref_begin, scaffolded_seqs_alns[j].ref_begin) <= std::min(scaffolded_seqs_alns[i].ref_end, scaffolded_seqs_alns[j].ref_end)) {
				scaffolding_failed = true;
				break;
			}
		}
		if (scaffolding_failed) break;
	}

	if (!scaffolding_failed) scaffolded_seqs_alns.swap(consensus_contigs_alns);
	return scaffolding_failed ? retained_assembled_sequences : scaffolded_sequences;
}

#endif // GUIDED_REFERENCE_ASSEMBLE_H
