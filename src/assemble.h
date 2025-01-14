#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_

#include <unordered_set>
#include <mutex>
#include <htslib/sam.h>

#include "sw_utils.h"
#include "types.h"
#include "dc_remapper.h"

std::mutex failed_assembly_mtx;

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

void correct_contig(std::string& contig, std::vector<std::string>& reads, std::vector<StripedSmithWaterman::Alignment>& alns, config_t& config,
	bool ambiguous_as_N = false) {

	std::vector<int> As(contig.length()), Cs(contig.length()), Gs(contig.length()), Ts(contig.length());
	for (int i = 0; i < reads.size(); i++) {
		StripedSmithWaterman::Alignment& aln = alns[i];
		if (accept(aln, config.min_clip_len)) {
			for (int j = aln.query_begin; j < aln.query_end; j++) {
				char c = reads[i][j];
				if (c == 'A') As[j-aln.query_begin+aln.ref_begin]++;
				else if (c == 'C') Cs[j-aln.query_begin+aln.ref_begin]++;
				else if (c == 'G') Gs[j-aln.query_begin+aln.ref_begin]++;
				else if (c == 'T') Ts[j-aln.query_begin+aln.ref_begin]++;
			}
		}
	}

	for (int i = 0; i < contig.length(); i++) {
		int max_freq = max(As[i], Cs[i], Gs[i], Ts[i]);
		if (max_freq == 0) continue;

		int count_max = (As[i] == max_freq) + (Cs[i] == max_freq) + (Gs[i] == max_freq) + (Ts[i] == max_freq);
		if (ambiguous_as_N && count_max > 1) contig[i] = 'N';
		else if (max_freq == As[i]) contig[i] = 'A';
		else if (max_freq == Cs[i]) contig[i] = 'C';
		else if (max_freq == Gs[i]) contig[i] = 'G';
		else if (max_freq == Ts[i]) contig[i] = 'T';
	}
}

void correct_contig(std::string& contig, std::vector<std::string>& reads, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config,
	bool ambiguous_as_N = false) {

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::vector<StripedSmithWaterman::Alignment> alns;
	for (std::string& read : reads) {
		harsh_aligner.Align(read.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
		alns.push_back(aln);
	}
	correct_contig(contig, reads, alns, config, ambiguous_as_N);
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
		build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 0, config.min_clip_len);

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

	const int TOO_MANY_READS = 1000;
	if (unstable_read_seqs.size() + left_stable_read_seqs.size() + right_stable_read_seqs.size() >= TOO_MANY_READS) return {"TOO_MANY_READS"};

	std::vector<std::string> assembled_sequences = assemble_reads(left_stable_read_seqs, unstable_read_seqs, right_stable_read_seqs,
			harsh_aligner, config, stats);

	return assembled_sequences;
}

sv_t* detect_de_novo_insertion(std::string& contig_name, chr_seqs_map_t& contigs,
		insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		std::ofstream& assembly_failed_no_seq, std::ofstream& assembly_failed_cycle_writer, std::ofstream& assembly_failed_too_many_reads_writer,
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
			chosen_ins->sample_info.alt_bp1.supp_pairs++;
			if (chosen_ins->disc_pairs_lf_maxmapq < read->core.qual) chosen_ins->disc_pairs_lf_maxmapq = read->core.qual;
			if (read->core.qual >= config.high_confidence_mapq) chosen_ins->sample_info.alt_bp1.supp_pairs_high_mapq++;
			chosen_ins->disc_pairs_lf_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
			bam1_t* d = bam_dup1(read);
			d->core.mpos = 0;
			d->core.mtid = d->core.tid;
			d->core.flag |= BAM_FMUNMAP;
			assembled_reads.push_back(d);
		}
	}
	chosen_ins->disc_pairs_lf_avg_nm /= std::max(1, chosen_ins->sample_info.alt_bp1.supp_pairs);
	for (bam1_t* read : l_cluster->cluster->reads) {
		std::string mate_seq = get_mate_seq(read, mateseqs);
		std::string mate_qual = get_mate_qual(read, matequals);
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(mate_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len, config.max_seq_error, mate_qual, stats.min_avg_base_qual)) {
			chosen_ins->sample_info.alt_bp2.supp_pairs++;
			if (chosen_ins->disc_pairs_rf_maxmapq < read->core.qual) chosen_ins->disc_pairs_rf_maxmapq = read->core.qual;
			if (read->core.qual >= config.high_confidence_mapq) chosen_ins->sample_info.alt_bp2.supp_pairs_high_mapq++;
			chosen_ins->disc_pairs_rf_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
			bam1_t* d = bam_dup1(read);
			d->core.mpos = 0;
			d->core.mtid = d->core.tid;
			d->core.flag |= BAM_FMUNMAP;
			assembled_reads.push_back(d);
		}
	}
	chosen_ins->disc_pairs_rf_avg_nm /= std::max(1, chosen_ins->sample_info.alt_bp2.supp_pairs);
	if (l_cluster->clip_consensus) {
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(l_cluster->clip_consensus->sequence.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.min_clip_len, config.max_seq_error)) {
			chosen_ins->lc_consensus = l_cluster->clip_consensus;
		}
	}

	// start and end of inserted sequence within the full assembled sequence
	int ins_seq_start = chosen_ins->left_anchor_aln->seq_len - chosen_ins->mh_len;
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
			chosen_ins->sample_info.alt_bp1.supp_pairs++;
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
			chosen_ins->sample_info.alt_bp2.supp_pairs++;
			assembled_reads.push_back(bam_dup1(read));
		}
	}

	for (bam1_t* read : assembled_reads) {
		bam_aux_update_str(read, "ID", chosen_ins->id.length(), chosen_ins->id.c_str());
	}

	chosen_ins->source = "DE_NOVO_ASSEMBLY";
	chosen_ins->mh_len = 0; // current value is just a temporary approximation, reset it. genotype will calculate it correctly
	return chosen_ins;
}

struct base_score_t {
    int freq = 0, qual = 0;
    char base;

    base_score_t(char base) : base(base) {}
};
bool operator < (const base_score_t& bs1, const base_score_t& bs2) {
    if (bs1.freq != bs2.freq) return bs1.freq < bs2.freq;
    return bs1.qual < bs2.qual;
}

std::string build_full_consensus_seq(std::vector<std::string>& seqs,
    std::vector<uint8_t*>& quals, std::vector<hts_pos_t>& read_start_offsets) {
    const int MAX_CONSENSUS_LEN = 100000;
    char consensus[MAX_CONSENSUS_LEN];

    std::vector<hts_pos_t> read_end_offsets;
    hts_pos_t consensus_len = 0;
    for (int i = 0; i < seqs.size(); i++) {
        hts_pos_t start_offset = read_start_offsets[i];
        hts_pos_t end_offset = start_offset + seqs[i].length() - 1;
        read_end_offsets.push_back(end_offset);
        consensus_len = std::max(consensus_len, end_offset+1);
    }

    int s = 0;
    for (int i = 0; i < consensus_len; i++) {
        while (s < seqs.size() && read_end_offsets[s] < i) s++;

        base_score_t base_scores[4] = { base_score_t('A'), base_score_t('C'), base_score_t('G'), base_score_t('T') };
        for (int j = s; j < seqs.size() && read_start_offsets[j] <= i; j++) {
            if (read_end_offsets[j] < i) continue;

            char nucl = seqs[j][i - read_start_offsets[j]];
            uint8_t qual = quals[j][i - read_start_offsets[j]];
            if (nucl == 'A') {
                base_scores[0].freq++;
                base_scores[0].qual += qual;
            } else if (nucl == 'C') {
                base_scores[1].freq++;
                base_scores[1].qual += qual;
            } else if (nucl == 'G') {
                base_scores[2].freq++;
                base_scores[2].qual += qual;
            } else if (nucl == 'T') {
                base_scores[3].freq++;
                base_scores[3].qual += qual;
            }
        }

        base_score_t best_base_score = max(base_scores[0], base_scores[1], base_scores[2], base_scores[3]);

        consensus[i] = best_base_score.base;
    }
    consensus[consensus_len] = '\0';

    return std::string(consensus);
}


#endif /* ASSEMBLE_H_ */
