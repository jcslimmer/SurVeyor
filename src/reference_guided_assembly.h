#ifndef GUIDED_REFERENCE_ASSEMBLE_H
#define GUIDED_REFERENCE_ASSEMBLE_H

#include <iostream>

#include "../libs/ssw_cpp.h"
#include "utils.h"
#include "sam_utils.h"
#include "dc_remapper.h"
#include "assemble.h"
#include "remapping.h"

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

std::string generate_consensus_sequences(std::string contig_name, chr_seqs_map_t& contigs, contig_map_t& contig_map, insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster, 
        std::vector<remap_info_t>& rc_remap_infos, std::vector<remap_info_t>& lc_remap_infos, region_t& best_region, bool is_rc, 
		bool& left_bp_precise, bool& right_bp_precise, std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, stats_t& stats) {

	left_bp_precise = false, right_bp_precise = false;

	if (best_region.score.remap_start >= best_region.score.remap_end) return "";

	insertion_cluster_t* refined_r_cluster = new insertion_cluster_t(new cluster_t());
	insertion_cluster_t* refined_l_cluster = new insertion_cluster_t(new cluster_t());
	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
		if (rc_remap_infos[i].accepted) {
			refined_r_cluster->add_stable_read(r_cluster->cluster->reads[i]);
		}
	}
	for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
		if (lc_remap_infos[i].accepted) {
			refined_l_cluster->add_stable_read(l_cluster->cluster->reads[i]);
		}
	}
	if (r_cluster->clip_consensus && rc_remap_infos.rbegin()->accepted) {
		refined_r_cluster->add_clip_cluster(r_cluster->clip_consensus);
	}
	if (l_cluster->clip_consensus && lc_remap_infos.rbegin()->accepted) {
		refined_l_cluster->add_clip_cluster(l_cluster->clip_consensus);
	}

	if (refined_r_cluster->empty() || refined_l_cluster->empty()) return "";
	if (refined_r_cluster->end-refined_r_cluster->start <= 0 || refined_l_cluster->end-refined_l_cluster->start <= 0) return "";

    // build reference-guide as left flanking + predicted transposed region + right flanking
	char* left_flanking = new char[refined_r_cluster->end-refined_r_cluster->start+2];
	strncpy(left_flanking, contigs.get_seq(contig_name)+refined_r_cluster->start, refined_r_cluster->end-refined_r_cluster->start+1);
	left_flanking[refined_r_cluster->end-refined_r_cluster->start+1] = '\0';
	to_uppercase(left_flanking);

	int ins_seq_start = best_region.start + best_region.score.remap_start,
		ins_seq_len = best_region.score.remap_end - best_region.score.remap_start + 1;
	char* pred_ins_seq = new char[ins_seq_len+1];
	strncpy(pred_ins_seq, contigs.get_seq(contig_map.get_name(best_region.contig_id))+ins_seq_start, ins_seq_len);
	pred_ins_seq[ins_seq_len] = '\0';
	if (is_rc) rc(pred_ins_seq);
	to_uppercase(pred_ins_seq);

	char* right_flanking = new char[refined_l_cluster->end-refined_l_cluster->start+2];
	strncpy(right_flanking, contigs.get_seq(contig_name)+refined_l_cluster->start, refined_l_cluster->end-refined_l_cluster->start+1);
	right_flanking[refined_l_cluster->end-refined_l_cluster->start+1] = '\0';
	to_uppercase(right_flanking);

	// NOTE: prefixes and suffixes may be duplicated - i.e., suffix of left flanking also appearing as prefix of predicted ins seq
	// This is a problem because it predicts an insertion where there is none. One such case is fragment LINE insertions.
	// There are small regions in LINE/L1P* that accumulate a lot of mutations compared to reference, and it may be predicted as a block-swap
	// (which is not an insertion as the length of alt allele is not changed)
	// However, because of the duplicated suffix/prefix, a small insertion is predicted
	std::string left_flanking_str = left_flanking, pred_ins_seq_str = pred_ins_seq, right_flanking_str = right_flanking;

	suffix_prefix_aln_t spa12 = aln_suffix_prefix(left_flanking_str, pred_ins_seq_str, 1, -4, 1.0, config.min_clip_len);
	int overlap12 = spa12.overlap;
	if (!overlap12) {
		StripedSmithWaterman::Alignment aln;
		StripedSmithWaterman::Filter filter;
		aligner.Align(pred_ins_seq_str.c_str(), left_flanking_str.c_str(), left_flanking_str.length(), filter, &aln, 0);
		if (aln.query_begin == 0 && aln.ref_end == left_flanking_str.length()-1) overlap12 = aln.query_end;
	}

	suffix_prefix_aln_t spa23 = aln_suffix_prefix(pred_ins_seq_str, right_flanking_str, 1, -4, 1.0, config.min_clip_len);
	int overlap23 = spa23.overlap;
	if (!overlap23) {
		StripedSmithWaterman::Alignment aln;
		StripedSmithWaterman::Filter filter;
		aligner.Align(right_flanking_str.c_str(), pred_ins_seq_str.c_str(), pred_ins_seq_str.length(), filter, &aln, 0);
		if (aln.query_begin == 0 && aln.ref_end == pred_ins_seq_str.length()-1) overlap23 = aln.query_end;
	}

	std::string full_junction_sequence = left_flanking_str;
	full_junction_sequence += pred_ins_seq_str.substr(overlap12);
	full_junction_sequence += right_flanking_str.substr(overlap23);

	delete[] left_flanking;
	delete[] pred_ins_seq;
	delete[] right_flanking;


	// assembled contigs
	std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
	std::vector<std::string> consensus_contigs = generate_reference_guided_consensus(full_junction_sequence, r_cluster, l_cluster, mateseqs,
			aligner, harsh_aligner, consensus_contigs_alns, config, stats);

	if (consensus_contigs.empty()) return full_junction_sequence;

	int l_bp_in_seq = refined_r_cluster->end-refined_r_cluster->start, r_bp_in_seq = full_junction_sequence.length()-(refined_l_cluster->end-refined_l_cluster->start);

	const int OFFSET = 50;
	for (StripedSmithWaterman::Alignment& aln : consensus_contigs_alns) {
		left_bp_precise |= (aln.ref_begin < l_bp_in_seq-OFFSET && aln.ref_end >= l_bp_in_seq+OFFSET);
		right_bp_precise |= (aln.ref_begin < r_bp_in_seq-OFFSET && aln.ref_end >= r_bp_in_seq+OFFSET);
	}

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
	for (int i = 0; i < consensus_contigs.size(); i++) {
		contigs_sorted_by_pos.push_back({consensus_contigs[i], consensus_contigs_alns[i]});
	}
	std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	// complement holes in contigs with reference
	std::string corrected_junction_seq = full_junction_sequence.substr(0, contigs_sorted_by_pos[0].second.ref_begin);
	corrected_junction_seq += contigs_sorted_by_pos[0].first;
	for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
		corrected_junction_seq += full_junction_sequence.substr(contigs_sorted_by_pos[i-1].second.ref_end+1, contigs_sorted_by_pos[i].second.ref_begin-contigs_sorted_by_pos[i-1].second.ref_end);
		corrected_junction_seq += contigs_sorted_by_pos[i].first;
	}
	corrected_junction_seq += full_junction_sequence.substr(contigs_sorted_by_pos[contigs_sorted_by_pos.size()-1].second.ref_end+1);

	StripedSmithWaterman::Alignment aln;
	StripedSmithWaterman::Filter filter;
	if (r_cluster->clip_consensus) {
		harsh_aligner.Align(r_cluster->clip_consensus->sequence.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error);
		rc_remap_infos.back() = remap_info_t(aln, accepted, config.min_clip_len);
	}
	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
		bam1_t* read = r_cluster->cluster->reads[i];
		std::string read_seq = get_sequence(read);
		std::string read_qual = get_qual_ascii(read);
		std::string mateseq = get_mate_seq(read, mateseqs);
		std::string matequal = get_mate_qual(read, matequals);
		rc(mateseq);
		matequal = std::string(matequal.rbegin(), matequal.rend());

		harsh_aligner.Align(read_seq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error, read_qual, stats.min_avg_base_qual);
		harsh_aligner.Align(mateseq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		accepted &= accept(aln, config.min_clip_len, config.max_seq_error, matequal, stats.min_avg_base_qual);
		rc_remap_infos[i] = remap_info_t(aln, accepted, config.min_clip_len);
	}
	for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
		bam1_t* read = l_cluster->cluster->reads[i];
		std::string read_seq = get_sequence(read);
		std::string read_qual = get_qual_ascii(read);
		std::string mateseq = get_mate_seq(read, mateseqs);
		std::string matequal = get_mate_qual(read, matequals);

		harsh_aligner.Align(read_seq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error, read_qual, stats.min_avg_base_qual);
		harsh_aligner.Align(mateseq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		accepted &= accept(aln, config.min_clip_len, config.max_seq_error, matequal, stats.min_avg_base_qual);
		lc_remap_infos[i] = remap_info_t(aln, accepted, config.min_clip_len);
	}
	if (l_cluster->clip_consensus) {
		harsh_aligner.Align(l_cluster->clip_consensus->sequence.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error);
		lc_remap_infos.back() = remap_info_t(aln, accepted, config.min_clip_len);
	}

	return corrected_junction_seq;
}


void update_read(bam1_t* read, region_t& chosen_region, remap_info_t& remap_info, bool region_is_rc) {
    read->core.mtid = chosen_region.original_bam_id;
    read->core.mpos = chosen_region.start + remap_info.start;
    if (region_is_rc == bam_is_rev(read)) {
        read->core.flag |= BAM_FMREVERSE; //sets flag to true
    } else {
        read->core.flag &= ~BAM_FMREVERSE; //sets flag to false
    }
    bam_aux_update_str(read, "MC", remap_info.cigar.length() + 1, remap_info.cigar.c_str());
    char ok = remap_info.accepted ? 'T' : 'F';
    bam_aux_append(read, "OK", 'A', 1, (const uint8_t*) &ok);
    bam_aux_append(read, "MS", 'i', 4, (const uint8_t*) &remap_info.score);
}

insertion_t* detect_reference_guided_assembly_insertion(std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::string& junction_seq,
		insertion_cluster_t* r_cluster, insertion_cluster_t* l_cluster, std::vector<remap_info_t>& ro_remap_infos, std::vector<remap_info_t>& lo_remap_infos,
		region_t& best_region, bool is_rc, std::vector<bam1_t*>& kept, bool left_bp_precise, bool right_bp_precise, StripedSmithWaterman::Aligner& aligner, config_t& config) {

	int remap_start = std::max(hts_pos_t(0), r_cluster->start-50);
	int remap_end = std::min(l_cluster->end+50, contig_len-1);
	if (remap_start >= remap_end) return NULL;

	 std::vector<sv_t*> insertions = detect_svs_from_junction(contig_name, contig_seq, junction_seq, 
                remap_start, remap_end, remap_start, remap_end, aligner, config.min_clip_len);

	if (insertions.empty() || insertions[0]->svtype() != "INS" || insertions[0]->svlen() < config.min_sv_size) return NULL;

	insertion_t* insertion = (insertion_t*) insertions[0];

	insertion_cluster_t* refined_r_cluster = new insertion_cluster_t(new cluster_t());
    for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
        bam1_t* r = r_cluster->cluster->reads[i];
        if (ro_remap_infos[i].accepted) {
            refined_r_cluster->add_stable_read(r);
        }
    }
	if (r_cluster->clip_consensus && ro_remap_infos.rbegin()->accepted) {
        refined_r_cluster->add_clip_cluster(r_cluster->clip_consensus);
    }

	insertion_cluster_t* refined_l_cluster = new insertion_cluster_t(new cluster_t());
    for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
        bam1_t* r = l_cluster->cluster->reads[i];
        if (lo_remap_infos[i].accepted) {
            refined_l_cluster->add_stable_read(r);
        }
    }
    if (l_cluster->clip_consensus && lo_remap_infos.rbegin()->accepted) {
        refined_l_cluster->add_clip_cluster(l_cluster->clip_consensus);
    }

	if (refined_r_cluster->empty() || refined_l_cluster->empty()) {
		delete refined_r_cluster;
		delete refined_l_cluster;
		return NULL;
	}

	// Get transposed sequence coverage
	std::vector<std::pair<int,int>> covered_segments;
	for (remap_info_t ri : ro_remap_infos) {
		if (ri.accepted) covered_segments.push_back({ri.start, ri.end});
	}
	for (remap_info_t ri : lo_remap_infos) {
		if (ri.accepted) covered_segments.push_back({ri.start, ri.end});
	}

	int ins_seq_start = insertion->left_anchor_aln->seq_len - insertion->prefix_mh_len;
	int ins_seq_end = ins_seq_start + insertion->ins_seq.length();

	int i = 0;
	std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.first < p2.first;});
	while (i < covered_segments.size() && covered_segments[i].second <= ins_seq_start) i++; // find left-most segment that either overlaps or is fully contained in the inserted sequence
	insertion->prefix_cov_start = 0, insertion->prefix_cov_end = 0;
	if (i < covered_segments.size()) {
		insertion->prefix_cov_start = std::max(ins_seq_start, covered_segments[i].first);
		insertion->prefix_cov_end = covered_segments[i].second;
		while (i < covered_segments.size() && covered_segments[i].first <= insertion->prefix_cov_end) {
			insertion->prefix_cov_end = std::max(insertion->prefix_cov_end, covered_segments[i].second);
			i++;
		}
		insertion->prefix_cov_start -= ins_seq_start;
		insertion->prefix_cov_end -= ins_seq_start;
		if (insertion->prefix_cov_end >= insertion->ins_seq.length()) insertion->prefix_cov_end = insertion->ins_seq.length()-1;
	}

	i = covered_segments.size()-1;
	std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.second < p2.second;});
	while (i >= 0 && covered_segments[i].first >= ins_seq_end) i--; // find right-most segment that either overlaps or is fully contained in the inserted sequence
	insertion->suffix_cov_start = 0, insertion->suffix_cov_end = 0;
	if (i >= 0) {
		insertion->suffix_cov_end = std::min(ins_seq_end-1, covered_segments[i].second);
		insertion->suffix_cov_start = covered_segments[i].first;
		while (i >= 0 && covered_segments[i].second >= insertion->suffix_cov_start) {
			insertion->suffix_cov_start = std::min(insertion->suffix_cov_start, covered_segments[i].first);
			i--;
		}
		insertion->suffix_cov_start -= ins_seq_start;
		if (insertion->suffix_cov_start < 0) insertion->suffix_cov_start = 0;
		insertion->suffix_cov_end -= ins_seq_start;
	}

	if (refined_r_cluster->clip_consensus) insertion->rc_consensus = refined_r_cluster->clip_consensus;
	if (refined_l_cluster->clip_consensus) insertion->lc_consensus = refined_l_cluster->clip_consensus;
	insertion->disc_pairs_lf = refined_r_cluster->cluster->reads.size();
	insertion->disc_pairs_rf = refined_l_cluster->cluster->reads.size();
	for (bam1_t* read : refined_r_cluster->cluster->reads) {
		if (read->core.qual >= config.high_confidence_mapq) insertion->disc_pairs_lf_high_mapq++;
	}
	for (bam1_t* read : refined_l_cluster->cluster->reads) {
		if (read->core.qual >= config.high_confidence_mapq) insertion->disc_pairs_rf_high_mapq++;
	}
	insertion->disc_pairs_lf_maxmapq = refined_r_cluster->cluster->max_mapq;
	insertion->disc_pairs_rf_maxmapq = refined_l_cluster->cluster->max_mapq;
	for (bam1_t* read : refined_r_cluster->cluster->reads) {
		insertion->disc_pairs_lf_avg_nm += get_nm(read);
	}
	insertion->disc_pairs_lf_avg_nm /= std::max(1, (int) refined_r_cluster->cluster->reads.size());
	for (bam1_t* read : refined_l_cluster->cluster->reads) {
		insertion->disc_pairs_rf_avg_nm += get_nm(read);
	}
	insertion->disc_pairs_rf_avg_nm /= std::max(1, (int) refined_l_cluster->cluster->reads.size());
	insertion->imprecise_bp = !left_bp_precise || !right_bp_precise;
	insertion->source = "REFERENCE_GUIDED_ASSEMBLY";

	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
		bam1_t* r = r_cluster->cluster->reads[i];
		update_read(r, best_region, ro_remap_infos[i], is_rc);
		kept.push_back(r);
	}
	for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
		bam1_t* r = l_cluster->cluster->reads[i];
		update_read(r, best_region, lo_remap_infos[i], is_rc);
		kept.push_back(r);
	}

	delete refined_r_cluster;
	delete refined_l_cluster;

	return insertion;
}

#endif // GUIDED_REFERENCE_ASSEMBLE_H
