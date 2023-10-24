#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <ctime>

#include "../libs/ssw_cpp.h"
#include "utils.h"
#include "vcf_utils.h"
#include "sam_utils.h"
#include "sw_utils.h"
#include "types.h"
#include "../libs/cptl_stl.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

std::mutex mtx;
config_t config;
std::string workspace;

chr_seqs_map_t contigs;

std::vector<sv_t*> svs;

int MAX_BP_DIST = 10;

struct clip_consensus_t {
    hts_pos_t pos;
    std::string seq;
    int fwd_clipped, rev_clipped;

    clip_consensus_t(hts_pos_t pos, std::string& seq, int fwd_clipped, int rev_clipped) :
    	pos(pos), seq(seq), fwd_clipped(fwd_clipped), rev_clipped(rev_clipped) {}
};

sv_t* detect_sv_from_junction(std::string& contig_name, std::string junction_seq, hts_pos_t ref_remap_start, hts_pos_t ref_remap_end,
                    StripedSmithWaterman::Aligner& aligner) {
    
    char* ref = contigs.get_seq(contig_name);
    hts_pos_t ref_remap_len = ref_remap_end - ref_remap_start;

    char ref_cstr[100000];
	for (int i = 0; i < ref_remap_len; i++) {
		ref_cstr[i] = toupper(ref[ref_remap_start+i]);
	} ref_cstr[ref_remap_len] = '\0';
    int* prefix_scores = smith_waterman_gotoh(ref_cstr, ref_remap_len, junction_seq.c_str(), junction_seq.length(), 1, -4, -6, -1);

    char ref_rev_cstr[100000];
    for (int i = 0; i < ref_remap_len; i++) {
        ref_rev_cstr[i] = ref_cstr[ref_remap_len-1-i];
    } ref_rev_cstr[ref_remap_len] = '\0';
    std::string junction_seq_rev = std::string(junction_seq.rbegin(), junction_seq.rend());
    int* suffix_scores = smith_waterman_gotoh(ref_rev_cstr, ref_remap_len, junction_seq_rev.c_str(), junction_seq_rev.length(), 1, -4, -6, -1);

    int max_score = 0, best_i = 0, best_j = 0;
	for (int i = config.min_clip_len; i < junction_seq.length()-config.min_clip_len; i++) {
        int prefix_score = prefix_scores[i-1]; // score of the best aln of [0..i-1]
        for (int j = i; j < junction_seq.length()-config.min_clip_len; j++) {
            int suffix_score = suffix_scores[junction_seq.length()-j-1]; // score of the best aln of [j..junction_seq.length()-1]
            // note that we want the score of the suffix of length junction_seq.length()-j, 
            // so we need to subtract 1 because suffix_scores[n] is the score of the best suffix of length n+1
            if (prefix_score + suffix_score > max_score) {
                max_score = prefix_score + suffix_score;
                best_i = i, best_j = j;
            }
        }
    }

    if (max_score == 0) return NULL;

    std::string left_part = junction_seq.substr(0, best_i);
    std::string middle_part = junction_seq.substr(best_i, best_j-best_i);
    std::string right_part = junction_seq.substr(best_j);

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment left_part_aln, right_part_aln;
    aligner.Align(left_part.c_str(), ref_cstr, ref_remap_len, StripedSmithWaterman::Filter(), &left_part_aln, 0);
    aligner.Align(right_part.c_str(), ref_cstr, ref_remap_len, StripedSmithWaterman::Filter(), &right_part_aln, 0);

    hts_pos_t left_bp = ref_remap_start + left_part_aln.ref_end;
    hts_pos_t right_bp = ref_remap_start + right_part_aln.ref_begin-1;
    
    if (left_bp > right_bp) { // there is microhomology in the inserted seq or it's a duplication
        int mh_len = left_bp - right_bp;
        std::pair<int, int> lp_suffix_score = find_aln_suffix_score(left_part_aln.cigar, mh_len, 1, -4, -6, -1);
        std::pair<int, int> rp_prefix_score = find_aln_prefix_score(right_part_aln.cigar, mh_len, 1, -4, -6, -1);
        
        if (right_part_aln.ref_end - left_part_aln.ref_begin < config.min_clip_len || 
            right_part_aln.ref_begin - left_part_aln.ref_begin < config.min_clip_len ||
            (lp_suffix_score.first == mh_len && rp_prefix_score.first == mh_len && middle_part.empty() &&
            !is_right_clipped(left_part_aln) && !is_left_clipped(right_part_aln))) { // it's a duplication
            duplication_t* sv = new duplication_t(contig_name, right_bp, left_bp, middle_part);
            sv->left_anchor = std::to_string(ref_remap_start + left_part_aln.ref_begin) + "-" + std::to_string(ref_remap_start + left_part_aln.ref_end); 
            sv->right_anchor = std::to_string(ref_remap_start + right_part_aln.ref_begin) + "-" + std::to_string(ref_remap_start + right_part_aln.ref_end);
            sv->left_anchor_cigar = left_part_aln.cigar_string;
            sv->right_anchor_cigar = right_part_aln.cigar_string;
            return sv;
        }

        // otherwise, it's an insertion
        std::string mh;
        if (lp_suffix_score.first > rp_prefix_score.first) { 
            // the suffix of the left part is more similar to the reference, hence we choose the prefix of the right part 
            //to add as part of the inserted sequence
            int right_bp_adjustment = 0;
            int query_mh_bases = query_prefix_len(right_part_aln.cigar, mh_len, right_bp_adjustment);
            mh = right_part.substr(0, query_mh_bases);
            right_bp = left_bp + right_bp_adjustment;
            middle_part = middle_part + mh;
        } else {
            // the prefix of the right part is more similar to the reference, hence we choose the suffix of the left part
            // to add as part of the inserted sequence
            int left_bp_adjustment = 0;
            int query_mh_bases = query_suffix_len(left_part_aln.cigar, mh_len, left_bp_adjustment);
            mh = left_part.substr(left_part.length() - query_mh_bases);
            left_bp = right_bp - left_bp_adjustment;
            middle_part = mh + middle_part;
        }
    }

    sv_t* sv = NULL;
    if (right_bp - left_bp > middle_part.length()) { // length of ALT < REF, deletion
        sv = new deletion_t(contig_name, left_bp, right_bp, middle_part);
    } else { // length of ALT > REF, insertion
        sv = new insertion_t(contig_name, left_bp, right_bp, 0, 0, middle_part);
    }
    sv->left_anchor = std::to_string(ref_remap_start + left_part_aln.ref_begin) + "-" + std::to_string(ref_remap_start + left_part_aln.ref_end); 
    sv->right_anchor = std::to_string(ref_remap_start + right_part_aln.ref_begin) + "-" + std::to_string(ref_remap_start + right_part_aln.ref_end);
    sv->left_anchor_cigar = left_part_aln.cigar_string;
    sv->right_anchor_cigar = right_part_aln.cigar_string;
    return sv;
}

sv_t* detect_sv(std::string& contig_name, clip_consensus_t& rc_consensus, clip_consensus_t& lc_consensus,
                           StripedSmithWaterman::Aligner& aligner) {
    suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus.seq, lc_consensus.seq, 1, -4, config.max_seq_error);
    if (spa.overlap < config.min_clip_len || is_homopolymer(lc_consensus.seq.c_str(), spa.overlap)) return NULL;

    int mm = spa.mismatches;
    int aln_len = spa.overlap;
    double mm_rate = double(mm)/aln_len;
    if (mm_rate > config.max_seq_error) return NULL;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln_rh, aln_lh;

    std::string consensus_junction_seq = rc_consensus.seq + lc_consensus.seq.substr(spa.overlap);
    hts_pos_t ref_remap_start = rc_consensus.pos - consensus_junction_seq.length(), 
              ref_remap_end = lc_consensus.pos + consensus_junction_seq.length();
    if (ref_remap_start < 0) ref_remap_start = 0;
    if (ref_remap_end > contigs.get_len(contig_name)) ref_remap_end = contigs.get_len(contig_name);

    sv_t* sv = detect_sv_from_junction(contig_name, consensus_junction_seq, ref_remap_start, ref_remap_end, aligner);
    if (sv == NULL) return NULL;
    sv->rc_fwd_reads = rc_consensus.fwd_clipped, sv->rc_rev_reads = rc_consensus.rev_clipped;
    sv->lc_fwd_reads = lc_consensus.fwd_clipped, sv->lc_rev_reads = lc_consensus.rev_clipped;
    sv->overlap = spa.overlap;
    return sv;
}

void call_insertions(int id, int contig_id, std::string contig_name) {
    mtx.lock();
    std::cout << "Calling insertions for " << contig_name << std::endl;
    mtx.unlock();

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);

    std::ifstream clip_fin(workspace + "/sr_consensuses/" + std::to_string(contig_id) + ".txt");
    std::vector<clip_consensus_t> rc_consensuses, lc_consensuses;
    std::string chr, dir, seq;
    hts_pos_t start, end, breakpoint;
    int fwd_clipped, rev_clipped;
    int max_mapq, lowq_clip_portion;
    hts_pos_t remap_boundary;
    while (clip_fin >> chr >> start >> end >> breakpoint >> dir >> seq >> 
        fwd_clipped >> rev_clipped >> max_mapq >> remap_boundary >> lowq_clip_portion) {
        if (dir == "R") {
            seq = seq.substr(0, seq.length()-lowq_clip_portion);
            rc_consensuses.push_back(clip_consensus_t(breakpoint, seq, fwd_clipped, rev_clipped));
        } else {
            seq = seq.substr(lowq_clip_portion);
            lc_consensuses.push_back(clip_consensus_t(breakpoint, seq, fwd_clipped, rev_clipped));
        }
    }

    int MAX_MH_LEN = 100;
    for (int i = 0; i < rc_consensuses.size(); i++) {
        clip_consensus_t& rc_consensus = rc_consensuses[i];
        for (int j = 0; j < lc_consensuses.size(); j++) {
            clip_consensus_t& lc_consensus = lc_consensuses[j];
            int mh_len = rc_consensus.pos - lc_consensus.pos;
            if (lc_consensus.pos-rc_consensus.pos <= MAX_BP_DIST && mh_len <= MAX_MH_LEN) {
                sv_t* sv = detect_sv(contig_name, rc_consensus, lc_consensus, aligner);
                if (sv == NULL) continue;

                mtx.lock();
                svs.push_back(sv);
                mtx.unlock();
            }
        }
    }
}

int main(int argc, char* argv[]) {

    std::string workdir = argv[1];
    workspace = workdir + "/workspace/";
    std::string reference_fname = argv[2];
    std::string sample_name = argv[3];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contigs.read_fasta_into_map(reference_fname);

    contig_map_t contig_map;
    contig_map.parse(workdir);
    config.parse(workdir + "/config.txt");

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool.push(call_insertions, contig_id, contig_name);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

    std::string out_vcf_fname = workdir + "/small_ins.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + out_vcf_fname + ".");
	}

	bcf_hdr_t* out_vcf_header = generate_vcf_header(contigs, sample_name, config, full_cmd_str);
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	std::sort(svs.begin(), svs.end(), [&out_vcf_header](sv_t* i1, sv_t* i2) {
		int contig_id1 = bcf_hdr_name2id(out_vcf_header, i1->chr.c_str());
		int contig_id2 = bcf_hdr_name2id(out_vcf_header, i2->chr.c_str());
		// negative because we want descending order
		int sc_score1 = -(i1->rc_reads()*i1->lc_reads()), sc_score2 = -(i2->rc_reads()*i2->lc_reads());
		int overlap1 = -i1->overlap, overlap2 = -i2->overlap;
		return std::tie(contig_id1, i1->start, i1->end, i1->ins_seq, sc_score2, overlap1) <
			   std::tie(contig_id2, i2->start, i2->end, i2->ins_seq, sc_score1, overlap2);
	});

    int del_id = 0, ins_id = 0, dup_id = 0;
    bcf1_t* bcf_entry = bcf_init();
    std::unordered_set<std::string> used_keys;
    for (sv_t* sv : svs) {
    	std::string key = sv->unique_key();
    	if (used_keys.count(key)) continue;
    	used_keys.insert(key);

        std::vector<std::string> filters;
        if (sv->svtype() == "DEL") {
            sv->id = sv->svtype() + "_" + std::to_string(del_id++);
            deletion_t* deletion = (deletion_t*) sv;
            del2bcf(out_vcf_header, bcf_entry, deletion, contigs.get_seq(deletion->chr), filters);
        } else if (sv->svtype() == "INS") {
            sv->id = sv->svtype() + "_" + std::to_string(ins_id++);
            insertion_t* insertion = (insertion_t*) sv;
            ins2bcf(out_vcf_header, bcf_entry, insertion, contigs.get_seq(insertion->chr), filters);
        } else if (sv->svtype() == "DUP") {
            sv->id = sv->svtype() + "_" + std::to_string(dup_id++);
            duplication_t* duplication = (duplication_t*) sv;
            dup2bcf(out_vcf_header, bcf_entry, duplication, contigs.get_seq(duplication->chr), filters);        
        }
        bcf_update_info_int32(out_vcf_header, bcf_entry, "OVERLAP", &sv->overlap, 1);
        bcf_update_info_string(out_vcf_header, bcf_entry, "ALGORITHM", "consensus_overlap");
        if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
            throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
        }
    }

    bcf_close(out_vcf_file);
}
