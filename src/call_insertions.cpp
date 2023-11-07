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

void call_insertions(int id, int contig_id, std::string contig_name) {
    mtx.lock();
    std::cout << "Calling insertions for " << contig_name << std::endl;
    mtx.unlock();

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);

    std::ifstream clip_fin(workspace + "/sr_consensuses/" + std::to_string(contig_id) + ".txt");
    std::vector<consensus_t> rc_consensuses, lc_consensuses;
    std::string chr, dir, seq;
    hts_pos_t start, end, breakpoint;
    int fwd_clipped, rev_clipped;
    int max_mapq, lowq_clip_portion;
    hts_pos_t remap_boundary;
    while (clip_fin >> chr >> start >> end >> breakpoint >> dir >> seq >> 
        fwd_clipped >> rev_clipped >> max_mapq >> remap_boundary >> lowq_clip_portion) {
        if (dir == "R") {
            seq = seq.substr(0, seq.length()-lowq_clip_portion);
            rc_consensuses.push_back(consensus_t(false, start, breakpoint, end, seq, fwd_clipped, rev_clipped, end-breakpoint, max_mapq, remap_boundary, lowq_clip_portion));
        } else {
            seq = seq.substr(lowq_clip_portion);
            lc_consensuses.push_back(consensus_t(true, start, breakpoint, end, seq, fwd_clipped, rev_clipped, breakpoint-start, max_mapq, remap_boundary, lowq_clip_portion));
        }
    }

    int MAX_MH_LEN = 100;
    for (int i = 0; i < rc_consensuses.size(); i++) {
        consensus_t& rc_consensus = rc_consensuses[i];
        for (int j = 0; j < lc_consensuses.size(); j++) {
            consensus_t& lc_consensus = lc_consensuses[j];
            int mh_len = rc_consensus.breakpoint - lc_consensus.breakpoint;
            if (lc_consensus.breakpoint-rc_consensus.breakpoint <= MAX_BP_DIST && mh_len <= MAX_MH_LEN) {
                sv_t* sv = detect_sv(contig_name, contigs.get_seq(contig_name), contigs.get_len(contig_name), rc_consensus, lc_consensus, aligner, config.min_clip_len, config.min_clip_len, config.max_seq_error);
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
