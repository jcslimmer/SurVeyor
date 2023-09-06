#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <algorithm>

#include "htslib/faidx.h"

#include "sam_utils.h"
#include "utils.h"
#include "../libs/cptl_stl.h"

std::string workdir;
std::string reference_fname;

config_t config;

std::mutex mtx;
std::mutex* mtx_contig;
std::vector<std::vector<std::string> > mate_seqs;

void categorize(int id, int contig_id, std::string contig_name, std::string bam_fname, std::string reference_fname, std::vector<hts_pos_t> rnd_positions) {
    mtx.lock();
    std::cerr << "Categorizing reads in " << contig_name << std::endl;
    mtx.unlock();

    std::sort(rnd_positions.begin(), rnd_positions.end());

    open_samFile_t* bam_file = open_samFile(bam_fname);
    if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
        throw "Failed to read reference " + reference_fname;
    }

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    if (iter == NULL) { // no reads
    	close_samFile(bam_file);
    	return;
    }

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (!is_primary(read)) continue;

        int64_t mq = get_mq(read);
        if (is_dc_pair(read) && (read->core.qual >= config.min_stable_mapq || mq >= config.min_stable_mapq || is_unmapped(read))) {
            if (bam_is_rev(read) && !is_right_clipped(read, config.min_clip_len)) {
                // if (!ldc_writer) ldc_writer = get_writer(workspace + "/L", std::to_string(contig_id) + ".noremap.bam", bam_file->header);
                // int ok = sam_write1(ldc_writer, bam_file->header, read);
                // if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
            } else if (!bam_is_rev(read) && !is_left_clipped(read, config.min_clip_len)) {
                // if (!rdc_writer) rdc_writer = get_writer(workspace + "/R", std::to_string(contig_id) + ".noremap.bam", bam_file->header);
                // int ok = sam_write1(rdc_writer, bam_file->header, read);
                // if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
            }
        }
    }
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    reference_fname = argv[3];

    config.parse(workdir + "/config.txt");

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());
	if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + reference_fname;
	}
	bam_hdr_t* header = bam_file->header;

    mtx_contig = new std::mutex[header->n_targets];
    mate_seqs.resize(header->n_targets);

    // read random positions
	std::string contig_name;
	std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
	hts_pos_t pos;
	std::unordered_map<std::string, std::vector<hts_pos_t> > rnd_pos_map;
	while (rnd_pos_fin >> contig_name >> pos) {
		rnd_pos_map[contig_name].push_back(pos);
	}

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (int contig_id = 0; contig_id < header->n_targets; contig_id++) {
        std::string contig_name = header->target_name[contig_id];
        std::future<void> future = thread_pool.push(categorize, contig_id, contig_name, bam_fname, reference_fname, rnd_pos_map[contig_name]);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
        }
    }
}