#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "htslib/faidx.h"

#include "sam_utils.h"
#include "utils.h"
#include "../libs/cptl_stl.h"

const int MIN_RND_POS = 1000;

std::string workdir;
std::string workspace;
std::string reference_fname;

config_t config;
stats_t stats;

std::mutex mtx;
std::mutex* mtx_contig;
std::vector<std::vector<std::string> > mate_seqs;

std::unordered_map<std::string, int> min_depth_by_contig, max_depth_by_contig, median_depth_by_contig;
std::vector<uint32_t> depths;
uint64_t qual_counts[256];
std::vector<std::vector<uint32_t> > dist_between_end_and_rnd;
std::vector<std::pair<uint64_t, uint32_t> > partial_sums;

std::vector<uint32_t> isize_counts;
std::vector<std::vector<uint32_t> > isizes_count_geq_i;

void categorize(int id, int contig_id, std::string contig_name, std::string bam_fname, std::string reference_fname, std::vector<hts_pos_t> rnd_positions) {
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

    samFile* sr_writer = NULL;
    samFile* hsr_writer = NULL;
    samFile* rdc_writer = NULL;
	samFile* ldc_writer = NULL;
    samFile* lp_writer = NULL;
    samFile* ow_writer = NULL;
    std::ofstream mateseqs_fout(workspace + "/sc_mateseqs/" + std::to_string(contig_id) + ".txt");

    int curr_pos = 0;
    std::vector<uint32_t> local_depths(rnd_positions.size());
    std::vector<std::vector<uint32_t> > rnd_positions_dist_between_end_and_rnd(rnd_positions.size());
    uint64_t contig_qual_counts[256];
	std::fill(contig_qual_counts, contig_qual_counts+256, uint64_t(0));
    std::vector<uint32_t> local_isize_counts(stats.max_is+1);
    std::vector<std::vector<uint32_t> > local_isize_dist_by_pos(rnd_positions.size());
    uint64_t sum_is = 0;
    uint32_t n_is = 0;

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (!is_primary(read)) continue;

        int64_t mq = get_mq(read);
        if (is_dc_pair(read)) {
            if (read->core.qual >= config.min_stable_mapq && read->core.qual >= mq && !is_unmapped(read)) { // stable end
                if (bam_is_rev(read) && !is_right_clipped(read, config.min_clip_len)) {
                    if (!ldc_writer) ldc_writer = open_writer(workspace + "/rev-stable/" + std::to_string(contig_id) + ".noremap.bam", bam_file->header);
                    int ok = sam_write1(ldc_writer, bam_file->header, read);
                    if (ok < 0) throw "Failed to write to " + std::string(ldc_writer->fn);
                } else if (!bam_is_rev(read) && !is_left_clipped(read, config.min_clip_len)) {
                    if (!rdc_writer) rdc_writer = open_writer(workspace + "/fwd-stable/" + std::to_string(contig_id) + ".noremap.bam", bam_file->header);
                    int ok = sam_write1(rdc_writer, bam_file->header, read);
                    if (ok < 0) throw "Failed to write to " + std::string(rdc_writer->fn);
                }
            }
            if (!is_mate_unmapped(read)) { // save read seq for remapping
                std::string qname = bam_get_qname(read), read_seq = get_sequence(read, true), qual_ascii = get_qual_ascii(read, true);
                if (is_samechr(read)) {
                    if (read->core.flag & BAM_FREAD1) qname += "_1";
                    else qname += "_2";
                }
                mtx_contig[read->core.mtid].lock();
                mate_seqs[read->core.mtid].push_back(qname + " " + read_seq + " " + qual_ascii + " " + std::to_string(read->core.qual));
                mtx_contig[read->core.mtid].unlock();
            }
        }
        if (is_samechr(read) && !is_samestr(read)) {
        	if (!bam_is_rev(read)) {
                if (read->core.isize > stats.max_is) {
                    if (!lp_writer) lp_writer = open_writer(workspace + "/long-pairs/" + std::to_string(contig_id) + ".bam", bam_file->header);

                    int ok = sam_write1(lp_writer, bam_file->header, read);
                    if (ok < 0) throw "Failed to write to " + std::string(lp_writer->fn);
                } else if (read->core.isize < 0) {
                    if (!ow_writer) ow_writer = open_writer(workspace + "/outward-pairs/" + std::to_string(contig_id) + ".bam", bam_file->header);

                    int ok = sam_write1(ow_writer, bam_file->header, read);
                    if (ok < 0) throw "Failed to write to " + std::string(ow_writer->fn);
                }
        	} else if (read->core.isize < -stats.max_is || read->core.isize > 0) {
        		mateseqs_fout << bam_get_qname(read) << " " << get_sequence(read) << "\n";
        	}
        }

        if (is_unmapped(read)) continue;

        while (curr_pos < rnd_positions.size() && read->core.pos > rnd_positions[curr_pos]) curr_pos++;

        // sample depth
		bool sampled = false;
		hts_pos_t read_endpos = bam_endpos(read);
		for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < read_endpos; i++) {
			if (read->core.pos <= rnd_positions[i] && rnd_positions[i] <= read_endpos) {
				sampled = true;
				local_depths[i]++;
			}
		}

        // sample pairs crossing
		hts_pos_t pair_startpos = read->core.pos + stats.read_len/2, pair_endpos = get_mate_endpos(read) - stats.read_len/2;
		if (read->core.isize > 0 && read->core.isize <= stats.max_is && is_samechr(read) && !is_samestr(read)
			&& !is_left_clipped(read, config.min_clip_len) && !is_right_clipped(read, config.min_clip_len)) {
			for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < pair_endpos; i++) {
				 if (pair_startpos <= rnd_positions[i] && rnd_positions[i] <= pair_endpos) {
					rnd_positions_dist_between_end_and_rnd[i].push_back(pair_endpos-rnd_positions[i]);
                    local_isize_counts[read->core.isize]++;
                    local_isize_dist_by_pos[curr_pos].push_back(read->core.isize);
                    sum_is += read->core.isize;
                    n_is++;
				}
			}
		}

        if (is_left_clipped(read, config.min_clip_len) || is_right_clipped(read, config.min_clip_len)) {
			if (!sr_writer) sr_writer = open_writer(workspace + "/clipped/" + std::to_string(contig_id) + ".bam", bam_file->header);

			int ok = sam_write1(sr_writer, bam_file->header, read);
			if (ok < 0) throw "Failed to write to " + std::string(sr_writer->fn);
		} else if (is_samechr(read) && is_hidden_split_read(read, config)) {
            if (!hsr_writer) hsr_writer = open_writer(workspace + "/hsr/" + std::to_string(contig_id) + ".bam", bam_file->header);

            int ok = sam_write1(hsr_writer, bam_file->header, read);
            if (ok < 0) throw "Failed to write to " + std::string(hsr_writer->fn);
            if (sampled) {
                contig_qual_counts[int(avg_qual(read)+0.5)]++;
            }
        } 
    }

    if (sr_writer) sam_close(sr_writer);
    if (rdc_writer) sam_close(rdc_writer);
    if (ldc_writer) sam_close(ldc_writer);
    if (hsr_writer) sam_close(hsr_writer);
    if (lp_writer) sam_close(lp_writer);
    if (ow_writer) sam_close(ow_writer);
    mateseqs_fout.close();

    local_depths.erase(std::remove(local_depths.begin(), local_depths.end(), 0), local_depths.end());
    mtx.lock();
    if (local_depths.size() >= MIN_RND_POS) {
		std::sort(local_depths.begin(), local_depths.end());
		min_depth_by_contig[contig_name] = local_depths[local_depths.size()/100];
		median_depth_by_contig[contig_name] = local_depths[local_depths.size()/2];
		max_depth_by_contig[contig_name] = local_depths[local_depths.size()-local_depths.size()/100];
    }
    depths.insert(depths.end(), local_depths.begin(), local_depths.end());

    partial_sums.push_back({sum_is, n_is});
	for (int i = 0; i < 256; i++) qual_counts[i] += contig_qual_counts[i];
	for (auto& v : rnd_positions_dist_between_end_and_rnd) {
		if (v.empty()) continue;
		dist_between_end_and_rnd.push_back(v);
	}

    for (int i = 0; i < rnd_positions.size(); i++) {
		std::vector<uint32_t> local_isizes_count_geq_i(stats.max_is+1); // position i contains the count of isizes that are geq than i
		for (uint32_t isize : local_isize_dist_by_pos[i]) local_isizes_count_geq_i[isize]++;
		for (int j = stats.max_is-1; j >= 0; j--) {
			local_isizes_count_geq_i[j] += local_isizes_count_geq_i[j+1];
		}
		if (local_isizes_count_geq_i[0] == 0) continue;
		for (int j = 0; j <= stats.max_is; j++) {
			isizes_count_geq_i[j].push_back(local_isizes_count_geq_i[j]);
		}
    }

    for (int i = 0; i <= stats.max_is; i++) {
    	isize_counts[i] += local_isize_counts[i];
    }
    mtx.unlock();

    bam_destroy1(read);
	hts_itr_destroy(iter);

    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];
    workdir = argv[2];
    workspace = workdir + "/workspace";
    reference_fname = argv[3];

    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contig_map_t contig_map(workdir);

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());
	if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + reference_fname;
	}

    mtx_contig = new std::mutex[contig_map.size()];
    mate_seqs.resize(contig_map.size());

    isize_counts.resize(stats.max_is+1);
    isizes_count_geq_i.resize(stats.max_is+1);

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
    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
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

    for (int i = 0; i < contig_map.size(); i++) {
    	if (mate_seqs[i].empty()) continue;
		std::string fname = std::to_string(i) + ".txt";
		std::ofstream mate_seqs_fout(workspace + "/mateseqs/" + fname);
		for (std::string& mate_seq : mate_seqs[i]) {
			mate_seqs_fout << mate_seq << std::endl;
		}
		mate_seqs_fout.close();
	}

    // get 1-percentile of base qualities
    uint64_t tot_quals = std::accumulate(qual_counts, qual_counts+256, uint64_t(0));
    uint64_t _1_perc_count = tot_quals/100;
    int min_avg_base_qual = 0;
    for (int i = 0; i < 256; i++) {
    	if (_1_perc_count <= qual_counts[i]) {
    		min_avg_base_qual = i;
    		break;
    	}
    	_1_perc_count -= qual_counts[i];
    }

    // TODO: double check the logic
    // for each insertion size, compute the estimate the minimum number of discordant pairs
    // in theory, it would be the number of discordant pairs *per breakpoint*, so we would need to require twice as many from the insertion
    // in practice, we do not because insertions often show a reduced number of pairs, perhaps due to overrepresentation of low complexity
    // subsequences in inserted sequences
    std::vector<std::vector<uint32_t>> pairs_crossing_dists(stats.max_is+1);
	for (int i = 0; i < dist_between_end_and_rnd.size(); i++) {
		std::vector<uint32_t> dist(stats.max_is+1);
		for (uint32_t val : dist_between_end_and_rnd[i]) {
			if (val <= stats.max_is) dist[val]++;
		}
		pairs_crossing_dists[0].push_back(dist[0]);
		for (int j = 1; j <= stats.max_is; j++) {
			dist[j] += dist[j-1];
			pairs_crossing_dists[j].push_back(dist[j]);
		}
	}

    std::ofstream mdpbs_fout(workdir + "/min_disc_pairs_by_size.txt");
	for (int i = 0; i <= stats.max_is; i++) {
		std::sort(pairs_crossing_dists[i].begin(), pairs_crossing_dists[i].end());
		mdpbs_fout << i << " " << pairs_crossing_dists[i][pairs_crossing_dists[i].size()/200] << std::endl;
	}
	mdpbs_fout.close();

    uint64_t sum_is = 0;
    uint32_t n_is = 0;
    for (std::pair<uint64_t, uint32_t>& p : partial_sums) {
        sum_is += p.first;
        n_is += p.second;
    }

    std::sort(depths.begin(), depths.end());

    std::ofstream stats_out(workdir + "/stats.txt", std::ios_base::app);
	std::sort(depths.begin(), depths.end());
	stats_out << "min_depth . " << depths[depths.size()/100] << std::endl;
	stats_out << "median_depth . " << depths[depths.size()/2] << std::endl;
	stats_out << "max_depth . " << depths[depths.size()-depths.size()/100] << std::endl;
    stats_out << "min_avg_base_qual . " << min_avg_base_qual << std::endl;
    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		if (!min_depth_by_contig.count(contig_name)) continue;
    	stats_out << "min_depth " << contig_name << " " << min_depth_by_contig[contig_name] << std::endl;
    	stats_out << "median_depth " << contig_name << " " << median_depth_by_contig[contig_name] << std::endl;
    	stats_out << "max_depth " << contig_name << " " << max_depth_by_contig[contig_name] << std::endl;
    }
    stats_out << "pop_avg_crossing_is . " << sum_is / n_is << std::endl;

    std::ofstream crossing_isizes_dist_fout(workdir + "/crossing_isizes.txt");
    for (int i = 0; i <= stats.max_is; i++) {
    	crossing_isizes_dist_fout << i << " " << isize_counts[i] << std::endl;
    }
    crossing_isizes_dist_fout.close();

    std::ofstream crossing_isizes_count_geq_i_fout(workdir + "/crossing_isizes_count_geq_i.txt");
    int mid = isizes_count_geq_i[0].size()/2;
	for (int i = 0; i <= stats.max_is; i++) {
		crossing_isizes_count_geq_i_fout << i << " ";
		std::sort(isizes_count_geq_i[i].begin(), isizes_count_geq_i[i].end());
		crossing_isizes_count_geq_i_fout << isizes_count_geq_i[i][mid] << std::endl;
	}
	crossing_isizes_dist_fout.close();

	stats_out << "min_disc_pairs . " << isizes_count_geq_i[0][isizes_count_geq_i[0].size()/100] << std::endl;
	stats_out.close();
}
