#include <string>

#include "../src/sw_utils.h"
#include "../src/sam_utils.h"
#include "../src/vcf_utils.h"
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "htslib/vcf.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;

StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

int main(int argc, char* argv[]) {

    std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];
    std::string workdir = argv[4];

    chr_seqs.read_fasta_into_map(reference_fname);

    contig_map_t contig_map;
    contig_map.load(workdir);

    config_t config;
    config.parse(workdir + "/config.txt");

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* b = bcf_init();

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln_before, aln_after;

    std::vector<sv_t*> svs;
    while (bcf_read(in_vcf_file, hdr, b) == 0) {
        sv_t* sv = bcf_to_sv(hdr, b);
        if (sv == NULL) continue;

        sv_t* new_dup = NULL;

        std::string ins_seq = sv->ins_seq;
        if (sv->svtype() == "INS" && !ins_seq.empty() && ins_seq.length() < 16000 && ins_seq.find('-') == std::string::npos) {
            std::string contig_name = sv->chr;
            hts_pos_t contig_len = chr_seqs.get_len(contig_name);

            bool skip = false;
            if (sv->start+1 < ins_seq.length()) skip = true;
            if (sv->start+ins_seq.length() >= contig_len) skip = true;

            if (!skip) {
                aligner.Align(ins_seq.c_str(), chr_seqs.get_seq(contig_name)+sv->start-ins_seq.length()+1, ins_seq.length(), filter, &aln_before, 0);
                aligner.Align(ins_seq.c_str(), chr_seqs.get_seq(contig_name)+sv->start+1, ins_seq.length(), filter, &aln_after, 0);

                std::vector<int> prefix_scores = compute_prefix_scores(aln_after.cigar, ins_seq.length(), 1, -4, -6, -1);
                std::vector<int> suffix_scores = compute_suffix_scores(aln_before.cigar, ins_seq.length(), 1, -4, -6, -1);

                int max_score = 0, best_i = 0;
                for (int i = 0; i <= ins_seq.length(); i++) {
                    if (prefix_scores[i] + suffix_scores[i] > max_score) {
                        max_score = prefix_scores[i] + suffix_scores[i];
                        best_i = i;
                    }
                }

                bool aln_alter_lc = get_left_clip_size(aln_after) >= config.min_clip_len;
                bool aln_before_rc = get_right_clip_size(aln_before) >= config.min_clip_len;
                int aln_len = 0;
                if (!aln_alter_lc) aln_len += aln_after.query_end+aln_after.query_begin+1;
                if (!aln_before_rc) aln_len += aln_before.query_end-aln_before.query_begin+1;
                if (aln_len >= ins_seq.length()) {
                    sv_t::anchor_aln_t* anchor_aln = new sv_t::anchor_aln_t(sv->start-ins_seq.length()+best_i, sv->start+best_i, ins_seq.length(), 0);
                    new_dup = new duplication_t(contig_name, anchor_aln->start, anchor_aln->end, "", NULL, NULL, anchor_aln, anchor_aln);
                    new_dup->id = sv->id + "_DUP";
                    new_dup->source = sv->source;
                }
            }
        }

        sv->n_gt = 0;
        delete[] sv->sample_info.gt;
        sv->sample_info.gt = NULL;
        svs.push_back(sv);
        if (new_dup != NULL) {
            new_dup->n_gt = 0;
            delete[] new_dup->sample_info.gt;
            new_dup->sample_info.gt = NULL;
            svs.push_back(new_dup);
        }
    }

    std::sort(svs.begin(), svs.end(), [&contig_map](sv_t* a, sv_t* b) { 
        size_t chr_a = contig_map.get_id(a->chr);
        size_t chr_b = contig_map.get_id(b->chr);
        return std::tie(chr_a, a->start) < std::tie(chr_b, b->start);
    });
    bcf_hdr_t* out_hdr = generate_vcf_header(chr_seqs, "", config, "");

    if (bcf_hdr_set_samples(out_hdr, NULL, 0) != 0) {
        throw std::runtime_error("Failed to unset samples in VCF header");
    }

    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "w");
    if (bcf_hdr_write(out_vcf_file, out_hdr) != 0) {
        throw std::runtime_error("Failed to write VCF header to " +  out_vcf_fname);
    }
    for (sv_t* sv : svs) {
        sv2bcf(out_hdr, b, sv, chr_seqs.get_seq(sv->chr), true);
        if (bcf_write(out_vcf_file, out_hdr, b) != 0) {
            throw std::runtime_error("Failed to write to " + out_vcf_fname);
        }
    }

    bcf_destroy(b);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(out_hdr);
    bcf_close(in_vcf_file);
    bcf_close(out_vcf_file);
}
