#include <string>
#include <algorithm>
#include <vector>
#include <memory>
#include <mutex>

#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "../src/sam_utils.h"
#include "../src/vcf_utils.h"
#include "../libs/cptl_stl.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;

std::vector<bcf1_t*> normalised_vcf_records;
std::mutex mtx;

void atomize_del(std::shared_ptr<sv_t> sv) {
	if (sv->ins_seq.empty()) return;

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	char* deleted_seq = new char[sv->end - sv->start + 1];
	strncpy(deleted_seq, chr_seq + sv->start+1, sv->end - sv->start);
	deleted_seq[sv->end - sv->start] = '\0';

	// if deleted sequence contains sv->ins_seq as a substring, split into two smaller deletions
	char* pos_ptr = strstr(deleted_seq, sv->ins_seq.c_str());
	if (pos_ptr) {
		// Create two new SVs for the split deletions
		hts_pos_t sv1_start = sv->start;
		hts_pos_t sv1_end = sv->start + (pos_ptr - deleted_seq);
		hts_pos_t sv2_start = sv->start + (pos_ptr - deleted_seq) + sv->ins_seq.length();
		hts_pos_t sv2_end = sv->end;
		if (sv1_end-sv1_start >= sv2_end-sv2_start) {
			sv->end = sv1_end;
			std::shared_ptr<sv_t> sv2 = std::make_shared<deletion_t>(sv->chr, sv2_start, sv2_end, "", nullptr, nullptr, nullptr, nullptr);
			sv->aux_indels.push_back(sv2);
		} else {
			sv->start = sv2_start;
			std::shared_ptr<sv_t> sv1 = std::make_shared<deletion_t>(sv->chr, sv1_start, sv1_end, "", nullptr, nullptr, nullptr, nullptr);
			sv->aux_indels.push_back(sv1);
		}
		sv->ins_seq = "";
	} else if (sv->ins_seq.length() == 1) { // DEL + SNP
		hts_pos_t snp_pos = sv->end;
		char ref_base = chr_seq[snp_pos];
		char alt_base = sv->ins_seq[0];
		snp_t snp(snp_pos, alt_base);
		sv->aux_snps.push_back(snp);

		sv->end--;
		sv->ins_seq = "";
	}
	delete[] deleted_seq;
}

void atomize(int id, std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL" && sv->svsize() <= 50) { 
		// only atomize small deletions
		// splitting large deletions can create problems to the current genotyping algorithm
		// especially when calculating features like discordant pairs or read depth
		atomize_del(sv);
	} else if (svtype == "INS") {
		// atomize_ins(sv);
	}
}

void simplify_del(std::shared_ptr<sv_t> sv) {

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	// try to shorten deletion if ins_seq
	int start_is = 0;
	while (start_is < sv->ins_seq.length() && sv->ins_seq[start_is] == chr_seq[sv->start+1]) {
		start_is++;
		sv->start++;
	}
	sv->ins_seq = sv->ins_seq.substr(start_is);

	int end_is = sv->ins_seq.length();
	while (end_is > 0 && sv->ins_seq[end_is-1] == chr_seq[sv->end]) {
		end_is--;
		sv->end--;
	}
	sv->ins_seq = sv->ins_seq.substr(0, end_is);
}

void simplify_ins(std::shared_ptr<sv_t> sv) {

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	// try to shorten insertion if start != end
	int start_is = 0;
	while (start_is < sv->ins_seq.length() && sv->start < sv->end && toupper(sv->ins_seq[start_is]) == toupper(chr_seq[sv->start+1])) {
		start_is++;
		sv->start++;
	}
	sv->ins_seq = sv->ins_seq.substr(start_is);

	int end_is = sv->ins_seq.length();
	while (end_is > 0 && sv->start < sv->end && toupper(sv->ins_seq[end_is-1]) == toupper(chr_seq[sv->end])) {
		end_is--;
		sv->end--;
	}
	sv->ins_seq = sv->ins_seq.substr(0, end_is);
}

void simplify(std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL") {
		simplify_del(sv);
	} else if (svtype == "INS") {
		simplify_ins(sv);
	}
}

void left_align_del(std::shared_ptr<sv_t> sv) {

	if (!sv->ins_seq.empty()) return;

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	for (snp_t snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	hts_pos_t limit = 0; // we cannot left-align past this position
	for (std::shared_ptr<sv_t> indel : sv->aux_indels) {
		if (indel->end <= sv->start) {
			limit = std::max(limit, indel->end); 
		}
	}

	while (sv->start > limit && chr_seq[sv->start] == chr_seq[sv->end]) {
		sv->start--;
		sv->end--;
	}

	for (snp_t snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	// delete snps that fall within the deleted region
	sv->aux_snps.erase(std::remove_if(sv->aux_snps.begin(), sv->aux_snps.end(),
		[sv](const snp_t& snp) { return snp.pos > sv->start && snp.pos <= sv->end; }), sv->aux_snps.end());
}

void left_align_dup(std::shared_ptr<sv_t> sv) {

	if (!sv->ins_seq.empty()) return;
	if (!sv->aux_indels.empty()) return; // this may be complicated, let us skip for now

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	hts_pos_t limit = -1;
	for (snp_t& snp : sv->aux_snps) {
		if (snp.pos <= sv->start) {
			limit = snp.pos; // note that snps are sorted by position
		}
	}

	int i = 0;
	while (sv->start > 0 && chr_seq[sv->start] == chr_seq[sv->end]) {
		sv->start--;
		sv->end--;
		i++;
	}
}

void left_align_ins(std::shared_ptr<sv_t> sv) {

	if (!sv->aux_indels.empty()) return; // this may be complicated, let us skip for now

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	for (snp_t snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	if (sv->start == sv->end) {
		while (sv->start > 0 && toupper(chr_seq[sv->start]) == toupper(sv->ins_seq[sv->ins_seq.length()-1])) {
			for (int i = sv->ins_seq.length()-1; i >= 1; i--) {
				sv->ins_seq[i] = sv->ins_seq[i-1];
			}
			sv->ins_seq[0] = toupper(chr_seq[sv->start]);
			sv->start--;
			sv->end--;
		}
	}

	for (snp_t snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	// delete snps that fall within the deleted region
	sv->aux_snps.erase(std::remove_if(sv->aux_snps.begin(), sv->aux_snps.end(),
		[sv](const snp_t& snp) { return snp.pos > sv->start && snp.pos <= sv->end; }), sv->aux_snps.end());
}


void left_align(std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL") {
		left_align_del(sv);
	} else if (svtype == "DUP") {
		left_align_dup(sv);
	} else if (svtype == "INS") {
		left_align_ins(sv);
	}
}

int main(int argc, char* argv[]) {

	std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];

	chr_seqs.read_fasta_into_map(reference_fname);

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* vcf_record = bcf_init();

	std::vector<std::shared_ptr<sv_t>> svs;
	while (bcf_read(in_vcf_file, hdr, vcf_record) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, vcf_record);
		if (sv == nullptr) {
			std::cout << "Ignoring SV of unsupported type: " << vcf_record->d.id << std::endl; 
			continue;
		}
		sv->vcf_entry = bcf_dup(vcf_record);
		simplify(sv);
		atomize(0, sv);
		left_align(sv);
		svs.push_back(sv);
	}

	for (const auto& sv : svs) {
		bcf1_t* vcf_record_norm = bcf_init();
		sv2bcf(hdr, vcf_record_norm, sv.get(), chr_seqs.get_seq(sv->chr));
		copy_all_fmt(hdr, sv->vcf_entry, vcf_record_norm);
		normalised_vcf_records.push_back(vcf_record_norm);
	}

	std::sort(normalised_vcf_records.begin(), normalised_vcf_records.end(),
			[](const bcf1_t* b1, const bcf1_t* b2) { return std::tie(b1->rid, b1->pos) < std::tie(b2->rid, b2->pos); });

	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write VCF header to " +  out_vcf_fname);
	}
	for (bcf1_t* vcf_record_norm : normalised_vcf_records) {
		if (bcf_write(out_vcf_file, hdr, vcf_record_norm) != 0) {
			throw std::runtime_error("Failed to write VCF record to " +  out_vcf_fname);
		}
	}

	hts_close(in_vcf_file);
	hts_close(out_vcf_file);

	tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
