#include <string>
#include <algorithm>

#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "../src/sam_utils.h"
#include "../src/vcf_utils.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;

void normalise_del(std::shared_ptr<sv_t> sv) {

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

	if (sv->ins_seq.empty()) {
		while (sv->start > 0 && chr_seq[sv->start] == chr_seq[sv->end]) {
			sv->start--;
			sv->end--;
		}
	}
}

void normalise_dup(std::shared_ptr<sv_t> sv) {

	if (!sv->ins_seq.empty()) return;

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	int i = 0;
	while (sv->start > 0 && chr_seq[sv->start] == chr_seq[sv->end]) {
		sv->start--;
		sv->end--;
		i++;
	}
}

void normalise_ins(std::shared_ptr<sv_t> sv) {

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
}


void normalise(std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL") {
		normalise_del(sv);
	} else if (svtype == "DUP") {
		normalise_dup(sv);
	} else if (svtype == "INS") {
		normalise_ins(sv);
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

	std::vector<bcf1_t*> normalised_vcf_records;
	while (bcf_read(in_vcf_file, hdr, vcf_record) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, vcf_record);
		normalise(sv);
		bcf1_t* vcf_record_norm = bcf_init();
		sv2bcf(hdr, vcf_record_norm, sv.get(), chr_seqs.get_seq(sv->chr));
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
