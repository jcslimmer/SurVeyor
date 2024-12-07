#include <unordered_map>

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "types.h"
#include "vcf_utils.h"

// if store_in_rc, the new 2SR entry will be in rc_sv, otherwise in lc_sv
void merge_1sr_with_1sr(bcf_hdr_t* hdr, sv_t* rc_sv, sv_t* lc_sv, bool store_in_rc) {
	if (store_in_rc) {
		rc_sv->lc_consensus = lc_sv->lc_consensus;
	} else {
		lc_sv->rc_consensus = rc_sv->rc_consensus;
	}
}

int priority(sv_t* sv) {
	if (sv->source == "DE_NOVO_ASSEMBLY" || sv->source == "REFERENCE_GUIDED_ASSEMBLY") {
		insertion_t* ins = dynamic_cast<insertion_t*>(sv);
		if (ins->imprecise || ins->incomplete_ins_seq()) return 7;
		else return 1;
	} else if (sv->source == "2SR") return 2;
	else if (sv->source == "HSR-SR" || sv->source == "SR-HSR") return 3;
	else if (sv->source == "2HSR") return 4;
	else if (sv->source == "1SR_LC" || sv->source == "1SR_RC") return 5;
	else if (sv->source == "1HSR_LC" || sv->source == "1HSR_RC") return 6;
	else if (sv->source == "DP") return 7;
	throw std::runtime_error("Unknown source: " + sv->source);
}

int main(int argc, char* argv[]) {

	std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];

	std::unordered_map<std::string, int> source_priorities;
	source_priorities["2SR"] = 1;
	source_priorities["HSR-SR"] = 2;
	source_priorities["SR-HSR"] = 2;
	source_priorities["2HSR"] = 3;
	source_priorities["1SR_LC"] = 4;
	source_priorities["1SR_RC"] = 4;
	source_priorities["1HSR_LC"] = 5;
	source_priorities["1HSR_RC"] = 5;
	source_priorities["DP"] = 6;

	std::unordered_map<std::string, sv_t*> sv_entries;
	std::vector<std::string> ids_sorted;

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	bcf_hdr_t* in_vcf_hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* b = bcf_init();
	while (bcf_read(in_vcf_file, in_vcf_hdr, b) == 0) {
		sv_t* sv = bcf_to_sv(in_vcf_hdr, b);

		std::string unique_id = sv->unique_key();
		if (sv_entries.count(unique_id)) {
			// if this is not PASS and there is an equivalent indel already stored, ignore
			if (sv->is_fail()) continue;
			else if (sv_entries[unique_id]->is_fail()) {
				// if stored indel is not PASS, but there is an equivalent PASS indel, store that instead
				sv_entries[unique_id] = sv;
				continue;
			}

			std::string src_prev = sv_entries[unique_id]->source;
			std::string src_curr = sv->source;
			if (source_priorities[src_prev] > source_priorities[src_curr]) {
				std::swap(sv_entries[unique_id], sv);
				std::swap(src_prev, src_curr);
			}
 
			if (src_prev == "2SR") continue; // 2SR has highest priority. If prev indel is 2SR, keep it
			else if (src_prev == "HSR-SR" && src_curr == "1SR_RC") { // HSR-SR + 1SR_RC -> 2SR
				merge_1sr_with_1sr(in_vcf_hdr, sv, sv_entries[unique_id], false);
				sv_entries[unique_id]->source = "2SR";
			} else if (src_prev == "SR-HSR" && src_curr == "1SR_LC") { // SR-HSR + 1SR_LC -> 2SR
				merge_1sr_with_1sr(in_vcf_hdr, sv_entries[unique_id], sv, true);
				sv_entries[unique_id]->source = "2SR";
			} else if (src_prev == "1SR_LC" && src_curr == "1SR_RC") { // 1SR_RC + 1SR_LC -> 2SR
				merge_1sr_with_1sr(in_vcf_hdr, sv, sv_entries[unique_id], false);
				sv_entries[unique_id]->source = "2SR";
			} else if (src_prev == "1SR_RC" && src_curr == "1SR_LC") { // 1SR_RC + 1SR_LC -> 2SR
				merge_1sr_with_1sr(in_vcf_hdr, sv_entries[unique_id], sv, true);
				sv_entries[unique_id]->source = "2SR";
			} else if (src_prev == "1SR_RC" && src_curr == "1HSR_LC") { // 1SR_RC + 1HSR_LC -> SR-HSR
				merge_1sr_with_1sr(in_vcf_hdr, sv_entries[unique_id], sv, true);
				sv_entries[unique_id]->source = "SR-HSR";
			} else if (src_prev == "1SR_LC" && src_curr == "1HSR_RC") { // 1HSR_RC + 1SR_LC -> HSR-SR
				merge_1sr_with_1sr(in_vcf_hdr, sv, sv_entries[unique_id], false);
				sv_entries[unique_id]->source = "HSR-SR";
			} else if (src_prev == "1HSR_RC" && src_curr == "1HSR_LC") { // 1HSR_RC + 1HSR_LC -> 2HSR
				merge_1sr_with_1sr(in_vcf_hdr, sv_entries[unique_id], sv, true);
				sv_entries[unique_id]->source = "2HSR";
			} else if (src_prev == "1HSR_LC" && src_curr == "1HSR_RC") { // 1HSR_RC + 1HSR_LC -> 2HSR
				merge_1sr_with_1sr(in_vcf_hdr, sv, sv_entries[unique_id], false);
				sv_entries[unique_id]->source = "2HSR";
			}

			// TODO: other possible cases: HSR-SR + SR-HSR -> 2SR, 2HSR + 1SR -> HSR-SR (or SR-HSR)
		} else {
			ids_sorted.push_back(unique_id);
			sv_entries[unique_id] = sv;
		}
	}

	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, in_vcf_hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	chr_seqs_map_t chr_seqs;
	chr_seqs.read_fasta_into_map(reference_fname);

	for (std::string& id : ids_sorted) {
		sv_t* sv = sv_entries[id];
		sv2bcf(in_vcf_hdr, b, sv, chr_seqs.get_seq(sv->chr));
		if (bcf_write(out_vcf_file, in_vcf_hdr, b) != 0) {
			throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
		}
	}

	bcf_close(out_vcf_file);
}
