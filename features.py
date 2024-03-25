import os, pysam
from collections import defaultdict
import numpy as np

class Features:

    source_to_int = { "2SR" : 0, "2HSR" : 1, "HSR-SR" : 2, "SR-HSR" : 3, "1SR_LC" : 4, "1SR_RC" : 5,
                 "1HSR_LC" : 6, "1HSR_RC" : 7, "DP" : 8, "DE_NOVO_ASSEMBLY" : 9, "REFERENCE_GUIDED_ASSEMBLY" : 10 }

    svtype_to_int = { "DEL" : 0, "DUP" : 1, "INS" : 2 }

    def get_value(info, key, default, norm_factor = 1.0):
        if key in info:
            v = info[key]
        else:
            v = default
        if isinstance(v, list) or isinstance(v, tuple):
            return [float(x)/norm_factor for x in v]
        else:
            return float(v)/norm_factor

    def record_to_features(record, stats):
        min_depth = get_stat(stats, 'min_depth', record.chrom)
        median_depth = get_stat(stats, 'median_depth', record.chrom)
        max_depth = get_stat(stats, 'max_depth', record.chrom)
        max_is = stats['max_is']['.']
        read_len = stats['read_len']['.']

        info = record.info
        svtype_str = record.info['SVTYPE']
        source_str = record.info['SOURCE']
        svinsseq = ""
        if 'SVINSSEQ' in info:
            svinsseq = info['SVINSSEQ']
        if 'SVLEN' in info:
            svlen = abs(float(info['SVLEN']))
        else:
            svlen = len(svinsseq) - (record.stop - record.pos)
        svinslen = Features.get_value(info, 'SVINSLEN', 0)
        if svinslen == 0 and svinsseq:
            svinslen = len(svinsseq)

        if svtype_str == "INS":
            prefix_mh_len_ratio = Features.get_value(info, 'PREFIX_MH_LEN', 0, max(1, svinslen))
            suffix_mh_len_ratio = Features.get_value(info, 'SUFFIX_MH_LEN', 0, max(1, svinslen))

        ins_seq_cov_prefix_start, ins_seq_cov_prefix_end = 0, 1
        ins_seq_cov_suffix_start, ins_seq_cov_suffix_end = 0, 1
        if source_str == "DE_NOVO_ASSEMBLY":
            if '-' in svinsseq:
                ins_seq_cov_prefix_end = svinsseq.index('-')/len(svinsseq)
                ins_seq_cov_suffix_start = ins_seq_cov_prefix_end
        elif source_str == "REFERENCE_GUIDED_ASSEMBLY":
            ins_seq_cov_prefix_start, ins_seq_cov_prefix_end = Features.get_value(info, 'INS_PREFIX_COV', [0, 0], len(svinsseq)-1)
            ins_seq_cov_suffix_start, ins_seq_cov_suffix_end = Features.get_value(info, 'INS_SUFFIX_COV', [0, 0], len(svinsseq)-1)

        split_reads = Features.get_value(info, 'SPLIT_READS', [0, 0])
        split_reads_ratio = [split_reads[0]/median_depth, split_reads[1]/median_depth]
        fwd_split_reads = Features.get_value(info, 'FWD_SPLIT_READS', [0, 0])
        fwd_split_reads_ratio = [fwd_split_reads[0]/max(1, split_reads[0]), fwd_split_reads[1]/max(1, split_reads[1])]
        rev_split_reads = Features.get_value(info, 'REV_SPLIT_READS', [0, 0])
        rev_split_reads_ratio = [rev_split_reads[0]/max(1, split_reads[0]), rev_split_reads[1]/max(1, split_reads[1])]
        overlap = Features.get_value(info, 'OVERLAP', 0, read_len)
        mismatch_rate = Features.get_value(info, 'MISMATCH_RATE', 0)

        rcc_ext_1sr_reads = Features.get_value(info, 'RCC_EXT_1SR_READS', [0, 0], median_depth*max_is)
        lcc_ext_1sr_reads = Features.get_value(info, 'LCC_EXT_1SR_READS', [0, 0], median_depth*max_is)
        rcc_hq_ext_1sr_reads = Features.get_value(info, 'RCC_HQ_EXT_1SR_READS', [0, 0], median_depth*max_is)
        lcc_hq_ext_1sr_reads = Features.get_value(info, 'LCC_HQ_EXT_1SR_READS', [0, 0], median_depth*max_is)

        full_junction_score = Features.get_value(info, 'FULL_JUNCTION_SCORE', 0)
        split_junction_score1 = [float(x) for x in info['SPLIT_JUNCTION_SCORE']]
        split_junction_score2 = [float(x) for x in info['SPLIT_JUNCTION_SCORE2']]
        split_junction_size = [float(x) for x in info['SPLIT_JUNCTION_SIZE']]
        if sum(split_junction_score1) == 0:
            full_to_split_junction_score_ratio = 0
            full_to_split_junction_score_diff = 0
            split2_to_split1_junction_score_ratio = [0, 0]
            split2_to_split1_junction_score_diff_ratio = [0, 0]
            split_to_size_ratio = [0, 0]
        else:
            full_to_split_junction_score_ratio = full_junction_score/sum(split_junction_score1)
            full_to_split_junction_score_diff = (full_junction_score-sum(split_junction_score1))
            split2_to_split1_junction_score_ratio = [split_junction_score2[0]/max(1, split_junction_score1[0]), split_junction_score2[1]/max(1, split_junction_score1[1])]
            s1 = max(split_junction_size[0]-split_junction_score1[0], 0.01)/max(split_junction_size[0]-split_junction_score2[0], 0.01)
            s2 = max(split_junction_size[1]-split_junction_score1[1], 0.01)/max(split_junction_size[1]-split_junction_score2[1], 0.01)
            split2_to_split1_junction_score_diff_ratio = [s1, s2]
            split_to_size_ratio = [split_junction_score1[0]/split_junction_size[0], split_junction_score1[1]/split_junction_size[1]]
        split_junction_size_ratio = [split_junction_size[0]/sum(split_junction_size), split_junction_size[1]/sum(split_junction_size)]

        max_mapq = Features.get_value(info, 'MAX_MAPQ', [0, 0])
        lb_diff, ub_diff = 0, 0
        if 'REMAP_LB' in info and 'REMAP_UB' in info:
            remap_lb = float(info['REMAP_LB'])
            remap_ub = float(info['REMAP_UB'])
            lb_diff = max(0, remap_lb-record.pos)
            ub_diff = max(0, record.stop-remap_ub)

        disc_pairs = Features.get_value(info, 'DISC_PAIRS', [0, 0])
        disc_pairs_scaled = [0, 0]
        if svtype_str == "DEL":
            min_is_to_become_disc = int(max(0, max_is-svlen))
            min_disc_pairs = stats['min_pairs_crossing_gap'][str(min_is_to_become_disc)]
            max_disc_pairs = stats['max_pairs_crossing_gap'][str(min_is_to_become_disc)]
            disc_pairs_scaled = [(d-min_disc_pairs)/(max_disc_pairs-min_disc_pairs) for d in disc_pairs]
        elif svtype_str == "INS" and source_str in ("DE_NOVO_ASSEMBLY", "REFERENCE_GUIDED_ASSEMBLY"):
            min_inslen = int(min(max_is, svinslen))
            min_disc_pairs = stats['min_disc_pairs_by_insertion_size'][str(min_inslen)]
            max_disc_pairs = stats['max_disc_pairs_by_insertion_size'][str(min_inslen)]
            disc_pairs_scaled = [(d-min_disc_pairs)/(max_disc_pairs-min_disc_pairs) for d in disc_pairs]
        else:
            disc_pairs_scaled = [d/median_depth for d in disc_pairs]
        disc_pairs_highmapq = Features.get_value(info, 'DISC_PAIRS_HIGHMAPQ', [0, 0], median_depth)
        disc_pairs_highmapq_ratio = [d/max(1, disc_pairs[i]) for i, d in enumerate(disc_pairs_highmapq)]
        disc_pairs_maxmapq = Features.get_value(info, 'DISC_PAIRS_MAXMAPQ', [0, 0])
        conc_pairs = Features.get_value(info, 'CONC_PAIRS', 0, median_depth)
        disc_avg_nm = Features.get_value(info, 'DISC_AVG_NM', [0, 0], read_len)
        disc_pair_surrounding = Features.get_value(info, 'DISC_PAIRS_SURROUNDING', [0, 0], median_depth) 

        disc_pairs = [d/median_depth for d in disc_pairs]
        ptn_ratio = [d/max(1, d+conc_pairs) for d in disc_pairs]
        ks_pval = max(0, Features.get_value(info, 'KS_PVAL', 0.0))
        max_size_diff = 0
        if 'MAX_SIZE' in info:
            max_size = float(info['MAX_SIZE'])
            max_size_diff = max(0, svlen - 2*max_size)

        median_depths = [float(x) for x in info['MEDIAN_DEPTHS']]
        median_depths_norm = [float(x)/median_depth for x in info['MEDIAN_DEPTHS']]
        median_depths_ratio = [median_depths[0]/max(1, median_depths[1]), median_depths[3]/max(1, median_depths[2])]
        median_depths_above_max = [max(0, x-max_depth)/median_depth for x in median_depths]
        median_depths_below_min = [max(0, min_depth-median_depths[0])/median_depth, 
                                   max(0, min_depth-median_depths[3])/median_depth] 
        if 'CLUSTER_DEPTHS' in info:
            cluster_depths_above_max = [max(0, float(x)-max_depth)/median_depth for x in info['CLUSTER_DEPTHS']]
        else:
            cluster_depths_above_max = [0, 0]

        if svtype_str == "INS":
            prefix_base_count = Features.get_value(info, 'PREFIX_BASE_COUNT', [0, 0, 0, 0])
            max_prefix_base_count_perc = max(prefix_base_count)/max(1, sum(prefix_base_count))
            prefix_base_count_perc = [x/max(1, sum(prefix_base_count)) for x in prefix_base_count]
            suffix_base_count = Features.get_value(info, 'SUFFIX_BASE_COUNT', [0, 0, 0, 0])
            max_suffix_base_count_perc = max(suffix_base_count)/max(1, sum(suffix_base_count))
            suffix_base_count_perc = [x/max(1, sum(suffix_base_count)) for x in suffix_base_count]

        feature_values = [svlen, svinslen]
        if source_str == "DE_NOVO_ASSEMBLY" or source_str == "REFERENCE_GUIDED_ASSEMBLY":
            feature_values += [ins_seq_cov_prefix_start, ins_seq_cov_prefix_end, ins_seq_cov_suffix_start, ins_seq_cov_suffix_end]
        feature_values += split_reads_ratio + fwd_split_reads_ratio + rev_split_reads_ratio + [sum(fwd_split_reads)/max(1, sum(split_reads)), sum(rev_split_reads)/max(1, sum(split_reads)), overlap, mismatch_rate]
        feature_values += rcc_ext_1sr_reads + lcc_ext_1sr_reads + rcc_hq_ext_1sr_reads + lcc_hq_ext_1sr_reads
        feature_values += [full_to_split_junction_score_ratio, full_to_split_junction_score_diff] + split2_to_split1_junction_score_ratio
        feature_values += split2_to_split1_junction_score_diff_ratio + split_to_size_ratio + split_junction_size_ratio + [max(split_junction_size_ratio), min(split_junction_size_ratio)]
        feature_values += max_mapq + [lb_diff, ub_diff]
        feature_values += disc_pairs_scaled + disc_pairs_highmapq_ratio + disc_pairs_maxmapq + [conc_pairs] 
        feature_values += disc_pair_surrounding + disc_avg_nm
        feature_values += ptn_ratio + [ks_pval, max_size_diff]
        feature_values += median_depths_norm + median_depths_ratio + median_depths_above_max + median_depths_below_min + cluster_depths_above_max
        if svtype_str == "INS":
            feature_values += [prefix_mh_len_ratio, suffix_mh_len_ratio]
            feature_values += prefix_base_count_perc + suffix_base_count_perc + [max_prefix_base_count_perc, max_suffix_base_count_perc]
        return feature_values

def read_false_ids(file_path, tolerate_no_fps = False):
    if not os.path.exists(file_path) and tolerate_no_fps:
        return set()
    with open(file_path, 'r') as file:
        false_ids = set([line.strip() for line in file])
    return false_ids

def get_stat(stats, stat_name, chrom):
    if chrom in stats[stat_name]:
        return stats[stat_name][chrom]
    return stats[stat_name]['.']

# Function to parse the VCF file and extract relevant features using pysam
def parse_vcf(vcf_fname, stats_fname, fp_fname, svtype, tolerate_no_fps = False):
    false_ids = read_false_ids(fp_fname, tolerate_no_fps=tolerate_no_fps)
    vcf_reader = pysam.VariantFile(vcf_fname)
    stats_reader = open(stats_fname, 'r')
    data_by_source = defaultdict(list)
    labels_by_source = defaultdict(list)
    variant_ids_by_source = defaultdict(list)

    # read the stats file and extract the relevant values
    stats = defaultdict(dict)
    for line in stats_reader:
        sl = line.strip().split()
        stats[sl[0]][sl[1]] = int(sl[2])

    for record in vcf_reader.fetch():
        if svtype != 'ALL' and record.info['SVTYPE'] != svtype:
            continue
        source = record.info['SVTYPE'] + "_" + record.info['SOURCE']
        feature_values = Features.record_to_features(record, stats)
        data_by_source[source].append(feature_values)
        labels_by_source[source].append(0 if record.id in false_ids else 1)
        variant_ids_by_source[source].append(record.id)
    for source in data_by_source:
        data_by_source[source] = np.array(data_by_source[source])
        labels_by_source[source] = np.array(labels_by_source[source])
    return data_by_source, labels_by_source, variant_ids_by_source
