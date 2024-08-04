import os, pysam
from collections import defaultdict
import numpy as np

class Features:

    shared_features_names = ['START_STOP_DIST', 'SVLEN', 'SVINSLEN',
                             'SV_REF_PREFIX_A_RATIO', 'SV_REF_PREFIX_C_RATIO', 'SV_REF_PREFIX_G_RATIO', 'SV_REF_PREFIX_T_RATIO', 'MAX_SV_REF_PREFIX_BASE_RATIO',
                             'SV_REF_SUFFIX_A_RATIO', 'SV_REF_SUFFIX_C_RATIO', 'SV_REF_SUFFIX_G_RATIO', 'SV_REF_SUFFIX_T_RATIO', 'MAX_SV_REF_SUFFIX_BASE_RATIO',
                             'LEFT_ANCHOR_A_RATIO', 'LEFT_ANCHOR_C_RATIO', 'LEFT_ANCHOR_G_RATIO', 'LEFT_ANCHOR_T_RATIO', 'MAX_LEFT_ANCHOR_BASE_RATIO',
                             'RIGHT_ANCHOR_A_RATIO', 'RIGHT_ANCHOR_C_RATIO', 'RIGHT_ANCHOR_G_RATIO', 'RIGHT_ANCHOR_T_RATIO', 'MAX_RIGHT_ANCHOR_BASE_RATIO',
                             'INS_PREFIX_A_RATIO', 'INS_PREFIX_C_RATIO', 'INS_PREFIX_G_RATIO', 'INS_PREFIX_T_RATIO', 'MAX_INS_PREFIX_BASE_COUNT_RATIO',
                             'INS_SUFFIX_A_RATIO', 'INS_SUFFIX_C_RATIO', 'INS_SUFFIX_G_RATIO', 'INS_SUFFIX_T_RATIO', 'MAX_INS_SUFFIX_BASE_COUNT_RATIO']
    
    denovo_features_names = ['SPLIT_READS_RATIO', 'SPLIT_READS_RATIO1', 'SPLIT_READS_RATIO2', 'FWD_SPLIT_READS_RATIO1', 'FWD_SPLIT_READS_RATIO2', 'REV_SPLIT_READS_RATIO1', 'REV_SPLIT_READS_RATIO2',
                      'FWD_SPLIT_READS_RATIO', 'REV_SPLIT_READS_RATIO', 'OVERLAP', 'MISMATCH_RATE', 'RCC_EXT_1SR_READS1', 'RCC_EXT_1SR_READS2', 'LCC_EXT_1SR_READS1', 'LCC_EXT_1SR_READS2',
                      'EXT_1SR_READS1', 'EXT_1SR_READS2', 'RCC_HQ_EXT_1SR_READS1', 'RCC_HQ_EXT_1SR_READS2', 'LCC_HQ_EXT_1SR_READS1', 'LCC_HQ_EXT_1SR_READS2', 'HQ_EXT_1SR_READS1', 'HQ_EXT_1SR_READS2', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                      'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                      'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                      'MAX_MAPQ1', 'MAX_MAPQ2', 'MAX_MAPQ', 'MAX_MAPQ_EXT1', 'MAX_MAPQ_EXT2', 'MAX_MAPQ_EXT',
                      'LB_DIFF', 'UB_DIFF', 'B_DIFF', 'DISC_PAIRS_SCALED1', 'DISC_PAIRS_SCALED2', 'DISC_PAIRS_HIGHMAPQ_RATIO1', 'DISC_PAIRS_HIGHMAPQ_RATIO2', 'DISC_PAIRS_MAXMAPQ1', 'DISC_PAIRS_MAXMAPQ2',
                      'CONC_PAIRS_SCALED', 'DISC_PAIRS_SURROUNDING1', 'DISC_PAIRS_SURROUNDING2', 'DISC_AVG_NM1', 'DISC_AVG_NM2', 'PTN_RATIO1', 'PTN_RATIO2', 'KS_PVAL', 'KS_PVAL_HIGHMQ', 'SIZE_NORM', 'SIZE_NORM_HIGHMQ',
                      'PREFIX_MH_LEN_RATIO', 'SUFFIX_MH_LEN_RATIO', 'INS_SEQ_COV_PREFIX_START', 'INS_SEQ_COV_PREFIX_END', 'INS_SEQ_COV_SUFFIX_START', 'INS_SEQ_COV_SUFFIX_END',
                      'MEDIAN_DEPTHS_NORM1', 'MEDIAN_DEPTHS_NORM2', 'MEDIAN_DEPTHS_NORM3', 'MEDIAN_DEPTHS_NORM4', 'MEDIAN_DEPTHS_RATIO1', 'MEDIAN_DEPTHS_RATIO2', 
                      'MEDIAN_DEPTHS_HIGHMQ_NORM1', 'MEDIAN_DEPTHS_HIGHMQ_NORM2', 'MEDIAN_DEPTHS_HIGHMQ_NORM3', 'MEDIAN_DEPTHS_HIGHMQ_NORM4', 'MEDIAN_DEPTHS_HIGHMQ_RATIO1', 'MEDIAN_DEPTHS_HIGHMQ_RATIO2',
                      'CLUSTER_DEPTHS_ABOVE_MAX1', 'CLUSTER_DEPTHS_ABOVE_MAX2']

    def get_denovo_model_name(record, max_is):
        svtype_str = record.info['SVTYPE']
        source_str = Features.get_string_value(record.info, 'SOURCE', "")
        split_reads = Features.get_number_value(record.info, 'SPLIT_READS', [0, 0])
        if svtype_str == "DEL":
            if split_reads[0] > 0 and split_reads[1] > 0:
                source_str = "2SR"
            elif split_reads[0] > 0 or split_reads[1] > 0:
                source_str = "1SR"
            else:
                source_str = "DP"
            if abs(Features.get_svlen(record)) >= max_is:
                source_str += "_LARGE"
        elif svtype_str == "DUP":
            if split_reads[0] > 0 and split_reads[1] > 0:
                source_str = "2SR"
            elif split_reads[0] > 0 or split_reads[1] > 0:
                source_str = "1SR"
        elif svtype_str == "INS":
            if source_str == "DE_NOVO_ASSEMBLY" or source_str == "REFERENCE_GUIDED_ASSEMBLY":
                pass
            elif split_reads[0] > 0 and split_reads[1] > 0:
                source_str = "2SR"
            elif split_reads[0] > 0 or split_reads[1] > 0:
                source_str = "1SR"
        return svtype_str + "_" + source_str

    def get_regt_model_name(record, max_is, read_len):
        svtype_str = record.info['SVTYPE']
        if svtype_str == "DEL" and abs(Features.get_svlen(record)) >= max_is:
                svtype_str += "_LARGE"
        elif svtype_str == "DUP" and record.stop-record.start > read_len-30:
            svtype_str += "_LARGE"
        if "IMPRECISE" in record.info:
            svtype_str += "_IMPRECISE"
        return svtype_str

    def get_denovo_feature_names(model_name):
        return Features.shared_features_names + Features.denovo_features_names

    regt_shared_features_names = \
            ['IP', 'AR1', 'ARC1', 'ARCF1', 'ARCR1', 'MAXARCD1', 'ARCAS1', 'ARC1MQ', 'ARC1HQ',
             'AR2', 'ARC2', 'ARCF2', 'ARCR2', 'MAXARCD2', 'ARCAS2', 'ARC2MQ', 'ARC2HQ',
             'RR1', 'RRC1', 'RR2', 'RRC2', 'ER', 
             'AR1_RATIO', 'AR2_RATIO', 'RR1_RATIO', 'RR2_RATIO', 'ARC1_RATIO', 'ARC2_RATIO', 'RRC1_RATIO', 'RRC2_RATIO',
             'AR1_OVER_RR1', 'RR1_OVER_AR1', 'AR2_OVER_RR2', 'RR2_OVER_AR2', 'ARC1_OVER_RRC1', 'RRC1_OVER_ARC1', 'ARC2_OVER_RRC2', 'RRC2_OVER_ARC2',
             'MDLF', 'MDSP', 'MDSF', 'MDRF', 'MDLC', 'MDRC', 'MDLFHQ', 'MDSPHQ', 'MDSFHQ', 'MDRFHQ',
             'MDSP_OVER_MDLF', 'MDSF_OVER_MDRF', 'MDLF_OVER_MDSP', 'MDRF_OVER_MDSF',
             'MDSP_OVER_MDLF_HQ', 'MDSF_OVER_MDRF_HQ', 'MDLF_OVER_MDSP_HQ', 'MDRF_OVER_MDSF_HQ', 'DPSL', 'DPSR', 'CP', 
             'AXR', 'AXRHQ', 'EXL', 'EXAS', 'EXRS', 'EXAS_EXRS_RATIO', 'EXAS_EXRS_DIFF']

    regt_stat_test_features_names = ['FMT_KSPVAL', 'FMT_KSPVAL_HQ', 'FMT_SIZE_NORM', 'FMT_SIZE_NORM_HQ']

    regt_dp_features_names = ['DP1', 'DP2', 'DP1_HQ_RATIO', 'DP2_HQ_RATIO', 'DP1MQ', 'DP2MQ', 'DPLANM', 'DPRANM', 'PTNR']

    def get_regt_feature_names(model_name):
        extra_feature_names = []
        if model_name in ["DEL", "DEL_IMPRECISE", "DUP"]:
            extra_feature_names = Features.regt_stat_test_features_names
        elif model_name in ["DEL_LARGE", "DEL_LARGE_IMPRECISE", "INS", "INS_IMPRECISE"]:
            extra_feature_names = Features.regt_dp_features_names
        return Features.shared_features_names + Features.regt_shared_features_names + extra_feature_names

    def get_number_value(info, key, default, norm_factor = 1.0):
        if key in info:
            v = info[key]
        else:
            v = default
        if isinstance(v, list) or isinstance(v, tuple):
            return [float(x)/norm_factor for x in v]
        else:
            return float(v)/norm_factor

    def get_string_value(info, key, default):
        if key in info:
            return info[key]
        else:
            return default

    def get_svlen(record):
        if 'SVLEN' not in record.info:
            svinsseq = ""
            if 'SVINSSEQ' in record.info:
                svinsseq = record.info['SVINSSEQ']
            return len(svinsseq) - (record.stop - record.pos)
        svlen = record.info['SVLEN']
        if isinstance(svlen, list) or isinstance(svlen, tuple):
            return svlen[0]
        else:
            return svlen

    def normalise(value, min, max):
        if max == min:
            return value - min
        return (value - min) / (max - min)

    def record_to_features(record, stats):
        min_depth = get_stat(stats, 'min_depth', record.chrom)
        median_depth = get_stat(stats, 'median_depth', record.chrom)
        max_depth = get_stat(stats, 'max_depth', record.chrom)
        max_is = stats['max_is']['.']
        read_len = stats['read_len']['.']
        
        denovo_model_name = Features.get_denovo_model_name(record, max_is)
        regt_model_name = Features.get_regt_model_name(record, max_is, read_len)
        features = dict()

        info = record.info
        svtype_str = record.info['SVTYPE']
        source_str = Features.get_string_value(info, 'SOURCE', "")
        features['START_STOP_DIST'] = record.stop - record.pos

        svinsseq = ""
        if 'SVINSSEQ' in info:
            svinsseq = info['SVINSSEQ']
        svlen = abs(Features.get_svlen(record))
        features['SVLEN'] = svlen
        
        svinslen = Features.get_number_value(info, 'SVINSLEN', 0)
        if svinslen == 0 and svinsseq:
            svinslen = len(svinsseq)
        features['SVINSLEN'] = svinslen

        features['INS_SEQ_COV_PREFIX_START'], features['INS_SEQ_COV_PREFIX_END'] = 0, 1
        features['INS_SEQ_COV_SUFFIX_START'], features['INS_SEQ_COV_SUFFIX_END'] = 0, 1
        if source_str == "DE_NOVO_ASSEMBLY":
            if '-' in svinsseq:
                features['INS_SEQ_COV_PREFIX_END'] = svinsseq.index('-')/len(svinsseq)
                features['INS_SEQ_COV_SUFFIX_START'] = features['INS_SEQ_COV_PREFIX_END']
        elif source_str == "REFERENCE_GUIDED_ASSEMBLY":
            features['INS_SEQ_COV_PREFIX_START'], features['INS_SEQ_COV_PREFIX_END'] = Features.get_number_value(info, 'INS_PREFIX_COV', [0, 0], len(svinsseq)-1)
            features['INS_SEQ_COV_SUFFIX_START'], features['INS_SEQ_COV_SUFFIX_END'] = Features.get_number_value(info, 'INS_SUFFIX_COV', [0, 0], len(svinsseq)-1)

        split_reads = Features.get_number_value(info, 'SPLIT_READS', [0, 0])
        features['SPLIT_READS_RATIO1'], features['SPLIT_READS_RATIO2'] = Features.normalise(split_reads[0], min_depth, max_depth), Features.normalise(split_reads[1], min_depth, max_depth)
        features['SPLIT_READS_RATIO'] = sum(split_reads)/median_depth
        fwd_split_reads = Features.get_number_value(info, 'FWD_SPLIT_READS', [0, 0])
        features['FWD_SPLIT_READS_RATIO1'], features['FWD_SPLIT_READS_RATIO2'] = fwd_split_reads[0]/max(1, split_reads[0]), fwd_split_reads[1]/max(1, split_reads[1])
        rev_split_reads = Features.get_number_value(info, 'REV_SPLIT_READS', [0, 0])
        features['REV_SPLIT_READS_RATIO1'], features['REV_SPLIT_READS_RATIO2'] = rev_split_reads[0]/max(1, split_reads[0]), rev_split_reads[1]/max(1, split_reads[1])
        features['FWD_SPLIT_READS_RATIO'], features['REV_SPLIT_READS_RATIO'] = sum(fwd_split_reads)/max(1, sum(split_reads)), sum(rev_split_reads)/max(1, sum(split_reads))
        features['OVERLAP'] = Features.get_number_value(info, 'OVERLAP', 0, read_len)
        features['MISMATCH_RATE'] = Features.get_number_value(info, 'MISMATCH_RATE', 0)

        features['RCC_EXT_1SR_READS1'], features['RCC_EXT_1SR_READS2'] = Features.get_number_value(info, 'RCC_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['LCC_EXT_1SR_READS1'], features['LCC_EXT_1SR_READS2'] = Features.get_number_value(info, 'LCC_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['EXT_1SR_READS1'] = features['RCC_EXT_1SR_READS1'] + features['LCC_EXT_1SR_READS1']
        features['EXT_1SR_READS2'] = features['RCC_EXT_1SR_READS2'] + features['LCC_EXT_1SR_READS2']
        features['RCC_HQ_EXT_1SR_READS1'], features['RCC_HQ_EXT_1SR_READS2'] = Features.get_number_value(info, 'RCC_HQ_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['LCC_HQ_EXT_1SR_READS1'], features['LCC_HQ_EXT_1SR_READS2'] = Features.get_number_value(info, 'LCC_HQ_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['HQ_EXT_1SR_READS1'] = features['RCC_HQ_EXT_1SR_READS1'] + features['LCC_HQ_EXT_1SR_READS1']
        features['HQ_EXT_1SR_READS2'] = features['RCC_HQ_EXT_1SR_READS2'] + features['LCC_HQ_EXT_1SR_READS2']

        full_junction_score = Features.get_number_value(info, 'FULL_JUNCTION_SCORE', 0)
        split_junction_score1 = Features.get_number_value(info, 'SPLIT_JUNCTION_SCORE', [0, 0])
        split_junction_score2 = Features.get_number_value(info, 'SPLIT_JUNCTION_SCORE2', [0, 0])
        split_junction_size = Features.get_number_value(info, 'SPLIT_JUNCTION_SIZE', [0, 0])
        if sum(split_junction_score1) == 0:
            features['FULL_TO_SPLIT_JUNCTION_SCORE_RATIO'] = 0
            features['FULL_TO_SPLIT_JUNCTION_SCORE_DIFF'] = 0
            features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1'], features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2'] = 0, 0
            features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1'], features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2'] = 0, 0
            features['SPLIT_TO_SIZE_RATIO1'], features['SPLIT_TO_SIZE_RATIO2'] = 0, 0
        else:
            features['FULL_TO_SPLIT_JUNCTION_SCORE_RATIO'] = full_junction_score/sum(split_junction_score1)
            features['FULL_TO_SPLIT_JUNCTION_SCORE_DIFF'] = full_junction_score-sum(split_junction_score1)
            features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1'], features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2'] = split_junction_score2[0]/max(1, split_junction_score1[0]), split_junction_score2[1]/max(1, split_junction_score1[1])
            s1 = (split_junction_size[0]-split_junction_score1[0])/max(split_junction_size[0]-split_junction_score2[0], 1)
            s2 = (split_junction_size[1]-split_junction_score1[1])/max(split_junction_size[1]-split_junction_score2[1], 1)
            features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1'], features['SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2'] = s1, s2
            features['SPLIT_TO_SIZE_RATIO1'], features['SPLIT_TO_SIZE_RATIO2'] = split_junction_score1[0]/split_junction_size[0], split_junction_score1[1]/split_junction_size[1]
        features['SPLIT_JUNCTION_SIZE_RATIO1'], features['SPLIT_JUNCTION_SIZE_RATIO2'] = split_junction_size[0]/max(1, sum(split_junction_size)), split_junction_size[1]/max(1, sum(split_junction_size))
        features['MAX_SPLIT_JUNCTION_SIZE_RATIO'] = max(features['SPLIT_JUNCTION_SIZE_RATIO1'], features['SPLIT_JUNCTION_SIZE_RATIO2'])
        features['MIN_SPLIT_JUNCTION_SIZE_RATIO'] = min(features['SPLIT_JUNCTION_SIZE_RATIO1'], features['SPLIT_JUNCTION_SIZE_RATIO2'])

        features['MAX_MAPQ1'], features['MAX_MAPQ2'] = Features.get_number_value(info, 'MAX_MAPQ', [0, 0])
        features['MAX_MAPQ'] = max(features['MAX_MAPQ1'], features['MAX_MAPQ2'])
        features['MAX_MAPQ_EXT1'], features['MAX_MAPQ_EXT2'] = Features.get_number_value(info, 'MAX_MAPQ_EXT', [0, 0])
        features['MAX_MAPQ_EXT'] = max(features['MAX_MAPQ_EXT1'], features['MAX_MAPQ_EXT2'])
        remap_lb = Features.get_number_value(info, 'REMAP_LB', record.pos)
        remap_ub = Features.get_number_value(info, 'REMAP_UB', record.stop)
        features['LB_DIFF'], features['UB_DIFF'] = max(0, remap_lb-record.pos), max(0, record.stop-remap_ub)
        features['B_DIFF'] = features['LB_DIFF'] + features['UB_DIFF']

        disc_pairs = Features.get_number_value(info, 'DISC_PAIRS', [0, 0])
        if svtype_str == "DEL":
            min_is_to_become_disc = int(max(0, max_is-svlen))
            min_disc_pairs = stats['min_pairs_crossing_gap'][str(min_is_to_become_disc)]
            max_disc_pairs = stats['max_pairs_crossing_gap'][str(min_is_to_become_disc)]
            disc_pairs_scaled = [Features.normalise(d, min_disc_pairs, max_disc_pairs) for d in disc_pairs]
        elif svtype_str == "INS" and source_str in ("DE_NOVO_ASSEMBLY", "REFERENCE_GUIDED_ASSEMBLY"):
            min_inslen = int(min(max_is, svinslen))
            min_disc_pairs = stats['min_disc_pairs_by_insertion_size'][str(min_inslen)]
            max_disc_pairs = stats['max_disc_pairs_by_insertion_size'][str(min_inslen)]
            disc_pairs_scaled = [Features.normalise(d, min_disc_pairs, max_disc_pairs) for d in disc_pairs]
        else:
            disc_pairs_scaled = [d/median_depth for d in disc_pairs]
        features['DISC_PAIRS_SCALED1'], features['DISC_PAIRS_SCALED2'] = disc_pairs_scaled
        disc_pairs_highmapq = Features.get_number_value(info, 'DISC_PAIRS_HIGHMAPQ', [0, 0], median_depth)
        features['DISC_PAIRS_HIGHMAPQ_RATIO1'], features['DISC_PAIRS_HIGHMAPQ_RATIO2'] = [d/max(1, disc_pairs[i]) for i, d in enumerate(disc_pairs_highmapq)]
        features['DISC_PAIRS_MAXMAPQ1'], features['DISC_PAIRS_MAXMAPQ2'] = Features.get_number_value(info, 'DISC_PAIRS_MAXMAPQ', [0, 0])

        conc_pairs = Features.get_number_value(info, 'CONC_PAIRS', 0, median_depth)
        min_pairs_crossing_point = stats['min_pairs_crossing_gap']["0"]
        max_pairs_crossing_point = stats['max_pairs_crossing_gap']["0"]
        features['CONC_PAIRS_SCALED'] = Features.normalise(conc_pairs, min_pairs_crossing_point, max_pairs_crossing_point)

        features['DISC_PAIRS_SURROUNDING1'], features['DISC_PAIRS_SURROUNDING2'] = Features.get_number_value(info, 'DISC_PAIRS_SURROUNDING', [0, 0], median_depth) 
        features['DISC_AVG_NM1'], features['DISC_AVG_NM2'] = Features.get_number_value(info, 'DISC_AVG_NM', [0, 0], read_len)

        disc_pairs = [d/median_depth for d in disc_pairs]
        features['PTN_RATIO1'], features['PTN_RATIO2'] = [d/max(1, d+conc_pairs) for d in disc_pairs]
        features['KS_PVAL'] = max(0, Features.get_number_value(info, 'KS_PVAL', 1.0))
        features['KS_PVAL_HIGHMQ'] = max(0, Features.get_number_value(info, 'KS_PVAL_HIGHMQ', 1.0))
        features['SIZE_NORM'] = 2
        features['SIZE_NORM_HIGHMQ'] = 2
        if 'MAX_SIZE' in info:
            min_size = float(info['MIN_SIZE'])
            max_size = float(info['MAX_SIZE'])
            features['SIZE_NORM'] = Features.normalise(svlen/2, min_size, max_size)
        if 'MAX_SIZE_HIGHMQ' in info:
            min_size = float(info['MIN_SIZE_HIGHMQ'])
            max_size = float(info['MAX_SIZE_HIGHMQ'])
            features['SIZE_NORM_HIGHMQ'] = Features.normalise(svlen/2, min_size, max_size)

        median_depths = Features.get_number_value(info, 'MEDIAN_DEPTHS', [0, 0, 0, 0])
        features['MEDIAN_DEPTHS_NORM1'], features['MEDIAN_DEPTHS_NORM2'], features['MEDIAN_DEPTHS_NORM3'], features['MEDIAN_DEPTHS_NORM4'] = [Features.normalise(x, min_depth, max_depth) for x in median_depths]
        features['MEDIAN_DEPTHS_RATIO1'], features['MEDIAN_DEPTHS_RATIO2'] = median_depths[0]/max(1, median_depths[1]), median_depths[3]/max(1, median_depths[2])

        median_depths_highmq = Features.get_number_value(info, 'MEDIAN_DEPTHS_HIGHMQ', [0, 0, 0, 0])
        features['MEDIAN_DEPTHS_HIGHMQ_NORM1'], features['MEDIAN_DEPTHS_HIGHMQ_NORM2'], features['MEDIAN_DEPTHS_HIGHMQ_NORM3'], features['MEDIAN_DEPTHS_HIGHMQ_NORM4'] = [Features.normalise(x, min_depth, max_depth) for x in median_depths_highmq]
        features['MEDIAN_DEPTHS_HIGHMQ_RATIO1'], features['MEDIAN_DEPTHS_HIGHMQ_RATIO2'] = median_depths_highmq[0]/max(1, median_depths_highmq[1]), median_depths_highmq[3]/max(1, median_depths_highmq[2])

        if 'CLUSTER_DEPTHS' in info:
            features['CLUSTER_DEPTHS_ABOVE_MAX1'], features['CLUSTER_DEPTHS_ABOVE_MAX2'] = [max(0, float(x)-max_depth)/median_depth for x in info['CLUSTER_DEPTHS']]
        else:
            features['CLUSTER_DEPTHS_ABOVE_MAX1'], features['CLUSTER_DEPTHS_ABOVE_MAX2'] = 0, 0

        features['PREFIX_MH_LEN_RATIO'] = Features.get_number_value(info, 'PREFIX_MH_LEN', 0, max(1, svinslen))
        features['SUFFIX_MH_LEN_RATIO'] = Features.get_number_value(info, 'SUFFIX_MH_LEN', 0, max(1, svinslen))

        left_anchor_base_count = Features.get_number_value(info, 'LEFT_ANCHOR_BASE_COUNT', [0, 0, 0, 0])
        left_anchor_base_count_ratio = [x/max(1, sum(left_anchor_base_count)) for x in left_anchor_base_count]
        features['MAX_LEFT_ANCHOR_BASE_RATIO'] = max(left_anchor_base_count_ratio)
        features['LEFT_ANCHOR_A_RATIO'], features['LEFT_ANCHOR_C_RATIO'], features['LEFT_ANCHOR_G_RATIO'], features['LEFT_ANCHOR_T_RATIO'] = left_anchor_base_count_ratio

        right_anchor_base_count = Features.get_number_value(info, 'RIGHT_ANCHOR_BASE_COUNT', [0, 0, 0, 0])
        right_anchor_base_count_ratio = [x/max(1, sum(right_anchor_base_count)) for x in right_anchor_base_count]
        features['MAX_RIGHT_ANCHOR_BASE_RATIO'] = max(right_anchor_base_count_ratio)
        features['RIGHT_ANCHOR_A_RATIO'], features['RIGHT_ANCHOR_C_RATIO'], features['RIGHT_ANCHOR_G_RATIO'], features['RIGHT_ANCHOR_T_RATIO'] = right_anchor_base_count_ratio

        sv_ref_prefix_base_count = Features.get_number_value(info, 'SV_REF_PREFIX_BASE_COUNT', [0, 0, 0, 0])
        sv_ref_prefix_base_count_ratio = [x/max(1, sum(sv_ref_prefix_base_count)) for x in sv_ref_prefix_base_count]
        features['MAX_SV_REF_PREFIX_BASE_RATIO'] = max(sv_ref_prefix_base_count_ratio)
        features['SV_REF_PREFIX_A_RATIO'], features['SV_REF_PREFIX_C_RATIO'], features['SV_REF_PREFIX_G_RATIO'], features['SV_REF_PREFIX_T_RATIO'] = sv_ref_prefix_base_count_ratio

        sv_ref_suffix_base_count = Features.get_number_value(info, 'SV_REF_SUFFIX_BASE_COUNT', [0, 0, 0, 0])
        sv_ref_suffix_base_count_ratio = [x/max(1, sum(sv_ref_suffix_base_count)) for x in sv_ref_suffix_base_count]
        features['MAX_SV_REF_SUFFIX_BASE_RATIO'] = max(sv_ref_suffix_base_count_ratio)
        features['SV_REF_SUFFIX_A_RATIO'], features['SV_REF_SUFFIX_C_RATIO'], features['SV_REF_SUFFIX_G_RATIO'], features['SV_REF_SUFFIX_T_RATIO'] = sv_ref_suffix_base_count_ratio

        ins_prefix_base_count = Features.get_number_value(info, 'INS_PREFIX_BASE_COUNT', [0, 0, 0, 0])
        ins_prefix_base_count_ratio = [x/max(1, sum(ins_prefix_base_count)) for x in ins_prefix_base_count]
        features['MAX_INS_PREFIX_BASE_COUNT_RATIO'] = max(ins_prefix_base_count_ratio)
        features['INS_PREFIX_A_RATIO'], features['INS_PREFIX_C_RATIO'], features['INS_PREFIX_G_RATIO'], features['INS_PREFIX_T_RATIO'] = ins_prefix_base_count_ratio
        
        ins_suffix_base_count = Features.get_number_value(info, 'INS_SUFFIX_BASE_COUNT', [0, 0, 0, 0])
        ins_suffix_base_count_ratio = [x/max(1, sum(ins_suffix_base_count)) for x in ins_suffix_base_count]
        features['MAX_INS_SUFFIX_BASE_COUNT_RATIO'] = max(ins_suffix_base_count_ratio)
        features['INS_SUFFIX_A_RATIO'], features['INS_SUFFIX_C_RATIO'], features['INS_SUFFIX_G_RATIO'], features['INS_SUFFIX_T_RATIO'] = ins_suffix_base_count_ratio

        features['IP'] = True if "IMPRECISE" in record.info else False

        ar1 = Features.get_number_value(record.samples[0], 'AR1', 0)
        arc1 = Features.get_number_value(record.samples[0], 'ARC1', 0)
        ar2 = Features.get_number_value(record.samples[0], 'AR2', 0)
        arc2 = Features.get_number_value(record.samples[0], 'ARC2', 0)
        rr1 = Features.get_number_value(record.samples[0], 'RR1', 0)
        rrc1 = Features.get_number_value(record.samples[0], 'RRC1', 0)
        rr2 = Features.get_number_value(record.samples[0], 'RR2', 0)
        rrc2 = Features.get_number_value(record.samples[0], 'RRC2', 0)
        er = Features.get_number_value(record.samples[0], 'ER', 0)

        features['AR1'] = Features.normalise(ar1, min_depth, max_depth)
        features['ARC1'] = Features.normalise(arc1, min_depth, max_depth)
        features['ARCF1'] = Features.get_number_value(record.samples[0], 'ARCF1', 0, max(1, arc1))
        features['ARCR1'] = Features.get_number_value(record.samples[0], 'ARCR1', 0, max(1, arc1))
        features['MAXARCD1'] = max(features['ARCF1'], features['ARCR1'])
        features['ARCAS1'] = Features.get_number_value(record.samples[0], 'ARCAS1', 0)
        features['ARC1MQ'] = Features.get_number_value(record.samples[0], 'ARC1MQ', 0)
        features['ARC1HQ'] = Features.get_number_value(record.samples[0], 'ARC1HQ', 0, max(1, arc1))

        features['AR2'] = Features.normalise(ar2, min_depth, max_depth)
        features['ARC2'] = Features.normalise(arc2, min_depth, max_depth)
        features['ARCF2'] = Features.get_number_value(record.samples[0], 'ARCF2', 0, max(1, arc2))
        features['ARCR2'] = Features.get_number_value(record.samples[0], 'ARCR2', 0, max(1, arc2))
        features['MAXARCD2'] = max(features['ARCF2'], features['ARCR2'])
        features['ARCAS2'] = Features.get_number_value(record.samples[0], 'ARCAS2', 0)
        features['ARC2MQ'] = Features.get_number_value(record.samples[0], 'ARC2MQ', 0)
        features['ARC2HQ'] = Features.get_number_value(record.samples[0], 'ARC2HQ', 0, max(1, arc2))

        features['RR1'] = Features.normalise(rr1, min_depth, max_depth)
        features['RRC1'] = Features.normalise(rrc1, min_depth, max_depth)
        features['RR2'] = Features.normalise(rr2, min_depth, max_depth)
        features['RRC2'] = Features.normalise(rrc2, min_depth, max_depth)
        features['ER'] = Features.normalise(er, min_depth, max_depth)
        features['AR1_RATIO'] = ar1/max(1, ar1+rr1+er)
        features['ARC1_RATIO'] = arc1/max(1, arc1+rrc1+er)
        features['AR2_RATIO'] = ar2/max(1, ar2+rr2+er)
        features['ARC2_RATIO'] = arc2/max(1, arc2+rrc2+er)
        features['RR1_RATIO'] = rr1/max(1, ar1+rr1+er)
        features['RRC1_RATIO'] = rrc1/max(1, arc1+rrc1+er)
        features['RR2_RATIO'] = rr2/max(1, ar2+rr2+er)
        features['RRC2_RATIO'] = rrc2/max(1, arc2+rrc2+er)
        features['AR1_OVER_RR1'] = ar1/max(1, ar1+rr1)
        features['ARC1_OVER_RRC1'] = arc1/max(1, arc1+rrc1)
        features['RR1_OVER_AR1'] = rr1/max(1, ar1+rr1)
        features['RRC1_OVER_ARC1'] = rrc1/max(1, arc1+rrc1)
        features['AR2_OVER_RR2'] = ar2/max(1, ar2+rr2)
        features['ARC2_OVER_RRC2'] = arc2/max(1, arc2+rrc2)
        features['RR2_OVER_AR2'] = rr2/max(1, ar2+rr2)
        features['RRC2_OVER_ARC2'] = rrc2/max(1, arc2+rrc2)

        mdlf = Features.get_number_value(record.samples[0], 'MDLF', 0)
        features['MDLF'] = Features.normalise(mdlf, min_depth, max_depth)
        mdsp = Features.get_number_value(record.samples[0], 'MDSP', 0)
        features['MDSP'] = Features.normalise(mdsp, min_depth, max_depth)
        mdsf = Features.get_number_value(record.samples[0], 'MDSF', 0)
        features['MDSF'] = Features.normalise(mdsf, min_depth, max_depth)
        mdrf = Features.get_number_value(record.samples[0], 'MDRF', 0)
        features['MDRF'] = Features.normalise(mdrf, min_depth, max_depth)
        mdlc = Features.get_number_value(record.samples[0], 'MDLC', 0)
        features['MDLC'] = Features.normalise(mdlc, min_depth, max_depth)
        mdrc = Features.get_number_value(record.samples[0], 'MDRC', 0)
        features['MDRC'] = Features.normalise(mdrc, min_depth, max_depth)

        mdlfhq = Features.get_number_value(record.samples[0], 'MDLFHQ', 0)
        features['MDLFHQ'] = Features.normalise(mdlfhq, min_depth, max_depth)
        mdsphq = Features.get_number_value(record.samples[0], 'MDSPHQ', 0)
        features['MDSPHQ'] = Features.normalise(mdsphq, min_depth, max_depth)
        mdsfhq = Features.get_number_value(record.samples[0], 'MDSFHQ', 0)
        features['MDSFHQ'] = Features.normalise(mdsfhq, min_depth, max_depth)
        mdrfhq = Features.get_number_value(record.samples[0], 'MDRFHQ', 0)
        features['MDRFHQ'] = Features.normalise(mdrfhq, min_depth, max_depth)

        features['MDSP_OVER_MDLF'] = mdsp/max(1, mdlf)
        features['MDLF_OVER_MDSP'] = mdlf/max(1, mdsp)
        features['MDSF_OVER_MDRF'] = mdsf/max(1, mdrf)
        features['MDRF_OVER_MDSF'] = mdrf/max(1, mdsf)

        features['MDSP_OVER_MDLF_HQ'] = mdsphq/max(1, mdlfhq)
        features['MDLF_OVER_MDSP_HQ'] = mdlfhq/max(1, mdsphq)
        features['MDSF_OVER_MDRF_HQ'] = mdsfhq/max(1, mdrfhq)
        features['MDRF_OVER_MDSF_HQ'] = mdrfhq/max(1, mdsfhq)

        features['FMT_KSPVAL'] = max(0, Features.get_number_value(record.samples[0], 'KSPVAL', 1.0))
        features['FMT_KSPVAL_HQ'] = max(0, Features.get_number_value(record.samples[0], 'KSPVALHQ', 1.0))
        features['FMT_SIZE_NORM'] = 2
        features['FMT_SIZE_NORM_HQ'] = 2
        if 'MAXSIZE' in record.samples[0]:
            min_size = float(record.samples[0]['MINSIZE'])
            max_size = float(record.samples[0]['MAXSIZE'])
            features['FMT_SIZE_NORM'] = Features.normalise(svlen/2, min_size, max_size)
        if 'MAXSIZEHQ' in record.samples[0]:
            min_size = float(record.samples[0]['MINSIZEHQ'])
            max_size = float(record.samples[0]['MAXSIZEHQ'])
            features['FMT_SIZE_NORM_HQ'] = Features.normalise(svlen/2, min_size, max_size)

        dp1 = Features.get_number_value(record.samples[0], 'DP1', 0)
        dp2 = Features.get_number_value(record.samples[0], 'DP2', 0)
        if svtype_str == "DEL":
            min_is_to_become_disc = int(max(0, max_is-svlen))
            min_disc_pairs = stats['min_pairs_crossing_gap'][str(min_is_to_become_disc)]
            max_disc_pairs = stats['max_pairs_crossing_gap'][str(min_is_to_become_disc)]
            dp1_scaled = Features.normalise(dp1, min_disc_pairs, max_disc_pairs)
            dp2_scaled = Features.normalise(dp2, min_disc_pairs, max_disc_pairs)
        elif svtype_str == "INS" and source_str in ("DE_NOVO_ASSEMBLY", "REFERENCE_GUIDED_ASSEMBLY"):
            min_inslen = int(min(max_is, svinslen))
            min_disc_pairs = stats['min_disc_pairs_by_insertion_size'][str(min_inslen)]
            max_disc_pairs = stats['max_disc_pairs_by_insertion_size'][str(min_inslen)]
            dp1_scaled = Features.normalise(dp1, min_disc_pairs, max_disc_pairs)
            dp2_scaled = Features.normalise(dp2, min_disc_pairs, max_disc_pairs)
        else:
            dp1_scaled = dp1/median_depth
            dp2_scaled = dp2/median_depth

        features['DP1'] = dp1_scaled
        features['DP2'] = dp2_scaled
        features['DP1_HQ_RATIO'] = Features.get_number_value(record.samples[0], 'DP1HQ', 0, max(1, dp1))
        features['DP2_HQ_RATIO'] = Features.get_number_value(record.samples[0], 'DP2HQ', 0, max(1, dp2))
        features['DP1MQ'] = Features.get_number_value(record.samples[0], 'DP1MQ', 0)
        features['DP2MQ'] = Features.get_number_value(record.samples[0], 'DP2MQ', 0)
        features['DPLANM'] = Features.get_number_value(record.samples[0], 'DPLANM', 0)
        features['DPRANM'] = Features.get_number_value(record.samples[0], 'DPRANM', 0)
        features['DPSL'] = Features.get_number_value(record.samples[0], 'DPSL', 0, median_depth)
        features['DPSR'] = Features.get_number_value(record.samples[0], 'DPSR', 0, median_depth)

        cp = Features.get_number_value(record.samples[0], 'CP', 0)
        features['CP'] = Features.normalise(cp, min_pairs_crossing_point, max_pairs_crossing_point)
        features['PTNR'] = dp1/max(1, cp)

        features['AXR'] = Features.get_number_value(record.samples[0], 'AXR', 0, median_depth*max_is)
        features['AXRHQ'] = Features.get_number_value(record.samples[0], 'AXRHQ', 0, median_depth*max_is)
        exl = Features.get_number_value(record.samples[0], 'EXL', 0, max_is)
        features['EXL'] = exl
        exas = Features.get_number_value(record.samples[0], 'EXAS', 0)
        features['EXAS'] = exas/max(1, exl)
        exrs = Features.get_number_value(record.samples[0], 'EXRS', 0)
        features['EXRS'] = exrs/max(1, exl)
        features['EXAS_EXRS_RATIO'] = exas/max(1, exrs)
        features['EXAS_EXRS_DIFF'] = exas-exrs

        denovo_feature_values, regt_feature_values = [], []
        for feature_name in Features.get_denovo_feature_names(denovo_model_name):
            denovo_feature_values.append(features[feature_name])
        for feature_name in Features.get_regt_feature_names(regt_model_name):
            regt_feature_values.append(features[feature_name])

        return denovo_feature_values, regt_feature_values

def select_gt(gt1, gt2):
    if gt1 == "./." and gt2 != "./.":
        return gt2
    if gt1 != "./." and gt2 == "./.":
        return gt1
    elif gt1 == gt2:
        return gt1
    else:
        return "./."

def read_gts(file_path, tolerate_no_gts = False):
    if not os.path.exists(file_path) and tolerate_no_gts:
        return defaultdict(lambda: "0/0")
    with open(file_path, 'r') as file:
        gts = defaultdict(lambda: "0/0")
        for line in file:
            id, gt = line.strip().split()
            if id not in gts:
                gts[id] = gt
            else:
                gts[id] = select_gt(gts[id], gt)
    return gts

def get_stat(stats, stat_name, chrom):
    if chrom in stats[stat_name]:
        return stats[stat_name][chrom]
    return stats[stat_name]['.']

# Function to parse the VCF file and extract relevant features using pysam
def parse_vcf(vcf_fname, stats_fname, fp_fname, svtype, tolerate_no_gts = False):
    gts = read_gts(fp_fname, tolerate_no_gts=tolerate_no_gts)
    vcf_reader = pysam.VariantFile(vcf_fname)
    stats_reader = open(stats_fname, 'r')
    denovo_features_by_source, regt_features_by_source = defaultdict(list), defaultdict(list)
    denovo_gts_by_source, regt_gts_by_source = defaultdict(list), defaultdict(list)
    denovo_variant_ids_by_source, regt_variant_ids_by_source = defaultdict(list), defaultdict(list)

    # read the stats file and extract the relevant values
    stats = defaultdict(dict)
    for line in stats_reader:
        sl = line.strip().split()
        stats[sl[0]][sl[1]] = int(sl[2])

    read_len = stats['read_len']['.']
    for record in vcf_reader.fetch():
        if svtype != 'ALL' and record.info['SVTYPE'] != svtype:
            continue

        denovo_model_name = Features.get_denovo_model_name(record, stats['max_is']['.'])
        regt_model_name = Features.get_regt_model_name(record, stats['max_is']['.'], stats['read_len']['.'])
        denovo_feature_values, regt_feature_values = Features.record_to_features(record, stats)
        denovo_features_by_source[denovo_model_name].append(denovo_feature_values)
        denovo_gts_by_source[denovo_model_name].append(gts[record.id])
        denovo_variant_ids_by_source[denovo_model_name].append(record.id)
        if 'TD' not in record.samples[0] and gts[record.id] != "./.": # if too deep or no genotype is available, skip the record
            regt_features_by_source[regt_model_name].append(regt_feature_values)
            regt_gts_by_source[regt_model_name].append(gts[record.id])
            regt_variant_ids_by_source[regt_model_name].append(record.id)
    for model_name in denovo_features_by_source:
        denovo_features_by_source[model_name] = np.array(denovo_features_by_source[model_name])
        denovo_gts_by_source[model_name] = np.array(denovo_gts_by_source[model_name])
        denovo_variant_ids_by_source[model_name] = np.array(denovo_variant_ids_by_source[model_name])
    for model_name in regt_features_by_source:
        regt_features_by_source[model_name] = np.array(regt_features_by_source[model_name])
        regt_gts_by_source[model_name] = np.array(regt_gts_by_source[model_name])
        regt_variant_ids_by_source[model_name] = np.array(regt_variant_ids_by_source[model_name])
    return denovo_features_by_source, regt_features_by_source, denovo_gts_by_source, regt_gts_by_source, \
        denovo_variant_ids_by_source, regt_variant_ids_by_source
