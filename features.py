import os, pysam
from collections import defaultdict
import numpy as np

class Features:

    svtype_to_int = { "DEL" : 0, "DUP" : 1, "INS" : 2 }

    shared_features_names = ['START_STOP_DIST', 'SVLEN', 'SVINSLEN',
                             'SV_REF_PREFIX_A_RATIO', 'SV_REF_PREFIX_C_RATIO', 'SV_REF_PREFIX_G_RATIO', 'SV_REF_PREFIX_T_RATIO', 'MAX_SV_REF_PREFIX_BASE_RATIO',
                             'SV_REF_SUFFIX_A_RATIO', 'SV_REF_SUFFIX_C_RATIO', 'SV_REF_SUFFIX_G_RATIO', 'SV_REF_SUFFIX_T_RATIO', 'MAX_SV_REF_SUFFIX_BASE_RATIO',
                             'LEFT_ANCHOR_A_RATIO', 'LEFT_ANCHOR_C_RATIO', 'LEFT_ANCHOR_G_RATIO', 'LEFT_ANCHOR_T_RATIO', 'MAX_LEFT_ANCHOR_BASE_RATIO',
                             'RIGHT_ANCHOR_A_RATIO', 'RIGHT_ANCHOR_C_RATIO', 'RIGHT_ANCHOR_G_RATIO', 'RIGHT_ANCHOR_T_RATIO', 'MAX_RIGHT_ANCHOR_BASE_RATIO',
                             'INS_PREFIX_A_RATIO', 'INS_PREFIX_C_RATIO', 'INS_PREFIX_G_RATIO', 'INS_PREFIX_T_RATIO', 
                             'INS_SUFFIX_A_RATIO', 'INS_SUFFIX_C_RATIO', 'INS_SUFFIX_G_RATIO', 'INS_SUFFIX_T_RATIO', 
                             'MAX_INS_PREFIX_BASE_COUNT_RATIO', 'MAX_INS_SUFFIX_BASE_COUNT_RATIO']

    features_names = ['INS_SEQ_COV_PREFIX_START', 'INS_SEQ_COV_PREFIX_END', 'INS_SEQ_COV_SUFFIX_START', 'INS_SEQ_COV_SUFFIX_END',
                      'SPLIT_READS_RATIO1', 'SPLIT_READS_RATIO2', 'FWD_SPLIT_READS_RATIO1', 'FWD_SPLIT_READS_RATIO2', 'REV_SPLIT_READS_RATIO1', 'REV_SPLIT_READS_RATIO2',
                      'FWD_SPLIT_READS_RATIO', 'REV_SPLIT_READS_RATIO', 'OVERLAP', 'MISMATCH_RATE', 'RCC_EXT_1SR_READS1', 'RCC_EXT_1SR_READS2', 'LCC_EXT_1SR_READS1', 'LCC_EXT_1SR_READS2',
                      'RCC_HQ_EXT_1SR_READS1', 'RCC_HQ_EXT_1SR_READS2', 'LCC_HQ_EXT_1SR_READS1', 'LCC_HQ_EXT_1SR_READS2', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                      'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                      'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                      'MAX_MAPQ1', 'MAX_MAPQ2', 'LB_DIFF', 'UB_DIFF', 'DISC_PAIRS_SCALED1', 'DISC_PAIRS_SCALED2', 'DISC_PAIRS_HIGHMAPQ_RATIO1', 'DISC_PAIRS_HIGHMAPQ_RATIO2', 'DISC_PAIRS_MAXMAPQ1', 'DISC_PAIRS_MAXMAPQ2',
                      'CONC_PAIRS_SCALED', 'DISC_PAIRS_SURROUNDING1', 'DISC_PAIRS_SURROUNDING2', 'DISC_AVG_NM1', 'DISC_AVG_NM2', 'PTN_RATIO1', 'PTN_RATIO2', 'KS_PVAL', 'MAX_SIZE_DIFF', 'MEDIAN_DEPTHS_NORM1', 'MEDIAN_DEPTHS_NORM2',
                      'MEDIAN_DEPTHS_NORM3', 'MEDIAN_DEPTHS_NORM4', 'MEDIAN_DEPTHS_RATIO1', 'MEDIAN_DEPTHS_RATIO2', 'MEDIAN_DEPTHS_ABOVE_MAX1', 'MEDIAN_DEPTHS_ABOVE_MAX2', 'MEDIAN_DEPTHS_ABOVE_MAX3', 'MEDIAN_DEPTHS_ABOVE_MAX4',
                      'MEDIAN_DEPTHS_BELOW_MIN1', 'MEDIAN_DEPTHS_BELOW_MIN2', 'CLUSTER_DEPTHS_ABOVE_MAX1', 'CLUSTER_DEPTHS_ABOVE_MAX2', 'PREFIX_MH_LEN_RATIO', 'SUFFIX_MH_LEN_RATIO']

    def get_denovo_model_name(record, max_is):
        svtype_str = record.info['SVTYPE']
        source_str = record.info['SOURCE']
        if svtype_str == "DEL":
            split_reads = Features.get_value(record.info, 'SPLIT_READS', [0, 0])
            if split_reads[0] > 0 and split_reads[1] > 0:
                source_str = "2SR"
            elif split_reads[0] > 0 or split_reads[1] > 0:
                source_str = "1SR"
            else:
                source_str = "DP"
            if abs(record.info['SVLEN']) >= max_is:
                source_str += "_LARGE"
        return svtype_str + "_" + source_str

    def get_denovo_feature_names(model_name):
        if model_name == "DEL_DP_LARGE":
            return Features.shared_features_names + \
                    ['OVERLAP', 'MISMATCH_RATE', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                    'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                    'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                    'DISC_PAIRS_SCALED1', 'DISC_PAIRS_SCALED2', 'DISC_PAIRS_HIGHMAPQ_RATIO1', 'DISC_PAIRS_HIGHMAPQ_RATIO2', 'DISC_PAIRS_MAXMAPQ1', 'DISC_PAIRS_MAXMAPQ2',
                    'CONC_PAIRS_SCALED', 'DISC_PAIRS_SURROUNDING1', 'DISC_PAIRS_SURROUNDING2', 'DISC_AVG_NM1', 'DISC_AVG_NM2', 'PTN_RATIO1', 'PTN_RATIO2', 'MEDIAN_DEPTHS_NORM1', 'MEDIAN_DEPTHS_NORM2',
                    'MEDIAN_DEPTHS_NORM3', 'MEDIAN_DEPTHS_NORM4', 'MEDIAN_DEPTHS_RATIO1', 'MEDIAN_DEPTHS_RATIO2', 'MEDIAN_DEPTHS_ABOVE_MAX1', 'MEDIAN_DEPTHS_ABOVE_MAX2', 'MEDIAN_DEPTHS_ABOVE_MAX3', 'MEDIAN_DEPTHS_ABOVE_MAX4',
                    'MEDIAN_DEPTHS_BELOW_MIN1', 'MEDIAN_DEPTHS_BELOW_MIN2', 'CLUSTER_DEPTHS_ABOVE_MAX1', 'CLUSTER_DEPTHS_ABOVE_MAX2']
        elif model_name == "DEL_DP":
            return Features.shared_features_names + \
                    ['OVERLAP', 'MISMATCH_RATE', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                     'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                     'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                     'DISC_PAIRS_SCALED1', 'DISC_PAIRS_SCALED2', 'DISC_PAIRS_HIGHMAPQ_RATIO1', 'DISC_PAIRS_HIGHMAPQ_RATIO2', 'DISC_PAIRS_MAXMAPQ1', 'DISC_PAIRS_MAXMAPQ2',
                     'DISC_PAIRS_SURROUNDING1', 'DISC_PAIRS_SURROUNDING2', 'DISC_AVG_NM1', 'DISC_AVG_NM2', 'PTN_RATIO1', 'PTN_RATIO2', 'KS_PVAL', 'MAX_SIZE_DIFF', 'MEDIAN_DEPTHS_NORM1', 'MEDIAN_DEPTHS_NORM2',
                     'MEDIAN_DEPTHS_NORM3', 'MEDIAN_DEPTHS_NORM4', 'MEDIAN_DEPTHS_RATIO1', 'MEDIAN_DEPTHS_RATIO2', 'MEDIAN_DEPTHS_ABOVE_MAX1', 'MEDIAN_DEPTHS_ABOVE_MAX2', 'MEDIAN_DEPTHS_ABOVE_MAX3', 'MEDIAN_DEPTHS_ABOVE_MAX4',
                     'MEDIAN_DEPTHS_BELOW_MIN1', 'MEDIAN_DEPTHS_BELOW_MIN2', 'CLUSTER_DEPTHS_ABOVE_MAX1', 'CLUSTER_DEPTHS_ABOVE_MAX2', 
                     'PREFIX_MH_LEN_RATIO', 'SUFFIX_MH_LEN_RATIO']
        elif model_name == "DEL_2SR":
            return Features.shared_features_names + \
                    ['SPLIT_READS_RATIO1', 'SPLIT_READS_RATIO2', 'FWD_SPLIT_READS_RATIO1', 'FWD_SPLIT_READS_RATIO2', 'REV_SPLIT_READS_RATIO1', 'REV_SPLIT_READS_RATIO2',
                      'FWD_SPLIT_READS_RATIO', 'REV_SPLIT_READS_RATIO', 'OVERLAP', 'MISMATCH_RATE', 'RCC_EXT_1SR_READS1', 'LCC_EXT_1SR_READS2',
                      'RCC_HQ_EXT_1SR_READS1', 'LCC_HQ_EXT_1SR_READS2', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                      'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                      'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                      'MAX_MAPQ1', 'MAX_MAPQ2', 'LB_DIFF', 'UB_DIFF', 'DISC_PAIRS_SCALED1', 'DISC_PAIRS_SCALED2', 'DISC_PAIRS_HIGHMAPQ_RATIO1', 'DISC_PAIRS_HIGHMAPQ_RATIO2', 'DISC_PAIRS_MAXMAPQ1', 'DISC_PAIRS_MAXMAPQ2',
                      'CONC_PAIRS_SCALED', 'DISC_PAIRS_SURROUNDING1', 'DISC_PAIRS_SURROUNDING2', 'DISC_AVG_NM1', 'DISC_AVG_NM2', 'PTN_RATIO1', 'PTN_RATIO2', 'KS_PVAL', 'MAX_SIZE_DIFF', 'MEDIAN_DEPTHS_NORM1', 'MEDIAN_DEPTHS_NORM2',
                      'MEDIAN_DEPTHS_NORM3', 'MEDIAN_DEPTHS_NORM4', 'MEDIAN_DEPTHS_RATIO1', 'MEDIAN_DEPTHS_RATIO2', 'MEDIAN_DEPTHS_ABOVE_MAX1', 'MEDIAN_DEPTHS_ABOVE_MAX2', 'MEDIAN_DEPTHS_ABOVE_MAX3', 'MEDIAN_DEPTHS_ABOVE_MAX4',
                      'MEDIAN_DEPTHS_BELOW_MIN1', 'MEDIAN_DEPTHS_BELOW_MIN2', 'CLUSTER_DEPTHS_ABOVE_MAX1', 'CLUSTER_DEPTHS_ABOVE_MAX2', 
                      'PREFIX_MH_LEN_RATIO', 'SUFFIX_MH_LEN_RATIO']
        elif model_name == "DEL_1SR":
            return Features.shared_features_names + \
                    ['SPLIT_READS_RATIO', 'FWD_SPLIT_READS_RATIO', 'REV_SPLIT_READS_RATIO', 'EXT_1SR_READS1', 'EXT_1SR_READS2',
                     'HQ_EXT_1SR_READS1', 'HQ_EXT_1SR_READS2', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                     'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                     'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                     'MAX_MAPQ', 'B_DIFF', 'DISC_PAIRS_SCALED1', 'DISC_PAIRS_SCALED2', 'DISC_PAIRS_HIGHMAPQ_RATIO1', 'DISC_PAIRS_HIGHMAPQ_RATIO2', 'DISC_PAIRS_MAXMAPQ1', 'DISC_PAIRS_MAXMAPQ2',
                     'CONC_PAIRS_SCALED', 'DISC_PAIRS_SURROUNDING1', 'DISC_PAIRS_SURROUNDING2', 'DISC_AVG_NM1', 'DISC_AVG_NM2', 'PTN_RATIO1', 'PTN_RATIO2', 'KS_PVAL', 'MAX_SIZE_DIFF', 'MEDIAN_DEPTHS_NORM1', 'MEDIAN_DEPTHS_NORM2',
                     'MEDIAN_DEPTHS_NORM3', 'MEDIAN_DEPTHS_NORM4', 'MEDIAN_DEPTHS_RATIO1', 'MEDIAN_DEPTHS_RATIO2', 'MEDIAN_DEPTHS_ABOVE_MAX1', 'MEDIAN_DEPTHS_ABOVE_MAX2', 'MEDIAN_DEPTHS_ABOVE_MAX3', 'MEDIAN_DEPTHS_ABOVE_MAX4',
                     'MEDIAN_DEPTHS_BELOW_MIN1', 'MEDIAN_DEPTHS_BELOW_MIN2', 'CLUSTER_DEPTHS_ABOVE_MAX1', 'CLUSTER_DEPTHS_ABOVE_MAX2']
        else:
            return Features.shared_features_names + Features.features_names
    
    def get_regt_feature_names(model_name):
        return Features.shared_features_names

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
        
        model_name = Features.get_denovo_model_name(record, max_is)
        features = dict()

        info = record.info
        svtype_str = record.info['SVTYPE']
        source_str = record.info['SOURCE']
        features['START_STOP_DIST'] = record.stop - record.pos

        svinsseq = ""
        if 'SVINSSEQ' in info:
            svinsseq = info['SVINSSEQ']
        if 'SVLEN' in info:
            svlen = abs(float(info['SVLEN']))
        else:
            svlen = len(svinsseq) - (record.stop - record.pos)
        features['SVLEN'] = svlen
        
        svinslen = Features.get_value(info, 'SVINSLEN', 0)
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
            features['INS_SEQ_COV_PREFIX_START'], features['INS_SEQ_COV_PREFIX_END'] = Features.get_value(info, 'INS_PREFIX_COV', [0, 0], len(svinsseq)-1)
            features['INS_SEQ_COV_SUFFIX_START'], features['INS_SEQ_COV_SUFFIX_END'] = Features.get_value(info, 'INS_SUFFIX_COV', [0, 0], len(svinsseq)-1)

        split_reads = Features.get_value(info, 'SPLIT_READS', [0, 0])
        features['SPLIT_READS_RATIO1'], features['SPLIT_READS_RATIO2'] = split_reads[0]/median_depth, split_reads[1]/median_depth
        features['SPLIT_READS_RATIO'] = sum(split_reads)/median_depth
        fwd_split_reads = Features.get_value(info, 'FWD_SPLIT_READS', [0, 0])
        features['FWD_SPLIT_READS_RATIO1'], features['FWD_SPLIT_READS_RATIO2'] = fwd_split_reads[0]/max(1, split_reads[0]), fwd_split_reads[1]/max(1, split_reads[1])
        rev_split_reads = Features.get_value(info, 'REV_SPLIT_READS', [0, 0])
        features['REV_SPLIT_READS_RATIO1'], features['REV_SPLIT_READS_RATIO2'] = rev_split_reads[0]/max(1, split_reads[0]), rev_split_reads[1]/max(1, split_reads[1])
        features['FWD_SPLIT_READS_RATIO'], features['REV_SPLIT_READS_RATIO'] = sum(fwd_split_reads)/max(1, sum(split_reads)), sum(rev_split_reads)/max(1, sum(split_reads))
        features['OVERLAP'] = Features.get_value(info, 'OVERLAP', 0, read_len)
        features['MISMATCH_RATE'] = Features.get_value(info, 'MISMATCH_RATE', 0)
        
        features['RCC_EXT_1SR_READS1'], features['RCC_EXT_1SR_READS2'] = Features.get_value(info, 'RCC_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['LCC_EXT_1SR_READS1'], features['LCC_EXT_1SR_READS2'] = Features.get_value(info, 'LCC_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['EXT_1SR_READS1'] = features['RCC_EXT_1SR_READS1'] + features['LCC_EXT_1SR_READS1']
        features['EXT_1SR_READS2'] = features['RCC_EXT_1SR_READS2'] + features['LCC_EXT_1SR_READS2']
        features['RCC_HQ_EXT_1SR_READS1'], features['RCC_HQ_EXT_1SR_READS2'] = Features.get_value(info, 'RCC_HQ_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['LCC_HQ_EXT_1SR_READS1'], features['LCC_HQ_EXT_1SR_READS2'] = Features.get_value(info, 'LCC_HQ_EXT_1SR_READS', [0, 0], median_depth*max_is)
        features['HQ_EXT_1SR_READS1'] = features['RCC_HQ_EXT_1SR_READS1'] + features['LCC_HQ_EXT_1SR_READS1']
        features['HQ_EXT_1SR_READS2'] = features['RCC_HQ_EXT_1SR_READS2'] + features['LCC_HQ_EXT_1SR_READS2']

        full_junction_score = Features.get_value(info, 'FULL_JUNCTION_SCORE', 0)
        split_junction_score1 = [float(x) for x in info['SPLIT_JUNCTION_SCORE']]
        split_junction_score2 = [float(x) for x in info['SPLIT_JUNCTION_SCORE2']]
        split_junction_size = [float(x) for x in info['SPLIT_JUNCTION_SIZE']]
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
        features['SPLIT_JUNCTION_SIZE_RATIO1'], features['SPLIT_JUNCTION_SIZE_RATIO2'] = split_junction_size[0]/sum(split_junction_size), split_junction_size[1]/sum(split_junction_size)
        features['MAX_SPLIT_JUNCTION_SIZE_RATIO'] = max(features['SPLIT_JUNCTION_SIZE_RATIO1'], features['SPLIT_JUNCTION_SIZE_RATIO2'])
        features['MIN_SPLIT_JUNCTION_SIZE_RATIO'] = min(features['SPLIT_JUNCTION_SIZE_RATIO1'], features['SPLIT_JUNCTION_SIZE_RATIO2'])

        features['MAX_MAPQ1'], features['MAX_MAPQ2'] = Features.get_value(info, 'MAX_MAPQ', [0, 0])
        features['MAX_MAPQ'] = max(features['MAX_MAPQ1'], features['MAX_MAPQ2'])
        remap_lb = Features.get_value(info, 'REMAP_LB', record.pos)
        remap_ub = Features.get_value(info, 'REMAP_UB', record.stop)
        features['LB_DIFF'], features['UB_DIFF'] = max(0, remap_lb-record.pos), max(0, record.stop-remap_ub)
        features['B_DIFF'] = features['LB_DIFF'] + features['UB_DIFF']

        disc_pairs = Features.get_value(info, 'DISC_PAIRS', [0, 0])
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
        features['DISC_PAIRS_SCALED1'], features['DISC_PAIRS_SCALED2'] = disc_pairs_scaled
        disc_pairs_highmapq = Features.get_value(info, 'DISC_PAIRS_HIGHMAPQ', [0, 0], median_depth)
        features['DISC_PAIRS_HIGHMAPQ_RATIO1'], features['DISC_PAIRS_HIGHMAPQ_RATIO2'] = [d/max(1, disc_pairs[i]) for i, d in enumerate(disc_pairs_highmapq)]
        features['DISC_PAIRS_MAXMAPQ1'], features['DISC_PAIRS_MAXMAPQ2'] = Features.get_value(info, 'DISC_PAIRS_MAXMAPQ', [0, 0])

        conc_pairs = Features.get_value(info, 'CONC_PAIRS', 0, median_depth)
        min_pairs_crossing_point = stats['min_pairs_crossing_gap']["0"]
        max_pairs_crossing_point = stats['max_pairs_crossing_gap']["0"]
        conc_pairs_scaled = (conc_pairs-min_pairs_crossing_point)/(max_pairs_crossing_point-min_pairs_crossing_point)
        features['CONC_PAIRS_SCALED'] = conc_pairs_scaled

        features['DISC_PAIRS_SURROUNDING1'], features['DISC_PAIRS_SURROUNDING2'] = Features.get_value(info, 'DISC_PAIRS_SURROUNDING', [0, 0], median_depth) 
        features['DISC_AVG_NM1'], features['DISC_AVG_NM2'] = Features.get_value(info, 'DISC_AVG_NM', [0, 0], read_len)

        disc_pairs = [d/median_depth for d in disc_pairs]
        features['PTN_RATIO1'], features['PTN_RATIO2'] = [d/max(1, d+conc_pairs) for d in disc_pairs]
        features['KS_PVAL'] = max(0, Features.get_value(info, 'KS_PVAL', 0.0))
        features['MAX_SIZE_DIFF'] = 0
        if 'MAX_SIZE' in info:
            max_size = float(info['MAX_SIZE'])
            features['MAX_SIZE_DIFF'] = max(0, svlen - 2*max_size)

        median_depths = [float(x) for x in info['MEDIAN_DEPTHS']]
        features['MEDIAN_DEPTHS_NORM1'], features['MEDIAN_DEPTHS_NORM2'], features['MEDIAN_DEPTHS_NORM3'], features['MEDIAN_DEPTHS_NORM4'] = [float(x)/median_depth for x in info['MEDIAN_DEPTHS']]
        features['MEDIAN_DEPTHS_RATIO1'], features['MEDIAN_DEPTHS_RATIO2'] = median_depths[0]/max(1, median_depths[1]), median_depths[3]/max(1, median_depths[2])
        features['MEDIAN_DEPTHS_ABOVE_MAX1'], features['MEDIAN_DEPTHS_ABOVE_MAX2'], features['MEDIAN_DEPTHS_ABOVE_MAX3'], features['MEDIAN_DEPTHS_ABOVE_MAX4'] = [max(0, x-max_depth)/median_depth for x in median_depths]
        features['MEDIAN_DEPTHS_BELOW_MIN1'], features['MEDIAN_DEPTHS_BELOW_MIN2'] = max(0, min_depth-median_depths[0])/median_depth, max(0, min_depth-median_depths[3])/median_depth
        if 'CLUSTER_DEPTHS' in info:
            features['CLUSTER_DEPTHS_ABOVE_MAX1'], features['CLUSTER_DEPTHS_ABOVE_MAX2'] = [max(0, float(x)-max_depth)/median_depth for x in info['CLUSTER_DEPTHS']]
        else:
            features['CLUSTER_DEPTHS_ABOVE_MAX1'], features['CLUSTER_DEPTHS_ABOVE_MAX2'] = 0, 0

        features['PREFIX_MH_LEN_RATIO'] = Features.get_value(info, 'PREFIX_MH_LEN', 0, max(1, svinslen))
        features['SUFFIX_MH_LEN_RATIO'] = Features.get_value(info, 'SUFFIX_MH_LEN', 0, max(1, svinslen))

        left_anchor_base_count = Features.get_value(info, 'LEFT_ANCHOR_BASE_COUNT', [0, 0, 0, 0])
        left_anchor_base_count_ratio = [x/max(1, sum(left_anchor_base_count)) for x in left_anchor_base_count]
        features['MAX_LEFT_ANCHOR_BASE_RATIO'] = max(left_anchor_base_count_ratio)
        features['LEFT_ANCHOR_A_RATIO'], features['LEFT_ANCHOR_C_RATIO'], features['LEFT_ANCHOR_G_RATIO'], features['LEFT_ANCHOR_T_RATIO'] = left_anchor_base_count_ratio

        right_anchor_base_count = Features.get_value(info, 'RIGHT_ANCHOR_BASE_COUNT', [0, 0, 0, 0])
        right_anchor_base_count_ratio = [x/max(1, sum(right_anchor_base_count)) for x in right_anchor_base_count]
        features['MAX_RIGHT_ANCHOR_BASE_RATIO'] = max(right_anchor_base_count_ratio)
        features['RIGHT_ANCHOR_A_RATIO'], features['RIGHT_ANCHOR_C_RATIO'], features['RIGHT_ANCHOR_G_RATIO'], features['RIGHT_ANCHOR_T_RATIO'] = right_anchor_base_count_ratio

        sv_ref_prefix_base_count = Features.get_value(info, 'SV_REF_PREFIX_BASE_COUNT', [0, 0, 0, 0])
        sv_ref_prefix_base_count_ratio = [x/max(1, sum(sv_ref_prefix_base_count)) for x in sv_ref_prefix_base_count]
        features['MAX_SV_REF_PREFIX_BASE_RATIO'] = max(sv_ref_prefix_base_count_ratio)
        features['SV_REF_PREFIX_A_RATIO'], features['SV_REF_PREFIX_C_RATIO'], features['SV_REF_PREFIX_G_RATIO'], features['SV_REF_PREFIX_T_RATIO'] = sv_ref_prefix_base_count_ratio

        sv_ref_suffix_base_count = Features.get_value(info, 'SV_REF_SUFFIX_BASE_COUNT', [0, 0, 0, 0])
        sv_ref_suffix_base_count_ratio = [x/max(1, sum(sv_ref_suffix_base_count)) for x in sv_ref_suffix_base_count]
        features['MAX_SV_REF_SUFFIX_BASE_RATIO'] = max(sv_ref_suffix_base_count_ratio)
        features['SV_REF_SUFFIX_A_RATIO'], features['SV_REF_SUFFIX_C_RATIO'], features['SV_REF_SUFFIX_G_RATIO'], features['SV_REF_SUFFIX_T_RATIO'] = sv_ref_suffix_base_count_ratio

        ins_prefix_base_count = Features.get_value(info, 'INS_PREFIX_BASE_COUNT', [0, 0, 0, 0])
        ins_prefix_base_count_ratio = [x/max(1, sum(ins_prefix_base_count)) for x in ins_prefix_base_count]
        features['MAX_INS_PREFIX_BASE_COUNT_RATIO'] = max(ins_prefix_base_count_ratio)
        features['INS_PREFIX_A_RATIO'], features['INS_PREFIX_C_RATIO'], features['INS_PREFIX_G_RATIO'], features['INS_PREFIX_T_RATIO'] = ins_prefix_base_count_ratio
        
        ins_suffix_base_count = Features.get_value(info, 'INS_SUFFIX_BASE_COUNT', [0, 0, 0, 0])
        ins_suffix_base_count_ratio = [x/max(1, sum(ins_suffix_base_count)) for x in ins_suffix_base_count]
        features['MAX_INS_SUFFIX_BASE_COUNT_RATIO'] = max(ins_suffix_base_count_ratio)
        features['INS_SUFFIX_A_RATIO'], features['INS_SUFFIX_C_RATIO'], features['INS_SUFFIX_G_RATIO'], features['INS_SUFFIX_T_RATIO'] = ins_suffix_base_count_ratio

        denovo_feature_values, regt_feature_values = [], []
        for feature_name in Features.get_denovo_feature_names(model_name):
            denovo_feature_values.append(features[feature_name])
        for feature_name in Features.get_regt_feature_names(model_name):
            regt_feature_values.append(features[feature_name])
        return denovo_feature_values, regt_feature_values

def read_gts(file_path, tolerate_no_gts = False):
    if not os.path.exists(file_path) and tolerate_no_gts:
        return defaultdict(lambda: "0/0")
    with open(file_path, 'r') as file:
        gts = defaultdict(lambda: "0/0")
        for line in file:
            sl = line.strip().split()
            gts[sl[0]] = sl[1]
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
    data_by_source = defaultdict(list)
    gts_by_source = defaultdict(list)
    variant_ids_by_source = defaultdict(list)

    # read the stats file and extract the relevant values
    stats = defaultdict(dict)
    for line in stats_reader:
        sl = line.strip().split()
        stats[sl[0]][sl[1]] = int(sl[2])

    for record in vcf_reader.fetch():
        if svtype != 'ALL' and record.info['SVTYPE'] != svtype:
            continue

        model_name = Features.get_denovo_model_name(record, stats['max_is']['.'])
        denovo_feature_values, regt_feature_values = Features.record_to_features(record, stats)
        data_by_source[model_name].append(denovo_feature_values)
        gts_by_source[model_name].append(gts[record.id])
        variant_ids_by_source[model_name].append(record.id)
    for model_name in data_by_source:
        data_by_source[model_name] = np.array(data_by_source[model_name])
        gts_by_source[model_name] = np.array(gts_by_source[model_name])
        variant_ids_by_source[model_name] = np.array(variant_ids_by_source[model_name])
    return data_by_source, gts_by_source, variant_ids_by_source
