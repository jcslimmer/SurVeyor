from __future__ import division
import os, pysam
from collections import defaultdict
import numpy as np

class Features:

    info_features_names = [ 'START_STOP_DIST', 'SVLEN', 'SVINSLEN', 
                            'SV_REF_PREFIX_A_RATIO', 'SV_REF_PREFIX_C_RATIO', 'SV_REF_PREFIX_G_RATIO', 'SV_REF_PREFIX_T_RATIO', 'MAX_SV_REF_PREFIX_BASE_RATIO',
                            'SV_REF_SUFFIX_A_RATIO', 'SV_REF_SUFFIX_C_RATIO', 'SV_REF_SUFFIX_G_RATIO', 'SV_REF_SUFFIX_T_RATIO', 'MAX_SV_REF_SUFFIX_BASE_RATIO',
                            'LEFT_ANCHOR_A_RATIO', 'LEFT_ANCHOR_C_RATIO', 'LEFT_ANCHOR_G_RATIO', 'LEFT_ANCHOR_T_RATIO', 'MAX_LEFT_ANCHOR_BASE_RATIO',
                            'RIGHT_ANCHOR_A_RATIO', 'RIGHT_ANCHOR_C_RATIO', 'RIGHT_ANCHOR_G_RATIO', 'RIGHT_ANCHOR_T_RATIO', 'MAX_RIGHT_ANCHOR_BASE_RATIO',
                            'INS_PREFIX_A_RATIO', 'INS_PREFIX_C_RATIO', 'INS_PREFIX_G_RATIO', 'INS_PREFIX_T_RATIO', 'MAX_INS_PREFIX_BASE_COUNT_RATIO',
                            'INS_SUFFIX_A_RATIO', 'INS_SUFFIX_C_RATIO', 'INS_SUFFIX_G_RATIO', 'INS_SUFFIX_T_RATIO', 'MAX_INS_SUFFIX_BASE_COUNT_RATIO',
                            'INS_SEQ_COV_PREFIX_LEN', 'INS_SEQ_COV_SUFFIX_LEN', 'MH_LEN', 'MH_LEN_RATIO']

    fmt_features_names = [  'ARCMQ',
                            'MDLF', 'MDSP', 'MDSF', 'MDRF', 'MDSP_OVER_MDLF', 'MDSF_OVER_MDRF', 'MDLF_OVER_MDSP', 'MDRF_OVER_MDSF', 
                            'MDLFHQ', 'MDSPHQ', 'MDSFHQ', 'MDRFHQ', 'MDSP_OVER_MDLF_HQ', 'MDSF_OVER_MDRF_HQ', 'MDLF_OVER_MDSP_HQ', 'MDRF_OVER_MDSF_HQ',
                            'MDLC', 'MDRC', 'MDLCHQ', 'MDRCHQ', 'DPSL', 'DPSR', 'DPSLHQ', 'DPSRHQ', 'CP1', 'CP2', 'CP3']

    denovo_features_names = ['IMPRECISE', 'SPLIT_READS_RATIO', 'SPLIT_READS_RATIO1', 'SPLIT_READS_RATIO2', 'FWD_SPLIT_READS_RATIO1', 'FWD_SPLIT_READS_RATIO2', 'REV_SPLIT_READS_RATIO1', 'REV_SPLIT_READS_RATIO2',
                      'FWD_SPLIT_READS_RATIO', 'REV_SPLIT_READS_RATIO', 'RCC_EXT_1SR_READS1', 'RCC_EXT_1SR_READS2', 'LCC_EXT_1SR_READS1', 'LCC_EXT_1SR_READS2',
                      'EXT_1SR_READS1', 'EXT_1SR_READS2', 'RCC_HQ_EXT_1SR_READS1', 'RCC_HQ_EXT_1SR_READS2', 'LCC_HQ_EXT_1SR_READS1', 'LCC_HQ_EXT_1SR_READS2', 'HQ_EXT_1SR_READS1', 'HQ_EXT_1SR_READS2', 'FULL_TO_SPLIT_JUNCTION_SCORE_RATIO', 'FULL_TO_SPLIT_JUNCTION_SCORE_DIFF',
                      'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_RATIO2', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO1', 'SPLIT2_TO_SPLIT1_JUNCTION_SCORE_DIFF_RATIO2',
                      'SPLIT_TO_SIZE_RATIO1', 'SPLIT_TO_SIZE_RATIO2', 'SPLIT_JUNCTION_SIZE_RATIO1', 'SPLIT_JUNCTION_SIZE_RATIO2', 'MAX_SPLIT_JUNCTION_SIZE_RATIO', 'MIN_SPLIT_JUNCTION_SIZE_RATIO',
                      'B_DIFF']

    regt_shared_features_names = \
            ['AR1', 'ARC1', 'ARCF1', 'ARCR1', 'MAXARCD1', 'ARCAS1', 'ARC1HQ',
             'AR2', 'ARC2', 'ARCF2', 'ARCR2', 'MAXARCD2', 'ARCAS2', 'ARC2HQ',
             'RR1', 'RRC1', 'RR2', 'RRC2', 'ER', 
             'AR1_RATIO', 'AR2_RATIO', 'RR1_RATIO', 'RR2_RATIO', 'ARC1_RATIO', 'ARC2_RATIO', 'RRC1_RATIO', 'RRC2_RATIO',
             'AR1_OVER_RR1', 'RR1_OVER_AR1', 'AR2_OVER_RR2', 'RR2_OVER_AR2', 'ARC1_OVER_RRC1', 'RRC1_OVER_ARC1', 'ARC2_OVER_RRC2', 'RRC2_OVER_ARC2',
             'AXR', 'AXRHQ', 'EXL', 'EXAS', 'EXRS', 'EXAS_EXRS_RATIO', 'EXAS_EXRS_DIFF']

    stat_test_features_names = ['KS_PVAL', 'SIZE_NORM']

    dp_features_names = ['DP1', 'DP2', 'DP1HQ', 'DP2HQ', 'DP1_HQ_RATIO', 'DP2_HQ_RATIO', 'DP1MQ', 'DP2MQ', 'DPLANM', 'DPRANM', 'PTNR1', 'PTNR2']

    def get_denovo_feature_names(model_name):
        return Features.info_features_names + Features.fmt_features_names + \
            Features.denovo_features_names + Features.stat_test_features_names + Features.dp_features_names

    def get_regt_feature_names(model_name):
        features_names = Features.info_features_names + Features.fmt_features_names + Features.regt_shared_features_names + Features.dp_features_names
        if model_name in ["DEL", "DEL_IMPRECISE", "DUP", "DUP_IMPRECISE"]:
            features_names += Features.stat_test_features_names
        return features_names

    def get_feature_names(model_name, denovo):
        if denovo:
            return Features.get_denovo_feature_names(model_name)
        else:
            return Features.get_regt_feature_names(model_name)

    def get_denovo_model_name(record, max_is):
        svtype_str = Features.get_svtype(record)
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
        svtype_str = Features.get_svtype(record)
        if svtype_str == "DEL" and abs(Features.get_svlen(record)) >= max_is:
                svtype_str += "_LARGE"
        elif svtype_str == "DUP" and record.stop-record.start > read_len-30:
            svtype_str += "_LARGE"

        if Features.get_number_value(record.samples[0], 'EXL', 0) == 0:
            svtype_str += "_IMPRECISE"
        return svtype_str
    
    def get_model_name(record, max_is, read_len, denovo):
        if denovo:
            return Features.get_denovo_model_name(record, max_is)
        else:
            return Features.get_regt_model_name(record, max_is, read_len)

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
        
    def generate_id(record):
        svinsseq = Features.get_svinsseq(record)
        return f"{record.chrom}:{record.pos}-{record.stop}:{Features.get_svtype(record)}:{Features.get_svlen(record)}:{hash(svinsseq)}"

    def get_svinsseq(record):
        if "<" not in record.alts[0]:
            return record.alts[0]
        elif 'SVINSSEQ' in record.info:
            return record.info['SVINSSEQ']
        return ""

    def get_svlen(record):
        if 'SVLEN' not in record.info:
            svinsseq = Features.get_svinsseq(record)
            return len(svinsseq) - (record.stop - record.pos)
        svlen = record.info['SVLEN']
        if isinstance(svlen, list) or isinstance(svlen, tuple):
            return svlen[0]
        else:
            return svlen
        
    def get_svtype(record):
        if isinstance(record.info['SVTYPE'], list) or isinstance(record.info['SVTYPE'], tuple):
            return record.info['SVTYPE'][0]
        return record.info['SVTYPE']

    def normalise(value, min, max):
        if max == min:
            return value - min
        return (value - min) / (max - min)

    def record_to_features(record, stats, denovo):
        min_depth = get_stat(stats, 'min_depth', record.chrom)
        median_depth = get_stat(stats, 'median_depth', record.chrom)
        max_depth = get_stat(stats, 'max_depth', record.chrom)
        max_is = stats['max_is']['.']
        read_len = stats['read_len']['.']
        model_name = Features.get_model_name(record, max_is, read_len, denovo)

        features = dict()
        info = record.info
        svtype_str = Features.get_svtype(record)
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

        features['INS_SEQ_COV_PREFIX_LEN'] = 1
        features['INS_SEQ_COV_SUFFIX_LEN'] = 1
        if '-' in svinsseq:
            i = svinsseq.index('-')
            features['INS_SEQ_COV_PREFIX_LEN'] = i/len(svinsseq)
            features['INS_SEQ_COV_SUFFIX_LEN'] = (len(svinsseq)-i)/len(svinsseq)

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

        remap_lb = Features.get_number_value(info, 'REMAP_LB', record.pos)
        remap_ub = Features.get_number_value(info, 'REMAP_UB', record.stop)
        lb_diff, ub_diff = max(0, remap_lb-record.pos), max(0, record.stop-remap_ub)
        features['B_DIFF'] = lb_diff + ub_diff

        features['MH_LEN'] = Features.get_number_value(info, 'MH_LEN', 0)
        features['MH_LEN_RATIO'] = features['MH_LEN']/abs(max(1, svlen))

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

        features['IMPRECISE'] = True if "IMPRECISE" in record.info else False

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
        features['ARC1HQ'] = Features.get_number_value(record.samples[0], 'ARC1HQ', 0, max(1, arc1))

        features['AR2'] = Features.normalise(ar2, min_depth, max_depth)
        features['ARC2'] = Features.normalise(arc2, min_depth, max_depth)
        features['ARCF2'] = Features.get_number_value(record.samples[0], 'ARCF2', 0, max(1, arc2))
        features['ARCR2'] = Features.get_number_value(record.samples[0], 'ARCR2', 0, max(1, arc2))
        features['MAXARCD2'] = max(features['ARCF2'], features['ARCR2'])
        features['ARCAS2'] = Features.get_number_value(record.samples[0], 'ARCAS2', 0)
        features['ARC2HQ'] = Features.get_number_value(record.samples[0], 'ARC2HQ', 0, max(1, arc2))

        arc1 = Features.get_number_value(record.samples[0], 'ARC1MQ', 0)
        arc2 = Features.get_number_value(record.samples[0], 'ARC2MQ', 0)
        features['ARCMQ'] = max(arc1, arc2)

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

        md = Features.get_number_value(record.samples[0], 'MD', [0, 0, 0, 0])
        features['MDLF'] = Features.normalise(md[0], min_depth, max_depth)
        features['MDSP'] = Features.normalise(md[1], min_depth, max_depth)
        features['MDSF'] = Features.normalise(md[2], min_depth, max_depth)
        features['MDRF'] = Features.normalise(md[3], min_depth, max_depth)
        features['MDSP_OVER_MDLF'] = md[1]/max(1, md[0])
        features['MDSF_OVER_MDRF'] = md[2]/max(1, md[3])
        features['MDLF_OVER_MDSP'] = md[0]/max(1, md[1])
        features['MDRF_OVER_MDSF'] = md[3]/max(1, md[2])

        mdhq = Features.get_number_value(record.samples[0], 'MDHQ', [0, 0, 0, 0])
        features['MDLFHQ'] = Features.normalise(mdhq[0], min_depth, max_depth)
        features['MDSPHQ'] = Features.normalise(mdhq[1], min_depth, max_depth)
        features['MDSFHQ'] = Features.normalise(mdhq[2], min_depth, max_depth)
        features['MDRFHQ'] = Features.normalise(mdhq[3], min_depth, max_depth)
        features['MDSP_OVER_MDLF_HQ'] = mdhq[1]/max(1, mdhq[0])
        features['MDSF_OVER_MDRF_HQ'] = mdhq[2]/max(1, mdhq[3])
        features['MDLF_OVER_MDSP_HQ'] = mdhq[0]/max(1, mdhq[1])
        features['MDRF_OVER_MDSF_HQ'] = mdhq[3]/max(1, mdhq[2])

        clmd = Features.get_number_value(record.samples[0], 'CLMD', [0, 0])
        features['MDLC'] = Features.normalise(clmd[0], min_depth, max_depth)
        features['MDRC'] = Features.normalise(clmd[1], min_depth, max_depth)

        clmdhq = Features.get_number_value(record.samples[0], 'CLMDHQ', [0, 0])
        features['MDLCHQ'] = Features.normalise(clmdhq[0], min_depth, max_depth)
        features['MDRCHQ'] = Features.normalise(clmdhq[1], min_depth, max_depth)

        features['KS_PVAL'] = max(0, Features.get_number_value(record.samples[0], 'KSPVAL', 1.0))
        features['SIZE_NORM'] = 2
        if 'MAXSIZE' in record.samples[0]:
            min_size = float(record.samples[0]['MINSIZE'])
            max_size = float(record.samples[0]['MAXSIZE'])
            features['SIZE_NORM'] = Features.normalise(svlen/2, min_size, max_size)

        dp1, dp2 = Features.get_number_value(record.samples[0], 'DP', 0)
        dp1hq, dp2hq = Features.get_number_value(record.samples[0], 'DPHQ', 0)
        if svtype_str == "DEL":
            min_is_to_become_disc = int(max(0, max_is-svlen))
            min_disc_pairs = stats['min_pairs_crossing_gap'][str(min_is_to_become_disc)]
            max_disc_pairs = stats['max_pairs_crossing_gap'][str(min_is_to_become_disc)]
            dp1_scaled = Features.normalise(dp1, min_disc_pairs, max_disc_pairs)
            dp2_scaled = Features.normalise(dp2, min_disc_pairs, max_disc_pairs)
            dp1hq_scaled = Features.normalise(dp1hq, min_disc_pairs, max_disc_pairs)
            dp2hq_scaled = Features.normalise(dp2hq, min_disc_pairs, max_disc_pairs)
        elif svtype_str == "INS" and source_str in ("DE_NOVO_ASSEMBLY", "REFERENCE_GUIDED_ASSEMBLY"):
            min_inslen = int(min(max_is, svinslen))
            min_disc_pairs = stats['min_disc_pairs_by_insertion_size'][str(min_inslen)]
            max_disc_pairs = stats['max_disc_pairs_by_insertion_size'][str(min_inslen)]
            dp1_scaled = Features.normalise(dp1, min_disc_pairs, max_disc_pairs)
            dp2_scaled = Features.normalise(dp2, min_disc_pairs, max_disc_pairs)
            dp1hq_scaled = Features.normalise(dp1hq, min_disc_pairs, max_disc_pairs)
            dp2hq_scaled = Features.normalise(dp2hq, min_disc_pairs, max_disc_pairs)
        else:
            dp1_scaled = dp1/median_depth
            dp2_scaled = dp2/median_depth
            dp1hq_scaled = dp1hq/median_depth
            dp2hq_scaled = dp2hq/median_depth

        features['DP1'] = dp1_scaled
        features['DP2'] = dp2_scaled
        features['DP1HQ'] = dp1hq_scaled
        features['DP2HQ'] = dp2hq_scaled
        features['DP1_HQ_RATIO'] = dp1hq/max(1, dp1)
        features['DP2_HQ_RATIO'] = dp2hq/max(1, dp2)
        features['DP1MQ'], features['DP2MQ'] = Features.get_number_value(record.samples[0], 'DPMQ', 0)
        features['DPLANM'], features['DPRANM'] = Features.get_number_value(record.samples[0], 'DPNM', 0, read_len)
        features['DPSL'], features['DPSR'] = Features.get_number_value(record.samples[0], 'DPS', [0, 0], median_depth)
        features['DPSLHQ'], features['DPSRHQ'] = Features.get_number_value(record.samples[0], 'DPSHQ', [0, 0], median_depth)

        min_pairs_crossing_point, max_pairs_crossing_point = stats['min_pairs_crossing_gap']["0"], stats['max_pairs_crossing_gap']["0"]
        cp = Features.get_number_value(record.samples[0], 'CP', [0, 0, 0])
        features['CP1'], features['CP2'], features['CP3'] = [Features.normalise(c, min_pairs_crossing_point, max_pairs_crossing_point) for c in cp]
        features['PTNR1'], features['PTNR2'] = dp1/max(1, dp1+cp[0]), dp2/max(1, dp2+cp[2])

        features['AXR'] = Features.get_number_value(record.samples[0], 'AXR', 0, median_depth*max_is)
        features['AXRHQ'] = Features.get_number_value(record.samples[0], 'AXRHQ', 0, median_depth*max_is)
        exl = Features.get_number_value(record.samples[0], 'EXL', 0)
        features['EXL'] = exl/max_is
        exas = Features.get_number_value(record.samples[0], 'EXAS', 0, max(1, exl))
        features['EXAS'] = exas
        exrs = Features.get_number_value(record.samples[0], 'EXRS', 0, max(1, exl))
        features['EXRS'] = exrs
        features['EXAS_EXRS_RATIO'] = exas/max(0.01, exrs)
        features['EXAS_EXRS_DIFF'] = exas-exrs

        feature_values = []
        for feature_name in Features.get_feature_names(model_name, denovo):
            feature_values.append(features[feature_name])
        return feature_values

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
def parse_vcf(vcf_fname, stats_fname, fp_fname, denovo, tolerate_no_gts = False):
    gts = read_gts(fp_fname, tolerate_no_gts=tolerate_no_gts)
    vcf_reader = pysam.VariantFile(vcf_fname)
    stats_reader = open(stats_fname, 'r')
    features_by_source, gts_by_source, variant_ids_by_source = defaultdict(list), defaultdict(list), defaultdict(list)

    # read the stats file and extract the relevant values
    stats = defaultdict(dict)
    for line in stats_reader:
        sl = line.strip().split()
        stats[sl[0]][sl[1]] = int(sl[2])

    for record in vcf_reader.fetch():
        record_svtype = Features.get_svtype(record)
        if record_svtype.startswith('INV'):
            continue

        model_name = Features.get_model_name(record, stats['max_is']['.'], stats['read_len']['.'], denovo)
        feature_values = Features.record_to_features(record, stats, denovo)
        if denovo or ('TD' not in record.samples[0] and gts[record.id] != "./."): # if too deep or no genotype is available, skip the record
            features_by_source[model_name].append(feature_values)
            gts_by_source[model_name].append(gts[record.id])
            variant_ids_by_source[model_name].append(Features.generate_id(record))

    for model_name in features_by_source:
        features_by_source[model_name] = np.array(features_by_source[model_name])
        gts_by_source[model_name] = np.array(gts_by_source[model_name])
        variant_ids_by_source[model_name] = np.array(variant_ids_by_source[model_name])
    
    return features_by_source, gts_by_source, variant_ids_by_source
