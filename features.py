from __future__ import division
import os, pysam
from collections import defaultdict
import numpy as np
import math

class Features:

    NAN = np.nan

    info_features_names = [ 'START_STOP_DIST', 'SVLEN', 'SVINSLEN', 
                            'SV_REF_PREFIX_A_RATIO', 'SV_REF_PREFIX_C_RATIO', 'SV_REF_PREFIX_G_RATIO', 'SV_REF_PREFIX_T_RATIO', 'MAX_SV_REF_PREFIX_BASE_RATIO',
                            'SV_REF_SUFFIX_A_RATIO', 'SV_REF_SUFFIX_C_RATIO', 'SV_REF_SUFFIX_G_RATIO', 'SV_REF_SUFFIX_T_RATIO', 'MAX_SV_REF_SUFFIX_BASE_RATIO',
                            'LEFT_ANCHOR_A_RATIO', 'LEFT_ANCHOR_C_RATIO', 'LEFT_ANCHOR_G_RATIO', 'LEFT_ANCHOR_T_RATIO', 'MAX_LEFT_ANCHOR_BASE_RATIO',
                            'RIGHT_ANCHOR_A_RATIO', 'RIGHT_ANCHOR_C_RATIO', 'RIGHT_ANCHOR_G_RATIO', 'RIGHT_ANCHOR_T_RATIO', 'MAX_RIGHT_ANCHOR_BASE_RATIO',
                            'INS_PREFIX_A_RATIO', 'INS_PREFIX_C_RATIO', 'INS_PREFIX_G_RATIO', 'INS_PREFIX_T_RATIO', 'MAX_INS_PREFIX_BASE_COUNT_RATIO',
                            'INS_SUFFIX_A_RATIO', 'INS_SUFFIX_C_RATIO', 'INS_SUFFIX_G_RATIO', 'INS_SUFFIX_T_RATIO', 'MAX_INS_SUFFIX_BASE_COUNT_RATIO',
                            'INS_SEQ_COV_PREFIX_LEN', 'INS_SEQ_COV_SUFFIX_LEN', 'MH_LEN', 'MH_LEN_RATIO']
    
    reads_features_names = ['AR1', 'AR1C', 'AR1CmQ', 'AR1CMQ', 'AR1CHQ', 'AR1C_HQ_RATIO',
                            'AR2', 'AR2C', 'AR2CmQ', 'AR2CMQ', 'AR2CHQ', 'AR2C_HQ_RATIO', 'MAXARCD',
                            'RR1', 'RR1C', 'RR1CmQ', 'RR1CMQ', 'RR1CHQ',
                            'RR2', 'RR2C', 'RR2CmQ', 'RR2CMQ', 'RR2CHQ', 'MAXRRCD',
                            'AR1_RR1_CAQ_Z_SCORE', 'AR2_RR2_CAQ_Z_SCORE', 'AR1_RR1_CAS_Z_SCORE', 'AR2_RR2_CAS_Z_SCORE', 
                            'AR1_OVER_RR1', 'AR2_OVER_RR2', 'AR1C_OVER_RR1C', 'AR2C_OVER_RR2C']

    fmt_features_names = [  'AXR1', 'AXR2', 'AXR1HQ', 'AXR2HQ',
                            'EXSS1_1', 'EXSS1_2', 'EXSS2_1', 'EXSS2_2',
                            'EXSS1_RATIO1', 'EXSS1_RATIO2', 'EXSS2_RATIO1', 'EXSS2_RATIO2',
                            'EXAS_EXRS_RATIO', 'EXAS_EXRS_DIFF',
                            'EXSCC1_1_IA_RATIO', 'EXSCC1_2_IA_RATIO', 'EXSCC2_1_IA_RATIO', 'EXSCC2_2_IA_RATIO',
                            'EXSSC1_1_IA_DIFF_RATIO', 'EXSSC1_2_IA_DIFF_RATIO', 'EXSSC2_1_IA_DIFF_RATIO', 'EXSSC2_2_IA_DIFF_RATIO',
                            'MEXL', 'mEXL', 'EXL',
                            'MDLF', 'MDSP', 'MDSF', 'MDRF', 
                            'MDSP_OVER_MDLF', 'MDSF_OVER_MDRF', 'MDLF_OVER_MDSP', 'MDRF_OVER_MDSF', 
                            'MDLFHQ', 'MDSPHQ', 'MDSFHQ', 'MDRFHQ', 
                            'MDSP_OVER_MDLF_HQ', 'MDSF_OVER_MDRF_HQ', 'MDLF_OVER_MDSP_HQ', 'MDRF_OVER_MDSF_HQ',
                            'MDLC', 'MDRC', 'MDLCHQ', 'MDRCHQ',
                            'TD']

    stat_test_features_names = ['KS_PVAL', 'SIZE_NORM']

    dp_features_names = ['ASP1', 'ASP1HQ_1', 'ASP1HQ_2', 'ASP1HQ_1_RATIO', 'ASP1HQ_2_RATIO',
                         'ASP2', 'ASP2HQ_1', 'ASP2HQ_2', 'ASP2HQ_1_RATIO', 'ASP2HQ_2_RATIO',
                         'ASP1_ASP2_RATIO',
                         'ASP1mQ_1', 'ASP1mQ_2', 'ASP1MQ_1', 'ASP1MQ_2', 'ASP1SPAN_1', 'ASP1SPAN_2',
                         'ASP2mQ_1', 'ASP2mQ_2', 'ASP2MQ_1', 'ASP2MQ_2', 'ASP2SPAN_1', 'ASP2SPAN_2',
                         'ASP1_OVER_RSP1', 'ASP2_OVER_RSP2',
                         'ASP1_RSP1_1_NM_Z_SCORE', 'ASP1_RSP1_2_NM_Z_SCORE', 'ASP2_RSP2_1_NM_Z_SCORE', 'ASP2_RSP2_2_NM_Z_SCORE',
                         'ASP1_RSP1_1_AQ_Z_SCORE', 'ASP1_RSP1_2_AQ_Z_SCORE', 'ASP2_RSP2_1_AQ_Z_SCORE', 'ASP2_RSP2_2_AQ_Z_SCORE',
                         'RSP1', 'RSP1HQ_1', 'RSP1HQ_2', 'RSP1HQ_1_RATIO', 'RSP1HQ_2_RATIO',
                         'RSP2', 'RSP2HQ_1', 'RSP2HQ_2', 'RSP2HQ_1_RATIO', 'RSP2HQ_2_RATIO',
                         'RSP1mQ_1', 'RSP1mQ_2', 'RSP1MQ_1', 'RSP1MQ_2', 
                         'RSP2mQ_1', 'RSP2mQ_2', 'RSP2MQ_1', 'RSP2MQ_2',
                         'SSP1', 'SSP1HQ_1', 'SSP1HQ_2', #'SSP1HQ_1_RATIO', 'SSP1HQ_2_RATIO',
                         'SSP2', 'SSP2HQ_1', 'SSP2HQ_2', #'SSP2HQ_1_RATIO', 'SSP2HQ_2_RATIO'
                         'SSP1mQ_1', 'SSP1mQ_2', 'SSP1MQ_1', 'SSP1MQ_2',
                         'SSP2mQ_1', 'SSP2mQ_2', 'SSP2MQ_1', 'SSP2MQ_2',
                         'SSP1_RSP1_1_NM_Z_SCORE', 'SSP1_RSP1_2_NM_Z_SCORE', 'SSP2_RSP2_1_NM_Z_SCORE', 'SSP2_RSP2_2_NM_Z_SCORE',
                         'SSP1_RSP1_1_AQ_Z_SCORE', 'SSP1_RSP1_2_AQ_Z_SCORE', 'SSP2_RSP2_1_AQ_Z_SCORE', 'SSP2_RSP2_2_AQ_Z_SCORE'
    ]

    def get_feature_names(model_name):
        return Features.info_features_names + Features.fmt_features_names + Features.reads_features_names + \
            Features.stat_test_features_names + Features.dp_features_names

    def get_model_name(record, max_is, read_len):
        svtype_str = Features.get_svtype(record)
        if svtype_str == "DEL" and abs(Features.get_svlen(record)) >= max_is:
            svtype_str += "_LARGE"
        elif svtype_str == "DUP" and record.stop-record.start > read_len-30:
            svtype_str += "_LARGE"

        if Features.get_number_value(record.samples[0], 'EXL', 0) == 0 or \
           Features.get_number_value(record.samples[0], 'EXL2', 1) == 0: # note that EXL2 should be PRESENT and 0
            svtype_str += "_IMPRECISE"
        return svtype_str

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
    
    def piecewise_normalise(value, minv, maxv):
        if value < minv:
            return value/minv * 0.25
        elif value >= maxv:
            return value/max(1, maxv) * 0.25 + 0.75
        else:
            return 0.25 + (value - minv) / (maxv - minv) * 0.75

    def calculate_z_score(mean1, stddev1, n1, mean2, stddev2, n2):
        if np.isnan(mean1) or np.isnan(mean2) or n1 == 0 or n2 == 0:
            return Features.NAN
        std_error = math.sqrt((stddev1**2 / n1) + (stddev2**2 / n2))
        if std_error == 0:
            std_error = 1
        z_score = (mean1 - mean2) / std_error
        return z_score

    def record_to_features(record, stats):
        min_depth = get_stat(stats, 'min_depth', record.chrom)
        median_depth = get_stat(stats, 'median_depth', record.chrom)
        max_depth = get_stat(stats, 'max_depth', record.chrom)
        max_is = stats['max_is']['.']
        read_len = stats['read_len']['.']
        min_pairs_crossing_point = stats['min_pairs_crossing_gap']["0"]
        max_pairs_crossing_point = stats['max_pairs_crossing_gap']["0"]
        model_name = Features.get_model_name(record, max_is, read_len)

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

        td = Features.get_number_value(record.samples[0], 'TD', 0)
        features['TD'] = td

        ar1 = Features.get_number_value(record.samples[0], 'AR1', 0)
        ar1c = Features.get_number_value(record.samples[0], 'AR1C', 0)
        ar1caq = Features.get_number_value(record.samples[0], 'AR1CAQ', Features.NAN)
        ar1csq = Features.get_number_value(record.samples[0], 'AR1CSQ', Features.NAN)
        ar1cas = Features.get_number_value(record.samples[0], 'AR1CAS', Features.NAN)
        ar1css = Features.get_number_value(record.samples[0], 'AR1CSS', Features.NAN)
        features['AR1'] = Features.piecewise_normalise(ar1, min_depth, max_depth)
        features['AR1C'] = Features.piecewise_normalise(ar1c, min_depth, max_depth)
        features['AR1CmQ'] = Features.get_number_value(record.samples[0], 'AR1CmQ', Features.NAN)
        features['AR1CMQ'] = Features.get_number_value(record.samples[0], 'AR1CMQ', Features.NAN)
        arc1hq = Features.get_number_value(record.samples[0], 'AR1CHQ', 0)
        features['AR1CHQ'] = Features.piecewise_normalise(arc1hq, min_depth, max_depth)
        features['AR1C_HQ_RATIO'] = arc1hq/max(1, ar1c)

        ar2 = Features.get_number_value(record.samples[0], 'AR2', 0)
        ar2c = Features.get_number_value(record.samples[0], 'AR2C', 0)
        ar2caq = Features.get_number_value(record.samples[0], 'AR2CAQ', Features.NAN)
        ar2csq = Features.get_number_value(record.samples[0], 'AR2CSQ', Features.NAN)
        ar2cas = Features.get_number_value(record.samples[0], 'AR2CAS', Features.NAN)
        ar2css = Features.get_number_value(record.samples[0], 'AR2CSS', Features.NAN)
        features['AR2'] = Features.piecewise_normalise(ar2, min_depth, max_depth)
        features['AR2C'] = Features.piecewise_normalise(ar2c, min_depth, max_depth)
        features['AR2CmQ'] = Features.get_number_value(record.samples[0], 'AR2CmQ', Features.NAN)
        features['AR2CMQ'] = Features.get_number_value(record.samples[0], 'AR2CMQ', Features.NAN)
        arc2hq = Features.get_number_value(record.samples[0], 'AR2CHQ', 0)
        features['AR2CHQ'] = Features.piecewise_normalise(arc2hq, min_depth, max_depth)
        features['AR2C_HQ_RATIO'] = arc2hq/max(1, ar2c)

        ar1cf = Features.get_number_value(record.samples[0], 'AR1CF', 0, max(1, ar1c))
        ar1cr = Features.get_number_value(record.samples[0], 'AR1CR', 0, max(1, ar1c))
        ar2cf = Features.get_number_value(record.samples[0], 'AR2CF', 0, max(1, ar2c))
        ar2cr = Features.get_number_value(record.samples[0], 'AR2CR', 0, max(1, ar2c))
        features['ARCF'] = ar1cf + ar2cf
        features['ARCR'] = ar1cr + ar2cr
        features['MAXARCD'] = max(features['ARCF'], features['ARCR'])

        rr1 = Features.get_number_value(record.samples[0], 'RR1', 0)
        rr1c = Features.get_number_value(record.samples[0], 'RR1C', 0)
        rr1caq = Features.get_number_value(record.samples[0], 'RR1CAQ', Features.NAN)
        rr1csq = Features.get_number_value(record.samples[0], 'RR1CSQ', Features.NAN)
        rr1cas = Features.get_number_value(record.samples[0], 'RR1CAS', Features.NAN)
        rr1css = Features.get_number_value(record.samples[0], 'RR1CSS', Features.NAN)
        features['RR1'] = Features.piecewise_normalise(rr1, min_depth, max_depth)
        features['RR1C'] = Features.piecewise_normalise(rr1c, min_depth, max_depth)
        features['RR1CmQ'] = Features.get_number_value(record.samples[0], 'RR1CmQ', Features.NAN)
        features['RR1CMQ'] = Features.get_number_value(record.samples[0], 'RR1CMQ', Features.NAN)
        features['RR1CHQ'] = Features.get_number_value(record.samples[0], 'RR1CHQ', 0, max(1, rr1c))

        rr2 = Features.get_number_value(record.samples[0], 'RR2', 0)
        rr2c = Features.get_number_value(record.samples[0], 'RR2C', 0)
        rr2caq = Features.get_number_value(record.samples[0], 'RR2CAQ', Features.NAN)
        rr2csq = Features.get_number_value(record.samples[0], 'RR2CSQ', Features.NAN)
        rr2cas = Features.get_number_value(record.samples[0], 'RR2CAS', Features.NAN)
        rr2css = Features.get_number_value(record.samples[0], 'RR2CSS', Features.NAN)
        features['RR2'] = Features.piecewise_normalise(rr2, min_depth, max_depth)
        features['RR2C'] = Features.piecewise_normalise(rr2c, min_depth, max_depth)
        features['RR2CmQ'] = Features.get_number_value(record.samples[0], 'RR2CmQ', Features.NAN)
        features['RR2CMQ'] = Features.get_number_value(record.samples[0], 'RR2CMQ', Features.NAN)
        features['RR2CHQ'] = Features.get_number_value(record.samples[0], 'RR2CHQ', 0, max(1, rr2c))

        rr1cf = Features.get_number_value(record.samples[0], 'RR1CF', 0, max(1, rr1c))
        rr1cr = Features.get_number_value(record.samples[0], 'RR1CR', 0, max(1, rr1c))
        rr2cf = Features.get_number_value(record.samples[0], 'RR2CF', 0, max(1, rr2c))
        rr2cr = Features.get_number_value(record.samples[0], 'RR2CR', 0, max(1, rr2c))
        features['RRCF'] = rr1cf + rr2cf
        features['RRCR'] = rr1cr + rr2cr
        features['MAXRRCD'] = max(features['RRCF'], features['RRCR'])
        
        features['AR1_RR1_CAQ_Z_SCORE'] = Features.calculate_z_score(ar1caq, ar1csq, ar1c, rr1caq, rr1csq, rr1c)
        features['AR1_RR1_CAS_Z_SCORE'] = Features.calculate_z_score(ar1cas, ar1css, ar1c, rr1cas, rr1css, rr1c)

        if 'AR2' not in record.samples[0]:
            ar2 = ar1
            ar2c = ar1c
            ar2caq = ar1caq
            ar2csq = ar1csq
            ar2cas = ar1cas
            ar2css = ar1css
        if 'RR2' not in record.samples[0]:
            rr2 = rr1
            rr2c = rr1c
            rr2caq = rr1caq
            rr2csq = rr1csq
            rr2cas = rr1cas
            rr2css = rr1css

        features['AR2_RR2_CAQ_Z_SCORE'] = Features.calculate_z_score(ar2caq, ar2csq, ar2c, rr2caq, rr2csq, rr2c)
        features['AR2_RR2_CAS_Z_SCORE'] = Features.calculate_z_score(ar2cas, ar2css, ar2c, rr2cas, rr2css, rr2c)

        features['AR1_OVER_RR1'] = ar1/max(1, ar1+rr1)
        features['AR1C_OVER_RR1C'] = ar1c/max(1, ar1c+rr1c)
        features['AR2_OVER_RR2'] = ar2/max(1, ar2+rr2)
        features['AR2C_OVER_RR2C'] = ar2c/max(1, ar2c+rr2c)

        md = Features.get_number_value(record.samples[0], 'MD', [0, 0, 0, 0])
        features['MDLF'] = Features.piecewise_normalise(md[0], min_depth, max_depth)
        features['MDSP'] = Features.piecewise_normalise(md[1], min_depth, max_depth)
        features['MDSF'] = Features.piecewise_normalise(md[2], min_depth, max_depth)
        features['MDRF'] = Features.piecewise_normalise(md[3], min_depth, max_depth)
        features['MDSP_OVER_MDLF'] = md[1]/max(1, md[0])
        features['MDSF_OVER_MDRF'] = md[2]/max(1, md[3])
        features['MDLF_OVER_MDSP'] = md[0]/max(1, md[1])
        features['MDRF_OVER_MDSF'] = md[3]/max(1, md[2])

        mdhq = Features.get_number_value(record.samples[0], 'MDHQ', [0, 0, 0, 0])
        features['MDLFHQ'] = Features.piecewise_normalise(mdhq[0], min_depth, max_depth)
        features['MDSPHQ'] = Features.piecewise_normalise(mdhq[1], min_depth, max_depth)
        features['MDSFHQ'] = Features.piecewise_normalise(mdhq[2], min_depth, max_depth)
        features['MDRFHQ'] = Features.piecewise_normalise(mdhq[3], min_depth, max_depth)
        features['MDSP_OVER_MDLF_HQ'] = mdhq[1]/max(1, mdhq[0])
        features['MDSF_OVER_MDRF_HQ'] = mdhq[2]/max(1, mdhq[3])
        features['MDLF_OVER_MDSP_HQ'] = mdhq[0]/max(1, mdhq[1])
        features['MDRF_OVER_MDSF_HQ'] = mdhq[3]/max(1, mdhq[2])

        clmd = Features.get_number_value(record.samples[0], 'CLMD', [0, 0])
        features['MDLC'] = Features.piecewise_normalise(clmd[0], min_depth, max_depth)
        features['MDRC'] = Features.piecewise_normalise(clmd[1], min_depth, max_depth)

        clmdhq = Features.get_number_value(record.samples[0], 'CLMDHQ', [0, 0])
        features['MDLCHQ'] = Features.piecewise_normalise(clmdhq[0], min_depth, max_depth)
        features['MDRCHQ'] = Features.piecewise_normalise(clmdhq[1], min_depth, max_depth)

        features['KS_PVAL'] = Features.get_number_value(record.samples[0], 'KSPVAL', Features.NAN)
        features['SIZE_NORM'] = Features.NAN
        if 'MAXSIZE' in record.samples[0]:
            min_size = float(record.samples[0]['MINSIZE'])
            max_size = float(record.samples[0]['MAXSIZE'])
            features['SIZE_NORM'] = Features.normalise(svlen/2, min_size, max_size)

        if svtype_str == "DEL":
            min_is_to_become_disc = int(max(0, max_is-svlen))
            min_disc_pairs = stats['min_pairs_crossing_gap'][str(min_is_to_become_disc)]
            max_disc_pairs = stats['max_pairs_crossing_gap'][str(min_is_to_become_disc)]
        elif svtype_str == "INS" and source_str in ("DE_NOVO_ASSEMBLY", "REFERENCE_GUIDED_ASSEMBLY"):
            min_inslen = int(min(max_is, svinslen))
            min_disc_pairs = stats['min_disc_pairs_by_insertion_size'][str(min_inslen)]
            max_disc_pairs = stats['max_disc_pairs_by_insertion_size'][str(min_inslen)]
        else:
            min_disc_pairs = min_pairs_crossing_point
            max_disc_pairs = max_pairs_crossing_point

        asp1 = Features.get_number_value(record.samples[0], 'ASP1', 0)
        asp1hq_1, asp1hq_2 = Features.get_number_value(record.samples[0], 'ASP1HQ', [0, 0])
        asp1nma_1, asp1nma_2 = Features.get_number_value(record.samples[0], 'ASP1NMA', [Features.NAN, Features.NAN])
        asp1nms_1, asp1nms_2 = Features.get_number_value(record.samples[0], 'ASP1NMS', [Features.NAN, Features.NAN])
        asp1aq_1, asp1aq_2 = Features.get_number_value(record.samples[0], 'ASP1AQ', [Features.NAN, Features.NAN])
        asp1sq_1, asp1sq_2 = Features.get_number_value(record.samples[0], 'ASP1SQ', [Features.NAN, Features.NAN])
        features['ASP1'] = Features.piecewise_normalise(asp1, min_disc_pairs, max_disc_pairs)
        features['ASP1HQ_1'] = Features.piecewise_normalise(asp1hq_1, min_disc_pairs, max_disc_pairs)
        features['ASP1HQ_2'] = Features.piecewise_normalise(asp1hq_2, min_disc_pairs, max_disc_pairs)
        features['ASP1HQ_1_RATIO'], features['ASP1HQ_2_RATIO'] = asp1hq_1/max(1, asp1), asp1hq_2/max(1, asp1)
        features['ASP1mQ_1'], features['ASP1mQ_2'] = Features.get_number_value(record.samples[0], 'ASP1mQ', [Features.NAN, Features.NAN])
        features['ASP1MQ_1'], features['ASP1MQ_2'] = Features.get_number_value(record.samples[0], 'ASP1MQ', [Features.NAN, Features.NAN])
        
        asp2 = Features.get_number_value(record.samples[0], 'ASP2', 0)
        asp2hq_1, asp2hq_2 = Features.get_number_value(record.samples[0], 'ASP2HQ', [0, 0])
        asp2nma_1, asp2nma_2 = Features.get_number_value(record.samples[0], 'ASP2NMA', [Features.NAN, Features.NAN])
        asp2nms_1, asp2nms_2 = Features.get_number_value(record.samples[0], 'ASP2NMS', [Features.NAN, Features.NAN])
        asp2aq_1, asp2aq_2 = Features.get_number_value(record.samples[0], 'ASP2AQ', [Features.NAN, Features.NAN])
        asp2sq_1, asp2sq_2 = Features.get_number_value(record.samples[0], 'ASP2SQ', [Features.NAN, Features.NAN])
        features['ASP2'] = Features.piecewise_normalise(asp2, min_disc_pairs, max_disc_pairs)
        features['ASP2HQ_1'] = Features.piecewise_normalise(asp2hq_1, min_disc_pairs, max_disc_pairs)
        features['ASP2HQ_2'] = Features.piecewise_normalise(asp2hq_2, min_disc_pairs, max_disc_pairs)
        features['ASP2HQ_1_RATIO'], features['ASP2HQ_2_RATIO'] = asp2hq_1/max(1, asp2), asp2hq_2/max(1, asp2)
        features['ASP2mQ_1'], features['ASP2mQ_2'] = Features.get_number_value(record.samples[0], 'ASP2mQ', [Features.NAN, Features.NAN])
        features['ASP2MQ_1'], features['ASP2MQ_2'] = Features.get_number_value(record.samples[0], 'ASP2MQ', [Features.NAN, Features.NAN])
        
        features['ASP1_ASP2_RATIO'] = max(asp1, asp2)/max(1, asp1+asp2)

        asp1span_1, asp1span_2 = Features.get_number_value(record.samples[0], 'ASP1SPAN', [0, 0])
        asp2span_1, asp2span_2 = Features.get_number_value(record.samples[0], 'ASP2SPAN', [0, 0])
        features['ASP1SPAN_1'], features['ASP2SPAN_2'] = asp1span_1/max_is, asp2span_2/max_is
        if svtype_str == "INS":
            features['ASP1SPAN_2'] = asp1span_2/max(1, max_is, svinslen)
            features['ASP2SPAN_1'] = asp2span_1/max(1, max_is, svinslen)
        else:
            features['ASP1SPAN_2'] = asp1span_2/max_is
            features['ASP2SPAN_1'] = asp2span_1/max_is

        rsp1 = Features.get_number_value(record.samples[0], 'RSP1', 0)
        rsp1hq_1, rsp1hq_2 = Features.get_number_value(record.samples[0], 'RSP1HQ', [0, 0])
        rsp1nma_1, rsp1nma_2 = Features.get_number_value(record.samples[0], 'RSP1NMA', [Features.NAN, Features.NAN])
        rsp1nms_1, rsp1nms_2 = Features.get_number_value(record.samples[0], 'RSP1NMS', [Features.NAN, Features.NAN])
        rsp1aq_1, rsp1aq_2 = Features.get_number_value(record.samples[0], 'RSP1AQ', [Features.NAN, Features.NAN])
        rsp1sq_1, rsp1sq_2 = Features.get_number_value(record.samples[0], 'RSP1SQ', [Features.NAN, Features.NAN])
        features['RSP1'] = Features.piecewise_normalise(rsp1, min_disc_pairs, max_disc_pairs)
        features['RSP1HQ_1']= Features.piecewise_normalise(rsp1hq_1, min_disc_pairs, max_disc_pairs)
        features['RSP1HQ_2'] = Features.piecewise_normalise(rsp1hq_2, min_disc_pairs, max_disc_pairs)
        features['RSP1HQ_1_RATIO'], features['RSP1HQ_2_RATIO'] = rsp1hq_1/max(1, rsp1), rsp1hq_2/max(1, rsp1)
        features['RSP1mQ_1'], features['RSP1mQ_2'] = Features.get_number_value(record.samples[0], 'RSP1mQ', [Features.NAN, Features.NAN])
        features['RSP1MQ_1'], features['RSP1MQ_2'] = Features.get_number_value(record.samples[0], 'RSP1MQ', [Features.NAN, Features.NAN])

        rsp2 = Features.get_number_value(record.samples[0], 'RSP2', 0)
        rsp2hq_1, rsp2hq_2 = Features.get_number_value(record.samples[0], 'RSP2HQ', [0, 0])
        rsp2nma_1, rsp2nma_2 = Features.get_number_value(record.samples[0], 'RSP2NMA', [Features.NAN, Features.NAN])
        rsp2nms_1, rsp2nms_2 = Features.get_number_value(record.samples[0], 'RSP2NMS', [Features.NAN, Features.NAN])
        rsp2aq_1, rsp2aq_2 = Features.get_number_value(record.samples[0], 'RSP2AQ', [Features.NAN, Features.NAN])
        rsp2sq_1, rsp2sq_2 = Features.get_number_value(record.samples[0], 'RSP2SQ', [Features.NAN, Features.NAN])
        features['RSP2'] = Features.piecewise_normalise(rsp2, min_disc_pairs, max_disc_pairs)
        features['RSP2HQ_1'] = Features.piecewise_normalise(rsp2hq_1, min_disc_pairs, max_disc_pairs)
        features['RSP2HQ_2'] = Features.piecewise_normalise(rsp2hq_2, min_disc_pairs, max_disc_pairs)
        features['RSP2HQ_1_RATIO'], features['RSP2HQ_2_RATIO'] = rsp2hq_1/max(1, rsp2), rsp2hq_2/max(1, rsp2)
        features['RSP2mQ_1'], features['RSP2mQ_2'] = Features.get_number_value(record.samples[0], 'RSP2mQ', [Features.NAN, Features.NAN])
        features['RSP2MQ_1'], features['RSP2MQ_2'] = Features.get_number_value(record.samples[0], 'RSP2MQ', [Features.NAN, Features.NAN])

        if 'ASP2' not in record.samples[0]:
            asp2 = asp1
            asp2nma_1 = asp1nma_1
            asp2nms_1 = asp1nms_1
            asp2nma_2 = asp1nma_2
            asp2nms_2 = asp1nms_2
            asp2aq_1 = asp1aq_1
            asp2aq_2 = asp1aq_2
            asp2sq_1 = asp1sq_1
            asp2sq_2 = asp1sq_2
        if 'RSP2' not in record.samples[0]:
            rsp2 = rsp1
            rsp2nma_1 = rsp1nma_1
            rsp2nms_1 = rsp1nms_1
            rsp2nma_2 = rsp1nma_2
            rsp2nms_2 = rsp1nms_2
            rsp2aq_1 = rsp1aq_1
            rsp2aq_2 = rsp1aq_2
            rsp2sq_1 = rsp1sq_1
            rsp2sq_2 = rsp1sq_2

        features['ASP1_OVER_RSP1'], features['ASP2_OVER_RSP2'] = asp1/max(1, asp1+rsp1), asp2/max(1, asp2+rsp2)

        features['ASP1_RSP1_1_NM_Z_SCORE'] = Features.calculate_z_score(asp1nma_1, asp1nms_1, asp1, rsp1nma_1, rsp1nms_1, rsp1)
        features['ASP1_RSP1_2_NM_Z_SCORE'] = Features.calculate_z_score(asp1nma_2, asp1nms_2, asp1, rsp1nma_2, rsp1nms_2, rsp1)
        features['ASP2_RSP2_1_NM_Z_SCORE'] = Features.calculate_z_score(asp2nma_1, asp2nms_1, asp2, rsp2nma_1, rsp2nms_1, rsp2)
        features['ASP2_RSP2_2_NM_Z_SCORE'] = Features.calculate_z_score(asp2nma_2, asp2nms_2, asp2, rsp2nma_2, rsp2nms_2, rsp2)

        features['ASP1_RSP1_1_AQ_Z_SCORE'] = Features.calculate_z_score(asp1aq_1, asp1sq_1, asp1, rsp1aq_1, rsp1sq_1, rsp1)
        features['ASP1_RSP1_2_AQ_Z_SCORE'] = Features.calculate_z_score(asp1aq_2, asp1sq_2, asp1, rsp1aq_2, rsp1sq_2, rsp1)
        features['ASP2_RSP2_1_AQ_Z_SCORE'] = Features.calculate_z_score(asp2aq_1, asp2sq_1, asp2, rsp2aq_1, rsp2sq_1, rsp2)
        features['ASP2_RSP2_2_AQ_Z_SCORE'] = Features.calculate_z_score(asp2aq_2, asp2sq_2, asp2, rsp2aq_2, rsp2sq_2, rsp2)

        ssp1 = Features.get_number_value(record.samples[0], 'SSP1', 0)
        ssp1hq_1, ssp1hq_2 = Features.get_number_value(record.samples[0], 'SSP1HQ', [0, 0])
        ssp1nma_1, ssp1nma_2 = Features.get_number_value(record.samples[0], 'SSP1NMA', [Features.NAN, Features.NAN])
        ssp1nms_1, ssp1nms_2 = Features.get_number_value(record.samples[0], 'SSP1NMS', [Features.NAN, Features.NAN])
        ssp1aq_1, ssp1aq_2 = Features.get_number_value(record.samples[0], 'SSP1AQ', [Features.NAN, Features.NAN])
        ssp1sq_1, ssp1sq_2 = Features.get_number_value(record.samples[0], 'SSP1SQ', [Features.NAN, Features.NAN])
        features['SSP1'] = Features.piecewise_normalise(ssp1, min_disc_pairs, max_disc_pairs)
        features['SSP1HQ_1'] = Features.piecewise_normalise(ssp1hq_1, min_disc_pairs, max_disc_pairs)
        features['SSP1HQ_2'] = Features.piecewise_normalise(ssp1hq_2, min_disc_pairs, max_disc_pairs)
        features['SSP1HQ_1_RATIO'], features['SSP1HQ_2_RATIO'] = ssp1hq_1/max(1, ssp1), ssp1hq_2/max(1, ssp1)
        features['SSP1mQ_1'], features['SSP1mQ_2'] = Features.get_number_value(record.samples[0], 'SSP1mQ', [Features.NAN, Features.NAN])
        features['SSP1MQ_1'], features['SSP1MQ_2'] = Features.get_number_value(record.samples[0], 'SSP1MQ', [Features.NAN, Features.NAN])

        ssp2 = Features.get_number_value(record.samples[0], 'SSP2', 0)
        ssp2hq_1, ssp2hq_2 = Features.get_number_value(record.samples[0], 'SSP2HQ', [0, 0])
        ssp2nma_1, ssp2nma_2 = Features.get_number_value(record.samples[0], 'SSP2NMA', [Features.NAN, Features.NAN])
        ssp2nms_1, ssp2nms_2 = Features.get_number_value(record.samples[0], 'SSP2NMS', [Features.NAN, Features.NAN])
        ssp2aq_1, ssp2aq_2 = Features.get_number_value(record.samples[0], 'SSP2AQ', [Features.NAN, Features.NAN])
        ssp2sq_1, ssp2sq_2 = Features.get_number_value(record.samples[0], 'SSP2SQ', [Features.NAN, Features.NAN])
        features['SSP2'] = Features.piecewise_normalise(ssp2, min_disc_pairs, max_disc_pairs)
        features['SSP2HQ_1'] = Features.piecewise_normalise(ssp2hq_1, min_disc_pairs, max_disc_pairs)
        features['SSP2HQ_2'] = Features.piecewise_normalise(ssp2hq_2, min_disc_pairs, max_disc_pairs)
        features['SSP2HQ_1_RATIO'], features['SSP2HQ_2_RATIO'] = ssp2hq_1/max(1, ssp2), ssp2hq_2/max(1, ssp2)
        features['SSP2mQ_1'], features['SSP2mQ_2'] = Features.get_number_value(record.samples[0], 'SSP2mQ', [Features.NAN, Features.NAN])
        features['SSP2MQ_1'], features['SSP2MQ_2'] = Features.get_number_value(record.samples[0], 'SSP2MQ', [Features.NAN, Features.NAN])

        if svtype_str in ["DEL", "INS"]:
            ssp1nma_2 = Features.NAN
            ssp2nma_1 = Features.NAN
        elif svtype_str == "DUP":
            ssp1nma_1 = Features.NAN
            ssp2nma_2 = Features.NAN

        features['SSP1_RSP1_1_NM_Z_SCORE'] = Features.calculate_z_score(ssp1nma_1, ssp1nms_1, ssp1, rsp1nma_1, rsp1nms_1, rsp1)
        features['SSP1_RSP1_2_NM_Z_SCORE'] = Features.calculate_z_score(ssp1nma_2, ssp1nms_2, ssp1, rsp1nma_2, rsp1nms_2, rsp1)
        features['SSP2_RSP2_1_NM_Z_SCORE'] = Features.calculate_z_score(ssp2nma_1, ssp2nms_1, ssp2, rsp2nma_1, rsp2nms_1, rsp2)
        features['SSP2_RSP2_2_NM_Z_SCORE'] = Features.calculate_z_score(ssp2nma_2, ssp2nms_2, ssp2, rsp2nma_2, rsp2nms_2, rsp2)

        features['SSP1_RSP1_1_AQ_Z_SCORE'] = Features.calculate_z_score(ssp1aq_1, ssp1sq_1, ssp1, rsp1aq_1, rsp1sq_1, rsp1)
        features['SSP1_RSP1_2_AQ_Z_SCORE'] = Features.calculate_z_score(ssp1aq_2, ssp1sq_2, ssp1, rsp1aq_2, rsp1sq_2, rsp1)
        features['SSP2_RSP2_1_AQ_Z_SCORE'] = Features.calculate_z_score(ssp2aq_1, ssp2sq_1, ssp2, rsp2aq_1, rsp2sq_1, rsp2)
        features['SSP2_RSP2_2_AQ_Z_SCORE'] = Features.calculate_z_score(ssp2aq_2, ssp2sq_2, ssp2, rsp2aq_2, rsp2sq_2, rsp2)

        axr1, axr2 = Features.get_number_value(record.samples[0], 'AXR', [0, 0], median_depth*max_is)
        axr1hq, axr2hq = Features.get_number_value(record.samples[0], 'AXRHQ', [0, 0], median_depth*max_is)
        features['AXR1'], features['AXR2'] = axr1, axr2
        features['AXR1HQ'], features['AXR2HQ'] = axr1hq, axr2hq
        
        exl1 = Features.get_number_value(record.samples[0], 'EXL', 0)
        exl2 = Features.get_number_value(record.samples[0], 'EXL2', 0)
        features['MEXL'] = max(exl1/max_is, exl2/max_is)
        features['mEXL'] = min(exl1/max_is, exl2/max_is)
        features['EXL'] = features['MEXL'] + features['mEXL']

        exas1 = Features.get_number_value(record.samples[0], 'EXAS', 0, max(1, exl1))
        exas2 = Features.get_number_value(record.samples[0], 'EXAS2', 0, max(1, exl2))

        exrs1 = Features.get_number_value(record.samples[0], 'EXRS', 0, max(1, exl1))
        exrs2 = Features.get_number_value(record.samples[0], 'EXRS2', 0, max(1, exl2))

        features['EXAS_EXRS_RATIO'] = 0 if exrs1+exrs2 == 0 else (exas1+exas2)/(exrs1+exrs2)
        features['EXAS_EXRS_DIFF'] = (exas1-exrs1) + (exas2-exrs2)

        exss1_1, exss1_2 = Features.get_number_value(record.samples[0], 'EXSS', [0, 0])
        features['EXSS1_1'] = exss1_1/max_is
        features['EXSS1_2'] = exss1_2/max_is
        features['EXSS1_RATIO1'] = exss1_1/max(1, exl1)
        features['EXSS1_RATIO2'] = exss1_2/max(1, exl1)

        exss2_1, exss2_2 = Features.get_number_value(record.samples[0], 'EXSS2', [0, 0])
        features['EXSS2_1'] = exss2_1/max_is
        features['EXSS2_2'] = exss2_2/max_is
        features['EXSS2_RATIO1'] = exss2_1/max(1, exl2)
        features['EXSS2_RATIO2'] = exss2_2/max(1, exl2)

        exssc1_1, exssc1_2 = Features.get_number_value(record.samples[0], 'EXSSC', [0, 0])
        exssc2_1, exssc2_2 = Features.get_number_value(record.samples[0], 'EXSSC2', [0, 0])
        exsscia1_1, exsscia1_2 = Features.get_number_value(record.samples[0], 'EXSSCIA', [0, 0])
        exsscia2_1, exsscia2_2 = Features.get_number_value(record.samples[0], 'EXSSC2IA', [0, 0])
        features['EXSCC1_1_IA_RATIO'], features['EXSCC1_2_IA_RATIO'] = exssc1_1/max(1, exsscia1_1), exssc1_2/max(1, exsscia1_2)
        features['EXSCC2_1_IA_RATIO'], features['EXSCC2_2_IA_RATIO'] = exssc2_1/max(1, exsscia2_1), exssc2_2/max(1, exsscia2_2)
        features['EXSSC1_1_IA_DIFF_RATIO'], features['EXSSC1_2_IA_DIFF_RATIO'] = (exss1_1-exssc1_1)/max(1, exss1_1-exsscia1_1), (exss1_2-exssc1_2)/max(1, exss1_2-exsscia1_2)
        features['EXSSC2_1_IA_DIFF_RATIO'], features['EXSSC2_2_IA_DIFF_RATIO'] = (exss2_1-exssc2_1)/max(1, exss2_1-exsscia2_1), (exss2_2-exssc2_2)/max(1, exss2_2-exsscia2_2)

        feature_values = []
        for feature_name in Features.get_feature_names(model_name):
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
def parse_vcf(vcf_fname, stats_fname, fp_fname, tolerate_no_gts = False):
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

        model_name = Features.get_model_name(record, stats['max_is']['.'], stats['read_len']['.'])
        feature_values = Features.record_to_features(record, stats)
        if gts[record.id] != "./.": # if too deep or no genotype is available, skip the record
            features_by_source[model_name].append(feature_values)
            gts_by_source[model_name].append(gts[record.id])
            variant_ids_by_source[model_name].append(Features.generate_id(record))

    for model_name in features_by_source:
        features_by_source[model_name] = np.array(features_by_source[model_name])
        gts_by_source[model_name] = np.array(gts_by_source[model_name])
        variant_ids_by_source[model_name] = np.array(variant_ids_by_source[model_name])
    
    return features_by_source, gts_by_source, variant_ids_by_source
