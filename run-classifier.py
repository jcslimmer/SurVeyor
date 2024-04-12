import pysam, argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict
import joblib
import features_gt as features

cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
cmd_parser.add_argument('in_vcf', help='Input VCF file.')
cmd_parser.add_argument('out_vcf', help='Output VCF file.')
cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
cmd_parser.add_argument('svtype', help='SV type to filter.', choices=['DEL', 'DUP', 'INS', 'ALL'])
cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
cmd_args = cmd_parser.parse_args()

test_data, test_labels, test_variant_ids = features.parse_vcf(cmd_args.in_vcf, cmd_args.stats, "XXX", cmd_args.svtype, tolerate_no_gts = True)

def write_vcf(vcf_reader, vcf_header, svid_to_gt, fname):
    vcf_writer = pysam.VariantFile(fname, 'w', header=vcf_header)
    for record in vcf_reader.fetch():
        record.info['HARD_FILTERS'] = ",".join(record.filter.keys())
        record.filter.clear()
        record.filter.add('PASS')
        record.samples[0]['GT'] = (svid_to_gt[record.id]//2, 1 if svid_to_gt[record.id] >= 1 else 0)
        vcf_writer.write(record)
    vcf_writer.close()

svid_to_gt = dict()
for model_name in test_data:
    classifier = joblib.load(cmd_args.model_dir + '/yes_or_no/' + model_name + '.model')
    predictions = classifier.predict(test_data[model_name])
    for i in range(len(predictions)):
        svid_to_gt[test_variant_ids[model_name][i]] = predictions[i]
    
    positive_data = test_data[model_name][predictions == 1]
    positive_variant_ids = test_variant_ids[model_name][predictions == 1]

    classifier = joblib.load(cmd_args.model_dir + '/gts/' + model_name + '.model')
    predictions = classifier.predict(positive_data)
    for i in range(len(predictions)):
        svid_to_gt[positive_variant_ids[i]] = predictions[i]

# write the predictions to a VCF file
vcf_reader = pysam.VariantFile(cmd_args.in_vcf)
header = vcf_reader.header
header.add_line('##INFO=<ID=HARD_FILTERS,Number=.,Type=String,Description="PASS or not according to hard filters.">')
write_vcf(vcf_reader, header, svid_to_gt, cmd_args.out_vcf)
