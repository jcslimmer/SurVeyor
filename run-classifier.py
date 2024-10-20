import pysam, argparse
import joblib
import features
import os

cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
cmd_parser.add_argument('in_vcf', help='Input VCF file.')
cmd_parser.add_argument('out_vcf', help='Output VCF file.')
cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
cmd_parser.add_argument('svtype', help='SV type to filter.', choices=['DEL', 'DUP', 'INS', 'ALL'])
cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
cmd_parser.add_argument('--denovo', action='store_true', help='Use denovo models.')
cmd_args = cmd_parser.parse_args()

def write_vcf(vcf_reader, vcf_header, svid_to_gt, fname):
    vcf_writer = pysam.VariantFile(fname, 'w', header=vcf_header)
    for record in vcf_reader.fetch():
        if features.Features.get_svtype(record) != "INV":
            record.info['HARD_FILTERS'] = ",".join(record.filter.keys())
            record.filter.clear()
            record.filter.add('PASS')
            if record.id in svid_to_gt:
                record.samples[0]['GT'] = (svid_to_gt[record.id]//2, 1 if svid_to_gt[record.id] >= 1 else 0)
                record.samples[0]['EPR'] = svid_to_prob[record.id]
            else:
                record.samples[0]['GT'] = (None, None)
            vcf_writer.write(record)
        elif "PASS" in record.filter:
            vcf_writer.write(record)
    vcf_writer.close()

test_denovo_data, test_regt_data, test_denovo_labels, test_regt_labels, test_denovo_variant_ids, test_regt_variant_ids = \
    features.parse_vcf(cmd_args.in_vcf, cmd_args.stats, "XXX", cmd_args.svtype, tolerate_no_gts = True)

if cmd_args.denovo:
    test_data = test_denovo_data
    test_variant_ids = test_denovo_variant_ids
else:
    test_data = test_regt_data
    test_variant_ids = test_regt_variant_ids

svid_to_gt = dict()
svid_to_prob = dict()
for model_name in test_data:
    if cmd_args.denovo:
        model_file = os.path.join(cmd_args.model_dir, "denovo", "yes_or_no", model_name + '.model')
    else:
        model_file = os.path.join(cmd_args.model_dir, "regt", "yes_or_no", model_name + '.model')
    
    if model_name.startswith("INV"):
        continue
    
    classifier = joblib.load(model_file)
    predictions = classifier.predict(test_data[model_name])
    probs = classifier.predict_proba(test_data[model_name])
    for i in range(len(predictions)):
        svid_to_gt[test_variant_ids[model_name][i]] = predictions[i]
        svid_to_prob[test_variant_ids[model_name][i]] = probs[i][1]

    if len(test_data[model_name]) == 0:
        continue

    positive_data = test_data[model_name][predictions == 1]
    positive_variant_ids = test_variant_ids[model_name][predictions == 1]

    if len(positive_data) == 0:
        continue

    if cmd_args.denovo:
        model_file = os.path.join(cmd_args.model_dir, "denovo", "gts", model_name + '.model')
    else:
        model_file = os.path.join(cmd_args.model_dir, "regt", "gts", model_name + '.model')
    classifier = joblib.load(model_file)
    predictions = classifier.predict(positive_data)
    for i in range(len(predictions)):
        svid_to_gt[positive_variant_ids[i]] = predictions[i]

# write the predictions to a VCF file
vcf_reader = pysam.VariantFile(cmd_args.in_vcf)
header = vcf_reader.header
header.add_line('##INFO=<ID=HARD_FILTERS,Number=.,Type=String,Description="PASS or not according to hard filters.">')
header.add_line('##FORMAT=<ID=EPR,Number=1,Type=Float,Description="Probability of the SV existing in the sample, according to the ML model.">')
write_vcf(vcf_reader, header, svid_to_gt, cmd_args.out_vcf)
