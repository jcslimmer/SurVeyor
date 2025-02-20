import pysam, argparse
import joblib
import features
import os

class Classifier:
    def write_vcf(vcf_reader, vcf_header, svid_to_gt, svid_to_prob, fname):
        vcf_writer = pysam.VariantFile(fname, 'w', header=vcf_header)
        for record in vcf_reader.fetch():
            if features.Features.get_svtype(record) != "INV":
                record.filter.clear()
                record_id = features.Features.generate_id(record)
                if record_id in svid_to_gt:
                    record.samples[0]['GT'] = (svid_to_gt[record_id]//2, 1 if svid_to_gt[record_id] >= 1 else 0)
                    record.samples[0]['EPR'] = float(svid_to_prob[record_id])
                else:
                    record.samples[0]['GT'] = (None, None)
                vcf_writer.write(record)
            else:
                vcf_writer.write(record)
        vcf_writer.close()

    def run_classifier(in_vcf, out_vcf, stats_fname, model_dir):
        test_data, _, test_variant_ids = \
            features.parse_vcf(in_vcf, stats_fname, "XXX", tolerate_no_gts = True)

        svid_to_gt = dict()
        svid_to_prob = dict()
        for model_name in test_data:
            model_file = os.path.join(model_dir, "yes_or_no", model_name + '.model')

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

            model_file = os.path.join(model_dir, "gts", model_name + '.model')
            classifier = joblib.load(model_file)
            predictions = classifier.predict(positive_data)
            for i in range(len(predictions)):
                svid_to_gt[positive_variant_ids[i]] = predictions[i] + 1

        # write the predictions to a VCF file
        vcf_reader = pysam.VariantFile(in_vcf)
        header = vcf_reader.header
        header.add_line('##FORMAT=<ID=EPR,Number=1,Type=Float,Description="Probability of the SV existing in the sample, according to the ML model.">')
        Classifier.write_vcf(vcf_reader, header, svid_to_gt, svid_to_prob, out_vcf)

if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
    cmd_parser.add_argument('in_vcf', help='Input VCF file.')
    cmd_parser.add_argument('out_vcf', help='Output VCF file.')
    cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
    cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
    cmd_args = cmd_parser.parse_args()
    Classifier.run_classifier(cmd_args.in_vcf, cmd_args.out_vcf, cmd_args.stats, cmd_args.model_dir)
