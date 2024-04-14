import argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import joblib
import features
import os

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('svtype', help='SV type to filter.', choices=['DEL', 'DUP', 'INS', 'ALL'])
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--n-trees', type=int, default=5000, help='Number of trees in the random forest.')
cmd_parser.add_argument('--model_name', default='ALL', help='Restrict to this model.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for training.')
cmd_args = cmd_parser.parse_args()

denovo_yes_or_no_output_dir = os.path.join(cmd_args.outdir, "denovo", "yes_or_no")
os.makedirs(denovo_yes_or_no_output_dir, exist_ok=True)

denovo_gts_output_dir = os.path.join(cmd_args.outdir, "denovo", "gts")
os.makedirs(denovo_gts_output_dir, exist_ok=True)

output_dir = os.path.join(cmd_args.outdir, "regt", "yes_or_no")
os.makedirs(output_dir, exist_ok=True)

output_dir = os.path.join(cmd_args.outdir, "regt", "gts")
os.makedirs(output_dir, exist_ok=True)

training_prefixes = cmd_args.training_prefixes.split(",")
denovo_training_data, regt_training_data, training_gts, variant_ids = None, None, None, None
for training_prefix in training_prefixes:
    vcf_denovo_training_data, vcf_regt_training_data, vcf_training_gts, vcf_variant_ids = features.parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats", 
                                                                                                             training_prefix + ".gts", cmd_args.svtype, tolerate_no_gts = False)
    if denovo_training_data is None:
        denovo_training_data = vcf_denovo_training_data
        regt_training_data = vcf_regt_training_data
        training_gts = vcf_training_gts
    else:
        for source in vcf_denovo_training_data:
            denovo_training_data[source] = np.concatenate((denovo_training_data[source], vcf_denovo_training_data[source]))
            regt_training_data[source] = np.concatenate((regt_training_data[source], vcf_regt_training_data[source]))
            training_gts[source] = np.concatenate((training_gts[source], vcf_training_gts[source]))

for model_name in denovo_training_data:
    if cmd_args.model_name != "ALL" and model_name != cmd_args.model_name:
        continue

    classifier = RandomForestClassifier(n_estimators=cmd_args.n_trees, max_depth=15, n_jobs=cmd_args.threads, random_state=42)

    # train yes or no classifier
    training_labels = np.array([0 if x == "0/0" else 1 for x in training_gts[model_name]])
    classifier.fit(denovo_training_data[model_name], training_labels)

    model_fname = os.path.join(denovo_yes_or_no_output_dir, model_name + ".model")
    joblib.dump(classifier, open(model_fname, 'wb'))

    # train GT classifier    
    positive_training_data = denovo_training_data[model_name][(training_gts[model_name] == "0/1") | (training_gts[model_name] == "1/1")]
    positive_training_labels = np.array([2 if x == "1/1" else 1 for x in training_gts[model_name] if x == "0/1" or x == "1/1"])
    classifier.fit(positive_training_data, positive_training_labels)
    
    model_fname = os.path.join(denovo_gts_output_dir, model_name + ".model")
    joblib.dump(classifier, open(model_fname, 'wb'))
