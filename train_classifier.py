import argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import joblib
import features
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--n-trees', type=int, default=5000, help='Number of trees in the random forest.')
cmd_parser.add_argument('--model_name', default='ALL', help='Restrict to this model.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for training.')
cmd_parser.add_argument('--denovo', action='store_true', help='Build denovo model.')
cmd_args = cmd_parser.parse_args()

denovo_yes_or_no_outdir = os.path.join(cmd_args.outdir, "denovo", "yes_or_no")
os.makedirs(denovo_yes_or_no_outdir, exist_ok=True)

denovo_gts_outdir = os.path.join(cmd_args.outdir, "denovo", "gts")
os.makedirs(denovo_gts_outdir, exist_ok=True)

regt_yes_or_no_outdir = os.path.join(cmd_args.outdir, "regt", "yes_or_no")
os.makedirs(regt_yes_or_no_outdir, exist_ok=True)

regt_gts_outdir = os.path.join(cmd_args.outdir, "regt", "gts")
os.makedirs(regt_gts_outdir, exist_ok=True)

training_data, training_gts = defaultdict(list), defaultdict(list)

def process_vcf(training_prefix):
    vcf_training_data, vcf_training_gts, _ = \
        features.parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats", training_prefix + ".gts", 
                           cmd_args.denovo, tolerate_no_gts = False)
    return vcf_training_data, vcf_training_gts

def merge_data(vcf_training_data, vcf_training_gts):
    global training_data, training_gts
    for source in vcf_training_data:
        if len(training_data[source]) == 0:
            training_data[source] = vcf_training_data[source]
            training_gts[source] = vcf_training_gts[source]
        else:
            training_data[source] = np.concatenate((training_data[source], vcf_training_data[source]))
            training_gts[source] = np.concatenate((training_gts[source], vcf_training_gts[source]))

training_prefixes = cmd_args.training_prefixes.split(",")
with ProcessPoolExecutor(max_workers=cmd_args.threads) as executor:
    future_to_prefix = {executor.submit(process_vcf, prefix): prefix for prefix in training_prefixes}
    for future in as_completed(future_to_prefix):
        vcf_training_data, vcf_training_gts = future.result()
        merge_data(vcf_training_data, vcf_training_gts)

classifier = RandomForestClassifier(n_estimators=cmd_args.n_trees, max_depth=15, n_jobs=cmd_args.threads, random_state=42)

yes_or_no_outdir = denovo_yes_or_no_outdir if cmd_args.denovo else regt_yes_or_no_outdir
gts_outdir = denovo_gts_outdir if cmd_args.denovo else regt_gts_outdir

for model_name in training_data:
    if cmd_args.model_name != "ALL" and model_name != cmd_args.model_name:
        continue
    
    training_labels = np.array([0 if x == "0/0" else 1 for x in training_gts[model_name]])
    classifier.fit(training_data[model_name], training_labels)

    # print feature importance to file
    if cmd_args.denovo:
        features_names = features.Features.get_denovo_feature_names(model_name)
    else:
        features_names = features.Features.get_regt_feature_names(model_name)
    importances = classifier.feature_importances_
    indices = np.argsort(importances)[::-1]
    
    with open(os.path.join(yes_or_no_outdir, model_name + ".importance.txt"), 'w') as f:
        for i in range(len(features_names)):
            f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

    model_fname = os.path.join(yes_or_no_outdir, model_name + ".model")
    joblib.dump(classifier, open(model_fname, 'wb'))

    positive_training_data = training_data[model_name][(training_gts[model_name] == "0/1") | (training_gts[model_name] == "1/1")]
    positive_training_labels = np.array([2 if x == "1/1" else 1 for x in training_gts[model_name] if x == "0/1" or x == "1/1"])
    classifier.fit(positive_training_data, positive_training_labels)

    model_fname = os.path.join(gts_outdir, model_name + ".model")
    joblib.dump(classifier, open(model_fname, 'wb'))
