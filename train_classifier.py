import argparse, time
import numpy as np
import features
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import xgboost as xgb

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--model_name', default='ALL', help='Restrict to this model.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for training.')
cmd_parser.add_argument('--cross-species', action='store_true', help='Use cross-species model.')
cmd_args = cmd_parser.parse_args()

training_data, training_gts = defaultdict(list), defaultdict(list)

def process_vcf(training_prefix):
    vcf_training_data, vcf_training_gts, _ = \
        features.parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats", training_prefix + ".gts", 
                           tolerate_no_gts = False)
    ins_to_dup_vcf_training_data, ins_to_dup_vcf_training_gts, _ = \
        features.parse_vcf(training_prefix + ".INS_TO_DUP.vcf.gz", training_prefix + ".stats",
                            training_prefix + ".gts", tolerate_no_gts = False)
    
    vcf_training_data["INS_TO_DUP"] = ins_to_dup_vcf_training_data["INS_TO_DUP"]
    vcf_training_gts["INS_TO_DUP"] = ins_to_dup_vcf_training_gts["INS_TO_DUP"]
    vcf_training_data["INS_TO_DUP_LARGE"] = ins_to_dup_vcf_training_data["INS_TO_DUP_LARGE"]
    vcf_training_gts["INS_TO_DUP_LARGE"] = ins_to_dup_vcf_training_gts["INS_TO_DUP_LARGE"]
    return vcf_training_data, vcf_training_gts

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.set_start_method("spawn", force=True)  # optional but recommended

    training_prefixes = cmd_args.training_prefixes.split(",")
    with ProcessPoolExecutor(max_workers=cmd_args.threads) as executor:
        future_to_prefix = {executor.submit(process_vcf, prefix): prefix for prefix in training_prefixes}
        for future in as_completed(future_to_prefix):
            vcf_training_data, vcf_training_gts = future.result()
            for model in vcf_training_data:
                training_data[model].append(vcf_training_data[model])
                training_gts[model].append(vcf_training_gts[model])

    for model in training_data:
        training_data[model] = np.concatenate(training_data[model])
        training_gts[model] = np.concatenate(training_gts[model])

    yes_or_no_outdir = os.path.join(cmd_args.outdir, "yes_or_no")
    os.makedirs(yes_or_no_outdir, exist_ok=True)

    gts_outdir = os.path.join(cmd_args.outdir, "gts")
    os.makedirs(gts_outdir, exist_ok=True)

    for model_name in training_data:
        if cmd_args.model_name != "ALL" and model_name != cmd_args.model_name:
            continue

        training_labels = np.array([0 if x == "0/0" else 1 for x in training_gts[model_name]])
        num_pos = np.sum(training_labels == 1)
        num_neg = np.sum(training_labels == 0)
        scale_pos_weight = max(num_neg/(num_pos*10), 1.0) if num_pos > 0 else 1.0
        print(f"Training model {model_name} with scale_pos_weight={scale_pos_weight}, num_pos={num_pos}, num_neg={num_neg}")

        if cmd_args.cross_species:
            classifier = xgb.XGBClassifier(n_estimators=50, max_depth=7, min_child_weight=42, learning_rate=0.1, n_jobs=cmd_args.threads, random_state=42, tree_method='hist')
        else:
            classifier = xgb.XGBClassifier(n_estimators=1000, scale_pos_weight=scale_pos_weight, max_depth=8, min_child_weight=16, learning_rate=0.1, n_jobs=cmd_args.threads, random_state=42, tree_method='hist')

        start_time = time.time()
        classifier.fit(training_data[model_name], training_labels)

        features_names = features.Features.get_feature_names(model_name)
        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]

        with open(os.path.join(yes_or_no_outdir, model_name + ".importance.txt"), 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

        model_fname = os.path.join(yes_or_no_outdir, model_name + ".ubj")
        classifier.save_model(model_fname)

        positive_training_data = training_data[model_name][(training_gts[model_name] == "0/1") | (training_gts[model_name] == "1/1")]
        positive_training_labels = np.array([1 if x == "1/1" else 0 for x in training_gts[model_name] if x == "0/1" or x == "1/1"])

        unique_labels = np.unique(positive_training_labels)
        if len(unique_labels) == 1:
            print(f"Only one label present in positive training data for model {model_name}.")
            first_entry = positive_training_data[0:1]
            first_label = positive_training_labels[0]
            opposite_label = 1 - first_label  
            positive_training_data = np.concatenate([positive_training_data, first_entry], axis=0)
            positive_training_labels = np.concatenate([positive_training_labels, np.array([opposite_label])])

        unique, counts = np.unique(positive_training_labels, return_counts=True)

        classifier.fit(positive_training_data, positive_training_labels)

        model_fname = os.path.join(gts_outdir, model_name + ".ubj")
        classifier.save_model(model_fname)

        end_time = time.time()
        print(f"Training for model {model_name} took {end_time - start_time} seconds")
