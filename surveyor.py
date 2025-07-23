import sys, os, argparse, pysam, timeit
from run_classifier import Classifier

VERSION = "0.8"

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

parser = argparse.ArgumentParser(description='SurVeyor, an SV caller.')
parser.add_argument('--version', action='version', version="SurVeyor v%s" % VERSION, help='Print version number.')

subparsers = parser.add_subparsers(dest='command', help='Commands', required=True)

common_parser = argparse.ArgumentParser(add_help=False)
common_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
common_parser.add_argument('--seed', type=int, default=0, help='Seed for random sampling of genomic positions.')
common_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                         'If not provided, the basename of the bam/cram file will be used,'
                                                         'up until the first \'.\'')
common_parser.add_argument('--min-sv-size', type=int, default=50, help='Min SV size.')
common_parser.add_argument('--min-clip-len', type=int, default=15, help='Min length for a clip to be used.')
common_parser.add_argument('--max-seq-error', type=float, default=0.04, help='Max sequencing error admissible on the platform used.')
common_parser.add_argument('--max-clipped-pos-dist', type=int, default=5, help='Max distance (in bp) for two clips to be considered '
                                                                   'representing the same breakpoint.')
common_parser.add_argument('--sampling-regions', help='File in BED format containing a list of regions to be used to estimate'
                                                   'statistics such as depth.')
common_parser.add_argument('--per-contig-stats', action='store_true',
                        help='Depth statistics are computed separately for each contig. Useful when one or more of the target contigs are expected to have '
                        'dramatically different depth than others. Otherwise, it is not recommended to use this option.')
common_parser.add_argument('--generate-training-data', action='store_true', help='Generate data needed to train a genotyping ML model.')
common_parser.add_argument('--tr-bed', help='BED file with tandem repetitive regions. If provided, it will be used for a more aggrestive duplicate removal.')

# SurVIndel2 specific arguments
common_parser.add_argument('--min_size_for_depth_filtering', type=int, default=1000, help='Minimum size for depth filtering.')
common_parser.add_argument('--min-diff-hsr', type=int, default=3, help='Minimum number of differences with the reference \
                        (considered as number of insertions, deletions and mismatches) for a read to be considered a hidden split read.')

# INSurVeyor specific arguments
common_parser.add_argument('--max-trans-size', type=int, default=10000, help='Maximum size of the transpositions which '
                                                                          'SurVeyor will predict when only one side is available.')
common_parser.add_argument('--min-stable-mapq', type=int, default=20, help='Minimum MAPQ for a stable read.')

call_parser = subparsers.add_parser('call', parents=[common_parser], help='Call SVs denovo.')
call_parser.add_argument('bam_file', help='Input bam file.')
call_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
call_parser.add_argument('reference', help='Reference genome in FASTA format.')
call_parser.add_argument('--ml-model', help='Path to the ML model to be used for filtering and genotyping.')

genotype_parser = subparsers.add_parser('genotype', parents=[common_parser], help='Genotype SVs.')
genotype_parser.add_argument('--use-call-info', action='store_true', help='Reuse info in the workdir stored by the call commands. Assumes the workdir is the same used by the call command, and no file has been deleted.')
genotype_parser.add_argument('--two-pass', action='store_true', help='Activate two-pass genotyping.')
genotype_parser.add_argument('in_vcf_file', help='Input VCF file.')
genotype_parser.add_argument('bam_file', help='Input bam file.')
genotype_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
genotype_parser.add_argument('reference', help='Reference genome in FASTA format.')
genotype_parser.add_argument('ml_model', help='Path to the ML model to be used for genotyping.')

cmd_args = parser.parse_args()

def run_cmd(cmd, error_msg=None):
    start_time = timeit.default_timer()
    print("Executing:", cmd)
    return_code = os.system(cmd)
    if return_code != 0:
        if error_msg:
            print(error_msg)
        else:
            print("Error executing:", cmd)
            print("Return code:", return_code)
        exit(1)
    elapsed = timeit.default_timer() - start_time
    print(cmd, "was run in %.2f seconds" % elapsed)

SURVEYOR_PATH = os.path.dirname(os.path.realpath(__file__))

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def use_call_info():
    if not hasattr(cmd_args, 'use_call_info') or not cmd_args.use_call_info:
        return False
    return True

if not use_call_info():
    with open(cmd_args.workdir + "/call_cmd.txt", "w") as full_cmd_out:
        print(" ".join(sys.argv[:]), file=full_cmd_out)

with open(cmd_args.workdir + "/config.txt", "w") as config_file:
    config_file.write("threads %d\n" % cmd_args.threads)
    config_file.write("seed %d\n" % cmd_args.seed)
    config_file.write("min_sv_size %s\n" % cmd_args.min_sv_size)
    config_file.write("min_clip_len %s\n" % cmd_args.min_clip_len)
    config_file.write("max_seq_error %s\n" % cmd_args.max_seq_error)
    config_file.write("max_clipped_pos_dist %d\n" % cmd_args.max_clipped_pos_dist)
    if cmd_args.sampling_regions:
        config_file.write("sampling_regions %s\n" % cmd_args.sampling_regions)
    config_file.write("per_contig_stats %d\n" % cmd_args.per_contig_stats)
    config_file.write("version %s\n" % VERSION)

    config_file.write("min_size_for_depth_filtering %s\n" % cmd_args.min_size_for_depth_filtering)
    config_file.write("min_diff_hsr %s\n" % cmd_args.min_diff_hsr)
    config_file.write("max_trans_size %d\n" % cmd_args.max_trans_size)
    config_file.write("min_stable_mapq %d\n" % cmd_args.min_stable_mapq)

mkdir(cmd_args.workdir + "/intermediate_results")

def reads_categorizer():
    mkdir(cmd_args.workdir)
    mkdir(cmd_args.workdir + "/workspace")
    mkdir(cmd_args.workdir + "/workspace/clipped")
    mkdir(cmd_args.workdir + "/workspace/hsr")
    mkdir(cmd_args.workdir + "/workspace/fwd-stable")
    mkdir(cmd_args.workdir + "/workspace/rev-stable")
    mkdir(cmd_args.workdir + "/workspace/long-pairs")
    mkdir(cmd_args.workdir + "/workspace/outward-pairs")
    mkdir(cmd_args.workdir + "/workspace/same-strand")
    mkdir(cmd_args.workdir + "/workspace/mateseqs")
    mkdir(cmd_args.workdir + "/workspace/sr_consensuses")
    mkdir(cmd_args.workdir + "/workspace/hsr_consensuses")

    with open("%s/contig_map" % cmd_args.workdir, "w") as contig_map:
        bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
        for k in bam_file.references:
            contig_map.write("%s\n" % (k))

    read_categorizer_cmd = SURVEYOR_PATH + "/bin/reads_categorizer %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference)
    run_cmd(read_categorizer_cmd)

def deduplicate_vcf(vcf_fname, deduped_vcf_fname):
    if cmd_args.tr_bed:
        compare_cmd = SURVEYOR_PATH + "/bin/compare %s %s -R %s -T %s > %s/compare.txt" % (vcf_fname, vcf_fname, cmd_args.reference, cmd_args.tr_bed, cmd_args.workdir)
    else:
        compare_cmd = SURVEYOR_PATH + "/bin/compare %s %s -R %s > %s/compare.txt" % (vcf_fname, vcf_fname, cmd_args.reference, cmd_args.workdir) 
    run_cmd(compare_cmd)

    with open(cmd_args.workdir + "/compare.txt") as compare_file:
        epr_vals, imprecise_vals = {}, {}
        with pysam.VariantFile(vcf_fname) as vcf:
            for record in vcf:
                epr_vals[record.id] = record.samples[0].get('EPR', 0)
                imprecise_vals[record.id] = 'IMPRECISE' in record.info

        removed_ids = set()
        for line in compare_file:
            id1, id2 = line.strip().split()[:2]
            if id1 >= id2 or id2 == "NONE": # each pair will be output twice, e.g. A B and B A. Let us process it once
                continue

            if imprecise_vals[id1] and not imprecise_vals[id2]:
                removed_ids.add(id1)
            elif imprecise_vals[id2] and not imprecise_vals[id1]:
                removed_ids.add(id2)
            elif epr_vals[id1] > epr_vals[id2]:
                removed_ids.add(id2)
            elif epr_vals[id2] > epr_vals[id1]:
                removed_ids.add(id1)

        vcf_deduped = deduped_vcf_fname
        with pysam.VariantFile(vcf_fname) as vcf, pysam.VariantFile(vcf_deduped, 'w', header=vcf.header) as out_vcf:
            for record in vcf:
                if record.id not in removed_ids:
                    out_vcf.write(record)

def separate_ins_to_dup(in_vcf_fname, ins_to_dup_vcf_fname, remaining_vcf_fname):
    with pysam.VariantFile(in_vcf_fname) as in_vcf, \
         pysam.VariantFile(ins_to_dup_vcf_fname, 'w', header=in_vcf.header) as ins_to_dup_vcf, \
         pysam.VariantFile(remaining_vcf_fname, 'w', header=in_vcf.header) as remaining_vcf:
        for record in in_vcf:
            if "INS_TO_DUP" in record.info:
                record.id = record.id[:-4] # remove the _DUP suffix
                ins_to_dup_vcf.write(record)
            else:
                remaining_vcf.write(record)

if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]

if cmd_args.command == 'call':

    if not cmd_args.ml_model and not cmd_args.generate_training_data:
        print("At least one of --ml-model or --generate-training-data must be provided.")
        exit(0)

    reads_categorizer()

    clip_consensus_builder_cmd = SURVEYOR_PATH + "/bin/clip_consensus_builder %s %s" % (cmd_args.workdir, cmd_args.reference)
    run_cmd(clip_consensus_builder_cmd)

    find_svs_from_sr_consensuses_cmd = SURVEYOR_PATH + "/bin/find_svs_from_sr_consensuses %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
    run_cmd(find_svs_from_sr_consensuses_cmd)

    merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/merge_identical_calls %s/intermediate_results/sr.vcf.gz %s/intermediate_results/sr.dedup.vcf.gz %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
    run_cmd(merge_identical_calls_cmd)

    dp_clusterer = SURVEYOR_PATH + "/bin/dp_clusterer %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
    run_cmd(dp_clusterer)

    ins_assembler_cmd = SURVEYOR_PATH + "/bin/insertions_assembler %s %s %s" % (cmd_args.workdir, cmd_args.reference, sample_name)
    run_cmd(ins_assembler_cmd)

    concat_cmd = SURVEYOR_PATH + "/bin/concat_vcf %s/intermediate_results/sr_dp.vcf.gz %s/intermediate_results/assembled_ins.vcf.gz %s/intermediate_results/out.vcf.gz" % (cmd_args.workdir, cmd_args.workdir, cmd_args.workdir)
    run_cmd(concat_cmd)

    normalise_cmd = SURVEYOR_PATH + "/bin/normalise %s/intermediate_results/out.vcf.gz %s/intermediate_results/out.norm.vcf.gz %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
    run_cmd(normalise_cmd)

    merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/merge_identical_calls %s/intermediate_results/out.norm.vcf.gz %s/intermediate_results/calls-raw.vcf.gz %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
    run_cmd(merge_identical_calls_cmd)

    insertions_to_duplications_cmd = SURVEYOR_PATH + "/bin/insertions_to_duplications %s/intermediate_results/calls-raw.vcf.gz %s/intermediate_results/calls-for-genotyping.vcf.gz %s %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference, cmd_args.workdir)
    run_cmd(insertions_to_duplications_cmd)
    
    genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s/intermediate_results/calls-for-genotyping.vcf.gz %s/intermediate_results/calls-with-fmt.vcf.gz %s %s %s %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.bam_file, cmd_args.reference, cmd_args.workdir, sample_name)
    run_cmd(genotype_cmd)
    
    if cmd_args.generate_training_data:
        separate_ins_to_dup(cmd_args.workdir + "/intermediate_results/calls-with-fmt.vcf.gz", cmd_args.workdir + "/training-data.INS_TO_DUP.vcf.gz", cmd_args.workdir + "/training-data.vcf.gz")

    if not cmd_args.ml_model:
        print("No model provided. Skipping filtering and genotyping.")
        exit(0)

    Classifier.run_classifier(cmd_args.workdir + "/intermediate_results/calls-with-fmt.vcf.gz", cmd_args.workdir + "/intermediate_results/calls-with-gt.vcf.gz", cmd_args.workdir + "/stats.txt", cmd_args.ml_model)

    reconcile_vcf_gt_cmd = SURVEYOR_PATH + "/bin/reconcile_vcf_gt %s %s %s %s" % (cmd_args.workdir + "/intermediate_results/calls-raw.vcf.gz", cmd_args.workdir + "/intermediate_results/calls-with-gt.vcf.gz", cmd_args.workdir + "/calls-genotyped.vcf.gz", sample_name)
    run_cmd(reconcile_vcf_gt_cmd)

    deduplicate_vcf(cmd_args.workdir + "/calls-genotyped.vcf.gz", cmd_args.workdir + "/calls-genotyped-deduped.vcf.gz")

elif cmd_args.command == 'genotype':

    check_duplicate_ids_cmd = SURVEYOR_PATH + "/bin/check_duplicate_ids %s" % cmd_args.in_vcf_file
    run_cmd(check_duplicate_ids_cmd, "Error: Duplicate IDs found in the input VCF file. Please remove duplicates before running the genotype command.")

    if not use_call_info():
        reads_categorizer()

    vcf_for_genotyping_fname = cmd_args.workdir + "/intermediate_results/calls-for-genotyping.vcf.gz"
    insertions_to_duplications_cmd = SURVEYOR_PATH + "/bin/insertions_to_duplications %s %s %s %s" % (cmd_args.in_vcf_file, vcf_for_genotyping_fname, cmd_args.reference, cmd_args.workdir)
    run_cmd(insertions_to_duplications_cmd)

    vcf_with_fmt_fname = cmd_args.workdir + "/intermediate_results/vcf_with_fmt.vcf.gz"
    genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s %s %s %s %s %s" % (vcf_for_genotyping_fname, vcf_with_fmt_fname, cmd_args.bam_file, cmd_args.reference, cmd_args.workdir, sample_name)
    run_cmd(genotype_cmd)

    if cmd_args.generate_training_data:
        separate_ins_to_dup(cmd_args.workdir + "/intermediate_results/vcf_with_fmt.vcf.gz", cmd_args.workdir + "/training-data.INS_TO_DUP.vcf.gz", cmd_args.workdir + "/training-data.vcf.gz")

    vcf_with_gt_fname = cmd_args.workdir + "/intermediate_results/vcf_with_gt.vcf.gz"
    Classifier.run_classifier(vcf_with_fmt_fname, vcf_with_gt_fname, cmd_args.workdir + "/stats.txt", cmd_args.ml_model)

    if cmd_args.two_pass:
        vcf_with_fmt_reassigned_fname = cmd_args.workdir + "/intermediate_results/vcf_with_fmt.reassigned.vcf.gz"
        genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s %s %s %s %s %s --reassign-evidence" % (vcf_with_gt_fname, vcf_with_fmt_reassigned_fname, cmd_args.bam_file, cmd_args.reference, cmd_args.workdir, sample_name)
        run_cmd(genotype_cmd)

        vcf_with_gt_fname = cmd_args.workdir + "/intermediate_results/vcf_with_gt.reassigned.vcf.gz"
        Classifier.run_classifier(vcf_with_fmt_reassigned_fname, vcf_with_gt_fname, cmd_args.workdir + "/stats.txt", cmd_args.ml_model)

    reconcile_out_vcf = cmd_args.workdir + "/intermediate_results/vcf_with_gt.reconciled.vcf.gz"
    reconcile_vcf_gt_cmd = SURVEYOR_PATH + "/bin/reconcile_vcf_gt %s %s %s %s" % (cmd_args.in_vcf_file, vcf_with_gt_fname, reconcile_out_vcf, sample_name)
    run_cmd(reconcile_vcf_gt_cmd)

    out_vcf_file = cmd_args.workdir + "/genotyped.vcf.gz"
    cp_cmd = "cp %s %s" % (reconcile_out_vcf, out_vcf_file)
    run_cmd(cp_cmd)

    out_deduped_vcf_file = cmd_args.workdir + "/genotyped.deduped.vcf.gz"
    deduplicate_vcf(out_vcf_file, out_deduped_vcf_file)

    out_resolved_vcf_file = cmd_args.workdir + "/genotyped.resolved.vcf.gz"
    resolve_incompatible_cmd = SURVEYOR_PATH + "/bin/resolve_incompatible_gts %s %s" % (out_vcf_file, out_resolved_vcf_file)
    run_cmd(resolve_incompatible_cmd)

    out_resolved_deduped_vcf_file = cmd_args.workdir + "/genotyped.resolved.deduped.vcf.gz"
    deduplicate_vcf(out_resolved_vcf_file, out_resolved_deduped_vcf_file)
