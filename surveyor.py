import sys, os, argparse, pysam, timeit

VERSION = "0.1"

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
call_parser.add_argument('--genotype', action='store_true', help='Genotype the SVs.')

genotype_parser = subparsers.add_parser('genotype', parents=[common_parser], help='Genotype SVs.')
genotype_parser.add_argument('in_vcf_file', help='Input VCF file.')
genotype_parser.add_argument('out_vcf_file', help='Output VCF file.')
genotype_parser.add_argument('bam_file', help='Input bam file.')
genotype_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
genotype_parser.add_argument('reference', help='Reference genome in FASTA format.')
genotype_parser.add_argument('--use-call-info', action='store_true', help='Reuse info in the workdir stored by the call commands. Assumes the workdir is the same used by the call command, and no file has been deleted.')

cmd_args = parser.parse_args()

def exec(cmd):
    start_time = timeit.default_timer()
    print("Executing:", cmd)
    if os.system(cmd) != 0:
        print("Error executing:", cmd)
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

def reads_categorizer():
    mkdir(cmd_args.workdir)
    mkdir(cmd_args.workdir + "/workspace")
    mkdir(cmd_args.workdir + "/workspace/clipped")
    mkdir(cmd_args.workdir + "/workspace/hsr")
    mkdir(cmd_args.workdir + "/workspace/fwd-stable")
    mkdir(cmd_args.workdir + "/workspace/rev-stable")
    mkdir(cmd_args.workdir + "/workspace/long-pairs")
    mkdir(cmd_args.workdir + "/workspace/outward-pairs")
    mkdir(cmd_args.workdir + "/workspace/mateseqs")
    mkdir(cmd_args.workdir + "/workspace/sc_mateseqs")
    mkdir(cmd_args.workdir + "/workspace/sr_consensuses")
    mkdir(cmd_args.workdir + "/workspace/hsr_consensuses")

    with open("%s/contig_map" % cmd_args.workdir, "w") as contig_map:
        bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
        for k in bam_file.references:
            contig_map.write("%s\n" % (k))

    read_categorizer_cmd = SURVEYOR_PATH + "/bin/reads_categorizer %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference)
    exec(read_categorizer_cmd)

if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]

if cmd_args.command == 'call':
    reads_categorizer()

    clip_consensus_builder_cmd = SURVEYOR_PATH + "/bin/clip_consensus_builder %s %s" % (cmd_args.workdir, cmd_args.reference)
    exec(clip_consensus_builder_cmd)


    find_svs_from_sr_consensuses_cmd = SURVEYOR_PATH + "/bin/find_svs_from_sr_consensuses %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
    exec(find_svs_from_sr_consensuses_cmd)

    merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/merge_identical_calls %s/sr.vcf.gz %s/sr.dedup.vcf.gz %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
    exec(merge_identical_calls_cmd)

    dp_clusterer = SURVEYOR_PATH + "/bin/dp_clusterer %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
    exec(dp_clusterer)

    ins_assembler_cmd = SURVEYOR_PATH + "/bin/insertions_assembler %s %s %s" % (cmd_args.workdir, cmd_args.reference, sample_name)
    exec(ins_assembler_cmd)

    concat_cmd = SURVEYOR_PATH + "/bin/concat_vcf %s/sr_dp.vcf.gz %s/assembled_ins.vcf.gz %s/out.vcf.gz" % (cmd_args.workdir, cmd_args.workdir, cmd_args.workdir)
    exec(concat_cmd)

    add_filtering_info_cmd = SURVEYOR_PATH + "/bin/add_filtering_info %s %s/out.vcf.gz %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.workdir, cmd_args.reference, sample_name)
    exec(add_filtering_info_cmd)

    normalise_cmd = SURVEYOR_PATH + "/bin/normalise %s/out.annotated.vcf.gz %s/out.annotated.norm.vcf.gz %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
    exec(normalise_cmd)

    merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/merge_identical_calls %s/out.annotated.norm.vcf.gz %s/out.annotated.norm.dedup.vcf.gz %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.reference)
    exec(merge_identical_calls_cmd)

    if cmd_args.genotype:
        genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s/out.annotated.norm.dedup.vcf.gz %s/out.annotated.norm.dedup.gt.vcf.gz %s %s %s %s" % (cmd_args.workdir, cmd_args.workdir, cmd_args.bam_file, cmd_args.reference, cmd_args.workdir, sample_name)
        exec(genotype_cmd)

elif cmd_args.command == 'genotype':
    if not use_call_info():
        reads_categorizer()
    
    genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s %s %s %s %s %s" % (cmd_args.in_vcf_file, cmd_args.out_vcf_file, cmd_args.bam_file, cmd_args.reference, cmd_args.workdir, sample_name)
    exec(genotype_cmd)
