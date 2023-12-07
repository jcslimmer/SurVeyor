import sys, os, argparse, pysam, pyfaidx, timeit
from random_pos_generator import RandomPositionGenerator
import numpy as np

VERSION = "0.1"

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

cmd_parser = argparse.ArgumentParser(description='SurVeyor, an SV caller.')
cmd_parser.add_argument('bam_file', help='Input bam file.')
cmd_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
cmd_parser.add_argument('reference', help='Reference genome in FASTA format.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--seed', type=int, default=0, help='Seed for random sampling of genomic positions.')
cmd_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                         'If not provided, the basename of the bam/cram file will be used,'
                                                         'up until the first \'.\'')
cmd_parser.add_argument('--min-sv-size', type=int, default=50, help='Min SV size.')
cmd_parser.add_argument('--min-clip-len', type=int, default=15, help='Min length for a clip to be used.')
cmd_parser.add_argument('--max-seq-error', type=float, default=0.04, help='Max sequencing error admissible on the platform used.')
cmd_parser.add_argument('--max-clipped-pos-dist', type=int, default=5, help='Max distance (in bp) for two clips to be considered '
                                                                   'representing the same breakpoint.')
cmd_parser.add_argument('--sampling-regions', help='File in BED format containing a list of regions to be used to estimate'
                                                   'statistics such as depth.')
cmd_parser.add_argument('--per-contig-stats', action='store_true',
                        help='Statistics are computed separately for each contig (experimental).')
cmd_parser.add_argument('--log', action='store_true', help='Activate in-depth logging.')
cmd_parser.add_argument('--version', action='version', version="SurVeyor v%s" % VERSION, help='Print version number.')

# SurVIndel2 specific arguments
cmd_parser.add_argument('--min_size_for_depth_filtering', type=int, default=1000, help='Minimum size for depth filtering.')
cmd_parser.add_argument('--min-diff-hsr', type=int, default=3, help='Minimum number of differences with the reference \
                        (considered as number of insertions, deletions and mismatches) for a read to be considered a hidden split read.')

# INSurVeyor specific arguments
cmd_parser.add_argument('--max-trans-size', type=int, default=10000, help='Maximum size of the transpositions which '
                                                                          'SurVeyor will predict when only one side is available.')
cmd_parser.add_argument('--min-stable-mapq', type=int, default=20, help='Minimum MAPQ for a stable read.')
cmd_args = cmd_parser.parse_args()

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

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

with open(cmd_args.workdir + "/full_cmd.txt", "w") as full_cmd_out:
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
    config_file.write("log %d\n" % cmd_args.log)
    config_file.write("version %s\n" % VERSION)

    # SurVIndel2 specific
    config_file.write("min_size_for_depth_filtering %s\n" % cmd_args.min_size_for_depth_filtering)
    config_file.write("min_diff_hsr %s\n" % cmd_args.min_diff_hsr)

    # INSurVeyor specific
    config_file.write("max_trans_size %d\n" % cmd_args.max_trans_size)
    config_file.write("min_stable_mapq %d\n" % cmd_args.min_stable_mapq)

bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)

# TODO: remove contig_map use
with open("%s/contig_map" % cmd_args.workdir, "w") as contig_map:
    for k in bam_file.references:
        contig_map.write("%s\n" % (k))

reference_fa = pyfaidx.Fasta(cmd_args.reference)
rand_pos_gen = RandomPositionGenerator(reference_fa, cmd_args.seed, cmd_args.sampling_regions)
random_positions = []
n_rand_pos = int(rand_pos_gen.reference_len/1000)
for i in range(1,n_rand_pos):
    random_positions.append(rand_pos_gen.next())
print("%d random positions generated." % n_rand_pos)

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)

read_len = 0
count = 0
for read in bam_file.fetch(until_eof=True):
    if count > MAX_READS: break
    if read.is_secondary or read.is_supplementary: continue
    read_len = max(read_len, read.query_length)
    count += 1

general_dist = []
rnd_i = 0
while rnd_i < len(random_positions) and len(general_dist) < GEN_DIST_SIZE:
    chr, pos = random_positions[rnd_i]
    rnd_i += 1

    if pos > len(reference_fa[chr])-10000:
        continue

    i = 0
    for read in bam_file.fetch(contig=chr, start=pos, end=pos+10000):
        if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
            0 < read.template_length < MAX_ACCEPTABLE_IS:
            if i > 100: break
            i += 1
            general_dist.append(read.template_length)
reference_fa.close()

mean_is = np.mean(general_dist)
stddev_is = np.std(general_dist)

general_dist = [x for x in general_dist if abs(x-mean_is) < 5*stddev_is]

mean_is = int(np.mean(general_dist))
lower_stddev_is = int(np.sqrt(np.mean([(mean_is-x)**2 for x in general_dist if x < mean_is])))
higher_stddev_is = int(np.sqrt(np.mean([(x-mean_is)**2 for x in general_dist if x > mean_is])))
min_is, max_is = mean_is-3*lower_stddev_is, mean_is+3.5*higher_stddev_is

with open(cmd_args.workdir + "/stats.txt", "w") as stats_file: 
    stats_file.write("read_len %d\n" % read_len)
    stats_file.write("min_is %d\n" % min_is)
    stats_file.write("avg_is %d\n" % mean_is)
    stats_file.write("max_is %d\n" % max_is)

survindel2_workdir = cmd_args.workdir + "/survindel2"
mkdir(survindel2_workdir)
insurveyor_workdir = cmd_args.workdir + "/insurveyor"
mkdir(insurveyor_workdir)

def exec(cmd):
    start_time = timeit.default_timer()
    print("Executing:", cmd)
    if os.system(cmd) != 0:
        print("Error executing:", cmd)
        exit(1)
    elapsed = timeit.default_timer() - start_time
    print(cmd, "was run in %.2f seconds" % elapsed)
    

SURVEYOR_PATH = os.path.dirname(os.path.realpath(__file__))

read_categorizer_cmd = SURVEYOR_PATH + "/bin/reads_categorizer %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference)
exec(read_categorizer_cmd)

clip_consensus_builder_cmd = SURVEYOR_PATH + "/bin/clip_consensus_builder %s %s" % (cmd_args.workdir, cmd_args.reference)
exec(clip_consensus_builder_cmd)

if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]

find_svs_from_sr_consensuses_cmd = SURVEYOR_PATH + "/bin/find_svs_from_sr_consensuses %s %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)
exec(find_svs_from_sr_consensuses_cmd)

def cp(src, dst):
    os.system("cp %s %s" % (src, dst))

cp(cmd_args.workdir + "/sr.vcf.gz", survindel2_workdir)

cp(cmd_args.workdir + "/config.txt", survindel2_workdir)
cp(cmd_args.workdir + "/config.txt", insurveyor_workdir)

def append(src, dst):
    os.system("cat %s >> %s" % (src, dst))

append(cmd_args.workdir + "/stats.txt", survindel2_workdir + "/config.txt")
append(cmd_args.workdir + "/stats.txt", insurveyor_workdir + "/config.txt")

cp(cmd_args.workdir + "/stats.txt", survindel2_workdir)
cp(cmd_args.workdir + "/stats.txt", insurveyor_workdir)

cp(cmd_args.workdir + "/contig_map", survindel2_workdir)
cp(cmd_args.workdir + "/contig_map", insurveyor_workdir)

cp(cmd_args.workdir + "/full_cmd.txt", survindel2_workdir)
cp(cmd_args.workdir + "/full_cmd.txt", insurveyor_workdir)

cp(cmd_args.workdir + "/crossing_isizes.txt", survindel2_workdir)
cp(cmd_args.workdir + "/crossing_isizes_count_geq_i.txt", survindel2_workdir)
cp(cmd_args.workdir + "/min_disc_pairs_by_size.txt", insurveyor_workdir)

mkdir(survindel2_workdir + "/workspace")
mkdir(insurveyor_workdir + "/workspace")

exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/long-pairs", survindel2_workdir + "/workspace/long-pairs"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/mateseqs", survindel2_workdir + "/workspace/mateseqs"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/sc_mateseqs", survindel2_workdir + "/workspace/sc_mateseqs"))

exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/sr_consensuses", survindel2_workdir + "/workspace/sr_consensuses"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/hsr_consensuses", survindel2_workdir + "/workspace/hsr_consensuses"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/sr_consensuses", insurveyor_workdir + "/workspace/sr_consensuses"))

exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/clipped", insurveyor_workdir + "/workspace/clipped"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/mateseqs", insurveyor_workdir + "/workspace/mateseqs"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/fwd-stable", insurveyor_workdir + "/workspace/R"))
exec("cp -r %s %s" % (cmd_args.workdir + "/workspace/rev-stable", insurveyor_workdir + "/workspace/L"))

add_sr_filtering_info_cmd = SURVEYOR_PATH + "/bin/survindel2_add_sr_filtering_info %s %s %s %s" % (cmd_args.bam_file, survindel2_workdir, cmd_args.reference, sample_name)
exec(add_sr_filtering_info_cmd)

normalise_cmd = SURVEYOR_PATH + "/bin/survindel2_normalise %s/sr.annotated.vcf.gz %s/sr.annotated.norm.vcf.gz %s" % (survindel2_workdir, survindel2_workdir, cmd_args.reference)
exec(normalise_cmd)

merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/survindel2_merge_identical_calls %s/sr.annotated.norm.vcf.gz %s/sr.norm.dedup.vcf.gz" % (survindel2_workdir, survindel2_workdir)
exec(merge_identical_calls_cmd)

dp_clusterer = SURVEYOR_PATH + "/bin/survindel2_dp_clusterer %s %s %s %s" % (cmd_args.bam_file, survindel2_workdir, cmd_args.reference, sample_name)
exec(dp_clusterer)

dc_remapper_cmd = SURVEYOR_PATH + "/bin/insurveyor_dc_remapper %s %s %s" % (insurveyor_workdir, cmd_args.reference, sample_name)
exec(dc_remapper_cmd)

add_filtering_info_cmd = SURVEYOR_PATH + "/bin/insurveyor_add_filtering_info %s %s %s" % (cmd_args.bam_file, insurveyor_workdir, cmd_args.reference)
exec(add_filtering_info_cmd)

filter_cmd = SURVEYOR_PATH + "/bin/insurveyor_filter %s %s 0.25" % (insurveyor_workdir, cmd_args.reference)
exec(filter_cmd)

concat_cmd = "bcftools concat -a %s/out.pass.vcf.gz %s/out.pass.vcf.gz -Oz -o %s/out.pass.vcf.gz" % (survindel2_workdir, insurveyor_workdir, cmd_args.workdir)
exec(concat_cmd)
