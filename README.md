# SurVeyor
An SV caller from paired-end NGS data.

## Installation

**Please download the source code from the latest release, since an appropriate trained model is provided.**
**We discourage cloning the repository, as it may not be compatible with the latest trained model.**

In order to compile the code, the following are required:
- A C and a C++ compiler are required. If the GCC suite is used, version 4.9.3 or above is required.
- CMake (2.8 or above)

Downloaded the latest release from https://github.com/Mesh89/SurVeyor/releases , uncompress it and enter the directory.
The following commands should be sufficient

```
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

If you are compiling on the same platform as where you will execute it, you can use -DNATIVE=ON to create faster executables
```
cmake -DCMAKE_BUILD_TYPE=Release -DNATIVE=ON . && make
```

Python 3 is necessary to run SurVeyor. Libraries NumPy (http://www.numpy.org/), PySam (https://github.com/pysam-developers/pysam) and xgboost (https://xgboost.readthedocs.io/en/stable/) are also necessary.

Please also download trained-model.zip from the release you are using, and uncompress it in a location of your convenience.

## Data preparation

SurVeyor needs a BAM/CRAM file, a (possibly empty) working directory and reference genome in FASTA format. For the genotype command, it will also require an input VCF.

The BAM/CRAM file must be coordinate-sorted and indexed. Furthermore, the MD, MC and the MQ tag must be present for all primary alignments, when applicable.

Recent versions of BWA MEM will automatically add all the necessary tags. If you used an old version, a different aligner, or if for any reason your BAM/CRAM does not have the required tags, you can use the samtools command fixmate (for MC and MQ) and callmd (for MD) to add them.

## Denovo calling of structural variants

The basic command to call structural variants 
```
python surveyor.py call --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA --ml-model TRAINED_MODEL
```
TRAINED_MODEL should be the location of the uncompressed trained-model.zip

For other parameters, please see the help with
```
python surveyor.py call -h
```

## Genotyping a set of structural variants

The basic command to genotype a set of structural variants in VCF format (IN_VCF) is
```
python surveyor.py genotype --threads N_THREADS IN_VCF BAM_FILE WORKDIR REFERENCE_FASTA TRAINED_MODEL
```
The genotyped variants will be reported in genotyped.vcf.gz (version with all of the original calls) and genotyped.deduped.vcf.gz (version with duplicate calls removed).

For other parameters, please see the help with
```
python surveyor.py genotype -h
```

## Citation

Currently, there is no manuscript for SurVeyor.

If you use deletions or tandem duplications, please cite: 
Rajaby, R., Sung, WK. SurVIndel2: improving copy number variant calling from next-generation sequencing using hidden split reads. Nat Commun 15, 10473 (2024). https://doi.org/10.1038/s41467-024-53087-7

If you use insertions, please cite: Rajaby, R., Liu, DX., Au, C.H. et al. 
INSurVeyor: improving insertion calling from short read sequencing data. Nat Commun 14, 3243 (2023). https://doi.org/10.1038/s41467-023-38870-2
