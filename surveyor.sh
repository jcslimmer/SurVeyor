BAM=$1
WORKDIR=$2
REF=$3
THREADS=$4

if [ ! -d $WORKDIR ]; then
    mkdir -p $WORKDIR
fi
if [ ! -d $WORKDIR/survindel2 ]; then
    mkdir -p $WORKDIR/survindel2
fi
if [ ! -d $WORKDIR/insurveyor ]; then
    mkdir -p $WORKDIR/insurveyor
fi

python3 survindel2.py --threads $THREADS $BAM $WORKDIR/survindel2 $REF
python3 insurveyor.py --threads $THREADS $BAM $WORKDIR/insurveyor $REF

bcftools concat -a $WORKDIR/survindel2/out.pass.vcf.gz $WORKDIR/insurveyor/out.pass.vcf.gz -Oz -o $WORKDIR/out.pass.vcf.gz
