
set -eu
set -o pipefail

vcf=$1
base=$(dirname $1)/$(basename $1 .vcf)

flt=$base.flt
mrna=$base.mrna
mat=$base.mat

if [[ ! -s $flt ]]; then
    ./src/input/synonymous.py -O data/refGene.pkl filter $pcoord > $flt
fi

if [[ ! -s $mrna ]]; then
    ./src/input/synonymous.py -O data/refGene.pkl annotate $flt > $mrna
fi

if [[ ! -s $mat ]]; then
    ./src/convert/mrna2mat.sh $mrna $base
fi

if [[ ! -s $base.input ]]; then
    ./src/input/standardize.py --class=0 $base.mat \
	> $base.input
fi

# Update control symlinks
