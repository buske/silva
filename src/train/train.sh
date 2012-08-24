
set -eu
set -o pipefail

pcoord=$1
control=$2
base=$(dirname $1)/$(basename $1 .pcoord)

flt=$base.flt
mrna=$base.mrna
mat=$base.mat

if [[ ! -s $flt ]]; then
    ./src/input/synonymous.py --protein-coords -O data/refGene.pkl filter $pcoord > $flt
fi

if [[ ! -s $mrna ]]; then
    ./src/input/synonymous.py -O data/refGene.pkl annotate $flt > $mrna
fi

if [[ ! -s $mat ]]; then
    ./src/convert/mrna2mat.sh $mrna $base
fi

if [[ ! -s $base.input ]]; then
    ./src/input/standardize.py \
	--class=1 --control=$control $base.mat > $base.input
fi

cp $base.input $base.merged.input
grep -v "^#" $control >> $base.merged.input



# Update control symlinks
