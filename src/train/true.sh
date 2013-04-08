#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 TRUE.pcoord controlbase OUTDIR

controlbase.mat and controlbase.input must exist

Creates TRUE.*, OUTDIR/*.model, and OUTDIR/*.input
EOF
    exit 1
}

if [[ $# -ne 3 ]]; then
    usage
fi

pcoord=$1
controlbase=$2
outdir=$3
basename=$(basename $pcoord .pcoord)
base=$(dirname $pcoord)/$basename


flt=$base.flt
mrna=$base.mrna
mat=$base.mat

if [[ ! -s $flt ]]; then
    echo "Filtering and annotating variants: $flt" >&2
    ./src/input/synonymous.py --protein-coords -O data/refGene.pkl filter $pcoord > $flt
fi

if [[ ! -s $mrna ]]; then
    echo "Creating mRNA annotations: $mrna" >&2
    ./src/input/synonymous.py -O data/refGene.pkl annotate $flt > $mrna
fi

if [[ ! -s $mat ]]; then
    echo "Creating feature matrix: $mat" >&2
    ./src/input/mrna2mat.sh $mrna $base
fi

mkdir -pv $outdir
norm=$outdir/$basename.input
merged=$outdir/$basename.merged.input
groups=$outdir/$basename.groups
model=$outdir/0.model

echo "Standardizing against $(basename $controlbase.mat) to prepare input: ${norm}" >&2
./src/input/standardize.py \
    --class=1 --control=$controlbase.mat $base.mat > $norm

cp $norm $merged
grep -v "^#" $controlbase.input >> $merged
echo "Creating unified training file: $merged" >&2

echo "Creating groups file: $groups" >&2
sed -e '1d' $flt | cut -f 6 > $groups

echo "Training model: $model" >&2
./src/models/forest/train $model $merged

echo "$0: successfully completed" >&2

# Lines below are useful for variance estimation:
#   One can train models on subsets of the data
#
#./src/train/split_data.py  $merged $outdir
#./src/train/train_models.sh $outdir{,/models} forest
