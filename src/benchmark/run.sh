#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 INPUT.mat TRUE.groups OUTBASE

Creates: 
         OUTBASE_SPLIT/, OUTBASE_LOO/

Existing files are preserved.
EOF
    exit 1
}

if [[ $# -ne 3 ]]; then
    usage
fi
input=$1
groups=$2
outbase=$3
pwd=$(cd $(dirname $0); pwd)

if [[ ! -e $input || ! -e $groups ]]; then
    echo -e "Error: missing input files\n" >&2
    usage
fi

mkdir -pv $(dirname "$outbase")


function split_data {
    local input="$1"
    local outbase="$2"
    echo "Creating split stratified train/test datasets..." >&2
    $pwd/split_data.py -N 50 $input ${outbase}_SPLIT
}
function split_loo {
    local input="$1"
    local outbase="$2"
    echo "Creating pseudo-LOO train/test datasets..." >&2
    $pwd/split_data.py -N 5 --infection --true-groups=$groups \
	$input ${outbase}_LOO
}

split_split $input $outbase
split_loo $input $outbase


for dir in ${outbase}{,_FS}_{LOO,SPLIT}; do
    echo "Running models to prioritize variants..." >&2
    if [[ -d $dir ]]; then
	$pwd/run_models.sh ${dir}{,/scored} forest
    fi
done


# Feature selection (via ARFF files and weka)
# $SILVA_PATH/src/convert/mat2arff.sh $input > $outbase.arff
# java -Xmx3800m weka.filters.supervised.attribute.AttributeSelection \
#     -E "weka.attributeSelection.ChiSquaredAttributeEval" \
#     -S "weka.attributeSelection.Ranker -T 1" \
#     < $outbase.arff > $outbase.fs.arff
# $SILVA_PATH/src/convert/arff2mat.sh < $outbase.fs.arff > $outbase.fs.input

# split_data $outbase.fs.input ${outbase}_FS
