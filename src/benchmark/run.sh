#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 INPUT TRUE.GROUPS OUTBASE

Creates: 
         OUTBASE_5050/, OUTBASE_LOO/

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
    echo "Creating 50/50 stratified train/test datasets..." >&2
    $pwd/split_data.py -N 25 $input ${outbase}_5050
    
    echo "Creating pseudo-LOO train/test datasets..." >&2
    $pwd/split_data.py -N 5 --infection --true-groups=$groups \
	$input ${outbase}_LOO
}

split_data $input $outbase

# Feature selection (via ARFF files and weka)
$SILVA_PATH/src/convert/mat2arff.sh $input > $outbase.arff
java -Xmx3800m weka.filters.supervised.attribute.AttributeSelection \
    -E "weka.attributeSelection.ChiSquaredAttributeEval" \
    -S "weka.attributeSelection.Ranker -T 1" \
    < $outbase.arff > $outbase.fs.arff
$SILVA_PATH/src/convert/arff2mat.sh < $outbase.fs.arff > $outbase.fs.input

split_data $outbase.fs.input ${outbase}_FS

for dir in ${outbase}{,_FS}_{LOO,5050}; do
    echo "Running models to prioritize variants..." >&2
    $pwd/run_models.sh ${dir}{,/scored}
done
