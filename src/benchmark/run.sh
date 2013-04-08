#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 INPUT CASE.groups OUTBASE

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


function split_split {
    local input="$1"
    local outbase="$2"
    echo "Creating split stratified train/test datasets..." >&2
    $pwd/split_data.py -N 50 --true-groups=$groups \
	$input ${outbase}_SPLIT
}
function split_loo {
    local input="$1"
    local outbase="$2"
    echo "Creating pseudo-LOO train/test datasets..." >&2
    $pwd/split_data.py -N 10 --infection --true-groups=$groups \
	$input ${outbase}_LOO
}
function split_fs {
    local input="$1"
    local outbase="$2"
    echo "Creating feature-subset train/test datasets..." >&2
    
    function fs_subset {
	local cols="$1"
	local name="$2"
	local type="${3:-split}"
	local fs_input="${outbase}_${name}.mat"
	cut --complement -f "$cols" $input > $fs_input
	if [[ $type == "loo" ]]; then
	    split_loo $fs_input ${outbase}_${name}_FS
	else
	    split_split $fs_input ${outbase}_${name}_FS
	fi
    }

    fs_subset 2 gerp
    fs_subset 3-4 "codon_usage"
    fs_subset 5-6,22-23 "sequence"
    fs_subset 7-10,18-21 "motif"
    fs_subset 11-17 "splice"
    fs_subset 24-27 "folding"
}


# split_loo $input $outbase
split_split $input $outbase
split_fs $input $outbase

for dir in ${outbase}_SPLIT; do # {LOO,SPLIT}; do
    echo "Running models to prioritize variants..." >&2
    if [[ -d $dir ]]; then
	$pwd/run_models.sh ${dir}{,/scored}
    fi
done

for dir in ${outbase}*_FS_*; do
    echo "Running models to prioritize variants..." >&2
    if [[ -d $dir ]]; then
	$pwd/run_models.sh ${dir}{,/scored} forest
    fi
done
