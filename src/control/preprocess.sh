#!/usr/bin/env bash

set -eu
set -o pipefail

export SYNORDER_PATH=${SYNORDER_PATH:-$(cd $(dirname $0)/..; pwd;)}

function usage {
    cat <<EOF
Usage: $0 TRUE.arff CONTROL.arff OUTBASE
EOF
    exit 1
}

if [[ $# -ne 3 ]]; then
    usage
fi

true=$1
control=$2
outbase=$3

if [[ ! -e $true || ! -e $control ]]; then
    echo "Error: could not find input files" >&2
    exit 1
fi

outdir=$(dirname $outbase)
mkdir -pv $outdir
tmpbase=$outdir/.$(basename $outbase)

suffix=merged.arff
tmp=$tmpbase.$suffix
out=$outbase.$suffix
if [[ -e $out ]]; then
    echo "Found: $out" >&2
else
    echo "Merging control files to: $out..." >&2
    java -Xmx1000m weka.core.Instances append $true $control > $tmp \
	&& mv $tmp $out
fi

in=$out
suffix=merged.std.arff
tmp=$tmpbase.$suffix
out=$outbase.$suffix
if [[ -e $out ]]; then
    echo "Found: $out" >&2
else
    echo "Standardizing to: $out..." >&2
    java -Xmx1800m weka.filters.unsupervised.attribute.Standardize \
        -i $in -o $tmp \
	&& mv $tmp $out
fi

in=$out
$SYNORDER_PATH/control/split_data.py $in $outbase

echo "$0: SUCCESS" >&2
