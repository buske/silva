#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 CONTROL.vcf

Creates CONTROL.mat and CONTROL.input files for training controls.
EOF
    exit 1
}

if [[ $# -ne 1 ]]; then
    usage
fi

vcf=$1
base=$(dirname $1)/$(basename $1 .vcf)

$SILVA_PATH/silva-preprocess $(dirname $vcf) $vcf || true
echo "Overwriting standardization: $base.input" >&2
./src/input/standardize.py --class=0 $base.mat \
    > $base.input

