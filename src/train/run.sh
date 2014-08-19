#!/usr/bin/env bash

set -eu
set -o pipefail

function usage() {
    cat <<EOF
Usage: $0 CONTROL.vcf CASE.pcoord OUTDIR

To be run from the root of SILVA_PATH, using complete file paths
EOF
    exit 1
}

if [[ $# -ne 3 ]]; then
    usage
fi

control=$1
case=$2
outdir=$3
tmpdir=/tmp/buske
mkdir -pv $tmpdir $outdir
log=$outdir/LOG

controlbase=$(dirname $control)/$(basename $control .vcf)

export SILVA_PATH=$(pwd)
export TMPDIR=$tmpdir

set -x
./src/train/control.sh $control 2>&1 | tee $log
./src/train/case.sh $case $controlbase $outdir 2>&1 | tee -a $log
set +x

echo "$0: SUCCESS" | tee -a $log >&2 

