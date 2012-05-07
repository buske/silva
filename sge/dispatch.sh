#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Please run this script from the tool's root folder.

Usage: $0 OUTROOT VCF...

Dispatches ./preprocess_vcf and ./run jobs to SGE cluster.
VCF files should have distinct names.
Creates subdirectory in OUTROOT for each VCF file.
EOF
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
elif [[ ! -e `pwd`/preprocess_vcf ]]; then
    echo -e "Could not find ./preprocess_vcf.\n" >&2
    usage
fi

outroot="$(readlink -m "$1")"
mkdir -pv "$outroot"
shift 1

for vcf in "$@"; do
    vcf="$(readlink -e "$vcf")"
    outbase="$(basename "$vcf" .vcf)"
    outdir="$outroot/$outbase"
    mkdir -pv "$outdir"

    id=$(uuidgen)
    qsub -b y -N "pre_$id" -l h_vmem=6G -l num_proc=4 \
	-q lunchQ -e ~/logs/ -o ~/logs/ \
	-S /bin/bash -v "PATH=$PATH" \
	"$(pwd)/preprocess_vcf $vcf $outdir"

    qsub -b y -N "run_$id" -hold_jid "pre_$id" -l h_vmem=2G \
	-q coffeeQ -e ~/logs/ -o ~/logs/ \
	-S /bin/bash -v "PATH=$PATH" \
	"$(pwd)/run $outdir > $outdir/RESULTS"
done