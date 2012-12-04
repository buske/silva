#!/usr/bin/env bash

set -eu
set -o pipefail

#==============================
# Set the values below as needed, according to your SGE configuration. 
# For example:
#
#   custom_args="-q lunchQ -e $HOME/logs/ -o $HOME/logs/"
#

# Used to run the 'silva-preprocess' script, which uses 4 threads and about 6GB of memory
preprocess_args="-l h_vmem=6G -l num_proc=6"
# Used to run the 'run' script, which uses 1 thread and about 2GB of memory
run_args="-l h_vmem=2G"
# Any additional parameters to include with both runs
custom_args=""
#==============================

function usage {
    cat <<EOF
Please run this script from the tool's root folder.

Usage: $0 OUTROOT VCF...

Dispatches ./silva-preprocess and ./silva-run jobs to SGE cluster.
VCF files should have distinct names.
Creates subdirectory in OUTROOT for each VCF file.
EOF
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
elif [[ ! -e `pwd`/silva-preprocess ]]; then
    echo -e "Could not find ./silva-preprocess.\n" >&2
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
    qsub -b y -N "silva_pre_$id" \
	$preprocess_args \
	$custom_args \
	-S /bin/bash -v "PATH=$PATH" \
	"$(pwd)/silva-preprocess $outdir $vcf"

    qsub -b y -N "silva_run_$id" -hold_jid "pre_$id" \
	$run_args \
	$custom_args \
	-S /bin/bash -v "PATH=$PATH" \
	"$(pwd)/silva-run $outdir > \$TMPDIR/RESULTS && mv -v \$TMPDIR/RESULTS $outdir/RESULTS"
done