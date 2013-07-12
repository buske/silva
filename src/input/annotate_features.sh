#!/usr/bin/env bash

set -eu
set -o pipefail

featuredir=$SILVA_PATH/src/features
datadir=$SILVA_PATH/data
MRNA_COL=11
CODON_COLS="4-5,8-10"

function usage {
    cat >&2 <<EOF
Usage: $0 MRNA OUTBASE

Creates OUTBASE.mat
EOF
    exit 1
}
if [[ $# -ne 2 ]]; then
    usage
fi
mrna="$1"
outdir="$(cd -P $(dirname "$2"); pwd)"
outbase="$outdir/$(basename "$2")"

temp="$TMPDIR/$(uuidgen)"

function cleanup {
    rm -f $temp.*
}
function die {
    echo "Error encountered!" >&2
    cleanup
    trap - INT TERM EXIT
    kill 0
    exit 1
}


trap die INT TERM EXIT

in=$temp.mrna-col
grep -v "^#" $mrna \
    | cut -f $MRNA_COL \
    > $in || exit 1

function run {
    local name=$1.col
    local dir=$2
    shift 2
    if [[ ! -e $outbase.$name ]]; then
	if [[ $dir != "." ]]; then
	    pushd $dir > /dev/null
	fi
	echo "Running: $name" >&2
	"$@" > $temp.$name \
	    && mv $temp.$name $outbase.$name \
	    && echo "Completed: $name" >&2
	if [[ $dir != "." ]]; then
	    popd > /dev/null
	fi
    else
	echo "Output file already exists: $outbase.$name" >&2
    fi
}

n_threads=0
n_cols=0
function setup {
    if [[ $n_threads -ge ${SILVA_N_THREADS:-1} ]]; then
	wait
	n_threads=0
    fi

    n_threads=$(expr $n_threads + 1);
    n_cols=$(expr $n_cols + 1); 
}

setup; run cpg     $featuredir/other   ./cpg.py     -q $in &
setup; run splice  $featuredir/other   ./splice.py  -q $in &
setup; run ese3    $featuredir/ese3    ./ese3.py    -q $in &
setup; run pesx    $featuredir/pesx    ./pesx.py    -q $in &
setup; run fas-ess $featuredir/fas-ess ./fas-ess.py -q $in &
setup; run maxent  $featuredir/maxent  ./maxent.py  -q $in &
setup; grep -v "^#" $mrna | cut -f $CODON_COLS \
    | run codon $featuredir/other ./codon_usage.py -q - &

setup; cut -f 1,2 $mrna \
    | run 0_gerp $featuredir/other ./gerp.py $datadir/gerp.refGene.table.gz $datadir/gerp.refGene.pkl &

if [[ -z "${EXCLUDE_RNA_FOLDING:-}" ]]; then
    setup; run unafold-50 $featuredir/unafold ./unafold.py -d 50 -q $in &
    setup; run vienna-50 $featuredir/vienna ./vienna.py -d 50 -q $in &
else
    echo "EXCLUDE_RNA_FOLDING=${EXCLUDE_RNA_FOLDING}... excluding RNA folding features." >&2
fi


# Wait for any others to finish
wait

found_n_cols=$(ls -1 $outbase.*.col | wc -l)
if [[ $found_n_cols -ne $n_cols ]]; then
    echo "Error: found only $found_n_cols / $n_cols column annotations." >&2
    exit 1
fi

paste $outbase.*.col \
    | perl -pe 's/\bna\b/0/g' \
    > $temp.mat \
    && mv $temp.mat ${outbase}.mat

cleanup
trap - INT TERM EXIT
