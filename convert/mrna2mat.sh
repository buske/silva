#!/usr/bin/env bash

set -eu
set -o pipefail

srcdir=$HOME/src/synonymous
tooldir=$HOME/src/tools/synonymous
datadir=/data/buske/synonymous
MRNA_COL=10
CODON_COLS="3-4,7-9"

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
out="$(readlink -f $2)"
outdir="$(dirname "$out")"

temp="/tmp/buske.$(uuidgen)"
touch $temp.touch || $temp="/data/buske/tmp/$(uuidgen)"
rm -f $temp.touch

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
    > $in

function run {
    local name=$1.col
    local dir=$2
    shift 2
    if [[ ! -e $out.$name ]]; then
	if [[ $dir != "." ]]; then
	    pushd $dir > /dev/null
	fi
	echo "Running $name" >&2
	"$@" > $temp.$name 2> ${out}.$name.log \
	    && mv $temp.$name $out.$name \
	    && echo "Completed $name" >&2
	if [[ $dir != "." ]]; then
	    popd > /dev/null
	fi
    else
	echo "Output file already exists: $out.$name" >&2
    fi
}

n_cols=0
n_cols=$(expr $n_cols + 1); run cpg     .                $srcdir/cpg.py -q $in &
n_cols=$(expr $n_cols + 1); run splice  .                $srcdir/splice.py -q $in &
n_cols=$(expr $n_cols + 1); run ese3    $tooldir/ese3    ./ese3.py      -q $in &
n_cols=$(expr $n_cols + 1); run fas-ess $tooldir/fas-ess ./fas-ess.py   -q $in &
n_cols=$(expr $n_cols + 1); run pesx    $tooldir/pesx    ./pesx.py      -q $in &
n_cols=$(expr $n_cols + 1); cat $mrna \
    | run 0_gerp . /home/buske/bin/genomedata-query -q $datadir/hg19.genomedata gerp &
wait

n_cols=$(expr $n_cols + 1); run maxent  $tooldir/maxent  ./maxent.py    -q $in &
n_cols=$(expr $n_cols + 1); run unafold-50 $tooldir/unafold ./unafold.py -d 50 -q $in &
n_cols=$(expr $n_cols + 1); run unafold-100 $tooldir/unafold ./unafold.py -d 100 -q $in &
#n_cols=$(expr $n_cols + 1); run unafold-500 $tooldir/unafold ./unafold.py -d 500 -q $in &

n_cols=$(expr $n_cols + 1); grep -v "^#" $mrna | cut -f $CODON_COLS \
    | run codon . $srcdir/codon_usage.py -q - &
wait

#    | perl -pe 's/\t(\t|$)/\tna\1/g' \
found_n_cols=$(ls -1 $out.*.col | wc -l)
if [[ $found_n_cols -ne $n_cols ]]; then
    echo "Error: found only $found_n_cols / $n_cols column annotations." >&2
    exit 1
fi

paste $out.*.col \
    | sed -e 's/na/0/g' \
    > $temp.mat \
    && mv $temp.mat ${out}.mat

trap - INT TERM EXIT

echo "$0: SUCCESS" >&2
