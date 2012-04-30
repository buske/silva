#!/usr/bin/env bash

set -eu
set -o pipefail

bin=$SYNORDER_PATH/input/synonymous.py
refgene=$SYNORDER_PATH/data/refGene.pkl

function usage {
    cat <<EOF
Usage: $0 FILE.avout [rare|common]

Creates FILE.mrna from rare (default) or common variants in FILE.avout.
EOF
    exit 1
}


if [[ $# -eq 2 ]]; then
    afsubset=$2
elif [[ $# -eq 1 ]]; then
    afsubset=rare
else
    usage
fi
input=$1
if [[ ! -e $input ]]; then
    echo "Missing input file: $input" >&2
fi
outbase=$(basename $input .avout)
outdir=$(dirname $input)

function synonymous {
    awk '$3 ~ /\<synonymous/'
}
function filter_rare {
    awk -F"\t" '$6 == "." || $6 < 0.01'
}
function filter_common {
    awk -F"\t" '$6 != "." && $6 < 0.9 && $6 > 0.1'
}

#function expand_txs {
#    awk -F"\t" '{OFS="\t"; split($4,txs,","); for (i in txs) {if (txs[i]) {$4 = txs[i]; print $0}}}'
#}

function expand_genes {
    awk -F"\t" '{OFS="\t"; split($2,genes,","); for (i in genes) {if (genes[i]) {split(genes[i],gene,";"); $2=gene[1]; print $0}}}'
}

function recolumn {
    awk -F"\t" '{OFS="\t"; split($4,x,":"); print $13, $14, $16, $17, $2;}'
}

function lookup_mrna {
    # Lookup mrna and filter to third-position variants
    $bin -O $refgene - \
	| awk -F"\t" '/^#/ || $9 == 2;'
}

id=$TMPDIR/$(uuidgen)

function cleanup {
    rm -f $id.*
}
function die {
    cleanup
    kill 0
    exit 1
}

trap die INT TERM EXIT

mrna=$outbase.mrna
if [[ -e $outdir/$mrna ]]; then
    echo "Already exists: $outdir/$mrna" >&2
else
    # Roughly 5GB each
    echo "Creating $mrna..." >&2
    cat $input | synonymous | filter_$afsubset | expand_genes \
	| recolumn | lookup_mrna > $id.$mrna 2> $outdir/$mrna.log \
	&& mv $id.$mrna $outdir/$mrna
fi

trap - INT TERM EXIT

cleanup
