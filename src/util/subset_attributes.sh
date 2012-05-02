#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 IN_TRAIN.ARFF IN_TEST.ARFF OUT_TRAIN.ARFF OUT_TEST.ARFF
EOF
    exit 1
}

if [[ $# -ne 4 ]]; then
    usage
fi
itrain=$1
itest=$2
otrain=$3
otest=$4

if [[ ! -e $itrain || ! -e $itest ]]; then
    echo "Error: missing input files" >&2
    exit 1
elif [[ -e $otrain && -e $otest ]]; then
    echo "Output files already exist" >&2
    exit 0
fi

id=$(uuidgen)
function cleanup {
    rm -f "$id."*
}
function die {
    cleanup
    exit 1
}

function bump_gerp {
    local f=$1
    gerp_col=$(grep -i "^@attribute" "$f" | cat -n | grep -i "gerp" | awk '{print $1}' || true)
    if [[ $gerp_col ]]; then
	if [[ $gerp_col -eq 1 ]]; then
	    cols="$gerp_col,2-last"
	else
	    cols="$gerp_col,1-$((gerp_col-1)),$((gerp_col+1))-last"
	fi
	java -Xmx1800m weka.filters.unsupervised.attribute.Reorder -R "$cols" -i "$f"
    else
	cat "$f"
    fi
}

trap die INT TERM EXIT

traintmp=$id.train.arff
testtmp=$id.test.arff
java -Xmx3800m weka.filters.supervised.attribute.AttributeSelection -E "weka.attributeSelection.ChiSquaredAttributeEval" -S "weka.attributeSelection.Ranker -T 1" -b -i "$itrain" -o "$traintmp" -r "$itest" -s "$testtmp"

bump_gerp "$traintmp" > "$otrain"
bump_gerp "$testtmp" > "$otest"

trap - INT TERM EXIT
cleanup

#java -Xmx1g weka.filters.supervised.attribute.AttributeSelection -E "weka.attributeSelection.CfsSubsetEval" -S "weka.attributeSelection.BestFirst -D 2" -i "$arff"
