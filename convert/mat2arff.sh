#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 RELATION CLASS MAT OUTBASE

Saves ARFF to OUTBASE.arff
EOF
    exit 1
}

if [[ $# -ne 4 ]]; then
    usage
fi
class=$2
mat=$3
outbase=$4

(
    echo -e "@RELATION $1\n"

    head -n 1 "$mat" \
	| tr -d "#" \
	| tr "\t" "\n" \
	| awk '{if ($0 ~ /\?$/) {print $0, "{0,1}"} else {print $0, "NUMERIC"}}' \
	| sed -e 's/^/@ATTRIBUTE /'
    echo "@ATTRIBUTE class {0,1}"
    
    echo -e "\n@DATA"
    
    tail -n +2 "$mat" \
	| tr "\t" "," \
	| sed -e 's/$/,'"$3"'/'
) > $outbase.arff.temp

mv $outbase.arff.temp $outbase.arff