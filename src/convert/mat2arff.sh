#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 INPUT > ARFF
EOF
    exit 1
}

if [[ $# -ne 1 ]]; then
    usage
fi

# Put class column last
function recol {
    awk -F"\t" '{OFS="\t"; for (i=2; i<=NF; i++) {printf "%s\t", $i} print $1}'
}

input="$1"
echo -e "@RELATION SILVA\n"

if [[ $(head -n 1 "$input" | cut -f 1 | tr -d "#") != "class" ]]; then
    echo "Error: expected first column to be class" >&2
    exit 1
fi

# Convert header line into header
#    | awk '{if ($0 ~ /\?$/) {print $0, "{0,1}"} else {print $0, "NUMERIC"}}' \
head -n 1 "$input" \
    | cut -f 2- \
    | tr -d "#" \
    | tr "\t" "\n" \
    | awk '{print $0, "NUMERIC"}' \
    | sed -e 's/^/@ATTRIBUTE /'
echo "@ATTRIBUTE class {0,1}"
    
echo -e "\n@DATA"
# Convert data
tail -n +2 "$input" \
    | recol \
    | tr "\t" ","
