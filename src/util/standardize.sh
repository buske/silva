#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 TRUE.arff FALSE.arff

Prints merged standardized data to stdout
EOF
    exit 1
}

if [[ $# -ne 2 ]]; then
    usage
fi
true=$1
false=$2
java -Xmx1800m weka.core.Instances append "$false" "$true" \
    | java -Xmx1800m weka.filters.unsupervised.attribute.Standardize


# id=$(uuidgen)
# function cleanup {
#     rm -f "$id"*
# }
# function die {
#     cleanup
#     exit 1
# }

# # Add ID's separately, merge, and standardize
# trap die INT TERM EXIT

# java -Xmx1g weka.filters.unsupervised.attribute.Standardize \
#     -b -i "$2" -o "$id.2.arff" -r "$1" -s "$id.1.arff"

# # Put temp2 first so relation is taken from false dataset, which will be varied
# java weka.core.Instances append "$id.2.arff" "$id.1.arff"

# trap - INT TERM EXIT

# cleanup