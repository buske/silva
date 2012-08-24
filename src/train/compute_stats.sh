#!/usr/bin/env bash

set -eu
set -o pipefail

echo "SILVA_CONTROL: $SILVA_CONTROL"
read -p "Press ENTER to continue..."

TMPDIR=${TMPDIR:-.}

if [[ ! -e $SILVA_CONTROL.stats ]]; then
    cat $SILVA_CONTROL.merged.arff \
	| java weka.filters.unsupervised.attribute.Standardize \
	| egrep '(^@|,0$)' \
	| $(dirname $0)/compute_stats.py \
	> $SILVA_CONTROL.stats.tmp \
	&& mv $SILVA_CONTROL.stats{.tmp,}
fi

test -e $SILVA_CONTROL.stats

echo "$0: SUCCESS" >&2