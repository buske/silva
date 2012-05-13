#!/usr/bin/env bash

set -eu
set -o pipefail

SILVA_PATH=${SILVA_PATH:-$(pwd)}
TMPDIR=${TMPDIR:-.}
sample=NA10851

if [[ ! -e $SILVA_PATH/data/control/$sample.stats ]]; then
    cat $SILVA_PATH/data/control/$sample.merged.arff \
	| java weka.filters.unsupervised.attribute.Standardize \
	| egrep '(^@|,0$)' \
	| $(dirname $0)/compute_stats.py \
	> $TMPDIR/$sample.stats \
	&& mv $TMPDIR/$sample.stats $SILVA_PATH/data/control/$sample.stats
fi

test -e $SILVA_PATH/data/control/$sample.stats

echo "$0: SUCCESS" >&2