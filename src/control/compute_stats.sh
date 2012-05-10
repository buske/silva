#!/usr/bin/env bash

set -eu
set -o pipefail

SYMPRI_PATH=${SYMPRI_PATH:-$(pwd)}
TMPDIR=${TMPDIR:-.}
sample=NA10851

if [[ ! -e $SYMPRI_PATH/data/control/$sample.stats ]]; then
    cat $SYMPRI_PATH/data/control/$sample.merged.arff \
	| java weka.filters.unsupervised.attribute.Standardize \
	| egrep '(^@|,0$)' \
	| $(dirname $0)/compute_stats.py \
	> $TMPDIR/$sample.stats \
	&& mv $TMPDIR/$sample.stats $SYMPRI_PATH/data/control/$sample.stats
fi

test -e $SYMPRI_PATH/data/control/$sample.stats

echo "$0: SUCCESS" >&2