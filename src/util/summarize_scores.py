#!/usr/bin/env python

"""
Usage: $0 FILE SCORED...

Reports rank mean and range for each line in file
"""


import sys

from collections import defaultdict
from numpy import array, column_stack, median, argsort
from scipy.stats.mstats import rankdata

args = sys.argv[1:]
if len(args) < 2:
    print __doc__
    sys.exit(1)

def read_ranks(file):
    scores = []
    with open(file) as ifp:
        scores = [int(line) for line in ifp]

    # Rank the negated scores
    return rankdata(-array(scores))

file = args[0]
scored = args[1:]
ranks_list = []
for filename in scored:
    ranks = read_ranks(filename)
    ranks_list.append(ranks)
    
ranks = column_stack(ranks_list)

#scores = ranks.shape[1]/(1/ranks).sum(axis=1)  # Harmonic mean rank
#scores = ranks.mean(axis=1)  # Mean rank
scores = (1/ranks).sum(axis=1)/ranks.shape[1]  # Mean reciprocal rank
order = argsort(-scores)

range_lows = ranks.min(axis=1)
range_highs = ranks.max(axis=1)
range_mids = median(ranks, axis=1)
lines = []
with open(file) as ifp:
    for line in ifp:
        line = line.strip()
        if line.startswith('#'):
            print "#%s" % '\t'.join(["score", "low,median,high", line.strip('#')])
        else:
            lines.append(line)

assert len(lines) == len(order)

for i in order:
    score, low, mid, high, line = scores[i], range_lows[i], range_mids[i], range_highs[i], lines[i]
    print "%.2f\t%d,%d,%d\t%s" % (score, low, mid, high, line)

