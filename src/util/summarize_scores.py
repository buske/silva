#!/usr/bin/env python

"""
Usage: $0 STATS MAT FILE SCORED...

STATS is the control stats file.
FILE has variants, one per line.
MAT is the matrix file for FILE.
Each SCORED file has a score, one per variant in FILE.

Reports rank mean and range for each line in FILE.
"""


import sys

from collections import defaultdict
from numpy import array, column_stack, median, argsort
from scipy.stats.mstats import rankdata

args = sys.argv[1:]
if len(args) < 4:
    print __doc__
    sys.exit(1)

def read_stats(filename):
    attrs = {}
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            attr, mean, stdev = line.split()
            attrs[attr] = (mean, stdev)
            
    return attrs

def read_mat(filename):
    mat = []
    features = []
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line: continue
            if line.startswith('#'):
                features = line.translate(None, '#').split('\t')
            else:
                values = [float(val) for val in line.split()]
                mat.append(values)
            
    return features, mat
    
def read_ranks(filename):
    scores = []
    with open(filename) as ifp:
        scores = [int(line) for line in ifp]

    # Rank the negated scores
    return rankdata(-array(scores))

statsfile, matfile, file = args[:3]
scored = args[3:]

stats = read_stats(statsfile)
features, mat = read_mat(matfile)

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
            print "#%s" % '\t'.join(["score", "low,median,high", "features", line.strip('#')])
        else:
            lines.append(line)

assert len(lines) == len(order)

for i in order:
    score, low, mid, high, line = scores[i], range_lows[i], range_mids[i], range_highs[i], lines[i]
    mat_row = mat[i]
    assert len(mat_row) == len(features)
    extreme_features = []
    for feature, value in zip(features, mat_row):
        mean, stdev = stats[feature]
        if value > mean + 1.5 * stdev:
            extreme_features.append(feature)

    extreme_features = ','.join(extreme_features) if extreme_features else '.'
    print "%.2f\t%d,%d,%d\t%s\t%s" % (score, low, mid, high, extreme_features, line)

