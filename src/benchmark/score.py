#!/usr/bin/env python

"""
Usage: $0 SCOREDIR

SCOREDIR/N[-REP].scored

Prints summary of rank data to stdout
"""


import sys
import os

from glob import glob
from collections import defaultdict
from numpy import array, concatenate
from scipy.stats.mstats import rankdata

args = sys.argv[1:]
if len(args) != 1:
    print __doc__
    sys.exit(1)
    
def read_scores(filename):
    pos = []
    neg = []
    with open(filename) as ifp:
        for line in ifp:
            line = line.split('#', 1)[0].strip()
            if not line: continue
            score, group = line.split()

            group = int(group)
            assert group == 0 or group == 1
            score = float(score)

            if group == 0:
                neg.append(score)
            else:
                pos.append(score)

    return array(pos), array(neg)

dir_ = args[0]
cum_scores = []
is_loo = None
for filename in glob("%s/*.scored" % dir_):
    pos, neg = read_scores(filename)
    scores = concatenate((pos, neg))
    ranks = rankdata(-scores)
    if len(pos) == 1:
        assert is_loo is None or is_loo
        is_loo = True
        # Compute LOO thresholds
        cum_scores.append(ranks[0])
    else:
        assert is_loo is None or not is_loo
        is_loo = False
        # Compute 50/50 r-precision
        r_precision = (ranks <= len(pos))[:len(pos)].sum()
        cum_scores.append(r_precision)

if is_loo:
    # Determine how many replicates were involved
    filenames = glob("%s/*.scored" % dir_)
    try:
        major = int(os.path.basename(filenames[0]).split('-')[0])
        reps = len(glob("%s/%d-*.scored" % (dir_, major)))
    except ValueError:
        reps = 1

    print >>sys.stderr, "Found %d reps per iter..." % reps
    ranks = array(cum_scores)
    print '#cutoff\tn\trecall'
    for cutoff in range(1, int(len(ranks) / reps) + 1):
        count = sum(ranks <= cutoff) / reps
        print '%d\t%.1f' % (cutoff, count)
else:
    print '#r_precision'
    for val in cum_scores:
        print '%d' % val
