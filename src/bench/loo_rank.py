#!/usr/bin/env python

"""
Usage: $0 SCOREDIR

SCOREDIR/N.scored
Prints summary of rank data to stdout
"""


import sys
import os

from glob import glob
from collections import defaultdict
from numpy import array
from scipy.stats.mstats import rankdata

args = sys.argv[1:]
if len(args) != 1:
    print __doc__
    sys.exit(1)
    
def rankfile(file):
    scores = []
    n_pos = 0
    with open(file) as ifp:
        for line in ifp:
            line = line.split('#', 1)[0].strip()
            if not line: continue
            score, group = line.split()

            group = int(group)
            assert group == 0 or group == 1
            score = float(score)
            scores.append(-score)  # negative, to allow proper ranking
            n_pos += (group == 1)

    assert n_pos == 1
    ranks = rankdata(array(scores))
    return ranks[0]

dir = args[0]
ranks = []
for filename in glob("%s/*.scored" % dir):
    ranks.append(rankfile(filename))

print "#top\tn"
for thresh in [5, 10, 25]:
    count = sum([bool(rank <= thresh) for rank in ranks])
    print '%d\t%d' % (thresh, count)
