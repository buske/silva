#!/usr/bin/env python

"""
Usage: $0 SCOREDIR

SCOREDIR/N-REP.scored
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
    r = 0  # number of positives
    with open(file) as ifp:
        for line in ifp:
            line = line.split('#', 1)[0].strip()
            if not line: continue
            score, group = line.split()

            group = int(group)
            assert group == 0 or group == 1
            score = float(score)
            scores.append(-score)  # negative, to allow proper ranking
            r += (group == 1)

    assert r == 1
    ranks = rankdata(array(scores))
    return ranks[0]

dir = args[0]
ranks = defaultdict(list)
for filename in glob("%s/*.scored" % dir):
    group = os.path.basename(filename).split('.')[0].split('-')[0]
    rank = rankfile(filename)
    ranks[group].append(rank)
    #print "%s\t%s" % (rank, filename)

groups = list(sorted(ranks))
#print 'thresh\t', '\t'.join(sorted(ranks))

for thresh in [25]:
    print '%d' % thresh,
    for group in groups:
        vals = ranks[group]
        count = sum([bool(rank <= thresh) for rank in vals])
        print "\t%d" % count,
        
    print ""

if False:
    print 'MRR:',
    for group in groups:
        vals = ranks[group]
        print "\t%.3f" % (sum([1/rank for rank in vals])/len(vals)),
    print ""
