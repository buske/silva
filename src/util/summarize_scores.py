#!/usr/bin/env python

"""
Usage: $0 INPUT FILE SCORED...

FILE has variants, one per line.
INPUT is the standardized input file for FILE.
Each SCORED file has a score, one per variant in FILE.

Reports ranked, annoated lines.
"""


import sys
import signal

from collections import defaultdict
from numpy import array, column_stack, median, argsort


# The following two rank functions are used in place of scipy's rankdata,
# in order to remove scipy as a dependency.
# http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
def rank_simple(vector):
    return sorted(range(len(vector)), key=vector.__getitem__)

def rankdata(a):
    n = len(a)
    ivec=rank_simple(a)
    svec=[a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
                
            sumranks = 0
            dupcount = 0
            
    return newarray


args = sys.argv[1:]
if len(args) < 3:
    print __doc__
    sys.exit(1)

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
    
def read_scores(filename):
    scores = []
    with open(filename) as ifp:
        scores = [float(line) for line in ifp]

    # Rank the negated scores
    return array(scores)

matfile, file = args[:2]
scored = args[2:]

features, mat = read_mat(matfile)

ranks_list = []
for filename in scored:
    scores = read_scores(filename)
    ranks_list.append(rankdata(-scores))
    
ranks = column_stack(ranks_list)

#scores = ranks.shape[1]/(1/ranks).sum(axis=1)  # Harmonic mean rank
#scores = ranks.mean(axis=1)  # Mean rank
#Use score if there's only one, else use MRR
if len(scored) > 1:
    scores = (1/ranks).sum(axis=1)/ranks.shape[1]  # Mean reciprocal rank

order = argsort(-scores)
ranks = rankdata(-scores)

lines = []
with open(file) as ifp:
    for line in ifp:
        line = line.strip()
        if line.startswith('#'):
            print "#%s" % '\t'.join(["rank", "score", 
                                     "features", line.strip('#')])
        else:
            lines.append(line)

assert len(lines) == len(order)

# Treat SIGPIPE as Unix would expect
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
for i in order:
    score, rank, line = scores[i], ranks[i], lines[i]
    mat_row = mat[i]
    assert len(mat_row) == len(features)
    extreme_features = []
    for feature, value in zip(features, mat_row):
        # Ignore class values
        if feature == 'class': continue
        if value > 3.0:  # 3 sigma
            extreme_features.append(feature)

    extreme_features = ','.join(extreme_features) if extreme_features else '.'
    rank = '%.1f' % rank
    rank = rank[:-2] if rank.endswith('.0') else rank
    print "%s\t%.3f\t%s\t%s" % (rank, score, extreme_features, line)

