#!/usr/bin/env python

"""
Usage: $0 CASE CONTROL

CASE, CONTROL are directories that SilVA has been run on.
Compares scores of CASE variants against all the CONTROL variants.
Expects .flt and .scored files in each directory.
"""


import sys
import signal

from numpy import array, concatenate
from glob import glob

args = sys.argv[1:]
if len(args) != 2:
    print __doc__
    sys.exit(1)


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

def read_scores(filename):
    scores = []
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if line.startswith('#'):
                pass
            elif line.startswith('na'):
                scores.append(float('nan'))
            else:
                scores.append(float(line))

    return array(scores)

def calc_score(case_score, control_scores):
    combined = concatenate([array([case_score]), control_scores])
    # Invert data so highest is ranked first
    ranks = rankdata(-combined)
    return ranks[0], case_score


# Parse arguments
def find_files(d):
    flt = glob("%s/*.flt" % d)
    scored = glob("%s/*.scored" % d)
    assert len(flt) == len(scored) == 1
    return flt[0], scored[0]

case_dir, control_dir = args

case, case_scored = find_files(case_dir)
control, control_scored = find_files(control_dir)

case_scores = read_scores(case_scored)
control_scores = read_scores(control_scored)

scores = []
for row in case_scores:
    scores.append(calc_score(row, control_scores))

lines = []
with open(case) as ifp:
    for line in ifp:
        line = line.strip()
        if line.startswith('#'):
            print "#%s" % '\t'.join(["rank", "score", line.strip('#')])
        else:
            lines.append(line)

assert len(lines) == len(scores)

# Treat SIGPIPE as Unix would expect
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
for i in range(len(scores)):
    rank, score = scores[i]
    print "%.1f\t%.3f\t%s" % (rank, score, lines[i])

