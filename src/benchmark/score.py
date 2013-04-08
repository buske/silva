#!/usr/bin/env python

"""
SCOREDIR/N[-REP].scored

Prints summary of rank data to stdout
"""


import sys
import os

from glob import glob
from collections import defaultdict
from numpy import array, concatenate
from scipy.stats.mstats import rankdata

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


def run_loo(root, print_all=False):
    cum_scores = defaultdict(list)
    for filename in glob("%s/*.scored" % root):
        id_ = os.path.basename(filename).split('.')[0]
        major, rep = id_.split('-')

        pos, neg = read_scores(filename)
        scores = concatenate((pos, neg))
        ranks = rankdata(-scores)
        assert len(pos) == 1
        cum_scores[rep].append(ranks[0])

    reps = len(cum_scores)
    assert len(set(map(len, cum_scores.values()))) == 1

    print >>sys.stderr, "Found %d reps per iter..." % reps
    ranks = array([cum_scores[rep] for rep in sorted(cum_scores)]).transpose()
    if print_all:
        print '#cutoff\trecall...'
        for cutoff in range(1, ranks.shape[0] + 1):
            counts = (ranks <= cutoff).sum(axis=0)
            print '%d\t%s' % (cutoff, '\t'.join(['%.1f' % count 
                                                 for count in counts]))
    else:
        print '#cutoff\trecall'
        for cutoff in range(1, ranks.shape[0] + 1):
            count = float((ranks <= cutoff).sum()) / reps
            print '%d\t%.1f' % (cutoff, count)


def run_split(root, print_all=False):
    cum_scores = []
    n_pos = None
    for filename in glob("%s/*.scored" % root):
        pos, neg = read_scores(filename)
        scores = concatenate((pos, neg))
        ranks = rankdata(-scores)
        if n_pos is None:
            n_pos = len(pos)
        else:
            assert n_pos == len(pos)

        pos_ranks = ranks[:n_pos]
        counts = [(pos_ranks <= i).sum() for i in range(1, len(ranks))]
        cum_scores.append(tuple(counts))

    cum_scores = array(cum_scores).transpose()
    if print_all:
        print '#cutoff\tcounts...'
        for i, row in enumerate(cum_scores):
            print '%d\t%s' % (i + 1, 
                              '\t'.join(['%d' % count for count in row]))
    else:
        means = cum_scores.mean(axis=1)
        print '#cutoff\trecall'
        for i, mean in enumerate(means):
            print '%d\t%.1f' % (i + 1, mean)


def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] SCOREDIR"
    description = __doc__.strip()
    
    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option("--all", default=False,
                      dest="print_all", action='store_true',
                      help="Print all data points instead of average"
                      " for SPLIT mode")

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    # Peek at first file to see if LOO or SPLIT
    for filename in glob("%s/*.scored" % args[0]):
        pos, neg = read_scores(filename)
        break

    if len(pos) == 1:
        run_loo(*args, **kwargs)
    else:
        run_split(*args, **kwargs)
        

if __name__ == '__main__':
    sys.exit(main())
