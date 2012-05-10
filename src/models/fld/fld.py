#!/usr/bin/env python

"""
Run Fisher's Linear Discriminant analysis on the given vector file
(or '-' to read from stdin)

Positive/negative training examples start with 1/0.
Positive/negative test examples start with 2/-1.
"""

# Author:  Orion Buske
# Date:    11 December 2010
# Revised: 11 April 2011
from __future__ import division, with_statement

import os
import sys
import cPickle

from numpy import loadtxt, cov, dot
from numpy.linalg import inv

sys.path.insert(0, os.path.expandvars("$SYMPRI_PATH/src/share"))
from sympri import maybe_gzip_open


__version__ = 4

TRUE = 1
FALSE = 0

def load_data(vector_filename, log=sys.stderr):
    print >>log, "Loading vector data from file: %s" % vector_filename
    with maybe_gzip_open(vector_filename) as ifp:
        data = loadtxt(ifp, dtype=float)

    # Pop solution column
    solutions = data[:, 0]
    data = data[:,1:]
    print >>log, "Loaded data with %d examples and %d features" % data.shape

    return solutions, data

class FLD(object):
    def __init__(self, filename=None, train=True):
        self.filename = filename
        self.W = None
        if not train and filename and not os.path.isfile(filename):
            raise IOError("Could not find model file: %s" % filename)
        else:
            if self.filename and os.path.isfile(self.filename):
                try:
                    self.load()
                except Exception, e:
                    print >>sys.stderr, "Loading failed due to Exception: %s" % e

    def is_trained(self):
        return self.W is not None

    def save(self):
        print >>sys.stderr, "Saving FLD to file: %s" % self.filename
        with open(self.filename, 'wb') as ofp:
            cPickle.dump(self.W, ofp, 2)

    def load(self):
        print >>sys.stderr, "Loading FLD from file: %s" % self.filename
        with open(self.filename, 'rb') as ifp:
            self.W = cPickle.load(ifp)
            print >>sys.stderr, "W: %s" % (["%.4f" % w for w in self.W])

    def train(self, solutions, data, log=sys.stderr, out=sys.stderr):
        def calc_W(S, D):
            # Compute feature covariance matrix
            pos = D[S == TRUE]
            neg = D[S == FALSE]
            mean_pos = pos.mean(axis=0)
            mean_neg = neg.mean(axis=0)
            Sw = cov(pos, rowvar=0) + cov(neg, rowvar=0)
            W = dot(inv(Sw), (mean_pos - mean_neg).transpose())
            return W
            
        self.W = calc_W(solutions, data)
        
    def rank(self, solutions, data, out=sys.stdout):
        assert self.W is not None
        
        scores = dot(data, self.W)
        if out:
            for score, truth in zip(scores, solutions):
                print >>out, "%7.4f\t%d" % (score, truth)

def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] TRAIN [TEST]"
    description = __doc__.strip()
    version = __version__
    
    parser = OptionParser(usage=usage, version=version,
                          description=description)
    parser.add_option("--train", dest="train",
                      default=False, action="store_true",
                      help="Train model on VECTOR file,"
                      " overwriting any existing saved model")
    parser.add_option("--model", dest="model", metavar="FILENAME",
                      default=None,
                      help="Load/Save FLD model to this file")
    options, args = parser.parse_args()

    if len(args) != 1 and len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)

    trainfile = args[0]
    testfile = args[1] if len(args) >= 2 else None

    train_solns, train_data = load_data(trainfile)

    fld = FLD(train=options.train, filename=options.model)
    
    if options.train or not fld.is_trained():
        fld.train(train_solns, train_data)
        if options.model:
            fld.save()
        
    # Rank test data
    if testfile:
        test_solns, test_data = load_data(testfile)
        fld.rank(test_solns, test_data)
        
if __name__ == '__main__':
    sys.exit(main())
