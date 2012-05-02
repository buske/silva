#!/usr/bin/env python

"""
Generates stratified train/test datasets from the given
ARFF file, with positive examples in the same gene never
trained and tested together. Class assumed to be last column.
"""

# Author: Orion Buske
# Date:   04 January 2012
from __future__ import division, with_statement

import os
import sys

from numpy import append, array, concatenate
from numpy.random import shuffle

from buske import maybe_gzip_open

def mkdir(dir):
    assert not os.path.isfile(dir)
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except OSError:
            pass

def read_arff(filename):
    header = []
    true_egs = []
    false_egs = []
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('@'):
                header.append(line)
            else:
                if line.endswith(',1'):
                    true_egs.append(line)
                elif line.endswith(',0'):
                    false_egs.append(line)
                else:
                    assert False, "Invalid line: " + line

    return header, true_egs, false_egs

def read_groups(filename):
    with maybe_gzip_open(filename) as ifp:
        return [line.strip() for line in ifp
                if not line.startswith('#')]

def write_examples(dir, filename, examples_list, header=[]):
    out_filename = os.path.join(dir, filename)
    if not os.path.exists(out_filename):
        with open(out_filename, 'w') as ofp:
            if header:
                print >>ofp, '\n'.join(header)

            for examples in examples_list:
                print >>ofp, '\n'.join(examples)

def make_infection_datasets(true_egs, false_egs, outdir,
                            true_groups=None, iterations=10, header=None,
                            test_all_negs=False):
    print >>sys.stderr, "Generating dataset:"
    n_false_train = int(0.5 * len(false_egs)) 
    for i in xrange(len(true_egs)):
        print >>sys.stderr, i

        held_out = true_egs[i:i+1]
        
        if true_groups:
            group = true_groups[i]
            true_train = [e for j, e in enumerate(true_egs)
                             if j != i and true_groups[j] != group]
        else:
            true_train = true_egs[:i] + true_egs[i+1:]

        shuffle(true_train)
        shuffle(false_egs)
        for train_frac in [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
            n_true_train = int(train_frac * len(true_egs))
            #n_false_train = int(train_frac * 0.5 * len(false_egs))

            true_train_sub = true_train[:n_true_train]
            false_train = false_egs[:n_false_train]
            false_test = false_egs[n_false_train:]

            write_examples(outdir, "%02d-%02d.train.arff" % (n_true_train, i),
                           [true_train_sub, false_train], header=header)
            write_examples(outdir, "%02d-%02d.test.arff" % (n_true_train, i),
                           [held_out, false_test], header=header)

def script(filename, outdir, group_file=None, infection=False,
           iterations=10, train_frac=0.5, test_all_negs=False, **kwargs):
    assert 0 <= train_frac <= 1
    
    mkdir(outdir)
    header, true_egs, false_egs = read_arff(filename)
    true_groups = read_groups(group_file) if group_file else None

    make_infection_datasets(true_egs, false_egs, outdir, iterations=iterations,
                            true_groups=true_groups, header=header,
                            test_all_negs=test_all_negs)
        
def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] ARFF OUTDIR"
    description = __doc__.strip()

    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option('-G', '--true-groups', metavar='FILE',
                      dest='group_file', default=None,
                      help="Expects file with a line for each true example."
                      " Examples with the same group will not be trained and"
                      " tested together in --infection datasets.")
    parser.add_option('-i', '--infection',
                      dest='infection', default=False, action='store_true',
                      help="Perform simulated infection experiments"
                      " (modified LOO), rather than splits.")
    parser.add_option('-t', '--train', metavar='FLOAT', type='float',
                      dest='train_frac', default=0.5,
                      help="Instead of 50/50 split, train on FLOAT fraction"
                      " and test on the rest (default: %default)")
    parser.add_option('-N', '--iterations', metavar='N', type='int',
                      dest='iterations', default=10, 
                      help="Number of repetitions of each experiment.")
    parser.add_option('--test-all-negatives',
                      dest='test_all_negs', default=False, action='store_true',
                      help="Include all negative examples, in order, in the"
                      " test data for an --infection experiment.")
    options, args = parser.parse_args()

    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
