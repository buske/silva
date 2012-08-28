#!/usr/bin/env python

"""
Generates stratified train/test datasets from the given
ARFF file, with positive examples in the same group never
trained and tested together. Class assumed to be last column.
"""

# Author: Orion Buske
# Date:   04 January 2012
from __future__ import division, with_statement

import os
import sys

from numpy.random import shuffle

assert os.getenv('SILVA_PATH') is not None, \
    "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/src/share'))
from silva import maybe_gzip_open

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

    last_attributes = [line.lower() for line in header if line.startswith('@')][-2:]
    assert last_attributes[-1] == '@data'
    assert last_attributes[-2].split()[1] == 'class', \
        "Class does not appear to be last attribute"

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


def make_split_datasets(true_egs, false_egs, outdir, true_groups=None,
                        header=None, iterations=50, test_frac=0.5):
    # Test set is decided first
    # Training set is the remaining elements that aren't in the same groups
    print >>sys.stderr, "Generating dataset:"
    n_true_test = int(test_frac * len(true_egs))
    n_false_test = int(test_frac * len(false_egs))
    true_is = range(len(true_egs))  # indices of true examples
    for i in xrange(iterations):
        print >>sys.stderr, i

        # Shuffle indices instead of examples so correspondence 
        # with groups is maintained
        shuffle(true_is)
        shuffle(false_egs)

        false_test = false_egs[:n_false_test]
        false_train = false_egs[n_false_test:]

        true_test_is = true_is[:n_true_test]
        true_test = [true_egs[j] for j in true_test_is]
        # true_train is all other true examples (not in same group)
        if true_groups:
            groups = set([true_groups[j] for j in true_test_is])
            true_train = [e for j, e in enumerate(true_egs)
                          if true_groups[j] not in groups]
        else:
            true_train = [true_egs[j] for j in true_is[n_true_test:]]

        write_examples(outdir, "%02d.train.arff" % i,
                       [true_train, false_train], header=header)
        write_examples(outdir, "%02d.test.arff" % i,
                       [true_test, false_test], header=header)

def make_infection_datasets(true_egs, false_egs, outdir, true_groups=None,
                            header=None, test_frac=0.5, iterations=1):
    print >>sys.stderr, "Generating dataset:"
    n_false_test = int(test_frac * len(false_egs))
    for i in xrange(len(true_egs)):
        print >>sys.stderr, i

        true_test = true_egs[i:i+1]  # held out example

        # true_train is all other true examples (not in same group)
        if true_groups:
            group = true_groups[i]
            true_train = [e for j, e in enumerate(true_egs)
                          if j != i and true_groups[j] != group]
        else:
            true_train = true_egs[:i] + true_egs[i+1:]

        for j in xrange(iterations):
            shuffle(false_egs)

            false_test = false_egs[:n_false_test]
            false_train = false_egs[n_false_test:]

            write_examples(outdir, "%02d-%02d.train.arff" % (i, j),
                           [true_train, false_train], header=header)
            write_examples(outdir, "%02d-%02d.test.arff" % (i, j),
                           [true_test, false_test], header=header)

def script(filename, outdir, group_file=None, infection=False,
           iterations=1, test_frac=0.5, **kwargs):
    assert 0 <= test_frac <= 1
    
    mkdir(outdir)
    header, true_egs, false_egs = read_arff(filename)

    true_groups = None
    if group_file:
        true_groups = read_groups(group_file)
        assert len(true_groups) == len(true_egs), \
            ("Found %d groups but %d true examples" % 
             (len(true_groups), len(true_egs)))

    if infection:
        make_infection_datasets(true_egs, false_egs, outdir, header=header,
                                test_frac=test_frac, true_groups=true_groups,
                                iterations=iterations)
    else:
        make_split_datasets(true_egs, false_egs, outdir, header=header,
                            true_groups=true_groups,
                            iterations=iterations, test_frac=test_frac)
        
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
    parser.add_option('-t', '--test', metavar='FLOAT', type='float',
                      dest='test_frac', default=0.5,
                      help="Instead of 50/50 split, test on FLOAT fraction"
                      " and train on the rest that are in different groups"
                      " (default: %default)")
    parser.add_option('-N', '--iterations', metavar='N', type='int',
                      dest='iterations', default=1,
                      help="Number of repetitions of each split experiment"
                      " (default: %default).")
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
