#!/usr/bin/env python

"""
Generates stratified train/test datasets from the given
ARFF files, with positive examples in the same gene never
trained and tested together. Class assumed to be last column.
"""

# Author: Orion Buske
# Date:   04 January 2012
from __future__ import division, with_statement

import os
import sys

from numpy import append, array, concatenate
from numpy.random import shuffle
from itertools import combinations

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
    if os.path.exists(out_filename):
        print >>sys.stderr, "Output file already exists: %s" % out_filename
    else:
        with open(out_filename, 'w') as ofp:
            if header:
                print >>ofp, '\n'.join(header)

            for examples in examples_list:
                print >>ofp, '\n'.join(examples)

def overlapping_subsets(n, d, *args):
    """Return list of overlapping subsets from n chunks of size N/d elts from each arg"""
    for offsets in combinations(range(d), n):
        yield [[elt for i, elt in enumerate(arg) if int(i % d) in offsets]
               for arg in args]
            
def make_datasets(true_control, false_control, outdir,
                  true_groups=None, header=None, n=2, d=4):
    """Make overlapping training datasets from n chunks of size N/d examples"""
    
    print >>sys.stderr, "Generating dataset:"
    shuffle(true_control)
    shuffle(false_control)
    for i, (true_train, false_train) in \
            enumerate(overlapping_subsets(n, d, true_control, false_control)):
        print >>sys.stderr, i
        write_examples(outdir, "%02d.train.arff" % i,
                       [true_train, false_control], header=header)

def script(controlfile, outdir, group_file=None, **kwargs):
    mkdir(outdir)
    header, true_control, false_control = read_arff(controlfile)
    
    true_groups = read_groups(group_file) if group_file else None
    make_datasets(true_control, false_control, outdir,
                  true_groups=true_groups, header=header)
        
def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] CONTROL.ARFF OUTDIR"
    description = __doc__.strip()

    parser = OptionParser(usage=usage,
                          description=description)
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
