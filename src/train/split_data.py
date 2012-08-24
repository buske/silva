#!/usr/bin/env python

"""
Generates stratified train/test datasets from the given
training data (MAT format). Class must be first column.
"""

# Author: Orion Buske
# Date:   04 January 2012
from __future__ import division, with_statement

import os
import sys

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

def read_mat(filename):
    true_egs = []
    false_egs = []
    header = ""
    with maybe_gzip_open(filename) as ifp:
        header = ifp.readline().strip()
        assert header.startswith('#')
        assert header.strip('#').split()[0] == 'class'

        for line in ifp:
            line = line.strip()
            cls = line.split(None, 1)[0]
            if cls == '1':
                true_egs.append(line)
            elif cls == '0':
                false_egs.append(line)
            else:
                assert False, "Invalid line: " + line

    return header, true_egs, false_egs

def write_examples(dir, filename, examples_list, header=None):
    out_filename = os.path.join(dir, filename)
    if os.path.exists(out_filename):
        print >>sys.stderr, "Output file already exists: %s" % out_filename
    else:
        with open(out_filename, 'w') as ofp:
            if header:
                print >>ofp, header

            for examples in examples_list:
                print >>ofp, '\n'.join(examples)

def overlapping_subsets(n, d, *args):
    """Return list of overlapping subsets from n chunks of size N/d elts from each arg"""
    for offsets in combinations(range(d), n):
        yield [[elt for i, elt in enumerate(arg) if int(i % d) in offsets]
               for arg in args]
            
def make_datasets(true_control, false_control, outdir,
                  header=None, n=2, d=3):
    """Make overlapping training datasets from n chunks of size N/d examples"""
    
    print >>sys.stderr, "Generating dataset:"
    shuffle(true_control)
    shuffle(false_control)
    for i, (true_train, false_train) in \
            enumerate(overlapping_subsets(n, d, true_control, false_control)):
        print >>sys.stderr, i
        write_examples(outdir, "%02d.train.input" % i,
                       [true_train, false_control], header=header)

def script(controlfile, outdir, **kwargs):
    mkdir(outdir)
    header, true_control, false_control = read_mat(controlfile)
    
    make_datasets(true_control, false_control, outdir, header=header)
        
def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] TRAINING.MAT OUTDIR"
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
