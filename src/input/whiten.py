#!/usr/bin/env python

"""
Print whitened (zero-mean, unit-stdev) FILE.mat, 
optionally based upon values in a control MAT file.
"""

# Author: Orion Buske
# Date:   04 January 2012
from __future__ import division, with_statement

import os
import sys

from numpy import array

assert os.getenv('SILVA_PATH') is not None, \
    "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/src/share'))
from silva import maybe_gzip_open


def read_examples(filename):
    lines = []
    header = []
    with maybe_gzip_open(filename) as ifp:
        header = ifp.readline().strip()
        assert header.startswith('#')
        for line in ifp:
            line = line.strip()
            if not line: continue

            assert not line.startswith('#')
            lines.append([float(val) for val in line.split()])

    data = array(lines, dtype=float)
    return header, data

def whiten_data(mat, control=None):
    data = mat
    if control is not None:
        data = control

    means = data.mean(axis=0)
    stds = data.std(axis=0)
    # Whiten mat according to mean/std in source array
    for col in xrange(mat.shape[1]):
        if stds[col] == 0:
            mat[:, col] = 0
        else:
            mat[:, col] = (mat[:, col] - means[col]) / stds[col]
    
def script(filename, control=None):
    header, mat = read_examples(filename)
    if control is not None:
        control_header, control = read_examples(control)
        assert header == control_header

    whiten_data(mat, control=control)

    print header
    for row in mat:
        print '\t'.join(['%.4f' % val for val in row])
        
def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] FILE.mat"
    description = __doc__.strip()

    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option('-c', '--control', metavar='MAT',
                      dest='control', default=None,
                      help="")
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
