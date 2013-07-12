#!/usr/bin/env python

"""
Print standardized (zero-mean, unit-stdev) FILE.mat, 
optionally based upon values in a control MAT file.
If the header includes a 'class' column, it will be 
left as is.
"""

# Author: Orion Buske
# Date:   04 January 2012
from __future__ import division, with_statement

import os
import sys

from numpy import array

assert os.getenv('SILVA_PATH') is not None, \
    "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/lib/python'))
from silva import maybe_gzip_open


def read_examples(filename):
    lines = []
    header = []
    ncols = None
    with maybe_gzip_open(filename) as ifp:
        header = ifp.readline().strip()
        assert header.startswith('#')
        header = header.replace('#', '')
        for line in ifp:
            line = line.strip()
            if not line: continue

            assert not line.startswith('#')
            tokens = [float(val) for val in line.split()]
            if ncols is None:
                ncols = len(tokens)
            else:
                assert ncols == len(tokens), \
                    "Found row in %s with %d columns (%d expected)" % (filename, len(tokens), ncols)

            lines.append(tokens)

    try:
        data = array(lines, dtype=float)
    except ValueError:
        print >>sys.stderr

    # Find class column
    cols = header.split()
    if cols[0] == 'class':
        class_col = 0
    else:
        class_col = None
        
    return header, class_col, data

def whiten_data(mat, control=None, class_col=None):
    data = mat
    if control is not None:
        data = control

    means = data.mean(axis=0)
    stds = data.std(axis=0)
    # Standardize mat according to mean/std in source array
    for col in xrange(mat.shape[1]):
        # Don't standardize the class column
        if col == class_col:
            continue

        if stds[col] == 0:
            mat[:, col] = 0
        else:
            mat[:, col] = (mat[:, col] - means[col]) / stds[col]
    
def script(filename, control=None, add_class=None):
    header, class_col, mat = read_examples(filename)
    if control is not None:
        control_header, control_class_col, control = read_examples(control)
        assert header == control_header

    whiten_data(mat, control=control, class_col=class_col)

    if add_class is not None:
        print '#class\t%s' % header
    else:
        print '#%s' % header

    for row in mat:
        values = []
        if add_class is not None:
            values.append('%d' % add_class)

        for i, val in enumerate(row):
            if i == class_col:
                values.append('%d' % val)
            else:
                values.append('%.4f' % val)

        print '\t'.join(values)
        
def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] FILE.mat"
    description = __doc__.strip()

    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option('-c', '--control', metavar='MAT',
                      dest='control', default=None,
                      help="Standardize according to data in MAT instead"
                      " of FILE.mat")
    parser.add_option('--class', metavar='INT', type='int',
                      dest='add_class', default=None,
                      help="If set, prepends a class column with the given"
                      " value")
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
