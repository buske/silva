#!/usr/bin/env python

"""
Usage: gerp.py (-|FILE) TABLE [TABLE.pkl] > GERP

Given input lines of the form: chrom, pos, ...
(with pos 1-indexed), prints the corresponding GERP
score. Reads table of scores from TABLE, with lines
of the form: score, chrom, pos, .... If TABLE.pkl is
specified, and exists, data is read from it instead.
If it does not exist, it is created from the data in
TABLE, allowing faster data retrieval in the future.
"""

from __future__ import with_statement, division

import os
import sys
import cPickle

from collections import defaultdict

assert os.getenv('SILVA_PATH') is not None, \
       "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/lib/python'))
from silva import maybe_gzip_open

args = sys.argv[1:]

if sys.stdin.isatty():
    print >>sys.stderr, __doc__
    sys.exit(1)

input_file, table_file = args[:2]
if len(args) == 2:
    optfile = None
elif len(args) == 3:
    optfile = args[2]
else:
    print >>sys.stderr, __doc__
    sys.exit(1)

def read_table(filename, optfile=None):
    data = defaultdict(dict)
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            value, chrom, pos = line.split()
            data[chrom.lstrip('chr')][int(pos)] = float(value)

    if optfile and not os.path.isfile(optfile):
        print >>sys.stderr, "Saving optimized table to:", optfile
        with maybe_gzip_open(optfile, 'wb') as ofp:
            cPickle.dump(data, ofp, cPickle.HIGHEST_PROTOCOL)

    return data


if optfile and os.path.isfile(optfile):
    print >>sys.stderr, "Loading optimized table from:", optfile
    with maybe_gzip_open(optfile, 'rb') as ifp:
        table = cPickle.load(ifp)
else:
    print >>sys.stderr, "Loading table from:", table_file
    table = read_table(table_file, optfile)

print '#GERP++'
for line in maybe_gzip_open(input_file):
    line = line.strip()
    if not line or line.startswith('#'): continue
    chrom, pos, rest = line.split(None, 2)
    try:
        value = '%.4f' % table[chrom.lstrip('chr')][int(pos)]
    except (IndexError, KeyError):
        value = 'na'
        
    print value
