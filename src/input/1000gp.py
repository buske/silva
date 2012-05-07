#!/usr/bin/env python

"""
Usage: 1000gp.py VCF [VCF.pkl] < POSITIONS > AF

Given input lines of the form: chrom, pos, alt, ...
(with pos 1-indexed), prints the corresponding
allele frequencies (or '.' if not found).
Reads allele frequencies from VCF file, with info fields
containing AC and AN annotations. If VCF.pkl is
specified, and exists, data is read from it instead.
If it does not exist, it is created from the data in
VCF, allowing faster data retrieval in the future.
"""

from __future__ import with_statement, division

import os
import sys
import cPickle

from collections import defaultdict

args = sys.argv[1:]

if sys.stdin.isatty():
    print >>sys.stderr, __doc__
    sys.exit(1)

if len(args) == 1:
    optfile = None
elif len(args) == 2:
    optfile = args[1]
else:
    print >>sys.stderr, __doc__
    sys.exit(1)
    
tablefile = args[0]

def read_table(filename, optfile=None):
    data = defaultdict(dict)
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            tokens = line.split()
            chrom, pos = tokens[:2]
            alt = tokens[4].split(',')[0]
            info = tokens[7]
            # Split AC and AN from INFO field
            fields = info.split(';')
            fields.sort()
            ac = an = None
            for field in fields:
                if field.startswith('AC='):
                    ac = int(field.split('=')[1])
                    if an is not None:
                        break
                elif field.startswith('AN='):
                    an = int(field.split('=')[1])
                    if ac is not None:
                        break

            assert ac is not None and an is not None, \
                   "Error: entry with AC and AN: %s" % line
            
            data[chrom.lstrip('chr')][(int(pos), alt)] = float(ac) / an

    if optfile and not os.path.isfile(optfile):
        print >>sys.stderr, "Saving optimized table to:", optfile
        with open(optfile, 'wb') as ofp:
            cPickle.dump(data, ofp, cPickle.HIGHEST_PROTOCOL)

    return data


if optfile and os.path.isfile(optfile):
    print >>sys.stderr, "Loading optimized table from:", optfile
    with open(optfile, 'rb') as ifp:
        table = cPickle.load(ifp)
else:
    print >>sys.stderr, "Loading table from:", tablefile
    table = read_table(tablefile, optfile)

print '#1000GP_AF'
for line in sys.stdin:
    line = line.strip()
    if not line or line.startswith('#'): continue
    chrom, pos, alt = line.split()
    try:
        value = '%.4f' % table[chrom.lstrip('chr')][(int(pos), alt)]
    except (IndexError, KeyError):
        value = '.'
        
    print value
