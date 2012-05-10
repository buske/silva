#!/usr/bin/env python

"""
Usage: $0 < CONTROL.ARFF > CONTROL.STATS

Report stats for each feature.
"""
from __future__ import with_statement, division

import sys

from numpy import array

if sys.stdin.isatty():
    print __doc__
    sys.exit(1)

attrs = []
sums = []
sums2 = []
N = 0
for line in sys.stdin:
    line = line.strip()
    if not line: continue
    
    if line.startswith('@'):
        if line.startswith('@attribute'):
            attr = line.split()[1]
            attrs.append(attr)
            sums.append(0)
            sums2.append(0)
        elif line.startswith('@data'):
            # Set up arrays
            sums = array(sums, dtype=float)
            sums2 = array(sums2, dtype=float)
    else:
        # Data line
        values = array([float(val) for val in line.split(',')])
        sums += values
        sums2 += values * values
        N += 1

for attr, sum, sum2 in zip(attrs, sums, sums2):
    # Don't print stats for class
    if attr == 'class': continue
    avg = sum / N
    print "%s\t%.4f\t%.4f" % (attr, avg, (N / (N - 1)) * ((1 / N) * sum2 - (avg * avg)))
