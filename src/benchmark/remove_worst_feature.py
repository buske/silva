#!/usr/bin/env python

"""
Usage: $0 in.input out.input MODEL...
"""

from __future__ import with_statement

import sys
import os

from collections import defaultdict
from subprocess import Popen, PIPE

if len(sys.argv) < 4:
    print >>sys.stderr, __doc__
    sys.exit(1)

inputfile = sys.argv[1]
outputfile = sys.argv[2]
modelfiles = sys.argv[3:]
mydir = os.path.dirname(sys.argv[0])
p = Popen([os.path.join(mydir, 'get_rf_importance.R')] + modelfiles,
          stdout=PIPE)
stdout, stderr = p.communicate()

found = defaultdict(int)
for line in stdout.splitlines():
    line = line.strip()
    if line:
        found[line] += 1

# Remove feature with most frequent worst importance
freq = [(count, line) for line, count in found.iteritems()]
freq.sort(reverse=True)
toremove = freq[0][1]

with open(inputfile) as ifp:
    header = ifp.readline()
    assert header.startswith('#')
    # Find column to remove
    cols = header.strip().translate(None, '#').split()
    # Figure out how R bastardized the column naming, so we
    # can remove the appropriate column
    trans = {}
    seen = set()
    for col in cols:
        prefix = col.rstrip('+-?')
        if col.endswith('?'):
            trans[prefix + '.'] = col
        elif col.startswith('f_pE'):
            if prefix in seen:
                trans[prefix + '..1'] = col
            else:
                trans[prefix + '.'] = col
                seen.add(prefix)
        else:
            assert prefix == col
            trans[col] = col

        
    print >>sys.stderr, "\nRemoving feature:"
    print trans[toremove]
    print >>sys.stderr, ""
    rmcol = cols.index(trans[toremove])

    # Write other columns to output
    cols.pop(rmcol)
    with open(outputfile, 'w') as ofp:
        print >>ofp, '#%s' % '\t'.join(cols)
        for line in ifp:
            fields = line.strip().split()
            fields.pop(rmcol)
            print >>ofp, '\t'.join(fields)
