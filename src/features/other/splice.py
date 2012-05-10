#!/usr/bin/env python

"""
Input sequences are read, one per line, of the form: GGAG|AAGA[C/G]TCG.
Prints splicing-related numbers. Mutation must be in exon.
"""

# Author: Orion Buske
# Date:   27 December 2011

from __future__ import division, with_statement

import os
import sys
import re

assert os.getenv('SYNORDER_PATH') is not None, \
       "Error: SYNORDER_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SYNORDER_PATH/src/share'))
from synorder import maybe_gzip_open, print_args

def script(filename, quiet=False, verbose=False, **kwargs):
    fields = ['premrna_f', 'mrna_f', 'splice_dist']
    print '#%s' % '\t'.join(fields)
    for line in maybe_gzip_open(filename):
        line = line.strip().upper()
        assert line.count('/') == 1

        pre, post = line.split('/')
        # Trim off mutation nucs and brackets
        pre = pre[:-2]  # e,g, '[A'
        post = post[2:]  # e.g. 'C]'

        pre_chunks = pre.split('|')
        post_chunks = post.split('|')
        # Assume mutation is in exon
        
        premrna_f = min(len(pre), len(post)) \
                       / (len(pre) + len(post) + 1)
        pre_cds = ''.join(pre_chunks[::2])
        post_cds = ''.join(post_chunks[::2])
        mrna_f = min(len(pre_cds), len(post_cds)) \
                    / (len(pre_cds) + len(post_cds) + 1)
        splice_dist = min(len(pre_chunks[-1]), len(post_chunks[0]))
        
        print '%.4f\t%.4f\t%d' % (premrna_f, mrna_f, splice_dist)

def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] (SEQ|-)"
    description = __doc__.strip()
    
    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option("-q", "--quiet", default=False,
                      dest="quiet", action='store_true',
                      help="Quiet output, suitable"
                      " for additional processing")
    parser.add_option("-v", "--verbose", default=False,
                      dest="verbose", action='store_true')
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    if not options.quiet:
        print_args(args, kwargs)
    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
