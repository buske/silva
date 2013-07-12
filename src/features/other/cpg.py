#!/usr/bin/env python

"""
Given a genomic DNA mutation, computes cpg-related features
including whether the mutation alters a CpG and the local
obs/exp CpG ratio.

Input mutation should be of the form:
ACACAGGGTTT[A/C]CAAACCGAGCG
"""

# Author: Orion Buske
# Date:   29 December 2011

from __future__ import division, with_statement

import os
import sys
import re

assert os.getenv('SILVA_PATH') is not None, \
           "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/lib/python'))
from silva import maybe_gzip_open, print_args

def iter_sequences(filename):
    seq_re = re.compile(r'([ACGT]*)\[([ACGT])/([ACGT])\]([ACGT]*)')
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            seq = line.strip().upper().translate(None, 'N')
            chunks = line.split('|')
            mut_seqs = [chunk for chunk in chunks if '/' in chunk]
            assert len(mut_seqs) == 1
            exon = mut_seqs[0]
            m = seq_re.match(exon)
            if m:
                pre, old, new, post = m.groups()
                yield pre, old, new, post
            else:
                print >>sys.stderr, "Error, invalid sequence: %s..." % seq[:20]
                yield None
                

def calc_cpg(seq):
    """Return obs/exp CpG
    
    Calculated using formula in:
    Gardiner-Garden M, Frommer M.
    CpG islands in vertebrate genomes.
    J. Mol. Biol. 1987 Jul 20;196(2):261-82.
    """
    denom = seq.count('C') * seq.count('G')
    if denom > 0:
        return seq.count('CG') * len(seq) / float(denom)
    else:
        return 0
                    
def script(filename, quiet=False, **kwargs):
    fields = ['CpG?', 'CpG_exon']
    print '#%s' % '\t'.join(fields)
    for entry in iter_sequences(filename):
        if entry is None:
            print '\t'.join(['na'] * len(fields))
            continue
        
        pre, old, new, post = entry
        # Was a CpG created or destroyed by the mutation?
        mut_cpg = bool((pre and  pre[-1] == 'C' and (old == 'G' or new == 'G')) or 
                       (post and post[0] == 'G' and (old == 'C' or new == 'C')))

        exon = pre+old+post
        cpg_exon = calc_cpg(exon)
        print '%d\t%.4f' % (mut_cpg, cpg_exon)

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
