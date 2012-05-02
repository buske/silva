#!/usr/bin/env python

"""
Uses fas-hex3.txt in same directory to identify putative exonic
splice enhancers (ESEs) based upon hexamer subsequences.
Input sequences are read, one per line, of the form: AAGA[C/G]TCG.
Input sequence should be full exon (or more), if possible (with splice
junctions marked with pipe characters).
"""

# Author: Orion Buske
# Date:   27 December 2011

from __future__ import division, with_statement

import os
import sys
import re

sys.path.insert(0, os.path.expandvars("$SYNORDER_PATH/src/share"))
from synorder import maybe_gzip_open, print_args

def hexamer_subsequences(hexs, seq):
    """Return dict: pos -> hexamer found in seq"""
    found = dict()
    motif_len = 6  # hexamers
    for offset in xrange(0, len(seq) - motif_len + 1):
        target_seq = seq[offset: offset + motif_len]
        if target_seq in hexs:
            found[offset] = target_seq

    return found

def score_mutation(hexs, old_seq, new_seq):
    """Return differences in hex hits between new and old sequences
    @return: (number of hexamers lost, number of hexamers gained)
    """
    # hexamer_subsequences returns dict of {pos: hexamer}
    old = hexamer_subsequences(hexs, old_seq)
    new = hexamer_subsequences(hexs, new_seq)
    old_pos = set(old)
    new_pos = set(new)
    return len(old_pos - new_pos), len(new_pos - old_pos)

def read_hexamers(filename):
    hexs = []
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if line:
                if len(line) == 6:
                    hexs.append(line)
                else:
                    print >>sys.stderr, "Found invalid hexamer: %s" % line
    return set(hexs)

def iter_sequences(filename):
    # Get exon
    seq_re = re.compile(r'([ACGT]*)\[([ACGT])/([ACGT])\]([ACGT]*)')
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            seq = line.strip().upper()
            mut_exons = [chunk for chunk in seq.split('|') if '/' in chunk]
            assert len(mut_exons) == 1
            exon = mut_exons[0]
            m = seq_re.search(exon)
            if m:
                pre, old, new, post = m.groups()
                yield pre, old, new, post
            else:
                print >>sys.stderr, "Error, invalid sequence: %s" % exon
                yield None
                    
def script(filename, hex_filename='fas-hex3.txt', quiet=False, **kwargs):
    hexs = read_hexamers(hex_filename)
    fields = ['f_ESS_lost', 'f_ESS_gained']
    if quiet:
        print '#%s' % '\t'.join(fields)
        NULL = '\t'.join(['na'] * len(fields))
    else:
        print '#n_initial n_lost n_gained'
        NULL = 'na na na'
        
    def safe_div(num, denom):
        if num + denom > 0:
            return '%.4f' % (num / (num + denom))
        else:
            return 'na'
                
    for entry in iter_sequences(filename):
        if entry is None:
            print NULL
            continue
        
        pre, nuc_old, nuc_new, post = entry
        old = pre + nuc_old + post
        short_old = pre[-7:] + nuc_old + post[:7]
        short_new = pre[-7:] + nuc_new + post[:7]
        
        n_old = len(hexamer_subsequences(hexs, old))
        n_lost, n_gained = score_mutation(hexs, short_old, short_new)

        if quiet:
            print '\t'.join([safe_div(n_lost, n_old),
                             safe_div(n_gained, n_old)])
        else:
            print '%d\t%d\t%d' % (n_old, n_lost, n_gained)

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
