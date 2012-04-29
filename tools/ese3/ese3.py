#!/usr/bin/env python

"""
Uses weights.txt in same directory to identify putative exonic
splice enhancers (ESEs). Input sequences are read, one per line.
Sequences should be at least 15 nucleotides long, with 7 bases 5',
the mutated nucleotide: [OLD/NEW], and 7 bases 3'. For example:
AAGAGGT[C/G]TCGTTTACGGAGG. Suggested input is the full exon.
"""

# Author: Orion Buske
# Date:   27 December 2011

from __future__ import division, with_statement

import os
import sys
import re

sys.path.insert(0, os.path.expandvars("$SYNORDER_PATH/lib/synorder.py"))
from synorder import maybe_gzip_open, print_args

PRE_LEN = 7
POST_LEN = 7

def read_weights(filename):
    with open(filename) as ifp:
        return eval(ifp.read())

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

def score_target(motif, target):
    assert len(target) == len(motif['A'])
    score = 0.0
    for i, nuc in enumerate(target):
        score += motif[nuc][i]
    return score

def score_sequence(matrix, seq):
    """Return dict: pos -> score for all motif hits above threshold"""
    motif = matrix['motif']
    threshold = matrix['threshold']
    motif_len = len(motif['A'])

    hits = {}
    for offset in xrange(0, len(seq) - motif_len + 1):
        target_seq = seq[offset: offset + motif_len]
        score = score_target(motif, target_seq)
        if score >= threshold:
            hits[offset] = score

    return hits

def score_mutation(matrix, old_seq, new_seq, verbose=False):
    """Return tuple (# motifs lost, # motifs gained)"""
    old = score_sequence(matrix, old_seq)
    new = score_sequence(matrix, new_seq)
    if verbose:
        print >>sys.stderr, "Old scores:", str(old)
        print >>sys.stderr, "New scores:", str(new)
        
    old_pos = set(old)
    new_pos = set(new)
    return len(old_pos - new_pos), len(new_pos - old_pos)
                    
def script(filename, weight_filename='weights.txt',
           quiet=False, verbose=False, **kwargs):
    weights = read_weights(weight_filename)
    fields = ['f_motifs_lost', 'f_motifs_gained']

    if quiet:
        print '#%s' % '\t'.join(fields)
        NULL = '\t'.join(['na'] * 2)
    else:
        print "#PSSM  n_motifs  motifs_lost  motifs_gained"
        NULL = '\t'.join(['na'] * 3)

    def safe_div(num, denom):
        if num + denom > 0:
            return '%.4f' % (num / (num + denom))
        else:
            return 'na'

    for entry in iter_sequences(filename):
        if entry is None:
            print NULL
            continue
        
        pre, mut_old, mut_new, post = entry
        old = pre + mut_old + post
        
        tot = tot_lost = tot_gained = 0
        for name in sorted(weights):
            pssm = weights[name]
            n_old = len(score_sequence(pssm, old))

            short_old = pre[-7:] + mut_old + post[:7]
            short_new = pre[-7:] + mut_new + post[:7]
            n_lost, n_gained = score_mutation(pssm, short_old, short_new,
                                              verbose=verbose)

            tot += n_old
            tot_lost += n_lost
            tot_gained += n_gained
            if not quiet:
                print "%s\t%d\t%d\t%d" % (name, tot, n_lost, n_gained)

        if quiet:
            # Only print totalled effects
            print '\t'.join([safe_div(tot_lost, tot),
                             safe_div(tot_gained, tot)])

def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] (SEQ|-)"
    description = __doc__.strip()
    
    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option("-q", "--quiet", default=False,
                      dest="quiet", action='store_true',
                      help="Print one line per sequence, suitable"
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
