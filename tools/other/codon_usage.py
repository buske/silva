#!/usr/bin/env python

"""
Given a genomic DNA mutation, compute codon usage features.

Input should be of the form:
<ref> <alt> <strand> <codon> <codon_offset>
For example, to represent AGC -> AGA:
G T - AGC 2
"""

# Author: Orion Buske
# Date:   04 January 2012

from __future__ import division, with_statement

import os
import sys

sys.path.insert(0, os.path.expandvars('$SYNORDER_PATH/lib'))
from synorder import maybe_gzip_open, print_args

COMPLEMENT = dict(zip('ACGT', 'TGCA'))
CODON_FREQ = {
    "TTT": 17.6, "TCT": 15.2, "TAT": 12.2, "TGT": 10.6,
    "TTC": 20.3, "TCC": 17.7, "TAC": 15.3, "TGC": 12.6,
    "TTA":  7.7, "TCA": 12.2, "TAA":  1.0, "TGA":  1.6,
    "TTG": 12.9, "TCG":  4.4, "TAG":  0.8, "TGG": 13.2,
    "CTT": 13.2, "CCT": 17.5, "CAT": 10.9, "CGT":  4.5,
    "CTC": 19.6, "CCC": 19.8, "CAC": 15.1, "CGC": 10.4,
    "CTA":  7.2, "CCA": 16.9, "CAA": 12.3, "CGA":  6.2,
    "CTG": 39.6, "CCG":  6.9, "CAG": 34.2, "CGG": 11.4,
    "ATT": 16.0, "ACT": 13.1, "AAT": 17.0, "AGT": 12.1,
    "ATC": 20.8, "ACC": 18.9, "AAC": 19.1, "AGC": 19.5,
    "ATA":  7.5, "ACA": 15.1, "AAA": 24.4, "AGA": 12.2,
    "ATG": 22.0, "ACG":  6.1, "AAG": 31.9, "AGG": 12.0,
    "GTT": 11.0, "GCT": 18.4, "GAT": 21.8, "GGT": 10.8,
    "GTC": 14.5, "GCC": 27.7, "GAC": 25.1, "GGC": 22.2,
    "GTA":  7.1, "GCA": 15.8, "GAA": 29.0, "GGA": 16.5,
    "GTG": 28.1, "GCG":  7.4, "GAG": 39.6, "GGG": 16.5}
STOP = '*'
AA_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': STOP, 'TAG': STOP,
    'TGT': 'C', 'TGC': 'C', 'TGA': STOP, 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
AA = set(AA_CODE.values())
AA_CODONS = dict([(aa, set([codon for codon, aa2 in AA_CODE.iteritems()
                              if aa2 == aa])) for aa in AA])
AA_FREQ = dict([(aa, sum([CODON_FREQ[codon] for codon in codons]))
                for aa, codons in AA_CODONS.iteritems()])
assert abs(sum(CODON_FREQ.values()) - 1000) < 1
assert abs(sum(AA_FREQ.values()) - 1000) < 1

def iter_lines(filename):
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            ref, alt, strand, codon, offset = line.split()[:5]
            assert strand in set(['+', '-', '.'])
            assert len(ref) == len(alt) == 1
            assert len(codon) == 3
            offset = int(offset)
            if strand == '-':
                ref = COMPLEMENT[ref]
                alt = COMPLEMENT[alt]
                
            assert codon[offset] == ref
            new_codon = codon[:offset] + alt + codon[offset+1:]

            yield codon, new_codon

def calc_delta_rscu(old_codon, new_codon):
    """Return change in relative synonymous codon usage (RSCU):
    RSCU = S*N_c/N_a, where:
    N_c: frequency of codon c
    N_a: frequency of amino acid a represented by codon c
    S: number of synonymous codons for a"""
    aa = AA_CODE[old_codon]
    aa2 = AA_CODE[new_codon]
    assert aa == aa2
    s = len(AA_CODONS[aa])
    n_c1 = CODON_FREQ[old_codon] / 1000
    n_c2 = CODON_FREQ[new_codon] / 1000
    n_a = AA_FREQ[aa] / 1000
    old_rscu = s * n_c1 / n_a
    new_rscu = s * n_c2 / n_a
    return old_rscu, new_rscu, abs(new_rscu - old_rscu)
                    
def script(filename, quiet=False, **kwargs):
    print '#%s' % '\t'.join(['mut_RSCU', 'delta_RSCU'])
    for codon, new_codon in iter_lines(filename):
        old, new, delta = calc_delta_rscu(codon, new_codon)
        print '%.4f\t%.4f' % (new, delta)

def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] (FILE|-)"
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
