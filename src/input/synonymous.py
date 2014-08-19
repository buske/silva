#!/usr/bin/env python

"""
VARIANTS is tab-delimited, with fields: chrom, pos, id, ref, alt, ...
(compatible with VCF format; id is not used)
If ACTION is 'filter', only synonymous variants are printed, and three columns
are added: 1) gene name, 2) transcript name, and 3) mutation information 
(for the longest transcript in which the variant is synonymous).
If ACTION is 'generate', VARIANTS should contain the 8 columns
outputted by first 'filter'ing. 
If ACTION is 'generate', either random (--random) or all (--all) synonymous
variants are generated and printed to stdout. With --random, variants are
gene-matched, and additional flags allow matching other dimensions.
"""

# Author: Orion Buske
# Date:   ...
from __future__ import division, with_statement

import os
import sys

from collections import defaultdict
from random import sample

assert os.getenv('SILVA_PATH') is not None, \
    "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/lib/python'))

from silva import maybe_gzip_open, COMPLEMENT, COMPLEMENT_TAB, AA_CODE, \
    get_transcript, get_transcript_from_protein, Transcript, get_genes

NUCS = set('ACGT')
# Create dict: CODON: set((frame, alt),...), # synonymous codon mutations
SYN_MUTATIONS = defaultdict(set)
for codon, aa in AA_CODE.iteritems():
    for pos in range(len(codon)):
        for alt in NUCS:
            if alt != codon[pos] and aa == AA_CODE[codon[:pos] + alt + codon[pos+1:]]:
                # Mutating codon[pos] to alt is a synonymous change
                SYN_MUTATIONS[codon].add((pos, alt))

def random_synonymous_site(tx, cpg=None, avoid_splice=False):
    """Return random synonymous site (cds offset, ref, alt)
    cpg: if True or False, returned site is matched to this
    avoid_splice: if True, no mutations within 3 bp of splice site will be chosen
    """
    splice_sites = []
    cum_sum = 0
    cds_intervals = reversed(tx._cds) if tx._strand == '-' else tx._cds
    for start, end in cds_intervals:
        splice_sites.append(cum_sum)
        cum_sum += (end - start)
        
    assert cum_sum == len(tx._mrna)

    matched_sites = set()
    syn_sites = set()
    for codon_offset in xrange(0, len(tx._mrna), 3):
        codon = tx._mrna[codon_offset:codon_offset + 3]
        for frame, alt in SYN_MUTATIONS[codon]:
            offset = codon_offset + frame
            ref = codon[frame]
            pre = tx._mrna[max(offset - 1, 0):offset]
            post = tx._mrna[offset + 1:offset + 2]
            has_cpg = bool((pre and pre[0] == 'C' and 
                            (ref == 'G' or alt == 'G')) or 
                           (post and post[0] == 'G' and 
                            (ref == 'C' or alt == 'C')))
            splice_dist = min([abs(s - offset) for s in splice_sites])
            syn_sites.add((offset, ref, alt))
            if cpg is None or cpg == has_cpg:
                if not avoid_splice or splice_dist > 3:
                    matched_sites.add((offset, ref, alt))

    if matched_sites:
        sites = matched_sites
    else:
        sites = syn_sites
        print >>sys.stderr, "Warning: no matched synonymous mutation possible"

    return sample(sites, 1)[0]

def random_controls(genes, filename, match_cpg=False, avoid_splice=False):
    fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx']
    print '#%s' % '\t'.join(fields)
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            line = line.rstrip()
            if not line or line.startswith('#'): continue

            tokens = line.split()
            chrom, pos, id, ref, alt, gene_id, tx_id = tokens[:7]
            chrom = chrom[3:] if chrom.startswith('chr') else chrom
            pos = int(pos)
            tx = get_transcript(genes, pos, ref, alt, gene_id, tx_id)
            if not tx:
                continue

            if match_cpg:
                offset = tx.project_to_premrna(pos)
                pre = tx.premrna()[offset-1:offset]
                post = tx.premrna()[offset+1:offset+2]
                tx_ref = ref
                tx_alt = alt
                if tx.strand() == '-':
                    tx_ref = ref.translate(COMPLEMENT_TAB)
                    tx_alt = alt.translate(COMPLEMENT_TAB)
                assert tx_ref == tx.premrna()[offset]
                cpg = bool((pre and pre[0] == 'C' and
                            (tx_ref == 'G' or tx_alt == 'G')) or
                           (post and post[0] == 'G' and
                            (tx_ref == 'C' or tx_alt == 'C')))
            else:
                cpg=None

            cds_offset, new_ref, new_alt = \
                random_synonymous_site(tx, cpg=cpg, avoid_splice=avoid_splice)
            new_pos = tx.project_from_cds(cds_offset)
            if tx.strand() == '-':
                new_ref = COMPLEMENT[new_ref]
                new_alt = COMPLEMENT[new_alt]

            print '\t'.join([chrom, str(new_pos), id, new_ref, new_alt, 
                             gene_id, tx.tx()] + tokens[7:])

def print_all_synonymous(genes):
    fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx']
    print '#%s' % '\t'.join(fields)
    
    for gene, txs in genes.iteritems():
        tx = max(txs)
        for offset in xrange(0, len(tx._mrna), 3):
            codon = tx._mrna[offset:offset + 3]
            for frame, alt in SYN_MUTATIONS[codon]:
                ref = codon[frame]
                if tx.strand() == '-':
                    ref = COMPLEMENT[ref]
                    alt = COMPLEMENT[alt]

                pos = tx.project_from_cds(offset+frame)
                print '\t'.join([tx.chrom(), str(pos), '.', ref, alt, 
                                 gene, tx.tx()])

def filter_variants(genes, filename, protein_coords=False):
    # Do chromosomal binning to efficiently lookup overlapping transcripts
    n_bins = 2048
    tx_locations = defaultdict(lambda: defaultdict(list))
    for gene, txs in genes.iteritems():
        for tx in txs:
            assert tx.gene() == gene
            start = tx._cds_start + 1
            end = tx._cds_end
            i = int(start / n_bins)
            j = int(end / n_bins)
            for bin in xrange(i, j+1):
                tx_locations[tx.chrom()][bin].append((start, end, tx))

    def find_overlapping_transcripts(chrom, pos):
        bin = int(pos / n_bins)
        txs = tx_locations[chrom][bin]
        return [tx for (start, end, tx) in txs if start <= pos <= end]
            
    fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx', 'mut']
    print '#%s' % '\t'.join(fields)
    n_total = 0
    n_kept = 0
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            line = line.rstrip()
            if not line or line.startswith('#'): continue
            n_total += 1

            tokens = line.split('\t')
            if protein_coords:
                gene, codon, aa, mut = tokens[:4]
                rest = tokens[1:]
                match = get_transcript_from_protein(genes, gene, codon, 
                                                    aa, mut)
                if match is None:
                    tx = None
                else:
                    (tx, chrom, pos, ref, alt) = match
                    id = '.'
            else:
                chrom, pos, id, ref, alts = tokens[:5]
                rest = tokens[5:]
                chrom = chrom[3:] if chrom.startswith('chr') else chrom
                alt = alts.split(',')[0]
                # Only process SNVs with one allele
                if len(ref) != 1 or len(alt) != 1:
                    continue

                pos = int(pos)
                txs = []
                for tx in find_overlapping_transcripts(chrom, pos):
                    try:
                        if tx.is_synonymous(pos, ref, alt):
                            txs.append(tx)
                    except AssertionError:
                        continue

                tx = max(txs) if txs else None  # Take longest valid transcript

            if not tx:
                continue

            mut_str = tx.cds_mutation(pos, ref, alt)

            n_kept += 1
            print '\t'.join([chrom, str(pos), id, ref, alt, tx.gene(), 
                             tx.tx(), mut_str] + rest)

        print >>sys.stderr, "Found %d synonymous variants (%d dropped)" % \
              (n_kept, n_total - n_kept)

# def annotate_variants(genes, filename):
#     fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx', 'strand',
#               'codon', 'frame', 'premrna']
#     print '#%s' % '\t'.join(fields)
#     with maybe_gzip_open(filename) as ifp:
#         for line in ifp:
#             line = line.rstrip()
#             if not line or line.startswith('#'): continue

#             tokens = line.split()
#             chrom, pos, id, ref, alt, gene_id, tx_id = tokens[:7]
#             chrom = chrom[3:] if chrom.startswith('chr') else chrom
#             pos = int(pos)

#             tx = get_transcript(genes, pos, ref, alt, gene_id, tx_id)
#             if not tx:
#                 continue

#             # Get codon, frame, and mrna
#             cds_offset = tx.project_to_cds(pos)
#             aa_pos = int(cds_offset / 3) + 1
#             codon = tx.get_codon(aa_pos)
#             frame = cds_offset % 3

#             mut_str = tx.mutation_str(pos, ref, alt)
#             print '\t'.join([chrom, str(pos), id, ref, alt, tx.gene(), tx.tx(),
#                               tx.strand(), codon, str(frame), mut_str] + tokens[7:])



def script(action, filename, protein_coords=False, genome_filename=None,
           all=False, random=False, match_cpg=False, avoid_splice=False,
           **kwargs):
    genes = get_genes(genome_filename=genome_filename, **kwargs)

    if action == 'generate':
        if random:
            random_controls(genes, filename, match_cpg=match_cpg, 
                            avoid_splice=avoid_splice)
        elif all:
            print_all_synonymous(genes)
        else:
            raise NotImplementedError()
    elif action == 'filter':
        filter_variants(genes, filename, protein_coords=protein_coords)
    else:
        raise NotImplementedError()

def run_tests():
    class seq(object):
        def __init__(self, s):
            self.value = s
        def tostring(self):
            return self.value
        def __getitem__(self, val):
            return seq(self.value[val])

    t1 = Transcript('name1', 'tx1', 'chr1', 1, 23, '+', 3, 22,
                    [1, 11], [10, 23], 
                    seq=seq('CAAATGCCCTATTCCCCCCTAATCCCC'))
    #print t1
    #print t1.__dict__
    assert t1.get_codon(1) == 'ATG'
    assert t1.get_codon(3) == 'TTT'
    assert t1.get_codon(6) == 'TAA'
    assert t1.premrna() == 'AAATGCCCTATTCCCCCCTAAT'
    assert t1.cds() == 'ATGCCCTTTCCCCCCTAA'
    try:
        t1.project_to_premrna(1) is None
    except:
        pass
    else:
        assert False
        
    try:
        assert t1.project_to_cds(3) is None
    except:
        pass
    else:
        assert False
        
    assert t1.project_to_premrna(2) == 0
    assert t1.project_to_cds(4) == 0
    assert t1.project_to_cds(10) == 6
    assert t1.project_from_cds(6) == 10, "Found: %s" % t1.project_from_cds(6)
    assert t1.project_to_cds(12) == 7
    for i in range(len(t1.cds())):
        r = t1.project_to_cds(t1.project_from_cds(i))
        assert r == i, "Inconsistent CDS translation for offset: %d" % i
    for i in range(len(t1)):
        r = t1.project_to_premrna(t1.project_from_premrna(i))
        assert r == i, "Inconsistent premrna translation for offset: %d" % i

    # MINUS STRAND!
    t1 = Transcript('name1', 'tx1', 'chr1', 1, 23, '-', 3, 22,
                    [1, 15], [14, 23], 
                    seq=seq('AAATTAGGGGGGAATAGGGCATTCCCCCCC'))
    #                         uueeeeeeeeeeeieeeeeeeu
    print t1
    print t1.__dict__
    for i in range(20):
        site = random_synonymous_site(t1, cpg=None, avoid_splice=True)
        print site

    assert t1.get_codon(1) == 'ATG'
    assert t1.get_codon(3) == 'TTT'
    assert t1.get_codon(6) == 'TAA'
    assert t1.premrna() == 'AATGCCCTATTCCCCCCTAATT', "Found: %s" % t1.premrna()
    assert t1.cds() == 'ATGCCCTTTCCCCCCTAA'
    try:
        t1.project_to_premrna(24) is None
    except:
        pass
    else:
        assert False
        
    try:
        assert t1.project_to_cds(23) is None
    except:
        pass
    else:
        assert False
        
    assert t1.project_to_premrna(23) == 0
    assert t1.project_to_cds(22) == 0, "Found: %s" % t1.project_to_cds(22)
    assert t1.project_to_cds(16) == 6
    assert t1.project_to_cds(14) == 7
    for i in range(len(t1.cds())):
        r = t1.project_to_cds(t1.project_from_cds(i))
        assert r == i, "Inconsistent CDS translation for offset: %d" % i
    for i in range(len(t1)):
        r = t1.project_to_premrna(t1.project_from_premrna(i))
        assert r == i, "Inconsistent premrna translation for offset: %d" % i

def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] ACTION (VARIANTS|-)"
    description = __doc__.strip()
    
    parser = OptionParser(usage=usage,
                          description=description)
    parser.add_option("-O", "--cache-genes", metavar="FILE",
                      dest="cache_filename", default=None,
                      help="Read/write parsed genes to speed up re-runs")
    parser.add_option("-g", "--genes", metavar="UCSC",
                      dest="gene_filename", default=None,
                      help="Read genes from UCSC file, unless cache already"
                      " exists and specified with -O")
    parser.add_option("-G", "--genome", metavar="FASTA",
                      dest="genome_filename", default=None,
                      help="Extract sequence data from FASTA file (same assembly as --genes)")
    parser.add_option("--protein-coords", action="store_true",
                      dest="protein_coords", default=False,
                      help="If ACTION is 'filter', VARIANTS file contains"
                      " protein coordinates, not chromosomal coordinates")
    parser.add_option("--all", action="store_true",
                      dest="all", default=False,
                      help="Print all synonymous variants.")
    parser.add_option("--random", action="store_true",
                      dest="random", default=False,
                      help="Generate random synonymous variants.")
    parser.add_option("--match-cpg", action="store_true",
                      dest="match_cpg", default=False,
                      help="Random variants will be CpG-matched.")
    parser.add_option("--avoid-splice", action="store_true",
                      dest="avoid_splice", default=False,
                      help="Random variants will not be within 3bp of"
                      " annotated splice sites.")
    parser.add_option("--test", action="store_true", dest="run_tests")
    options, args = parser.parse_args()

    if options.run_tests:
        sys.exit(run_tests())
    
    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
