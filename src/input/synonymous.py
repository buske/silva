#!/usr/bin/env python

"""
VARIANTS is tab-delimited, with fields: chrom, pos, id, ref, alt, ...
(compatible with VCF format; id is not used)
If ACTION is 'filter', only synonymous variants are printed, and two columns
are added: gene name and transcript name (for the longest transcript in
which the variant is synonymous).
If ACTION is 'annotate' or 'generate', VARIANTS should contain the 7 columns
outputted by first 'filter'ing. 
If ACTION is 'annotate', additional columns are added: strand, codon, 
codon_offset, and pre-mRNA sequence.
If ACTION is 'generate', either random (--random) or all (--all) synonymous
variants are generated and printed to stdout. With --random, variants are
gene-matched, and additional flags allow matching other dimensions.
"""

# Author: Orion Buske
# Date:   ...
from __future__ import division, with_statement

import os
import sys
import cPickle

from collections import defaultdict
from string import maketrans
from random import sample

assert os.getenv('SILVA_PATH') is not None, \
    "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars('$SILVA_PATH/src/share'))
from silva import maybe_gzip_open

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
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

NUCS = set('ACGT')
# Create dict: CODON: set((frame, alt),...), # synonymous codon mutations
SYN_MUTATIONS = defaultdict(set)
for codon, aa in AA_CODE.iteritems():
    for pos in range(len(codon)):
        for alt in NUCS:
            if alt != codon[pos] and aa == AA_CODE[codon[:pos] + alt + codon[pos+1:]]:
                # Mutating codon[pos] to alt is a synonymous change
                SYN_MUTATIONS[codon].add((pos, alt))

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
COMPLEMENT_TAB = maketrans('ACGT', 'TGCA')

class Transcript(object):
    def __init__(self, gene, tx, chrom, tx_start, tx_end, strand, 
                 cds_start, cds_end, exon_starts, exon_ends, seq=None):
        """
        starts and ends: 0-indexed half-open

        Use cases:
        1) Have gene name and protein coords, need to identify transcript
        -- Need to verify that mutation matches mRNA
        -- Need to look up codons
        2) Have gene and transcript, need to look up mRNA
        3) Need to represent mRNA
        """
        try:
            self._valid = False
            self._chrom = chrom
            self._strand = strand
            self._gene = gene
            self._tx = tx
            
            self._tx_start = int(tx_start)
            self._tx_end = int(tx_end)
            self._tx_length = self._tx_end - self._tx_start
            
            self._cds = []
            self._cds_start = int(cds_start)
            self._cds_end = int(cds_end)
            self._mrna = None
            self._premrna = None

            assert self._cds_start != self._cds_end, "Zero-length CDS"
            assert self._tx_start <= self._cds_start <= \
                   self._cds_end <= self._tx_end, "CDS outside transcript"
            assert self._strand == '+' or self._strand == '-', \
                   "Invalid strand: %s" % self._strand
            
            exon_starts = [int(x) for x in exon_starts]
            exon_ends = [int(x) for x in exon_ends]
            assert len(exon_starts) == len(exon_ends)
            self._exons = zip(exon_starts, exon_ends)
            
            #print >>sys.stderr, "TX:  0 - %d" % (self._tx_length)
            #print >>sys.stderr, "CDS: %d - %d" % (self._cds_start - self._tx_start, self._cds_end - self._tx_start)
            for exon_start, exon_end in self._exons:
                #print >>sys.stderr, "Exon: %d - %d" % (b_start, b_start + b_len)
                cds_start = max(exon_start, self._cds_start)
                cds_end = min(exon_end, self._cds_end)
                if cds_start < cds_end:
                    self._cds.append((cds_start, cds_end))
                    
            assert self._cds, "Empty CDS"
            self._cds_length = sum([end - start for start, end in self._cds])
            assert self._cds_length % 3 == 0, "CDS length not a multiple of 3"

            if seq is not None:
                self.load_seq(seq)
            
            self._valid = True
        except AssertionError, e:
            print >>sys.stderr, "Skipping transcript: %s: %s" % (self, e)
            #print >>sys.stderr, self.__dict__
            #raise
            pass

    def __str__(self):
        s = "%s/%s, %s:%d-%d (%s)" % (self._gene, self._tx, self._chrom,
                                    self._tx_start, self._tx_end, self._strand)
        if self._mrna is not None:
            s = "%s (%s...%s)" % (s, self._mrna[:6], self._mrna[-6:])
        return s

    def __repr__(self):
        return "<Transcript: %s>" % (self)

    def valid(self):
        return self._valid
    
    def gene(self):
        return self._gene

    def tx(self):
        return self._tx

    def __len__(self):
        return self._tx_length

    def chrom(self):
        return self._chrom

    def strand(self):
        return self._strand

    def cds(self):
        return self._mrna

    def premrna(self):
        return self._premrna

    def _project_to_intervals(self, pos, intervals):
        """Return interval offset (0-indexed), given genome pos (1-indexed)
        intervals are 0-indexed, half-open
        """
        offset = 0
        if self._strand == '+':
            for start, end in intervals:
                if start < pos <= end:
                    return offset + (pos - (start + 1))
                elif pos < start + 1:
                    return None
                else:
                    offset += (end - start)
        else:
            for start, end in reversed(intervals):
                if start < pos <= end:
                    return offset + (end - pos)
                elif pos > end:
                    return None
                else:
                    offset += (end - start)
        
    def project_to_premrna(self, pos):
        """Return premrna offset (0-indexed), given genome pos (1-indexed)"""
        assert self._tx_start < pos <= self._tx_end
        if self._strand == '+':
            return pos - self._tx_start - 1  # -1 makes 0-indexed
        else:
            return self._tx_end - pos
        
    def project_to_cds(self, pos):
        """Return cds offset (0-indexed), given genome pos (1-indexed)"""
        assert self._cds_start < pos <= self._cds_end, \
               "Error, pos (%r) outside of CDS: [%r, %r)" % \
               (pos, self._cds_start, self._cds_end)
        return self._project_to_intervals(pos, self._cds)

    def _project_from_intervals(self, offset, intervals):
        """
        Return genome position (1-indexed), given interval offset (0-indexed)
        and 0-indexed, half-open intervals
        """
        if self._strand == '+':
            for start, end in intervals:
                size = end - start
                if offset < 0:
                    return None
                elif offset < size:
                    return start + offset + 1
                else:
                    offset -= size
        else:
            for start, end in reversed(intervals):
                size = end - start
                if offset < 0:
                    return None
                elif offset < size:
                    return end - offset
                else:
                    offset -= size

    def project_from_premrna(self, offset):
        """Return genome position (1-indexed), given premrna offset (0-indexed)"""
        assert 0 <= offset < self._tx_length
        if self._strand == '+':
            return self._tx_start + offset + 1  # +1 makes 1-indexed
        else:
            return self._tx_end - offset  # half-open tx_end makes 1-indexed

    def project_from_cds(self, offset):
        """Return genome position (1-indexed), given cds offset (0-indexed)"""
        assert 0 <= offset < self._cds_length
        return self._project_from_intervals(offset, self._cds)
    
    def load_seq(self, seq):
        assert seq is not None
        stop_codons = set(['TAA', 'TAG', 'TGA'])

        premrna = seq[self._tx_start:self._tx_end].tostring().upper()

        seqs = []
        for start, end in self._cds:
            seqs.append(seq[start:end].tostring().upper())

        if self._strand == '-':
            premrna = premrna.translate(COMPLEMENT_TAB)[::-1]
            seqs = [s.translate(COMPLEMENT_TAB)[::-1] for s in reversed(seqs)]
            
        self._mrna = ''.join(seqs)
        self._premrna = premrna
        
        assert len(self._mrna) == self._cds_length
        assert self._mrna[-3:] in stop_codons, \
               "Transcript ends with : %s" % self._mrna[-3:]
    
    def get_codon(self, aa):
        """Get codon at amino acid position (1-indexed)"""
        assert self._mrna is not None
        start = (aa - 1) * 3
        return self._mrna[start:start+3]

    def mutation_str(self, pos, ref, alt):
        """Given genomic pos (1-indexed), ref nuc and alt nuc, return
        mrna string in 'standard' form: ACT|ACGCACA[G/T]|ACAACA

        Splice sites are marked with pipes, mutation is in brackets
        """
        premrna = self._premrna
        assert premrna

        # Make nucs tx strand
        if self._strand == '-':
            alt = COMPLEMENT[alt]
            ref = COMPLEMENT[ref]
        
        mut_offset = self.project_to_premrna(pos)
        assert premrna[mut_offset] == ref, \
               "Found %s insead of %s as ref at pos %d in %s" % \
               (ref, premrna[mut_offset], pos, self)

        # Set of all edges between exons and introns
        edges = set([0, self._tx_length])
        edges.update([start - self._tx_start for start, end in self._exons])
        edges.update([end - self._tx_start for start, end in self._exons])
        
        if self._strand == '-':
            edges = [self._tx_length - e for e in edges]
            
        edges = list(sorted(edges))
        
        seqs = []
        for start, end in zip(edges[:-1], edges[1:]):
            seqs.append(premrna[start:end])

        # find sequence block with mutation
        for i, seq in enumerate(seqs):
            if mut_offset < len(seq):
                # insert mutation
                assert seq[mut_offset] == ref
                seqs[i] = '%s[%s/%s]%s' % (seq[:mut_offset], ref, alt, seq[mut_offset+1:])
                break
            else:
                mut_offset -= len(seq)

        # add splice markers
        return '|'.join(seqs)

    def is_synonymous(self, pos, ref, alt):
        """Is ref -> alt at pos (1-indexed, genomic) a synonymous change?"""
        try:
            cds_pos = self.project_to_cds(pos)
            assert cds_pos is not None
        except AssertionError:
            return False
        
        # Make nucs tx strand
        if self._strand == '-':
            alt = COMPLEMENT[alt]
            ref = COMPLEMENT[ref]

        assert self._mrna[cds_pos] == ref, \
               "Reference mismatch: %s:%s (%s)  %s != %s" % \
               (self._chrom, pos, self._strand, self._mrna[cds_pos], ref)
        
        frame = cds_pos % 3
        codon_start = cds_pos - frame
        old_codon = self._mrna[codon_start:codon_start+3]
        new_codon = old_codon[:frame] + alt + old_codon[frame+1:]
        return AA_CODE[old_codon] == AA_CODE[new_codon]

def iter_ucsc_genes(filename):
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            if line.startswith('#'): continue
            tokens = line.strip().split()
            bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, \
                 exonCount, exonStarts, exonEnds, score, name2 = tokens[:13]
                 
            chrom = chrom[3:] if chrom.startswith('chr') else chrom
            exonStarts = exonStarts.strip(',').split(',')
            exonEnds = exonEnds.strip(',').split(',')

            yield {'chrom': chrom,
                   'tx_start': txStart,
                   'tx_end': txEnd,
                   'strand': strand,
                   'cds_start': cdsStart,
                   'cds_end': cdsEnd,
                   'exon_starts': exonStarts,
                   'exon_ends': exonEnds,
                   'gene': name2,
                   'tx': name}
    

def get_genes(gene_filename=None, cache_filename=None,
              genomedata=None, **kwargs):
    """Loads (potentially cached) dict: gene_name -> set(genes)

    If not cached, genomedata path expected to provide sequence data
    """
    assert gene_filename and genomedata or cache_filename
    if cache_filename is not None and os.path.isfile(cache_filename):
        print >>sys.stderr, "Loading genes from pickled file: %s" % cache_filename
        with maybe_gzip_open(cache_filename) as ifp:
            genes = cPickle.load(ifp)
    else:
        from genomedata import Genome

        genes = defaultdict(set)
        cur_chrom = None
        cur_seq = None
        with Genome(genomedata) as genome:
            for entry in iter_ucsc_genes(gene_filename):
                chrom = entry['chrom']
                if chrom != cur_chrom:
                    try:
                        chromosome = genome['chr%s' % chrom]
                    except KeyError:
                        continue
                    print >>sys.stderr, "Cached sequence for chr%s." % chrom

                    cur_seq = chromosome.seq[0:chromosome.end]
                    cur_chrom = chrom

                # Substitute id with gene name
                entry['seq'] = cur_seq
                t = Transcript(**entry)

                if t.valid():
                    genes[entry['gene']].add(t)

        genes = dict(genes)  # remove defaultdict
        if cache_filename:
            print >>sys.stderr, "Saving genes to pickled file: %s" % cache_filename
            with open(cache_filename, 'wb') as ofp:
                cPickle.dump(genes, ofp, cPickle.HIGHEST_PROTOCOL)

    return genes

def get_transcript_from_protein(genes, gene, aa_pos, aa, mutation, *args, **kwargs):
    """Find longest transcripts matching protein coordinates

    gene: gene name (e.g. 'ABCB1')
    aa_pos: amino acid position (e.g. 1145), 1-indexed
    aa: amino acid code (e.g. 'I')
    mutation: nucleotide from, to (e.g. 'C>T')

    returns (tx, chrom, pos, ref, alt)
    """
    aa_pos = int(aa_pos)
    txs = genes.get(gene, [])
    matches = []
    for tx in txs:
        codon = tx.get_codon(aa_pos)
        nuc_from, nuc_to = mutation.split('>')
        if not codon or AA_CODE[codon] != aa: continue
        
        # Verify synonymous mutation is unambiguous
        frame = None
        for i in range(len(codon)):
            if codon[i] == nuc_from and \
                   AA_CODE[codon[:i]+nuc_to+codon[i+1:]] == aa:
                if frame is None:
                    frame = i
                else:
                    frame = None  # ambiguous
                    break
                
        if frame is not None: # unambiguous match
            ref = nuc_from if tx.strand() == '+' else COMPLEMENT[nuc_from]
            alt = nuc_to if tx.strand() == '+' else COMPLEMENT[nuc_to]
            cds_offset = (aa_pos - 1) * 3 + frame
            pos = tx.project_from_cds(cds_offset)
            matches.append((tx, tx.chrom(), pos, ref, alt))
            
    if matches:
        return max(matches)  # longest transcript
    else:
        print >>sys.stderr, "No match found for %s, %s, %s, %s" % (gene, aa_pos, aa, mutation)
        print >>sys.stderr, "Found %d transcripts" % len(txs)
        for tx in txs:
            codon = tx.get_codon(aa_pos)
            nuc_from, nuc_to = mutation.split('>')
            print >>sys.stderr, "%s: codon: %s -> %s" % (tx.tx(), codon, AA_CODE.get(codon, ''))


def get_transcript(genes, pos, ref, alt, gene_id, tx_id):
    try:
        txs = genes[gene_id]
    except KeyError:
        print >>sys.stderr, "Missing entry for gene: %s" % gene_id
        return

    txs = [x for x in txs if x.tx() == tx_id]
    if len(txs) == 0:
        print >>sys.stderr, "Missing entry for transcript: %s" % tx_id
        return
    else:
        # Possibly multiple transcripts with same accession
        txs = filter(lambda tx: tx.is_synonymous(pos, ref, alt), txs)

    if len(txs) == 0:
        print >>sys.stderr, "Mutation g.%d%s>%s does not appear to be" \
              "synonymous in %s" % (pos, ref, alt, tx_id)
        return
    elif len(txs) == 1:
        return txs[0]
    else:
        print >>sys.stderr, "Mutation g.%d%s>%s is not uniquely" \
              "synonymous among copies of %s" % (pos, ref, alt, tx_id)
        return

def random_synonymous_site(tx, cpg=None, avoid_splice=False):
    """Return random synonymous site (cds offset, ref, alt)
    
    cpg: if True or False, returned site is matched to this
    avoid_splice: if True, no mutations within 3 bp of splice site will be chosen
    """
    splice_sites = []
    cum_sum = 0
    for start, end in tx._exons:
        splice_sites.append(cum_sum)
        cum_sum += (end - start)
        
    matched_sites = set()
    syn_sites = set()
    for offset in xrange(0, len(tx._mrna), 3):
        codon = tx._mrna[offset:offset + 3]
        for frame, alt in SYN_MUTATIONS[codon]:
            ref = codon[frame]
            has_cpg = bool('CG' in tx._mrna[offset+frame-1:offset+frame+2])
            splice_dist = min([abs(s - offset + frame) for s in splice_sites])
            syn_sites.add((offset+frame, ref, alt))
            if cpg is None or cpg == has_cpg:
                if not avoid_splice or splice_dist > 3:
                    matched_sites.add((offset+frame, ref, alt))

    if matched_sites:
        sites = matched_sites
    else:
        sites = syn_sites
        print >>sys.stderr, "Warning: no matched synonymous mutation possible"
        
    return sample(sites, 1)[0]

def random_controls(genes, filename, match_cpg=False, avoid_splice=False):
    fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx']
    print '#%s' % '\t'.join(fields)
    for line in maybe_gzip_open(filename):
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
            cpg_seq = tx.premrna()[offset-1:offset+2]
            if tx.strand() == '-':
                cpg_seq = cpg_seq.translate(COMPLEMENT_TAB)

            cpg = 'CG' in cpg_seq
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
            
    fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx']
    print '#%s' % '\t'.join(fields)
    n_total = 0
    n_kept = 0
    for line in maybe_gzip_open(filename):
        line = line.rstrip()
        if not line or line.startswith('#'): continue
        n_total += 1
        
        tokens = line.split()
        if protein_coords:
            match = get_transcript_from_protein(genes, *tokens)
            (tx, chrom, pos, ref, alt) = match
            id = '.'
        else:
            chrom, pos, id, ref, alts = tokens[:5]
            chrom = chrom[3:] if chrom.startswith('chr') else chrom
            alt = alts.split(',')[0]
            # Only process SNVs
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

        n_kept += 1
        print '\t'.join([chrom, str(pos), id, ref, alt, tx.gene(), tx.tx()] + tokens[5:])

    print >>sys.stderr, "Found %d synonymous variants (%d dropped)" % \
          (n_kept, n_total - n_kept)

def annotate_variants(genes, filename):
    fields = ['chrom', 'pos', 'id', 'ref', 'alt', 'gene', 'tx', 'strand',
              'codon', 'frame', 'premrna']
    print '#%s' % '\t'.join(fields)
    for line in maybe_gzip_open(filename):
        line = line.rstrip()
        if not line or line.startswith('#'): continue

        tokens = line.split()
        chrom, pos, id, ref, alt, gene_id, tx_id = tokens[:7]
        chrom = chrom[3:] if chrom.startswith('chr') else chrom
        pos = int(pos)

        tx = get_transcript(genes, pos, ref, alt, gene_id, tx_id)
        if not tx:
            continue

        # Get codon, frame, and mrna
        cds_offset = tx.project_to_cds(pos)
        aa_pos = int(cds_offset / 3) + 1
        codon = tx.get_codon(aa_pos)
        frame = cds_offset % 3

        mut_str = tx.mutation_str(pos, ref, alt)
        print '\t'.join([chrom, str(pos), id, ref, alt, tx.gene(), tx.tx(),
                          tx.strand(), codon, str(frame), mut_str] + tokens[7:])



def script(action, filename, protein_coords=False, genomedata=None,
           all=False, random=False, match_cpg=False, avoid_splice=False,
           **kwargs):
    genes = get_genes(genomedata=genomedata, **kwargs)


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
    elif action == 'annotate':
        annotate_variants(genes, filename)
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
    print t1
    print t1.__dict__
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
    
    print t1
    print t1.__dict__
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
    parser.add_option("-G", "--genomedata", metavar="GD",
                      dest="genomedata", default=None,
                      help="Extract sequence/track data from Genomedata")
    parser.add_option("--protein-coords", action="store_true",
                      dest="protein_coords", default=False,
                      help="If ACTION is 'filter', VARIANTS file contains"
                      "protein coordinates, not chromosomal coordinates")
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
