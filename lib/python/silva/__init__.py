from __future__ import with_statement, division

import os
import sys
import cPickle

from datetime import datetime
from gzip import open as _gzip_open
from contextlib import closing
from collections import defaultdict
from string import maketrans

from twobitreader import TwoBitFile as Genome


def maybe_gzip_open(filename, *args, **kwargs):
    if filename.endswith('.gz'):
        return closing(_gzip_open(filename, *args, **kwargs))
    elif filename == '-':
        return sys.stdin
    else:
        return open(filename, *args, **kwargs)

def print_args(args, kwargs, out=sys.stdout):
    params = [("command", ' '.join(sys.argv)),
              ("time", datetime.today().strftime("%Y-%m-%d %H:%M:%S"))]
    params.append(("options", repr(kwargs)))
    params.append(("arguments", repr(args)))
    for pair in params:
        print >>out, "##%s=%s" % pair


STOP = '*'
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
COMPLEMENT_TAB = maketrans('ACGT', 'TGCA')
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

class Transcript(object):
    def __init__(self, gene, tx, chrom, tx_start, tx_end, strand, 
                 cds_start, cds_end, exon_starts, exon_ends, seq=None):
        """
        starts and ends: 0-indexed half-open
        """
        self._valid = False
        self._chrom = chrom
        self._strand = strand
        self._gene = gene
        self._tx = tx

        self._tx_start = int(tx_start)
        self._tx_end = int(tx_end)
        self._tx_length = self._tx_end - self._tx_start

        self._cds = []  # List of pairs of genomic coordinates
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
        assert len(exon_starts) == len(exon_ends), \
            "Different number of exon starts and ends"
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

    def start(self):
        """Return start genomic coordinate (1-indexed)"""
        return self._tx_start + 1

    def end(self):
        """Return end genomic coordinate (1-indexed)"""
        return self._tx_end

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
        assert seq is not None, "seq is None"
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

        assert len(self._mrna) == self._cds_length, "mRNA length != CDS length"
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

    def cds_mutation(self, pos, ref, alt):
        """Given genomic pos (1-indexed), ref nuc and alt nuc, return
        cds mutation string in 'standard' form: c.74A>T
        """
        try:
            cds_offset = self.project_to_cds(pos)
            assert cds_offset is not None
        except AssertionError:
            return False
        
        # Make nucs tx strand
        if self._strand == '-':
            alt = COMPLEMENT[alt]
            ref = COMPLEMENT[ref]

        assert self._mrna[cds_offset] == ref, \
               "Reference mismatch: %s:%s (%s)  %s != %s" % \
               (self._chrom, pos, self._strand, self._mrna[cds_offset], ref)
        
        return "c.%d%s>%s" % (cds_offset + 1, ref, alt)

    def get_exon_bounds(self, pos):
        """Given genomic pos (1-indexed), return [left, right] offsets
        of exon containing position, such that exon sequence can be
        retrieved with premrna[left:right+1]

        Return None if pos doesn't overlap exon
        """
        pos = pos - 1
        for exon_start, exon_end in self._exons:
            if exon_start <= pos < exon_end:
                # position is within this exon
                # Get premrna offsets of exon boundaries
                pre_left = self.project_to_premrna(exon_start + 1)
                pre_right = self.project_to_premrna(exon_end)
                if self._strand == '-':
                    pre_left, pre_right = pre_right, pre_left

                assert pre_left < pre_right
                return pre_left, pre_right
        
        assert False, "Error, pos (%r) outside of premrna: [%d, %d)" % \
            (pos, self._tx_start, self._tx_end)

    def mrna_exon(self, pos, left=0, right=0):
        """Given genomic pos (1-indexed), return premrna sequence of exon,
        with buffer on left and right

        Return None if pos doesn't overlap exon
        """
        # return exon sequence, potentially buffered
        pre_left, pre_right = self.get_exon_bounds(pos)
        pre_left = max(pre_left - left, 0)
        pre_right = min(pre_right + right + 1, len(self._premrna))
        return self._premrna[pre_left:pre_right]

    def mrna_window(self, pos, left=0, right=0):
        """Given genomic pos (1-indexed), return premrna sequence around
        pos with buffer on left and right, and index of pos into sequence
        """
        premrna_offset = self.project_to_premrna(pos)
        pre_left = max(premrna_offset - left, 0)
        pre_right = min(premrna_offset + right + 1, len(self._premrna))
        return self._premrna[pre_left:pre_right], premrna_offset - pre_left

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
              genome_filename=None, **kwargs):
    """Loads (potentially cached) dict: gene_name -> set(genes)

    If not cached, genome_filename FASTA expected to provide sequence data
    """
    assert gene_filename and genome_filename or cache_filename
    if cache_filename is not None and os.path.isfile(cache_filename):
        print >>sys.stderr, "Loading genes from pickled file: %s" % cache_filename
        with maybe_gzip_open(cache_filename) as ifp:
            genes = cPickle.load(ifp)
    else:
        genome = Genome(genome_filename)

        genes = defaultdict(set)
        missed_chroms = set()
        n_zero_len = 0
        for entry in iter_ucsc_genes(gene_filename):
            chrom = entry['chrom']
            if not chrom.startswith('chr'):
                chrom = 'chr%s' % chrom
                
            if chrom not in genome:
                if chrom not in missed_chroms:
                    print >>sys.stderr, "Could not find sequence for %s" \
                        " in: %s" % (chrom, genome_filename)
                    missed_chroms.add(chrom)
                
                continue

            # Substitute id with gene name
            entry['seq'] = genome[chrom]
            try:
                t = Transcript(**entry)
            except AssertionError, e:
                if "Zero-length CDS" in str(e):
                    n_zero_len += 1
                else:
                    print >>sys.stderr, "Skipping transcript: %s: %s" \
                        % (entry['gene'], e)
                continue

            if t.valid():
                genes[entry['gene']].add(t)

        if n_zero_len:
            print >>sys.stderr, "Skipped %d transcripts with zero-length CDS" \
                " annotations" % n_zero_len

        if missed_chroms:
            print >>sys.stderr, "Missing sequences with gene annotations: %s" \
                % ', '.join(sorted(missed_chroms))
            
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

    # If many hits, first try removing those on alternative haplotypes
    if len(matches) > 1:
        matches_no_haps = [m for m in matches if '_' not in m[1]]
        if matches_no_haps:
            matches = matches_no_haps
            
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


def iter_mutation_seqs(filename, genes, left=0, right=0):
    """Reads VCF-like annotated SilVA variants from file and returns
    old, new, and exonic sequences around each line's variant
    left/right bp buffer returned around variant and exon sequence
    Yields tuple of (exon, old, new) sequences for each line.
    """
    global COMPLEMENT

    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            if line.startswith('#'): continue

            tokens = line.strip().split('\t')
            chrom, pos = tokens[:2]
            pos = int(pos)
            ref, alt = tokens[3:5]
            gene, tx_id = tokens[5:7]

            tx = None
            for t in genes[gene]:
                if t.tx() == tx_id and chrom == t.chrom() and t.start() <= pos <= t.end():
                    tx = t

            assert tx

            exon = tx.mrna_exon(pos, left=left, right=right)
            assert exon

            if tx.strand() == '-':
                ref, alt = COMPLEMENT[ref], COMPLEMENT[alt]
            
            old, offset = tx.mrna_window(pos, left=left, right=right)
            assert old[offset] == ref
            new = old[:offset] + alt + old[offset+1:]
            assert len(old) == len(new)

            yield exon, old, new
