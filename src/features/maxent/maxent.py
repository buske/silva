#!/usr/bin/env python

"""
Input sequences are read, one per line, of the form: AAGA[C/G]TCG.
If possible, the sequence should include the entire exon (or more,
with pipes marking all splice junctions).
"""

# Author: Orion Buske
# Date:   27 December 2011

from __future__ import division, with_statement

import os
import sys
import re

from datetime import datetime
from itertools import chain, imap
from subprocess import Popen, PIPE

assert os.getenv('SILVA_PATH') is not None, \
           "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars("$SILVA_PATH/src/share"))
from silva import maybe_gzip_open, print_args

MAXENT_PATH = os.path.expandvars('$SILVA_PATH/tools/maxent')
GOOD_SCORE = 2

class Seq(object):
    seq_re=re.compile('([ACGT]*)\[([ACGT])/([ACGT])\]([ACGT]*)')
    gt_re=re.compile('GT')
    ag_re=re.compile('AG')
    before = {5: 3, 3: 20}
    after = {5: 6, 3: 3}
    
    def __init__(self, seq):
        """Store a sequence and associated splicing info

        Given full sequence, find known and canonical splice sites (indices)
        and relevant subsequence (str)
        """
        assert seq.count('/') == 1
        seq = seq.strip().upper()
        chunks = seq.split('|')
        mut_offset = None
        for i, chunk in enumerate(chunks):
            if '/' in chunk:
                if i == 0:
                    pre = ''
                else:
                    pre_merge = ''.join(chunks[:i])
                    pre = pre_merge[-self.before[3]:]

                if i == len(chunks) - 1:
                    post = ''
                else:
                    post_merge = ''.join(chunks[i+1:])
                    post = post_merge[:self.after[5]]
                    
                m = self.seq_re.match(chunk)
                assert m
                head, nuc_old, nuc_new, tail = m.groups()
                mut_offset = len(pre) + len(head)
                
                break

        assert mut_offset is not None
        self._old_seq = pre + head + nuc_old + tail + post
        self._new_seq = pre + head + nuc_new + tail + post
        # Given splice junctions
        known = {}
        known[5] = len(self._old_seq) - len(post) if post else None
        known[3] = len(pre) if pre else None

        #print pre, head, nuc_old, "/", nuc_new, tail, post
        junctions = {}
        for side in [3, 5]:
            #print "\nSIDE", side
            # Add known sites, if they are there
            junctions[side] = set()
            if known[side]:
                junctions[side].add(known[side])
            # Make sure we stay within the bounds of the required site lengths
            start = max(len(pre), self.before[side])
            length = len(self._old_seq) - start - max(len(post), self.after[side])
            # Find all canonical sites in old and mutated exon
            junctions[side].update(self.find_canonical(self._old_seq, side, start=start,
                                                       length=length))
            # Add any junctions created with mutation
            junctions[side].update(self.find_canonical(self._new_seq, side,
                                                       start=start,
                                                       length=length))

        self.known = known
        self.junctions = junctions

    def iter_seqs(self, side):
        return chain(self.get_sites(side, self._old_seq, self.junctions[side]),
                     self.get_sites(side, self._new_seq, self.junctions[side]))

    def get_sites(self, side, seq, indices):
        return [seq[i - self.before[side]:i + self.after[side]] for i in indices]

    def find_canonical(self, seq, side, start=0, length=None):
        """Find all canonical splice sites of a particular side (5 or 3')

        Searches within seq[start:length] (end if length is None)
        Return list of indices into seq
        """
        assert side == 3 or side == 5
        end = None if length is None else start + length
        subseq = seq[slice(start, end)]
        #print "subseq", subseq
        if side == 3:
            # CCCCAGexon: + 2 to offset to splice junction
            junctions = [m.start() + start + 2 for m in self.ag_re.finditer(subseq)]
        else: # side == 5
            # exonGTCCCC
            junctions = [m.start() + start for m in self.gt_re.finditer(subseq)]

        #print "junctions", junctions
        return junctions

    def score(self, side, scores, min_score=GOOD_SCORE):
        """Return mutation scores using site scores dict: seq -> float(score)"""
        assert side == 3 or side == 5
        # self.known[side]: int
        # self.junctions[side]: set(int)
        
        # Do things in consistent order
        junctions = list(sorted(self.junctions[side]))
        known = self.known[side]

        # Lookup full sequences again
        old_sites = self.get_sites(side, self._old_seq, junctions)
        new_sites = self.get_sites(side, self._new_seq, junctions)
        assert len(junctions) == len(old_sites) == len(new_sites)
        # Lookup scores from sequences, and associate with junctions
        # (because of consistent order)
        try:
            junction_scores = dict([(i, (scores[old], scores[new])) for i, old, new
                                    in zip(junctions, old_sites, new_sites)])
        except:
            print side
            print junctions
            print old_sites
            print new_sites
            raise
        
        known_score = junction_scores[known][0] if known else None
        
        # Clean out any that aren't strong enough to be interesting
        for i, (old, new) in junction_scores.items():
            if abs(old - new) < 0.1 or (old < min_score and new < min_score):
                junction_scores.pop(i)

        # Aggregate in various ways
        old_scores = [(old, i) for i, (old, new) in junction_scores.iteritems()]
        old_scores.sort(reverse=True)
        new_scores = [(new, i) for i, (old, new) in junction_scores.iteritems()]
        new_scores.sort(reverse=True)
        junction_deltas = [(abs(new - old), i) for i, (old, new) in junction_scores.iteritems()]
        junction_gains = [(new - old, i) for i, (old, new) in junction_scores.iteritems()
                           if new > min_score]
        junction_losses = [(old - new, i) for i, (old, new) in junction_scores.iteritems()
                           if old > min_score]

        #print junction_scores
        #print junction_deltas
        #print junction_gains
        #print junction_losses

        def safe_max_value(scores, index=0):
            return max(scores)[index] if scores else None

        def safe_max(scores):
            return max(scores) if scores else None

        max_val = safe_max([safe_max_value(old_scores), safe_max_value(new_scores)])
#         try:
#             old1 = old_scores[0]
#             old2 = old_scores[1]
#             old_diff = (old1[0] - old2[0])
#             same_new_diff = (junction_scores[old1[1]][1] - junction_scores[old2[1]][1])
#             top_diff_diff = old_diff - same_new_diff
#         except IndexError:
#             top_diff_diff = None
            
        max_diff = safe_max(junction_deltas)
        old_max = safe_max(old_scores)
        new_max = safe_max(new_scores)
        is_max_cryptic = bool(new_max and new_max[1] != known)
        # If one of the maxes disappears or if the max changes to a different site
        did_max_change = bool((bool(old_max) + bool(new_max) == 1) or
                              (old_max and new_max and old_max[1] != new_max[1]))
        did_known_change_most = bool(max_diff and max_diff[1] == known)

        return max_val, safe_max_value(junction_deltas), \
               safe_max_value(junction_gains), safe_max_value(junction_losses), \
               did_max_change, is_max_cryptic, did_known_change_most

def parse_maxent_output(output):
    """Return dict: seq -> float(score)"""
    scores = {}
    for line in output.split('\n'):
        # Parse output lines of the form: <seq> <score>
        line = line.strip()
        if line:
            try:
                seq, score = line.split()
                scores[seq] = float(score)
            except:
                print >>sys.stderr, "Error parsing line: %r" % line
                raise
            
    return scores

def score_sites(side, sites):
    assert side == 5 or side == 3
    script = 'score%d.pl' % side
    p = Popen(['/usr/bin/perl', script, '-'],
              stdin=PIPE, stdout=PIPE, close_fds=True,
              cwd=MAXENT_PATH)

    stdin = '\n'.join(sites)
    stdout = p.communicate(stdin)[0]
    scores = parse_maxent_output(stdout)
    
    return scores


def print_row(scores):
    strs = []
    for val in scores:
        if val is None:
            strs.append('na')
        elif isinstance(val, bool) or isinstance(val, int):
            strs.append('%d' % val)
        else:
            try:
                strs.append('%.4f' % val)
            except TypeError:
                strs.append(str(val))

    print '\t'.join(strs)

def script(filename, quiet=False, verbose=False, **kwargs):
    fields = ['max', 'max_diff', 'gained', 'lost', 'max_changed?', 'max_is_cryptic?', 'known_changed_most?']
    print '#%s' % '\t'.join(fields)
    NULL = [None] * len(fields)

    seqs = []
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if line:
                seqs.append(Seq(line))
                
    sites = {}
    scores = {}
    # Accumulate sites for each side
    for side in [3, 5]:
        sites[side] = set()
        for s in seqs:
            sites[side].update(s.iter_seqs(side))

        # Score ALL sites at once!
        scores[side] = score_sites(side, sites[side])

    # Print stats for each object, given queried scores
    for s in seqs:
        # Compute field-wise max of rows for 3' and 5'
        max_row = imap(max, s.score(3, scores[3]), s.score(5, scores[5]))
        print_row(max_row)

def run_tests():
    s1 = Seq('c|ccaaaaaaaaaaaaaaaaaaAG|TTTGGGGGGGGGGG[T/A]TT|GTaaaaccccc|cc')
    assert set(s1.iter_seqs(3)) == set(['AAAAAAAAAAAAAAAAAAAGTTT'])
    assert set(s1.iter_seqs(5)) == set(['GGGGTTTGT',
                                        'GGGGATTGT',
                                        'TTTGTAAAA',
                                        'ATTGTAAAA'])

    s2 = Seq('c|ccaaaaaaaaaaaaaaaaaaAG|TTTGGGGGGG[T/A]GGGGGTTT|GTaaaaccccc|cc')
    assert set(s2.iter_seqs(3)) == set(['AAAAAAAAAAAAAAAAAAAGTTT',
                                        'AAAAAAAGTTTGGGGGGGTGGGG',
                                        'AAAAAAAGTTTGGGGGGGAGGGG'])
    assert set(s2.iter_seqs(5)) == set(['TTTGTAAAA',
                                        'GGGGTGGGG',
                                        'GGGGAGGGG',
                                        'GGGGTTTGT'])

    s3 = Seq('TT[T/A]GGGGGTTT|GTaaaaccccc|cc')
    assert set(s3.iter_seqs(3)) == set([])
    assert set(s3.iter_seqs(5)) == set(['TTTGTAAAA',
                                        'GGGGTTTGT'])

    seqs = [s1, s2]

    sites = {}
    for side in [3, 5]:
        sites[side] = set()
        for s in seqs:
            sites[side].update(s.iter_seqs(side))

        #print sites[side]
        scores = score_sites(side, sites[side])
        #print scores
        for s in seqs:
            print_row(s.score(side, scores, min_score=-30))

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
    parser.add_option("--test", default=False,
                      dest="test", action='store_true')
    options, args = parser.parse_args()

    if options.test:
        sys.exit(run_tests())

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
