#!/usr/bin/env python

"""
Input sequences are read, one per line, of the form: AAGA[C/G]TCG
"""

# Author: Orion Buske
# Date:   27 December 2011

from __future__ import division, with_statement

import os
import sys
import re

from subprocess import Popen, PIPE
from math import log10

assert os.getenv('SILVA_PATH') is not None, \
       "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars("$SILVA_PATH/lib/python"))
from silva import maybe_gzip_open, print_args

BIN = os.path.expandvars("$SILVA_PATH/tools/vienna/install/bin/RNAfold")
assert os.path.isfile(BIN), \
    "Error: missing file: %s" % BIN

def score_sequence(*seqs):
    """Scores one or more sequences"""
    p = Popen([BIN, '-p1', '-d2', '--noPS'],
              stdin=PIPE, stdout=PIPE, close_fds=True)
    input = '\n'.join(seqs)
    output = p.communicate(input)[0]
    
    output = output.splitlines()
    #print '\n'.join(output)
    re2 = re.compile(r'.*\[\s*([\d.-]+)\]')
    re4 = re.compile(r'.*ensemble diversity ([\d.-]+)')
    results = []
    for i, line in enumerate(output):
        if i % 5 == 4:
            m = re4.match(line)
            assert m
            results.append(float(m.group(1)))
        
    assert len(results) == len(seqs)
    return results

def iter_sequences(filename, domain, **kwargs):
    def get_mut_seqs(seq):
        pre, post = seq.split('/')
        pre, old = pre.split('[')
        new, post = post.split(']')
        pre_len = min(len(pre), domain)
        post_len = min(len(post), domain)
        # If too close to one end of sequence, accomodate
        if pre_len < domain:
            post_len = min(len(post), 2*domain - pre_len)
        if post_len < domain:
            pre_len = min(len(pre), 2*domain - post_len)

        pre = pre[-pre_len:]
        post = post[:post_len]
        assert len(pre) + len(post) == 2 * domain
        return pre + old + post, pre + new + post

    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            try:
                seq = line.strip().upper()
                premrna = seq.replace('|', '')
                postmrna = ''.join(seq.split('|')[::2])
                yield get_mut_seqs(premrna), get_mut_seqs(postmrna)
            except (ValueError, AssertionError):
                print >>sys.stderr, "Error parsing sequence: skipping"
                yield None

def script(filename, quiet=False, domain=None, **kwargs):
    fields = ['pvar_pre', 'pvar_post']
    if domain is not None:
        fields = ['%s_%d' % (field, domain) for field in fields]
        
    print "#%s" % '\t'.join(fields)
    seqs = []
    for entry in iter_sequences(filename, domain=domain, **kwargs):
        if entry is None:
            print '\t'.join(['na'] * len(fields))
            continue
        
        seqs.extend(entry[0])
        seqs.extend(entry[1])
        
    scores = score_sequence(*seqs)

    def safe_f(new, old):
        if old == 0:
            return 'na'
        else:
            return '%.4f' % -log10(new / old)
        
    for pre_old, pre_new, post_old, post_new in \
            zip(scores[::4], scores[1::4], scores[2::4], scores[3::4]):
        print '\t'.join([safe_f(pre_new, pre_old), safe_f(post_new, post_old)])

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
    parser.add_option("-d", "--domain", metavar="WIDTH",
                      dest="domain", type="int", default=None,
                      help="Limit analysis to within WIDTH bases"
                      " on either side of the mutation [default: None]")
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")
    elif options.domain is None:
        parser.error("Must specify domain width")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    if not options.quiet:
        print_args(args, kwargs)
    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
