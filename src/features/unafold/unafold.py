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

from math import log10
from subprocess import Popen, PIPE

assert os.getenv('SILVA_PATH') is not None, \
       "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars("$SILVA_PATH/lib/python"))
from silva import maybe_gzip_open, print_args

BIN = os.path.expandvars("$SILVA_PATH/tools/unafold/src/hybrid-ss-min")
DATA = os.path.expandvars("$SILVA_PATH/tools/unafold/data")

def score_sequence(*seqs):
    """Scores one or more sequences"""
    p = Popen([BIN, '--stream', '-E'],
              stdin=PIPE, stdout=PIPE, close_fds=True,
              env={'UNAFOLDDAT': DATA})
    input = ''.join(['>\n%s\n' % seq for seq in seqs])
    output = p.communicate(input)[0]
    output = [float(val) for val in output.split('\n') if val]
    assert len(output) == len(seqs)
    return output

def iter_sequences(filename, domain=None, **kwargs):
    def get_mut_seqs(seq):
        pre, post = seq.split('/')
        pre, old = pre.split('[')
        new, post = post.split(']')

        if domain:
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
            seq = line.strip().upper()
            try:
                premrna = seq.replace('|', '')
                postmrna = ''.join(seq.split('|')[::2])
                yield get_mut_seqs(premrna), get_mut_seqs(postmrna)
            except (ValueError, AssertionError):
                print >>sys.stderr, "Error, invalid sequence: %s" % seq
                yield None

def script(filename, quiet=False, domain=None, **kwargs):
    fields = ['pdG_pre', 'pdG_post']
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
        try:
            return '%.4f' % -log10(new / old)
        except ValueError:
            return 'na'
        
    for dg_pre_old, dg_pre_new, dg_post_old, dg_post_new in \
            zip(scores[::4], scores[1::4], scores[2::4], scores[3::4]):
        print '\t'.join([safe_f(dg_pre_new, dg_pre_old), safe_f(dg_post_new, dg_post_old)])

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

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)

    if not options.quiet:
        print_args(args, kwargs)
    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
