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

assert os.getenv('SILVA_PATH') is not None, \
       "Error: SILVA_PATH is unset."
sys.path.insert(0, os.path.expandvars("$SILVA_PATH/src/share"))
from silva import maybe_gzip_open, print_args

assert os.getenv('UNAFOLD_BIN') is not None, \
       "Error: UNAFOLD_BIN environment variable must be set"
assert os.getenv('UNAFOLDDAT') is not None, \
       "Error: UNAFOLDDAT environment variable must be set"
HYBRID_SS_MIN = os.path.expandvars("$UNAFOLD_BIN")
UNAFOLDDAT = os.path.expandvars("$UNAFOLDDAT")

def score_sequence(*seqs):
    """Scores one or more sequences"""
    p = Popen([HYBRID_SS_MIN, '--stream', '-E'],
              stdin=PIPE, stdout=PIPE, close_fds=True,
              env={'UNAFOLDDAT': UNAFOLDDAT})
    input = ''.join(['>\n%s\n' % seq for seq in seqs])
    output = p.communicate(input)[0]
    output = [float(val) for val in output.split('\n') if val]
    assert len(output) == len(seqs)
    return output

def iter_sequences(filename, domain=None, **kwargs):
    # A bit of buffer above and below in case there are up to one
    # splicing marks: '|'
    if domain is None:
        width = '+'
    else:
        width = '{,%d}' % domain
        
    seq_re = re.compile(r'([ACGT]%(width)s)\[([ACGT])/([ACGT])\]([ACGT]%(width)s)'
                        % {'width': width})
    def groups_to_old_and_new(g):
        pre, old, new, post = g
        return pre + old + post, pre + new + post
    
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            seq = line.strip().upper()
            premrna = seq.replace('|', '')
            postmrna = ''.join(seq.split('|')[::2])
            
            if premrna and postmrna:
                prem = seq_re.search(premrna)
                postm = seq_re.search(postmrna)
                if prem and postm:
                    prem_seqs = groups_to_old_and_new(prem.groups())
                    postm_seqs = groups_to_old_and_new(postm.groups())
                    yield prem_seqs + postm_seqs
                else:
                    print >>sys.stderr, "Error, invalid sequence: %s" % seq
                    yield None

def script(filename, quiet=False, domain=None, **kwargs):
    if domain is not None:
        fields = ['f_premrna_delta_G_%d' % domain,
                 'f_postmrna_delta_G_%d' % domain]
    else:
        fields = ['f_premrna_delta_G', 'f_postmrna_delta_G']
        
    print "#%s" % '\t'.join(fields)
    seqs = []
    for entry in iter_sequences(filename, domain=domain, **kwargs):
        if entry is None:
            print '\t'.join('na' * len(fields))
            continue
        
        seqs.extend(entry)
        
    scores = score_sequence(*seqs)

    def safe_f(new, old):
        if old == 0:
            return 'na'
        else:
            return '%.4f' % abs((new - old) / old)
        
    for dg_pre_old, dg_pre_new, dg_post_old, dg_post_new in \
            zip(scores[::4], scores[1::4], scores[2::4], scores[3::4]):
        print '\t'.join([safe_f(dg_pre_new, dg_pre_old),
                         safe_f(dg_post_new, dg_post_old)])

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
