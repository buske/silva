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
sys.path.insert(0, os.path.expandvars("$SILVA_PATH/src/share"))
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
    cur_result = []
    for i, line in enumerate(output):
        if i % 5 == 0:
            cur_result = []
        if i % 5 == 2:
            m = re2.match(line)
            assert m
            cur_result.append(float(m.group(1)))
            results.append(cur_result)
        elif i % 5 == 4:
            m = re4.match(line)
            assert m
            cur_result.append(float(m.group(1)))
        
    assert len(results) == len(seqs)
    return results

def iter_sequences(filename, domain, **kwargs):
    with maybe_gzip_open(filename) as ifp:
        for line in ifp:
            try:
                seq = line.strip().upper()
                premrna = seq.replace('|', '')
                pre, post = premrna.split('/')
                pre, old = pre.split('[')
                new, post = post.split(']')
                pre = pre[-domain:]
                post = post[:domain]

                assert len(pre) == len(post) == domain
                yield pre + old + post, pre + new + post
            except:
                print >>sys.stderr, "Error parsing sequence: skipping"
                yield None

def script(filename, quiet=False, domain=None, **kwargs):
    fields = ['ddG_pre', 'dvar_pre'][1:]
    if domain is not None:
        fields = ['%s_%d' % (field, domain) for field in fields]
        
    print "#%s" % '\t'.join(fields)
    seqs = []
    for entry in iter_sequences(filename, domain=domain, **kwargs):
        if entry is None:
            print '\t'.join('na' * len(fields))
            continue
        
        seqs.extend(entry)
        
    scores = score_sequence(*seqs)

    def safe_f(args):
        new, old = args
        if old == 0:
            return 'na'
        else:
            return '%.4f' % log10(new / old)
        
    for pre_old, pre_new in zip(scores[::2], scores[1::2]):
        print '\t'.join(map(safe_f, zip(pre_new, pre_old)[1:]))

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
