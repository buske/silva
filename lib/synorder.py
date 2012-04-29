from __future__ import with_statement, division

import sys

from datetime import datetime
from gzip import open as _gzip_open
from contextlib import closing


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
