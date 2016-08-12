#! /usr/bin/env python

from itertools import chain

import utils
from utils import GTFLine

def main(args):
    gtf = [GTFLine(l).asdict() for l in utils.tab_line_gen(args.infile)]
    allkeys = set(chain.from_iterable(d.keys() for d in gtf))
    columns = GTFLine.GTFCOLS + GTFLine.ATTRORDER
    columns += [k for k in sorted(allkeys) if k not in columns]
    print >>args.outfile, '\t'.join(columns)
    for d in gtf:
        print >>args.outfile, '\t'.join(str(d[c]) if c in d else '' for c in columns)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Create tab-delimited table from GTF')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help="Input GTF file")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output text file")
    main(parser.parse_args())
