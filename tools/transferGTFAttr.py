#! /usr/bin/env python

import utils
from utils import GTFLine

def main(args):
    combined_gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]
    for g in combined_gtf:
        if args.fromAttr in g.attr:
            g.attr[args.toAttr] = g.attr[args.fromAttr]
        print >>args.outfile, g

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Sort GTF file')
    parser.add_argument('fromAttr', 
                        help="Attribute to transfer from")
    parser.add_argument('toAttr', 
                        help="Attribute to transfer to")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help="Input GTF file")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF file")
    main(parser.parse_args())
