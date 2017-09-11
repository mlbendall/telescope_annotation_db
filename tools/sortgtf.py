#! /usr/bin/env python

import utils

def main(args):
    if args.chroms:
        chroms = [l.strip('\n').split('\t')[0] for l in args.chroms]
    else:
        chroms = None
    
    for l in utils.sort_gtf(utils.tab_line_gen(args.infile), chroms):
        print >>args.outfile, '\t'.join(l)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Sort GTF file')
    parser.add_argument('--chroms', type=argparse.FileType('rU'),
                        help="File with chromosome names. If tab-delimited, names must be in first column")                        
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help="Input GTF file")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF file")
    main(parser.parse_args())
