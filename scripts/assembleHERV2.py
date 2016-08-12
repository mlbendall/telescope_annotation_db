#! /usr/bin/env python

from collections import defaultdict
import itertools
from utils import GTFLine, tab_line_gen
import utils

def locus_category(_locus):
    ''' Determine the category of locus
        Check whether locus begins and/or ends with LTR.
    '''
    if _locus[0].strand == '+':
        _locus.sort(key=lambda x:x.start)
    else:
        _locus.sort(key=lambda x:x.end, reverse=True)
    # Determine label
    regions = utils.simplify_list(a.attr['geneRegion'] for a in _locus)
    if regions[0] != 'ltr' and regions[-1] != 'ltr':
        return 'internal'
    elif regions[0] == 'ltr' and regions[-1] == 'ltr':
        return 'prototype'
    elif regions[0] == 'ltr' or regions[-1] == 'ltr':
        return 'oneside'
    else:
        return 'unknown'

def main(args):
    # Load the internal GTF
    combined_gtf = [GTFLine(l) for l in tab_line_gen(args.internalGTF)]
    combined_gtf = [g for g in combined_gtf if not g.source=='span_internal']
    
    # Load the left GTF
    for l in tab_line_gen(args.ltrGTF):
        l_b = GTFLine(l[10:])
        l_b.attr['locus'] = GTFLine(l[:9]).attr['locus']
        combined_gtf.append(l_b)

    byloc = defaultdict(list)
    for g in combined_gtf:
        byloc[g.attr['locus']].append(g)
    
    #for locid,ilocus in itertools.groupby(combined_gtf, key=lambda x:x.attr['locus']):
    for locid,locus in byloc.iteritems(): 
        locus = utils.remove_dups(locus)
        locus = utils.adjust_overlaps(locus)
        category = locus_category(locus)
        for a in locus:
            a.source = category
        spanning = utils.create_spanning(locus)
        spanning.attr = {'category': category,
                         'locus': locid }
        print >>args.outfile, '### %s ###' % locid
        print >>args.outfile, spanning
        print >>args.outfile, '\n'.join(str(_) for _ in sorted(locus,key=lambda x:x.start))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Assemble HERV loci from internal and LTR annotations',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('internalGTF', type=argparse.FileType('rU'),
                        help="Merged GTF for internal regions")
    parser.add_argument('ltrGTF', type=argparse.FileType('rU'),
                        help="Overlapping LTR regions")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())
