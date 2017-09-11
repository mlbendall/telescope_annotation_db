#! /usr/bin/env python

from collections import defaultdict, Counter
import itertools
import utils
from utils import GTFLine

def main(args):
    # Filehandle for rejected loci
    if args.reject_gtf is None:
        import os
        args.reject_gtf = open(os.devnull,'w')

    # Filtering parameters
    min_internal_pct = args.min_internal_pct * 100
    min_internal_bases = args.min_internal_bases
    
    # Load the internal GTF
    combined_gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]
    model_lengths = utils.guess_rmsk_model_lengths(combined_gtf)
    
    loccounts = Counter()
    byloc = defaultdict(list)
    for g in combined_gtf:
        byloc[g.attr['locus']].append(g)
    
    if min_internal_pct > 0:
        print >>sys.stderr, "Removing loci matching less than %d percent of internal model..." % int(min_internal_pct)
    if min_internal_bases > 0:
        print >>sys.stderr, "Removing loci matching less than %d internal bases..." % min_internal_bases
    
    rejectflag = False
    for locid,locus in byloc.iteritems():
        spn = utils.get_span(locus)
        locus  = [a for a in locus if a != spn] # not a.feature.startswith('span')]
        category = spn.attr['category']
        model_pct = spn.attr['model_pct']
        model_cov = spn.attr['model_cov']
        
        if model_pct >= min_internal_pct and model_cov >= min_internal_bases:
            loccounts[category] += 1
            outh = args.outfile
        else:
            if not rejectflag:
                print >>sys.stderr, 'Removed loci:'
                print >>sys.stderr, '%-18s%-6s%-6s%s' % ('locus','bp','pct','category')
                rejectflag = True
            loccounts['rejected'] += 1
            print >>sys.stderr, '%-18s%-6d%-6.1f%s' % (locid, model_cov, model_pct, category)
            outh = args.reject_gtf
        
        print >>outh, '### %s ###' % locid
        print >>outh, spn
        print >>outh, '\n'.join(str(_) for _ in sorted(locus,key=lambda x:x.start))            
    
    if not rejectflag:
        print >>sys.stderr, 'All passed filter.'

    print >>sys.stderr, 'Summary:'
    for cat in ['internal','prototype','oneside','unusual','rejected']:
        print >>sys.stderr, '%s%d' % (cat.ljust(20), loccounts[cat])
    
if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Filter HERV loci according to match length',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--min_internal_bases', type=int, default=0,
                        help="Minimum number of bases matching internal model")
    parser.add_argument('--min_internal_pct', type=float, default=0.,
                        help="Minimum percentage of internal model matched")
    parser.add_argument('--reject_gtf', type=argparse.FileType('w'),
                        help="File to output loci failing filters")
    parser.add_argument('infile', type=argparse.FileType('rU'),
                        help="GTF with HERV loci")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())