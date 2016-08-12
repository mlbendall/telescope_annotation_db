#! /usr/bin/env python

import utils
from utils import GTFLine
from collections import defaultdict


def main(args):
    ### Read the GTF file ################################################################
    gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]
    
    bychrom = defaultdict(lambda:{'+':list(), '-':list()})
    for g in gtf:
        bychrom[g.chrom][g.strand].append(g)

    merged_hits = []
    for cchrom, strands in bychrom.iteritems():
        # Plus strand
        if len(strands['+']):
            strands['+'].sort(key=lambda x:x.start)
            cur = [ strands['+'][0] ]
            for g1 in strands['+'][1:]:
                g0 = cur[-1]
                # Genomic distance between hits
                gdist = g1.start - g0.end            
                if gdist <= args.shortdist:
                    domerge = True
                else:
                    domerge = g0.attr['repLeft'] < g1.attr['repLeft']
                    domerge &= gdist < args.longdist
                if domerge:
                    cur.append(g1)
                else:
                    merged_hits.append(cur)
                    cur = [ g1 ]
            merged_hits.append(cur)

        # Minus strand
        if len(strands['-']):
            strands['-'].sort(key=lambda x:x.end, reverse=True)
            cur = [ strands['-'][0] ]
            for g1 in strands['-'][1:]:
                g0 = cur[-1]
                # Genomic distance between hits
                gdist = g0.start - g1.end
                if gdist <= args.shortdist:
                    domerge = True
                else:
                    domerge = g0.attr['repStart'] < g1.attr['repStart']
                    domerge &= gdist < args.longdist
                if domerge:
                    cur.append(g1)
                else:
                    merged_hits.append(cur)
                    cur = [ g1 ]
            merged_hits.append(cur)
    
    for i,cur in enumerate(merged_hits):
        locid = '%s_%04d' % (args.prefix, i+1)
        spanning = utils.create_spanning(cur)
        spanning.attr['locus'] = locid
        for g in cur:
            g.attr['locus'] = locid
        
        print >>args.outfile, '### %s ###' % locid 
        print >>args.outfile, spanning
        print >>args.outfile, '\n'.join(str(_) for _ in sorted(cur,key=lambda x:x.start))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Merge hits in GTF file belonging to the same model',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--shortdist', type=int, default=10, 
                        help='''A short distance. Annotations with genomic distance of
                                <shortdist> or less are merged automatically.''')    
    parser.add_argument('--longdist', type=int, default=10000, 
                        help='''An extreme distance. Annotations with genomic distance of
                                <longdist> or less are checked for consistency. If the
                                model coordinates agree, annotations are merged.''')
    parser.add_argument('--prefix',
                        help='''Prefix for locus''')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help="Input GTF")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())
