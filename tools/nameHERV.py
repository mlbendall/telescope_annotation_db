#! /usr/bin/env python
from subprocess import Popen, PIPE
from string import letters as someletters
from collections import defaultdict, Counter

import utils
from utils import GTFLine

def main(args):
    combined_gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]
    if args.prefix is None:
        prefix = Counter(g.attr['locus'].split('_')[0] for g in combined_gtf).most_common()[0][0]
    else:
        prefix = args.prefix

    namemap = {}
    if args.cytoband:
        byband = defaultdict(list)
        p1 = Popen('bedtools intersect -wo -a - -b %s' % args.cytoband, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        o,e = p1.communicate(input='\n'.join(str(g) for g in combined_gtf if g.feature.startswith('span')))
        for l in utils.tab_line_gen(o.strip('\n').split('\n')):
            g1 = GTFLine(l[:9])
            g2 = GTFLine(l[9:-1])
            band = '%s%s' % (g2.chrom.strip('chr'),g2.attr['gene_id'])
            byband[band].append(g1)

        for band, locs in byband.iteritems():
            if len(locs) == 1:
                namemap[locs[0].attr['locus']] = '%s_%s' % (prefix, band)
            else:
                locs.sort(key=lambda x:x.start)
                for i,loc in enumerate(locs):
                    namemap[loc.attr['locus']] = '%s_%s%s' % (prefix, band, someletters[i])
    else:
        for g in combined_gtf:
            namemap[g.attr['locus']] = g.attr['locus']
    
    for g in combined_gtf:
        g.attr['oLocus'] = g.attr['locus']
        g.attr['locus'] = namemap[g.attr['locus']]
        print g

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Filter HERV loci according to match length',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--cytoband',
                        help="File name containing cytoband information.")
    parser.add_argument('--prefix',
                        help="Prefix for locus names")
    parser.add_argument('infile', type=argparse.FileType('rU'),
                        help="GTF with HERV loci")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())