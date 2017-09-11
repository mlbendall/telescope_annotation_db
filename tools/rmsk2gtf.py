#! /usr/bin/env python


def convert(infile, outfile, gene_region, source, chroms):
    from collections import defaultdict
    import math
    
    alllines = defaultdict(list)  
    lines = (l.strip('\n').split('\t') for l in infile)
    header = lines.next()

    # Find column indices
    cidx = header.index('genoName')
    sidx = header.index('genoStart')
    eidx = header.index('genoEnd')
    stidx = header.index('strand')
    
    # Organize original lines by chromosome
    for l in lines:
        alllines[ l[cidx] ].append(l)
    
    # Chroms was not provided
    if chroms is None:
        chroms = sorted(alllines.keys())

    # Sorting and formatting lines
    counter = 1
    for cchrom in chroms:
        if cchrom not in alllines: continue 
        for l in sorted(alllines[cchrom],key=lambda x:int(x[sidx])):
            # Add 1 to start because GTF is 1-based
            spos = '%d' % (int(l[sidx]) + 1)
            # Don't change end because GTF is inclusive
            epos = l[eidx]
            assert int(spos)<int(epos), "ERROR start >= end %s" % l      
            score = l[header.index('swScore')]
            
            # Attributes
            attrs = {'repStart': l[header.index('repStart')],
                     'repEnd': l[header.index('repEnd')],
                     'repLeft': l[header.index('repLeft')],
                     'id': '%s_%d' % (l[header.index('repName')], counter),
                     'repName': l[header.index('repName')],
                     'repClass': l[header.index('repClass')],
                     'repFamily': l[header.index('repFamily')],
                     }
            if gene_region is not None:
                attrs['geneRegion'] = gene_region
            attr = ' '.join('%s "%s";' % (k,v) for k,v in attrs.iteritems())
            print >>outfile, '\t'.join([cchrom,source,'exon',spos,epos,score,l[stidx],'.',attr])
            counter += 1

def main(parser):
    args = parser.parse_args()
    if args.chroms:
        chroms = [l.strip('\n').split('\t')[0] for l in args.chroms]
    else:
        chroms = None
    
    convert(args.infile, args.outfile, args.gene_region, args.source, chroms)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Convert UCSC RepeatMasker table to GTF')
    parser.add_argument('--source', default='rmsk',
                        help="Value for source column of GTF file")
    parser.add_argument('--chroms', type=argparse.FileType('rU'),
                        help="File with chromosome names. If tab-delimited, names must be in first column")    
    parser.add_argument('--gene_region',
                        help='Name for gene region, i.e. internal, ltr, 5end, orf2')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help="Input UCSC RepeatMasker table")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF file")
    main(parser)
