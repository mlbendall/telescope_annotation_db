#! /usr/bin/env python

from collections import defaultdict
import utils
from utils import GTFLine

def calculate_internal_coverage(_locus):
    ''' Calculate number of internal model bases covered in locus
        This is the number of "query" bases represented, not reference bases.
    '''
    if _locus[0].strand == '+':
        return sum(a.attr['repEnd'] - a.attr['repStart'] for a in _locus if a.attr['geneRegion'] == 'internal')
    else:
        return sum(a.attr['repEnd'] - a.attr['repLeft'] for a in _locus if a.attr['geneRegion'] == 'internal')

def get_internal_model(_locus):
    ''' Return the name of the internal model for locus '''
    im = set(a.attr['repName'] for a in _locus if a.attr['geneRegion']=='internal')
    assert len(im) == 1
    return im.pop()

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
    combined_gtf = [GTFLine(l) for l in utils.tab_line_gen(args.internalGTF)]
    combined_gtf = [g for g in combined_gtf if not g.feature.startswith('span')]
    
    # Load the LTR GTF
    for l in utils.tab_line_gen(args.ltrGTF):
        l_b = GTFLine(l[-10:])
        l_b.attr['locus'] = GTFLine(l[:9]).attr['locus']
        combined_gtf.append(l_b)

    model_lengths = utils.guess_rmsk_model_lengths(combined_gtf)
    
    # Organize by locus    
    byloc = defaultdict(list)
    for g in combined_gtf:
        byloc[g.attr['locus']].append(g)
    
    for locid,locus in byloc.iteritems(): 
        # Remove duplicate annotations
        locus = utils.remove_dups(locus)
        # Adjust overlaps
        locus = utils.adjust_overlaps(locus)
        # Determine category
        category = locus_category(locus)
        # Add information to all annotations
        strand = set([a.strand for a in locus])
        if len(strand) == 1 and '-' in locus:
            locus.sort(key=lambda x:x.end, reverse=True)
        else:
            locus.sort(key=lambda x:x.start)
        for i,a in enumerate(locus):
            a.source = category
            a.attr['exon_number'] = i+1
        
        # Create spanning annotation
        spanning = utils.create_spanning(locus)        
        # Calculate model coverage and percent
        internal_model = get_internal_model(locus)
        model_cov = calculate_internal_coverage(locus)
        model_pct = min(100, (float(model_cov) / model_lengths[internal_model])*100)
        spanning.attr = {'locus': locid ,
                         'category': category,
                         'model_cov': model_cov,
                         'model_pct': '%.1f' % model_pct,
                         'exons':len(locus)
                         }
        
        print >>args.outfile, '### %s ###' % locid
        print >>args.outfile, spanning
        print >>args.outfile, '\n'.join(str(_) for _ in sorted(locus,key=lambda x:x.start))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Assemble HERV loci from internal and LTR annotations')
    parser.add_argument('internalGTF', type=argparse.FileType('rU'),
                        help="Merged hits for internal HERV regions as GTF file.")
    parser.add_argument('ltrGTF', type=argparse.FileType('rU'),
                        help='''LTR hits that overlap (flanked) internal regions. This
                                should be output from bedtools intersect using the -wo 
                                option, so overlapping records are on the same line.''')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())
