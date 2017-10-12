#! /usr/bin/env python
import sys
from collections import defaultdict, Counter
import utils
from utils import GTFLine

def correct_rmsk_model_coordinates(gtflines, modellen):
    ''' Fix incorrect model coordinates using model length '''
    for g in gtflines:
        if g.strand == '+':
            trueend = modellen[g.attr['repName']] + g.attr['repLeft']
            if trueend != g.attr['repEnd']:
                replen = g.attr['repEnd'] - g.attr['repStart']
                g.attr['repEnd'] = trueend
                g.attr['repStart'] = trueend - replen
        else:
            trueend = modellen[g.attr['repName']] + g.attr['repStart']
            if trueend != g.attr['repEnd']:
                replen = g.attr['repEnd'] - g.attr['repLeft']
                g.attr['repEnd'] = trueend
                g.attr['repLeft'] = trueend - replen

def main(args):
    ### Read the GTF file ################################################################
    gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]
    
    ### Correct model coordinates ########################################################
    # The repStart, repEnd, and repLeft attributes downloaded from the UCSC rmsk database
    # does not always give the same model length. Here we guess what the correct model
    # length is then correct each record    
    mlen = utils.guess_rmsk_model_lengths(gtf)
    print >>sys.stderr, 'Model lengths:'
    print >>sys.stderr, '\n'.join('%s%d' % (k.ljust(16), mlen[k]) for k in sorted(mlen.keys()))
    correct_rmsk_model_coordinates(gtf,mlen)
    # Check that model coordinates are correct
    for g in gtf:
        if g.strand == '+':
            trueend = mlen[g.attr['repName']] + g.attr['repLeft']
        else:
            trueend = mlen[g.attr['repName']] + g.attr['repStart']
        assert trueend == g.attr['repEnd']
        
        print >>args.outfile, g

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Correct model coordinates of GTF converted from UCSC RepeatMasker tables')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help="Input GTF with UCSC RepeatMasker records")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF file")
    main(parser.parse_args())
