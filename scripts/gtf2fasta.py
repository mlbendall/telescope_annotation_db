#! /usr/bin/env python
import sys
from collections import defaultdict
from Bio import SeqIO

import utils
from utils import GTFLine

cmpl = {x:y for x,y in zip('ACGTacgt','TGCAtgca')}
def revcomp(s):
    return ''.join(cmpl.get(_,_) for _ in s[::-1])

def uppercase_substring(s, b, e):
    return s[:b] + s[b:e].upper() + s[e:]

def merge_interval_list(ivs, dist=0):
    """ Merge intervals

    Args:
        ivs (list): List of intervals. Each interval is represented by a tuple of
            integers (start, end) where end > start.
        dist (int): Distance between intervals to be merged. Setting dist=1 will merge
            adjacent intervals

    Returns:
        list: Merged list of intervals

    Examples:
        >>> merge_interval_list([])
        []
        >>> merge_interval_list([(1,10)])
        [(1, 10)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)])
        [(1, 3), (4, 9), (10, 14)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)], dist=1)
        [(1, 14)]
    """
    if len(ivs)<= 1: return ivs
    ivs.sort(key=lambda x:x[0])
    ret = [ivs[0]]
    for iv in ivs[1:]:
        if iv[0] - ret[-1][1] > dist:
            ret.append(iv)
        else:
           ret[-1] = (ret[-1][0], max(iv[1],ret[-1][1]))
    return ret

def wrap_fasta(s, ll=100):
    return '\n'.join(s[i:(i+ll)] for i in range(0,len(s),ll))

def main(args):
    # Load the annotation
    byloc = defaultdict(list)
    gtf = (GTFLine(l) for l in utils.tab_line_gen(args.gtf))
    for g in gtf:
        byloc[g.attr['locus']].append(g)
    
    print >>sys.stderr, "Loaded GTF"
    
    # Load the sequence
    seqs = {s.id:s for s in SeqIO.parse(args.ref, 'fasta')}
    print >>sys.stderr, "Loaded Seqs"
    
    # Spanning GTf
    span_gtf = (GTFLine(l) for l in utils.tab_line_gen(args.span))
    for i,sg in enumerate(span_gtf):
        # Get the full spanning sequence (lowercase)
        seqstr = str(seqs[sg.chrom][sg.start:(sg.end+1)].seq).lower()
        # Get relative locations of "exons"
        rel = [(g.start-sg.start, g.end-sg.start+1) for g in byloc[sg.attr['locus']]]
        # Merge exons if they are close together
        rel = merge_interval_list(rel, 10)
        # Change exons to uppercase
        for b,e in rel:
            seqstr = uppercase_substring(seqstr, b, e)
        # Reverse complement if it is on minus strand
        if sg.strand == '-':
            seqstr = revcomp(seqstr)
        print >>args.outfile, '>%s\n%s' % (sg.attr['locus'], wrap_fasta(seqstr))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Create FASTA from GTF.')
    parser.add_argument('gtf', type=argparse.FileType('rU'),
                        help="Input GTF file")
    parser.add_argument('span', type=argparse.FileType('rU'),
                        help="Input GTF span file")
    parser.add_argument('ref', type=argparse.FileType('rU'),
                        help="Reference FASTA file")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output text file")
    main(parser.parse_args())