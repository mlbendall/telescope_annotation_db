#! /bin/bash
import sys
import re
from collections import defaultdict, Counter
from itertools import groupby
from subprocess import Popen,PIPE

def tab_line_gen(infile):
    ''' Returns a generator for tab-delmited file '''
    return (l.strip('\n').split('\t') for l in infile if not l.startswith('#'))

def simplify_list(l):
    """ Collapse consecutive repeats in list """
    return [k for k, g in groupby(l)]

def sort_gtf(liter, chroms=None):
    ''' Sort GTF file '''
    alllines = defaultdict(list)  
    for l in liter:
        alllines[l[0]].append(l)
    
    if chroms is None:
        chroms = sorted(alllines.keys())
    
    for cchrom in chroms:
        if cchrom not in alllines: continue 
        for l in sorted(alllines[cchrom],key=lambda x:int(x[3])):
            yield l


# GTFCOLS = ['chrom','source','feature','start','end','score','strand','frame']

class GTFLine:
    # Columns in a GTF line
    GTFCOLS = ['chrom','source','feature','start','end','score','strand','frame']
    # Attributes that should be at the front of attribute string
    ATTRORDER = ['gene_id', 'transcript_id', 'locus', 'repName']
    
    def __init__(self,row):
        for c,v in zip(self.GTFCOLS, row[:8]):
            try:
                setattr(self,c,int(v))
            except ValueError:
                setattr(self,c,v)
        _attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',row[8]))
        self.attr = {}
        for k,v in _attrd.iteritems():
            try:
                self.attr[k] = int(v)
            except ValueError:
                try:
                    self.attr[k] = float(v)
                except ValueError:
                    self.attr[k] = v
    
    def fmt_attrs(self):
        ret = ['%s "%s";' % (k,self.attr[k]) for k in self.ATTRORDER if k in self.attr]
        ret += ['%s "%s";' % (k,v) for k,v in self.attr.iteritems() if k not in self.ATTRORDER]
        return ' '.join(ret)
    
    def fmt(self):
        _attrs = self.fmt_attrs() # ' '.join('%s "%s";' % (k,v) for k,v in self.attr.iteritems())
        return [str(getattr(self,c)) for c in self.GTFCOLS] + [_attrs]
    
    def sequence(self, gdict):
        if self.strand == '+':
            return gdict[self.chrom][self.start:self.end]
        else:
            return gdict[self.chrom][self.start:self.end].reverse_complement()
    
    def asdict(self):
        ret = {k:v for k,v in self.attr.iteritems()}
        ret.update({k:getattr(self,k) for k in self.GTFCOLS})
        return ret
    
    def __str__(self):
        return '\t'.join(self.fmt())

### Manipulate locus #####################################################################

def remove_dups(_locus):
    ''' Remove annotations that are duplicated within locus '''
    tmp = {a.attr['id']:a for a in _locus}
    return tmp.values()

def create_spanning(_locus, feat='span'):
    ''' Return GTF line that spans all annotations in _locus '''
    spanning = GTFLine(_locus[0].fmt())
    spanning.feature = feat
    spanning.score = '.'
    spanning.attr = {}
    spanning.start = min(a.start for a in _locus)
    spanning.end = max(a.end for a in _locus)
    return spanning

def adjust_overlaps(_locus, strand=None):
    ''' Adjust annotations within locus to eliminate overlaps
        For overlapping annotations on the plus strand, the start position of the second
        annotation is adjusted so that it begins after the first annotation ends.
        The same is true for the minus strand, except we are adjusting the end position
        and annotations are considered from right to left.
    '''
    if strand is None: strand = _locus[0].strand
    if strand == '+':
        _locus.sort(key=lambda x:x.start)
        for i,p1 in enumerate(_locus[1:]):
            p0 = _locus[i]
            gdist = p1.start - p0.end
            if gdist <= 0:
                p1.start = p1.start - gdist + 1
                p1.attr['repStart'] = p1.attr['repStart'] - gdist
    else:
        _locus.sort(key=lambda x:x.end, reverse=True)
        for i,p1 in enumerate(_locus[1:]):
            p0 = _locus[i]
            gdist = p0.start - p1.end
            if gdist <= 0:
                p1.end = p1.end + gdist - 1
                p1.attr['repLeft'] = p1.attr['repLeft'] - gdist
    
    return _locus

def guess_rmsk_model_lengths(gtf):
    ''' Calculate the length of models from attributes '''
    ret = defaultdict(list)
    for g in gtf:
        if 'repName' in g.attr:
            sz = (g.attr['repEnd'] - g.attr['repLeft']) if g.strand == '+' else (g.attr['repEnd'] - g.attr['repStart'])
            ret[g.attr['repName']].append(sz)
    return {k:Counter(v).most_common()[0][0] for k,v in ret.iteritems()}

def get_span(_locus):
    spn = [a for a in _locus if a.feature.startswith('span')]
    assert len(spn) == 1
    return spn[0]

def add_exon_number(_locus):
    strand = set([a.strand for a in _locus])
    if len(strand) == 1 and '-' in strand:
        tmp = sorted(_locus, key=lambda x:x.end, reverse=True)

def raw_input_stderr(*args):
    sys.stdout = sys.stderr
    x = raw_input(*args)
    sys.stdout = sys.__stdout__
    return x    

def find_overlaps(gtf):
    ''' Find overlapping annotations within GTF '''
    # Cluster the "span" lines using bedtools cluster
    p1 = Popen('bedtools cluster -i -', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    o,e = p1.communicate(input='\n'.join(str(g) for g in gtf))
    # Parse bedtools output
    _clust = defaultdict(list)
    for l in tab_line_gen(o.strip('\n').split('\n')):
        _clust[l[-1]].append(GTFLine(l[:9]))
    return _clust

def groupstr(k,cgroup):
    ret = "### Overlap group %s ###\n" % k
    for spn in cgroup:
        ret += '### %s\t%d\t%s\n' % (spn.attr['locus'], spn.attr['model_cov'], spn.attr['category'])
    for spn in cgroup:
        ret += '%s\n' % spn
    return ret

PROMPT_HELP = '''Options are:
    ignore                    - Do nothing, do not change annotations
    reject loc1[,loc2,...]    - Remove locs from annotations
    diff loc1-loc2            - Shorten loc1 to eliminate overlap with loc2
    merge loc1+loc2           - Combine loc1 and loc2
'''    

def prompt_cmd():
    ''' Get command from user input '''
    operation = 'begin'
    while operation.lower() not in ['ignore','reject','diff','merge']:
        curprompt = ''
        if operation not in ['begin','?']:
            curprompt += 'INPUT IS INVALID'.ljust(20)
        else:
            curprompt += ' '.ljust(20)
        curprompt += '***Action to take (? for help): '
        z = raw_input_stderr(curprompt).strip()
        if z == '?':
            print >>sys.stderr, PROMPT_HELP
        inputcmd = z.strip().split()
        if len(inputcmd) == 1 or len(inputcmd) == 2:
            operation = inputcmd[0]
        else:
            operation = ''
    
    assert operation.lower() in ['ignore','reject','diff','merge']
    return inputcmd
