#! /usr/bin/env python

import os, sys
from collections import defaultdict, Counter
import json
from itertools import chain

import utils
from utils import GTFLine
from assembleHERV import calculate_internal_coverage, get_internal_model, locus_category

def resolve_reject(cmd, cgroup, fulllocs):
    ''' Resolve conflict by rejecting (removing) loci
        If no argument is provided to 'reject', the locus with the longest total model
        coverage is retained and all others are removed. Otherwise, the loci listed are
        removed and all others are retained.
    '''
    if len(cmd) == 1:
        tmp = sorted(cgroup, key=lambda x:x.attr['model_cov'], reverse=True)
        rem = [a.attr['locus'] for a in tmp[1:]]
    else:
        rem = cmd[1].split(',')
        if not all(r in fulllocs for r in rem):
            sys.exit('ERROR: Command does not match with conflicts found:\n%s' % ' '.join(cmd))
    for r in rem: del fulllocs[r]
    return fulllocs

def resolve_merge(cmd, cgroup, fulllocs):
    ''' Resolve conflict by merging loci
        Take the union of annotations from conflicting loci (without spanning annotations),
        create a new spanning annotation. If no argument is provided, default is to merge
        all loci. If an argument is provided, it must be delimited by '+' and loci that 
        are not listed are removed.
    '''
    if len(cmd) == 1:
        mer = [a.attr['locus'] for a in cgroup]
    else:
        mer = cmd[1].split('+')
        if not all(r in fulllocs for r in mer):
            sys.exit('ERROR: Command does not match with conflicts found:\n%s' % ' '.join(cmd))

    # has_inversion = False    
    newloc_id = mer[0]
    newloc = [a for a in chain.from_iterable(fulllocs[m] for m in mer) if not a.feature.startswith('span')]
    newstrand = Counter(a.strand for a in newloc).most_common()[0][0]    
    # Calculate model information
    model_lengths = utils.guess_rmsk_model_lengths(newloc)
    internal_model = get_internal_model(newloc)
    # Remove duplicates, adjust overlaps
    newloc = utils.remove_dups(newloc)
    newloc = utils.adjust_overlaps(newloc, newstrand)
    # Adjust attributes for annotations
    if newstrand == '+':
        newloc.sort(key=lambda x:x.start)
    else:
        newloc.sort(key=lambda x:x.end, reverse=True)

    newcat = 'merged*' #  % locus_category(newloc)
    for i,a in enumerate(newloc):
        if a.strand != newstrand: a.strand = '.' # has_inversion = True
        a.source = newcat #'merged'
        a.attr.update({'locus':newloc_id, 'exon_number':(i+1)})
    
    # Create spanning annotation
    spn = utils.create_spanning(newloc)
    spn.strand = newstrand
    model_cov = calculate_internal_coverage(newloc)
    model_pct = min(100, (float(model_cov) / model_lengths[internal_model])*100)
    spn.attr = {'locus': newloc_id,
                'category': newcat, # 'merged',
                'model_cov': model_cov,
                'model_pct': '%.1f' % model_pct,                
                # 'repName': internal_model,
                'exons': len(newloc),
                'conflict': 'merged',
               }
    # if has_inversion:
    #     spn.attr['inversion'] = 'true'
    return { newloc_id : [spn] + newloc }


def resolve_diff(cmd, cgroup, fulllocs):
    ''' Resolve conflict by subtracting loci
        The overlapping portion is removed from A and the remaining portion of A is 
        reported. B is reported without changes.
    '''
    assert len(cgroup) == 2, "Only two loci allowed for diff"
    if len(cmd) == 1:
        tmp = sorted(cgroup, key=lambda x:x.attr['model_cov'], reverse=True)
        dA = tmp[0]
        idA = dA.attr['locus']        
        dB = tmp[1]
        idB = dB.attr['locus']
    else:
        idA,idB = cmd[1].split('-')
        dA = [c for c in cgroup if c.attr['locus'] == idA]
        dB = [c for c in cgroup if c.attr['locus'] == idB]
        if not len(dA) == len(dB) and len(dA) == 1:
            sys.exit('ERROR: Command does not match with conflicts found:\n%s' % ' '.join(cmd))
        dA = dA[0]
        dB = dB[0]

    newloc = sorted(fulllocs[idA], key=lambda x:x.start)
    newloc = [loc for loc in newloc if not loc.feature.startswith('span')]
    model_lengths = utils.guess_rmsk_model_lengths(newloc)
    internal_model = get_internal_model(newloc)
    
    trimside = None
    # We are trimming dA. Find where dA overlaps dB
    if dA.start < dA.end <= dB.start or dA.end > dA.start >= dB.end:
        if dA.start - dB.end == 1:
            trimside = 'left'
        elif dB.start - dA.end == 1:
            trimside = 'right'
        else:
            assert False, "No overlap between %s and %s" % (idA,idB)
    elif dB.start <= dA.start < dA.end <= dB.end:
        # dA contained in dB. Remove dA
        print >>sys.stderr, "%s contained within %s. Removing %s." % (idA, idB, idA)
        return { idB:fulllocs[idB] }
    elif dA.start <= dB.start < dB.end <= dA.end:
        # dB contained in dA
        trimside = 'left' if (dB.end - dA.start) < (dA.end - dB.start) else 'right'
    elif dB.start <= dA.start <= dB.end:
        trimside = 'left'
    elif dB.start <= dA.end <= dB.end:
        trimside = 'right'        

    assert trimside is not None    

    if trimside == 'left':
        print >>sys.stderr, "\tTrimming from left side of %s" % idA
        for loc in newloc:
            gdist = loc.start - dB.end
            if gdist <= 0:
                loc.start = dB.end + 2
                if loc.strand == '+':
                    loc.attr['repStart'] -= gdist
                else:
                    loc.attr['repEnd'] += gdist
        newloc = [loc for loc in newloc if loc.end - loc.start > 5]
    elif trimside == 'right':
        print >>sys.stderr, "\tTrimming from right side of %s" % idA
        for loc in newloc:
            gdist = dB.start - loc.end
            if gdist <= 0:
                loc.end = dB.start - 2
                if loc.strand == '+':
                    loc.attr['repEnd'] += gdist
                else:
                    loc.attr['repLeft'] -= gdist
        newloc = [loc for loc in newloc if loc.end - loc.start > 5]
    
    newcat = '%s*' % locus_category(newloc)
    for loc in newloc:
        loc.source = newcat

    # Create spanning annotation
    spn = utils.create_spanning(newloc)
    model_cov = calculate_internal_coverage(newloc)
    model_pct = min(100, (float(model_cov) / model_lengths[internal_model])*100)
    spn.attr = {'locus': idA,
                'category': newcat,
                'model_cov': model_cov,
                'model_pct': '%.1f' % model_pct,
                # 'repName': internal_model,
                'exons': len(newloc),
                'conflict': 'diff',                
               }
    
    fulllocs[idA] = [spn] + newloc
    return fulllocs

def resolve(cmd, cgroup, fulllocs):
    if cmd[0] == 'ignore':
        return fulllocs
    elif cmd[0] == 'reject':
        return resolve_reject(cmd, cgroup, fulllocs)
    elif cmd[0] == 'diff':
        return resolve_diff(cmd, cgroup, fulllocs)
    elif cmd[0] == 'merge':
        return resolve_merge(cmd, cgroup, fulllocs)

def main(args):
    combined_gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]

    # Find overlaps within the span features
    overlap_groups = utils.find_overlaps([g for g in combined_gtf if g.feature.startswith('span')])
    overlap_groups = {k:v for k,v in overlap_groups.iteritems() if len(v) > 1}
    print >>sys.stderr, "Found %d groups with conflict." % len(overlap_groups)    
    
    # Resolve commands
    if args.resolve_file is not None:
        resolve_cmds = json.load(args.resolve_file)
    else:
        if args.resolve is not None:
            # Resolve commands were provided as command-line argument
            resolve_cmds = args.resolve
        else:
            # Prompt user to enter resolve commands
            resolve_cmds = {}
            for groupid in sorted(overlap_groups.keys(), key=lambda x:int(x)):
                ogroup = overlap_groups[groupid]
                print >>sys.stderr, utils.groupstr(groupid, ogroup)                        
                concmd = utils.prompt_cmd()
                resolve_cmds[groupid] = concmd
        
        print >>sys.stderr, 'Commands for resolving conflicts (JSON):\n'
        print >>sys.stderr, json.dumps(resolve_cmds)

    # Resolve the conflicts
    byloc = defaultdict(list)
    for g in combined_gtf:
        byloc[g.attr['locus']].append(g)
    
    for groupid, ogroup in overlap_groups.iteritems():
        print >>sys.stderr, utils.groupstr(groupid, ogroup)    
        locids = [a.attr['locus'] for a in ogroup]
        # Temporarily remove the loci
        fulllocs = {}
        for locid in locids:
            print >>sys.stderr, '\tRemoving %s' % locid      
            fulllocs[locid] = byloc.pop(locid)
        assert groupid in resolve_cmds
        newlocs = resolve(resolve_cmds[groupid], ogroup, fulllocs)
        for newlocid, newlocus in newlocs.iteritems():
            print >>sys.stderr, '\tInserting %s' % newlocid
            # print >>sys.stderr, '\n'.join(str(_) for _ in newlocus)
        byloc.update(newlocs)

    for locid,locus in byloc.iteritems():
        print >>args.outfile, '### %s ###' % locid
        print >>args.outfile, '\n'.join(str(_) for _ in sorted(locus,key=lambda x:x.start))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Polish HERV loci',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--resolve', type=json.loads,
                        help='''String with resolve commands, in JSON format''')
    parser.add_argument('--resolve_file', type=argparse.FileType('rU'),
                        help='''File containing resolve commands, in JSON format''')    
    parser.add_argument('infile', type=argparse.FileType('rU'),
                        help="GTF with HERV loci")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())
