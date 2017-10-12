#! /usr/bin/env python

import os, sys
# from collections import defaultdict, Counter
from collections import OrderedDict
import json

import utils
from utils import GTFLine
from IGV import *

# from polishHERVLoci import find_conflicts, groupstr, prompt_cmd


def get_display_str(_locus, padding=0):
    spn = utils.get_span(_locus)
    return '%s:%d-%d' % (spn.chrom, spn.start - padding, spn.end + padding)

def main(args):
    # Setup snapshot directory
    if args.snapshot:
        snapdir = os.path.abspath(args.snapshot_dir)
        print >>sys.stderr, "Snapshot directory: %s" % snapdir
        if not os.path.exists(snapdir):
            os.makedirs(snapdir)
    
    # Create IGV object
    igv = IGV()
    igv.setSleepInterval(1)
    if args.clear:
        igv.new()
    
    # Load reference
    igv.genome(args.refbuild)

    # Load annotation
    input_path = os.path.abspath(args.infile.name)
    igv.load(input_path)

    # Load compare annotations
    if args.compare_tracks:
        for fn in args.compare_tracks.split(','):
            fp = os.path.abspath(fn)
            if os.path.isfile(fp):
                igv.load(fp)

    # Expand display
    igv.expand()
    if not args.interactive:
        igv.setSleepInterval(args.sleep_interval)
    
    # Load the GTF
    combined_gtf = [GTFLine(l) for l in utils.tab_line_gen(args.infile)]
    
    if args.inspect == 'all':
        byloc = OrderedDict()
        for g in combined_gtf:
            if g.attr['locus'] not in byloc:
                byloc[g.attr['locus']] = list()
            byloc[g.attr['locus']].append(g)
           
        for locid,locus in byloc.iteritems():
            spn = utils.get_span(locus)
            category = spn.attr['category'] if 'category' in spn.attr else None
            ds = get_display_str(locus)
            ds_pad = get_display_str(locus, args.padding)
            
            if not args.interactive:
                print >>sys.stderr, '%s%s%s ' % (locid.ljust(22), category.ljust(12), ds)            

            # Go to locus
            igv.goto(ds_pad)
            
            # Snapshot 
            if args.snapshot:
                if not os.path.exists(os.path.join(snapdir, category.strip('*'))):
                    os.mkdir(os.path.join(snapdir, category.strip('*')))
                igv.snapshotDirectory(os.path.join(snapdir, category.strip('*')))
                igv.snapshot('%s.png' % locid)

            # Pause if interactive
            if args.interactive:
                z = utils.raw_input_stderr('%s%s%s ' % (locid.ljust(22), category.ljust(12), ds)).strip()
        
    elif args.inspect == 'overlap':
        # Find overlaps within the span features
        overlap_groups = utils.find_overlaps([g for g in combined_gtf if g.feature.startswith('span')])
        overlap_groups = {k:v for k,v in overlap_groups.iteritems() if len(v) > 1}
        print >>sys.stderr, "Found %d groups with conflict." % len(overlap_groups)

        if args.interactive:
            resolve_cmds = {}
        
        for groupid in sorted(overlap_groups.keys(), key=lambda x:int(x)):
            ogroup = overlap_groups[groupid]
            spn = utils.create_spanning(ogroup)
            ds = get_display_str([spn])
            ds_pad = get_display_str([spn], args.padding)
            
            # Print information to stderr            
            print >>sys.stderr, utils.groupstr(groupid, ogroup)
            
            # Go to locus
            igv.goto(ds_pad)
            
            # Snapshot
            if args.snapshot:
                if not os.path.exists(os.path.join(snapdir, 'overlap')):
                    os.mkdir(os.path.join(snapdir, 'overlap'))
                igv.snapshotDirectory(os.path.join(snapdir, 'overlap'))
                igv.snapshot('overlap%s.png' % groupid)
                        
            if args.interactive:
                concmd = utils.prompt_cmd()
                resolve_cmds[groupid] = concmd
        
        if args.interactive:
            if len(resolve_cmds) > 0:
                print >>sys.stderr, 'Commands for resolving conflicts (JSON):\n'
                print >>sys.stderr, json.dumps(resolve_cmds)
                print >>sys.stderr

            if args.resolve_file is not None:
                print >>args.resolve_file, json.dumps(resolve_cmds)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Use IGV to inspect HERV loci. ',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inspect', default='all', choices=['all','overlap',],
                         help='''Select the mode in which the program will run. If 
                                 --inspect is "all", every locus in the GTF will be
                                 shown. If --inspect is "overlap", bedtools is used to 
                                 cluster overlapping loci; clusters with 2 or more loci 
                                 are shown.''')
    parser.add_argument('--interactive', action='store_true',
                        help='''Enables interactive features. If --inspect is "all", the
                                program will pause after displaying each locus and wait
                                for the user to press enter. If --inspect is "overlap",
                                the user may enter commands about how to resolve the 
                                conflict; all commands are printed as JSON and can be
                                provided to a polishing script.''')
    parser.add_argument('--sleep_interval', type=int, default=1,
                        help='''Sleep interval (in milliseconds). Only effective when
                                --interactive is False.''')                                
    parser.add_argument('--resolve_file', type=argparse.FileType('w'),
                        help='''Resolve commands are printed to --resolve_file (in
                                addition to stderr).''')
    parser.add_argument('--compare_tracks',
                        help='''Comma-separated list of other tracks to be loaded into
                                IGV. Must be a format that is recognized by IGV, see 
                                http://software.broadinstitute.org/software/igv/FileFormats.
                                Useful for visually comparing annotations.''')
    parser.add_argument('--snapshot', action='store_true',
                        help="Turn on snapshots.")
    parser.add_argument('--snapshot_dir', default='./snapshots',
                        help="Directory to store snapshots.")
    parser.add_argument('--clear', action='store_true',
                        help="Clear current IGV session.")    
    parser.add_argument('--padding', type=int, default=1000,
                        help="Padding to show around locus.")    
    parser.add_argument('--refbuild', default='hg19',
                        help="Reference genome build")
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help="GTF to be inspected")      
    main(parser.parse_args())
