# HERV_rmsk.hg19.v1

The annotation is constructed using RepeatMasker hits for 23 HERV families. See `families.tsv` for table of families, internal models, and LTR models.

## Quick Start:

### HERV annotation: [transcripts.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg19.v1/transcripts.gtf)

<sup> This annotation includes only genomic regions that are matched by RepeatMasker. Thus, a single transcript may be discontinuous and represented by multiple rows. We recommend using this annotation for running *Telescope*. </sup>

### HERV loci: [genes.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg19.v1/genes.gtf)

<sup> Annotates entire genomic region spanning HERV locus. A single transcriptional unit is represented by one row, and may contain regions that are not matched by RepeatMasker.</sup>


#### Load annotation into UCSC genome browser:

1.  Navigate to [UCSC Genome Browser on Human hg19 Assembly](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19).
2. Go to **My Data > Custom tracks**.
3. Paste this into the first box:


    ```
    track name='Telescope.v1' description='Telescope Annotation'
    https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg19.v1/transcripts.gtf
    ```

4. Click "Submit" to load the track into your browser session. 


## How to build the annotation

The main scripts for building this annotation are `build_family.sh` and `build.sh`. Additional scripts used by `build_family.sh` are contained within the `scripts/` and `python/` directories.

### Setup

```bash
ln -s ../../refs/hg19.chrom_sizes.txt chrom.sizes
ln -s ../../refs/hg19.cytoband.gtf cytoband.gtf
```

```bash
mkdir -p compare
ln -s ../../../refs/hg19.subramanianT1.gtf compare/hg19.subramanianT1.gtf
ln -s ../../../refs/hg19.subramanianT2.gtf compare/hg19.subramanianT2.gtf
ln -s ../../../refs/hg19.grandiS1.gtf compare/hg19.grandiS1.gtf
```

### Build the annotation

The `build.sh` script parses the family table, calls `build_family.sh` for each family, and combines all family annotations into a single genome-wide HERV annotation. The `resolve.*.json` files are used to resolve conflicts.

```bash
/bin/bash build.sh
```

## Legacy annotation

Below is the README information that appeared with version 1 of this annotation.

> > In [RepBase](http://www.girinst.org/repbase/), each family is represented by one HERV internal model and one or more corresponding LTR models. This pipeline identifies HERV proviruses by grouping nearby regions that match to a single family. These regions are checked for agreement with a _relaxed_ expectation of what a provirus looks like (meaning that a provirus could be missing one or both of its LTRs). Since there is no clear rule about what qualifies as a HERV provirus, there is a great deal of manual/visual curation involved in creating an annotation; we have sought to automate this process where possible. The final output of this pipeline is a GTF annotation of HERV proviruses that can be used by other programs such as *Telescope*.
> > 
> > ## Fastest
> > 
> > ### The complete annotations are saved in this repository and are ready to use.
> > 
> > + `HERV_rmsk.gtf`: Hits identified by RepeatMasker are annotated as exon features, while regions not identified as belonging to the same family are not annotated. Records belonging to the same provirus are linked by the __transcript_id__, __gene_id__, and __locus__ attributes. This file is suitable for use in *Telescope*, since this ensures that only reads matching the model are  assigned to a locus.
> > + `HERV_rmsk.spanning.gtf`: The spanning regions of all proviruses are annotated as span features. The  __transcript_id__, __gene_id__, and __locus__ attributes can be linked back to the exon file. This is useful for extracting sequences and identifying other elements overlapping HERV loci. 
> > + `HERV_rmsk.table.tsv`: The genomic location and attributes for all proviruses in an easily parsable format.
> > 
> > ##### Difference between spanning and non-spanning annotation:
> > 
> > ![gtf_example1](/docs/gtf_example1.png?raw=true "GTF example 1")
> > 
> > ## Fast
> > 
> > ### Build annotations using provided parameters.
> > 
> > This is primarily to illustrate the process we used to obtain the above annotations. The `build_family.sh` script performs the steps (described below) needed to build annotations for a single family, given the internal model, LTR models, and flank size. The flank size was chosen in a way that maximizes the number of prototype loci without overextending. For families with overlapping loci, this also relies on the `resolve.[FAM].json` files (included in the repository) that encode operations performed to resolve conflicts.
> > 
> > Our entire pipeline can be run with the `build.sh` script. This will build each family using the parameters in `families.tsv`, combine the annotations for each family, and produce `HERV_rmsk.gtf`, `HERV_rmsk.spanning.gtf`, and `HERV_rmsk.table.tsv`.
> > 
> > ```
> > ./build.sh
> > ```
> > 
> > ## Details
> > 
> > The pipeline we used to generate this annotation is implemented in `build_family.sh`. This can be used as-is to reproduce our annotation, or can be modified to accommodate different assumptions. All the scripts used by `build_family.sh` have usage information to provide details about the command-line arguments available.
> > 
> > `build_family.sh` performs 9 steps to produce a complete annotation:
> > 
> > ### Step 1: Download RepeatMasker tracks from UCSC
> > 
> > Obtain RepeatMasker annotations by querying UCSC's mysql server. For example:
> > 
> > ```bash
> > mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A  \
> >     -e "SELECT * FROM rmsk WHERE repName = '[modelname]';"
> > ```
> > 
> > For each family, download the rows corresponding to the internal model and each of the LTR models and save these as text files.
> > 
> > ### Step 2: Convert RepeatMasker tables to GTF
> > 
> > All scripts used in this pipeline operate on GTF files. Convert the tables to GTF using `rmsk2gtf.py`, then use `fixRmskCoords.py` to correct errors in the RepeatMasker model coordinates.
> > 
> > ### Step 3: Merge internal hits
> > 
> > The internal region of a provirus is most often composed of several distinct hits identified by RepeatMasker. This is because various parts of the provirus do not match the RepeatMasker models due to accumulated mutations, indels, or other differences. The goal here is to merge groups of adjacent hits to identify a continuous region that may represent the internal region.
> > 
> > `mergeHits.py` reads the internal hits (from GTF file) and creates spanning records for hits that may belong to the same provirus. The records are linked by the __locus__ attribute. For example, these two annotations belong to the same model and are only two bases apart:
> > 
> > ```
> > chr1 rmsk exon 80394070 80395995 13998 + . repName "HERVL-int"; repStart "61"; repEnd "2018";
> > chr1 rmsk exon 80395997 80399424 26731 + . repName "HERVL-int"; repStart "2225"; repEnd "5654";
> > ```
> > 
> > `mergeHits.py` automatically merges hits that are within a short distance of each other (default is 10). For hits separated by longer distances, `mergeHits.py` checks whether the model coordinates are sequential; this happens when there are divergent regions or insertions. These sequential hits are merged if they are on the same strand and the distance between is not extreme (default is 10000, the approximate size of a provirus). The resulting output includes the original annotations and a spanning annotation that spans the internal region identified:
> > 
> > ```
> > chr1 rmsk span 80394070 80399424 . + . locus "HERVL_1206";
> > chr1 rmsk exon 80394070 80395995 13998 + . locus "HERVL_1206"; repName "HERVL-int"; repStart "61"; ...
> > chr1 rmsk exon 80395997 80399424 26731 + . locus "HERVL_1206"; repName "HERVL-int"; repStart "2225"; ...
> > ```
> > 
> > ### Step 4: Find LTRs flanking internal regions
> > 
> > LTR regions associated with each provirus are identified by intersecting the internal regions with the LTR annotations. Since LTRs are usually found on either side of the internal genes, we expand the merged regions (from step 3) before intersecting with all LTR files:
> > 
> > ```bash
> > grep 'span' internal.merged.gtf | \
> >     bedtools slop -g [faidx] -b [flank_size] -i -| \
> >     bedtools intersect -wo -s -a - -b LTR.*.gtf
> > ```
> > 
> > The resulting file contains two GTF records per line. The first GTF record is the internal region, and the second is the LTR annotations that intersect it.
> > 
> > ### Step 5: Assemble HERV proviruses
> > 
> > `assembleHERV.py` combines the internal regions identified in step 3 with the flanking LTRs identified in step 4 to produce assembled HERV proviruses. We also determine whether the provirus has LTRs on both sides ("prototype"), one side ("oneside"), or neither ("internal"). Also, the number of bases matching the model is estimated from the __repStart__, __repEnd__, and __repLeft__ attributes.
> > 
> > ### Step 6: Filter short loci
> > 
> > Loci that are very short do not provide useful information to the analysis. `filterHERVLoci.py` provides basic filtering of loci based on the number of bases matching the internal model or the percentage of the internal model that is matched.
> > 
> > ### Step 7: Resolve conflicting loci
> > 
> > After filtering, we are left with a rough estimate of the genomic locations for a HERV family. However, visual inspection of these loci is necessary to confirm whether the loci identified truly resemble HERV proviruses. We have automated this process using a script, `inspectHERVLoci.py`, to control [IGV](http://software.broadinstitute.org/software/igv/), a genomic viewer developed by Broad Institute. There are several modes of operation available with this script, including interactive usage or snapshot. This is fully described by the program usage: `inspectHERVLoci.py --help`.
> > 
> > If overlapping loci are present in the annotation, these can be resolved with `polishHERVLoci.py`. This script reads in the filtered annotation and performs operations on the conflicting loci to produce a polished annotation. The operation commands are provided by the user and can be specified in several ways. The first method is to run `polishHERVLoci.py` interactively and type in the desired operations when prompted. This is a good solution when IGV is not available on your system. Alternatively, the user can provide a JSON object (string or file) that contains the desired commands. We do not recommend that the user manually create the JSON object; instead, the object is created by running `inspectHERVLoci.py` with the `--inspect overlap` and `--interactive` options set. This allows the user to decide how to resolve conflicts while visually inspecting overlapping loci.
> > 
> >  Possible operations include:
> > 
> > + __ignore__                    : Do nothing, do not change annotations
> > + __reject loc1[,loc2,...]__    : Remove locs from annotations
> > + __diff loc1-loc2__            : Shorten loc1 to eliminate overlap with loc2
> > + __merge loc1+loc2__           : Combine loc1 and loc2, (and loc3, etc.)
> > 
> > The JSON object is a hash mapping the overlap group number to the command(s) for that group:
> > 
> > ```json
> > {"15": ["diff", "HERVW_0262-HERVW_0266"]}
> > ```
> > 
> > __NOTE:__ The `resolve.[FAM].json` files included in this repository will not work if there are any changes to the filtered annotation, since any change to the annotation may result in different overlap groupings. We intend these files to be used only to recreate our annotation.
> > 
> > ### Step 8: Rename loci by genomic location
> > 
> > To this point, the locus names have been used mainly to identify groups of annotations belonging to the same locus, and are not very useful. `nameHERV.py` renames the loci according to their family prefix and cytogenetic band. Bands containing multiple loci have a lowercase letter appended to the name in sequential order.
> > 
> > ### Step 9 (optional): Visually inspect resulting annotation
> 
> `inspectHERVLoci.py` can be used to visualize or snapshot the final annotations.
