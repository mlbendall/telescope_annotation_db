# HERV_rmsk.hg19.v2

The annotation is constructed using RepeatMasker hits for 60 HERV families. See [`families.tsv`](./families.tsv) for table of families, internal models, and LTR models.

## Quick Start:

### HERV annotation: [transcripts.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg19.v2/transcripts.gtf)

<sup> This annotation includes only genomic regions that are matched by RepeatMasker. Thus, a single transcript may be discontinuous and represented by multiple rows. We recommend using this annotation for running *Telescope*. </sup>

### HERV loci: [genes.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg19.v2/genes.gtf)

<sup> Annotates entire genomic region spanning HERV locus. A single transcriptional unit is represented by one row, and may contain regions that are not matched by RepeatMasker.</sup>


#### Load annotation into UCSC genome browser:

1.  Navigate to [UCSC Genome Browser on Human hg19 Assembly](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19).
2. Go to **My Data > Custom tracks**.
3. Paste this into the first box:


    ```
    track name='Telescope.v2' description='Telescope Annotation'
    https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg19.v2/transcripts.gtf
    ```

4. Click "Submit" to load the track into your browser session. 



------


## Building the annotation

In [RepBase](http://www.girinst.org/repbase/), each family is represented by one HERV internal model and one or more corresponding LTR models. This pipeline identifies HERV proviruses by grouping nearby regions that match to a single family. These regions are checked for agreement with a _relaxed_ expectation of what a provirus looks like (meaning that a provirus could be missing one or both of its LTRs). Since there is no clear rule about what qualifies as a HERV provirus, there is manual/visual curation involved in creating an annotation; we have sought to automate this process where possible. The final output of this pipeline is a GTF annotation of HERV proviruses that can be used by other programs such as *Telescope*.


### Setup

The software **[telebuilder](https://github.com/mlbendall/telebuilder)** constructs annotations for Telescope. The `buildERV` script requires information about which genome build will be used (hg19, hg38, ...), chromosome sizes, cytogenetic band annotations, and other GTFs to compare with.

By default, `buildERV` will look for the genome build ID in a file called `build.txt`, so we create this in the current directory:

```bash
echo 'hg19' > build.txt
```

`buildERV` will also look for the chromosome sizes in a file called `chrom.sizes` and the cytogenetic band annotation in a file called `cytoband.gtf`. Create these links here:

```bash
ln -s ../../refs/hg19.chrom_sizes.txt chrom.sizes
ln -s ../../refs/hg19.cytoband.gtf cytoband.gtf
```

Finally, `buildERV` will search the `compare` directory for other annotations used for comparison purposes. Create links to these as well:

```bash
mkdir -p compare
ln -s ../../../refs/hg19.subramanianT1.gtf compare/hg19.subramanianT1.gtf
ln -s ../../../refs/hg19.subramanianT2.gtf compare/hg19.subramanianT2.gtf
ln -s ../../../refs/hg19.grandiS1.gtf compare/hg19.grandiS1.gtf
```

### Build HERV family annotations

Here we use the `buildERV` script to construct GTF annotations for all HERV families.

**Usage:** `buildERV --auto --save_intermediate [fam] [intmodel] [ltrmodel]`

where: 

+ `fam` family name
+ `intmodel` internal model name
+ `ltrmodel` a comma separated list of LTR model names
+ `--auto` perform automatic conflict resolution
+ `--save_intermediate` save intermediate GTF files for additional analysis (see below)

I have compiled the internal and LTR model names for 60 HERV families in [families.tsv](families.tsv). All families can be built by running the `buildERV` script on each row of the table:

```bash
IFS=$'\t' grep -v '^#' families.tsv | while read n im lm ; do
    mkdir -p $n
    buildERV --auto --save_intermediate $n $im $lm 2>&1 | tee $n/build.log    
done
```

### Analyze LTR assignments

The `analyzeERV` script examines internal HERV locations and identifies all the possible LTR families flanking these elements. This can help to determine whether you have selected the complete set of LTRs for a given internal model. NOTE: The script will download the LTR models from the RepeatMasker track for each chromosome, which will be saved in `all_ltr/`.

**Usage:** `analyzeERV [internal_gtf]`

where `internal_gtf` is a GTF file with internal annotations. It is recommended to cluster these annotations so that internal loci are fully extended. We use the `merged.gtf` files for each family for this analysis. 

```bash
IFS=$'\t' grep -v '^#' families.tsv | while read n im lm ; do
    analyzeERV $n/merged.gtf 2>&1 | tee $n/analyze.log
done
```

### Combine GTFs

Combine annotations for all families:

```bash
IFS=$'\t' grep -v '^#' families.tsv | while read n im lm ; do
    cat $n/$n.gtf
done | gtftools sortclust > transcripts.gtf
```

Create an annotation with only the spanning gene coordinates:

```bash
IFS=$'\t' grep -v '^#' families.tsv | while read n im lm ; do
    perl -lane 'print if $F[2]=~/gene/' $n/$n.gtf
done | gtftools sort > genes.gtf
```


### Additional steps (optional)

#### Compare this annotation to the previous

```python
from telebuilder.utils.gtfutils import *
v1 = list(read_gtf_file('HERV_rmsk.hg19.v1/genes.gtf'))
v2 = list(read_gtf_file('HERV_rmsk.hg19.v2/genes.gtf'))

isect = intersect_gtf(v1, [v2,], False)
for gA, gBs in isect:
    locA = gA.attr['locus']
    if not gBs:
        print "Not found: {}".format(gA)
    else:
        locBs = set(gB.attr['locus'] for i,gB in gBs)
        if locA not in locBs:
            print "locus in A: {}\tlocus in B: {}".format(locA, ','.join(sorted(locBs)))
```
