# [L1Base.hg38.v1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgct_customText=track%20name=L1Base.hg38.v1%20description=%27L1Base.hg38.v1%27%0Ahttps://github.com/mlbendall/telescope_annotation_db/raw/master/builds/L1Base.hg38.v1/transcripts.gtf)


Annotations of putatively active LINE-1 elements from [L1Base](http://l1base.charite.de/).

## Quick Start:

### L1 annotation: [transcripts.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/L1Base.hg38.v1/transcripts.gtf)

<sup>We recommend using this annotation for running *Telescope*. </sup>

#### Load annotation into UCSC genome browser:

Click here: [L1Base.hg38.v1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgct_customText=track%20name=L1Base.hg38.v1%20description=%27L1Base.hg38.v1%27%0Ahttps://github.com/mlbendall/telescope_annotation_db/raw/master/builds/L1Base.hg38.v1/transcripts.gtf)

-----



## Building the annotation

### Setup

The software **[telebuilder](https://github.com/mlbendall/telebuilder)** constructs annotations for Telescope. The `buildL1` script requires information about which genome build will be used, chromosome sizes, cytogenetic band annotations, and other GTFs to compare with.

By default, `buildL1` will look for the genome build ID in a file called `build.txt`, so we create this in the current directory:

```bash
echo 'hg38' > build.txt
```

`buildL1` will also look for the chromosome sizes in a file called `chrom.sizes` and the cytogenetic band annotation in a file called `cytoband.gtf`. Create these links here:

```bash
ln -s ../../refs/hg38.chrom_sizes.txt chrom.sizes
ln -s ../../refs/hg38.cytoband.gtf cytoband.gtf
```

### Download the BED files from L1Base

```bash
mkdir -p l1base_tracks && cd l1base_tracks
wget http://l1base.charite.de/BED/hsflil1_8438.bed
wget http://l1base.charite.de/BED/hsorf2l1_8438.bed
wget http://l1base.charite.de/BED/hsflnil1_8438_rm.bed
cd ..
```

### Build annotations

`buildL1` converts the L1Base annotations to GTF, removes redundant annotations, and names loci according to chromosome band.

```bash
buildL1
```

### Combine GTFs

Combine all GTFs into one annotation file:

```bash
cat L1FLI.gtf L1ORF2.gtf L1FLnI.gtf | gtftools sort > transcripts.gtf
```

### Locus-level summary

Data in TSV format, convenient for loading in R

```
gtftools tsv --feat exon transcripts.gtf | gzip > genes.tsv.gz
```
