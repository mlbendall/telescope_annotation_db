# L1Base.hg19.v1

Annotations of putatively active LINE-1 elements from [L1Base](http://l1base.charite.de/).

## Quick Start:

### L1 annotation: [transcripts.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/L1Base.hg19.v1/transcripts.gtf)

<sup>We recommend using this annotation for running *Telescope*. </sup>

#### Load annotation into UCSC genome browser:

1.  Navigate to [UCSC Genome Browser on Human hg38 Assembly](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19).
2. Go to **My Data > Custom tracks**.
3. Paste this into the first box:


    ```
    track name='L1Base.v1' description='LINE1 Annotation'
    https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/L1Base.hg19.v1/transcripts.gtf
    ```

4. Click "Submit" to load the track into your browser session. 

-----

## Building the annotation

### Setup

The software **[telebuilder](https://github.com/mlbendall/telebuilder)** constructs annotations for Telescope. The `buildL1` script requires information about which genome build will be used, chromosome sizes, cytogenetic band annotations, and other GTFs to compare with.

By default, `buildL1` will look for the genome build ID in a file called `build.txt`, so we create this in the current directory:

```bash
echo 'hg19' > build.txt
```

`buildL1` will also look for the chromosome sizes in a file called `chrom.sizes` and the cytogenetic band annotation in a file called `cytoband.gtf`. Create these links here:

```bash
ln -s ../../refs/hg19.chrom_sizes.txt chrom.sizes
ln -s ../../refs/hg19.cytoband.gtf cytoband.gtf
```

### Download the BED files from L1Base.

**IMPORTANT**: L1Base annotations use hg38; hg19 annotations are not available. We will use liftOver to convert these to the correct coordinates.

```bash
mkdir -p l1base_tracks && cd l1base_tracks
wget http://l1base.charite.de/BED/hsflil1_8438.bed
wget http://l1base.charite.de/BED/hsorf2l1_8438.bed
wget http://l1base.charite.de/BED/hsflnil1_8438_rm.bed
```

### Download liftOver

Download liftOver depending on operating system:

```bash
# OSX:
wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/liftOver
# Linux:
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
```

### Download hg38 to hg19 chain file

```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

### Use liftOver to convert coordinates

```bash
for f in *8438*.bed; do
    echo "Processing $f"
    ./liftOver <(grep -v '^track' $f) hg38ToHg19.over.chain.gz ${f%%_*}_hg19.bed unmapped.bed
    wc -l unmapped.bed
done

# clean up
rm *8438*.bed
rm unmapped.bed
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
