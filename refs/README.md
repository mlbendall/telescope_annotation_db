# Reference Annotations

How to build files related to specific reference genomes

```
# mkdir -p ucsc
# rsync -a -P rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ./ucsc
```

## hg19

Download the chromosome list with sizes from UCSC.

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
```

Filter the chromosome list so that only complete chromosomes are present.

```python
sizes = dict(l.strip('\n').split('\t') for l in open('hg19.chrom.sizes', 'rU'))
with open('hg19.chrom_sizes.txt','w') as outh:
    for c in ['chr%d' % _ for _ in range(1,23)] + ['chrX', 'chrY']:
        print >>outh, '%s\t%s' % (c, sizes[c])

```

Download locations of the cytogenetic bands from UCSC.

Sort the GTF according to our order above.

```bash
../tools/sortgtf.py --chroms hg19.chrom_sizes.txt < hg19_cyto.gtf > hg19.cytoband.gtf

```

Clean up:

```
rm -f hg19_cyto.gtf hg19.chrom.sizes
```

## hg38

Download the chromosome list with sizes from UCSC.

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

Sort the chromosome list so that complete chromosomes come first, then alternate haplotypes, then random sequences, then unplaced contigs.

```python
import re
sizes = dict(l.strip('\n').split('\t') for l in open('hg38.chrom.sizes', 'rU'))
sortord = map(str, range(1,23)) + ['X', 'Y', 'Un', 'M' ]
sortord = {'chr%s' %s : i for i,s in enumerate(sortord)}
rxs = [
    re.compile('chr\d?[XY\d]$'),
    re.compile('_alt$'),
    re.compile('_random$'),
    re.compile('')
]

with open('hg38.chrom_sizes.txt', 'w') as outh:
    for rx in rxs:
        ks = sorted(k for k in sizes.keys() if rx.search(k))
        ks.sort(key=lambda x:sortord[x.split('_')[0]])
        for c in ks:
            print >>outh, '%s\t%s' % (c, sizes.pop(c))

```

Download locations of the cytogenetic bands from UCSC.

![Table browser](../docs/table_browser_download.png)

Sort the table.

```
../tools/sortgtf.py --chroms hg38.chrom_sizes.txt < hg38_cyto.gtf > hg38.cytoband.gtf
```


Clean up:

```
rm -f hg38_cyto.gtf hg38.chrom.sizes
```


## Other annotations

Included are other annotations of transposable elements from the literature.


#### HERVK(HML-2)

+ **hg19.subramanianT1.gtf**
+ **hg19.subramanianT2.gtf** 


The locations of HERVK(HML-2) loci are described in tables 1 and 2 from:

> Subramanian RP, Wildschutte JH, Russo C, Coffin JM. 2011. Identification, characterization, and comparative genomic distribution of the HERV-K (HML-2) group of human endogenous retroviruses. Retrovirology 8: 90.

#### HERVW

+ **hg19.grandiS1.gtf**

The locations of HERVW loci are described in table S1 from:

>Grandi N, Cadeddu M, Blomberg J, Tramontano E. 2016. Contribution of type W human endogenous retroviruses to the human genome: characterization of HERV-W proviral insertions and processed pseudogenes. Retrovirology 13: 67.


#### liftOver

Convert hg19 annotations to hg38 coordinates using liftOver:

```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
chmod a+x liftOver
./liftOver -gff hg19.subramanianT1.gtf hg19ToHg38.over.chain.gz hg38.subramanianT1.gtf unMapped
./liftOver -gff hg19.subramanianT2.gtf hg19ToHg38.over.chain.gz hg38.subramanianT2.gtf unMapped
./liftOver -gff hg19.grandiS1.gtf hg19ToHg38.over.chain.gz hg38.grandiS1.gtf unMapped
```