# [retro.hg38.v1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgct_customText=track%20name=retro.hg38.v1%20description=%27retro.hg38.v1%27%0Ahttps://github.com/mlbendall/telescope_annotation_db/raw/master/builds/retro.hg38.v1/transcripts.gtf)

Combines L1Base and HERV_rmsk annotations

-----

## Building the annotation


### 1. Obtain HERV and L1 annotations

See details for [HERV_rmsk.hg38.v2](../HERV_rmsk.hg38.v2/README.md) and [L1Base.hg38.v1](../L1Base.hg38.v1/README.md).


### 2. Combine, filter and sort annotations

```bash
cat ../HERV_rmsk.hg38.v2/transcripts.gtf ../L1Base.hg38.v1/transcripts.gtf |\
 grep -v '^#' |\
 perl -lane 'print if $F[2]=~/exon/' |\
 gtftools sort --chrom_sizes ../../refs/hg38.chrom_sizes.txt > transcripts.gtf

```

### Locus-level summary

Data in TSV format, convenient for loading in R

```python
import gzip

chrom_order = [l.split('\t')[0] for l in open('../../refs/hg38.chrom_sizes.txt', 'rU')]
chrom_order = {v:i for i,v in enumerate(chrom_order)}

hlines = (l.strip('\n').split('\t') for l in gzip.open('../HERV_rmsk.hg38.v2/genes.tsv.gz', 'rb'))
hhead = next(hlines)
lod = [dict(zip(hhead, hl)) for hl in hlines]

llines = (l.strip('\n').split('\t') for l in gzip.open('../L1Base.hg38.v1/genes.tsv.gz', 'rb'))
lhead = next(llines)
lod += [dict(zip(lhead, ll)) for ll in llines]

lod.sort(key=lambda x: int(x['start']))
lod.sort(key=lambda x: chrom_order[x['chrom']])

header = hhead + [f for f in lhead if f not in hhead]

with gzip.open('genes.tsv.gz', 'wb') as outh:
    print >>outh, '\t'.join(header)
    for d in lod:
        print >>outh, '\t'.join(d.get(f, '.') for f in header)
```