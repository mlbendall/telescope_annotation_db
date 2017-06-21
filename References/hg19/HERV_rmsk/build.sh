#! /bin/bash

# Include the repository scripts in the PATH
which sortgtf.py
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    # Path to repository root
    RROOT=$(git rev-parse --show-toplevel)
    export PATH=$RROOT/scripts:$PATH
    export PYTHONPATH=$RROOT/python:$PYTHONPATH
fi

# This file is used to specify the chromosome order, used in sorting
CHROM=../hg19.chrom.sizes

# Build annotations for each family, auto mode
# resolve.[FAM].json must be present for families with overlapping annotations
grep -v '^#' families.tsv | sed "s/$(printf '\t')/ /g" | while read line; do
    . ./build_family.sh $line auto
done

# Generate the "exon" annotation (no "span" entries)
[[ -e HERV_rmsk.gtf ]] && mv HERV_rmsk.gtf HERV_rmsk.gtf.bak
grep -v '^#' families.tsv | cut -f1 | while read fam; do 
    perl -lane 'print if $F[2]!~/^span/' $fam/$fam.gtf
done | sortgtf.py --chroms $CHROM > HERV_rmsk.gtf

# Generate the spanning annotation
[[ -e HERV_rmsk.spanning.gtf ]] && mv HERV_rmsk.spanning.gtf HERV_rmsk.spanning.gtf.bak
grep -v '^#' families.tsv | cut -f1 | while read fam; do 
    perl -lane 'print unless $F[2]!~/^span/' $fam/$fam.gtf
done | sortgtf.py --chroms $CHROM > HERV_rmsk.spanning.gtf

# Generate the table
[[ -e HERV_rmsk.table.tsv.gz ]] && mv HERV_rmsk.table.tsv.gz HERV_rmsk.table.tsv.gz.bak
gtf2tsv.py HERV_rmsk.spanning.gtf | gzip > HERV_rmsk.table.tsv.gz
