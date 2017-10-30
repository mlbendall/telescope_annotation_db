#! /bin/bash

# Include the scripts in the PATH
which sortgtf.py
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $0)/scripts:$PATH
    export PYTHONPATH=$(dirname $0)/python:$PYTHONPATH
fi

# This file is used to specify the chromosome order, used in sorting
CHROM=chrom.sizes

# Build annotations for each family, auto mode
# resolve.[FAM].json must be present for families with overlapping annotations
sed "s/$(printf '\t')/ /g" families.tsv | while read line; do
    ./build_family.sh $line auto
done

# Generate the "exon" annotation (no "span" entries)
cut -f1 families.tsv | while read fam; do 
    perl -lane 'print if $F[2]!~/^span/' $fam/$fam.gtf
done | sortgtf.py --chroms $CHROM > transcripts.gtf

# Generate the spanning annotation
cut -f1 families.tsv | while read fam; do 
    perl -lane 'print unless $F[2]!~/^span/' $fam/$fam.gtf
done | sortgtf.py --chroms $CHROM | sed 's/span/gene/' > genes.gtf
