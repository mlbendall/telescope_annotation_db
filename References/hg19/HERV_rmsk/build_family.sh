#! /bin/bash

### Model specific variables #############################################################
[[ -n $1 ]] && FAM="$1" || exit 1
[[ -n $2 ]] && intmodel="$2" || exit 1
[[ -n $3 ]] && models=$(sed 's/,/ /g' <<<"$3") || exit 1
[[ -n $4 ]] && FLANKSZ="$4" || exit 1
[[ -n $5 ]] && AUTO=true || AUTO=false

# Only allow AUTO if resolve file is found
# [[ ! -e resolve.$FAM.json ]] && AUTO=false


echo "Family:         $FAM"
echo "Internal Model: $intmodel"
echo "LTR models:     $models"
echo "Flank size:     $FLANKSZ"
echo "Automatic:      $AUTO"

mkdir -p $FAM

### Set environment variables ############################################################
which assembleHERV.py
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    # Path to repository root
    RROOT=$(git rev-parse --show-toplevel)
    export PATH=$RROOT/scripts:$PATH
    export PYTHONPATH=$RROOT/python:$PYTHONPATH
fi

# Chromosome index
CHROM=../hg19.chrom.sizes

### Step 1: Download RepeatMasker tracks from UCSC #######################################
# Internal model
[[ ! -e $FAM/$intmodel.txt ]] && mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
    -e "SELECT * FROM rmsk WHERE repName = '$intmodel';" > $FAM/$intmodel.txt

# LTR
for model in $models; do
    echo "Query: $model"
    [[ ! -e $FAM/LTR.$model.txt ]] && mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > $FAM/LTR.$model.txt
done

wc -l $FAM/*.txt
     
### Step 2: Convert RepeatMasker tables to GTF ###########################################
rmsk2gtf.py --chroms $CHROM --gene_region internal $FAM/$intmodel.txt | \
    fixRmskCoords.py | \
    sortgtf.py --chrom $CHROM > $FAM/$intmodel.gtf

for model in $models; do
    rmsk2gtf.py --chroms $CHROM --gene_region ltr $FAM/LTR.$model.txt | \
        fixRmskCoords.py | \
        sortgtf.py --chrom $CHROM > $FAM/LTR.$model.gtf
done

### Step 3: Merge internal hits ##########################################################
mergeHits.py --prefix $FAM --longdist 10000 $FAM/$intmodel.gtf | \
    sortgtf.py --chrom $CHROM > $FAM/internal.merged.gtf

### Step 4: Find LTRs flanking internal regions ##########################################
# Use bedtools to identify LTR annotations that flank and overlap internal annotations
grep 'span' $FAM/internal.merged.gtf | \
    bedtools slop -g $CHROM -b $FLANKSZ -i -| \
    bedtools intersect -wo -s -a - -b $FAM/LTR.*.gtf > $FAM/ltr.intersect.gtf

### Step 5: Assemble HERV proviruses #####################################################
assembleHERV.py $FAM/internal.merged.gtf $FAM/ltr.intersect.gtf | \
    sortgtf.py --chrom $CHROM > $FAM/assembled.gtf

### Step 6: Filter short loci ############################################################
# Filter loci that do not meet minimum percent of internal model coverage
filterHERVLoci.py --min_internal_pct 0.1 --reject_gtf $FAM/rejected.gtf $FAM/assembled.gtf | \
    sortgtf.py --chrom $CHROM > $FAM/filtered.gtf

### Step 7: Resolve conflicting loci #####################################################
cmptrks=$(sed 's/ //g' <<<$(ls -m alternate/*.gtf*))
# Create transcript_id attribute for IGV visualization
transferGTFAttr.py locus transcript_id $FAM/filtered.gtf $FAM/tmp.igv.gtf

overlaps=$(perl -lane 'print unless $F[2]!~/^span/' $FAM/filtered.gtf | bedtools cluster -i -| cut -f10 | uniq -d | wc -l)

if [[ $overlaps -eq 0 ]]; then
    echo "No overlaps found"
    cp $FAM/filtered.gtf $FAM/polished.gtf
else
    echo "$overlaps overlaps found"
    if [[ $AUTO == false ]]; then
        # Launch interactive
        echo "Launching interactive inspection"
        inspectHERVLoci.py --clear \
            --inspect overlap \
            --interactive \
            --resolve_file resolve.$FAM.json \
            --compare_tracks $cmptrks \
            $FAM/tmp.igv.gtf
    fi
    # The resolve file must exist to continue
    [[ ! -e resolve.$FAM.json ]] && echo "$overlaps overlaps found but resolve.$FAM.json does not exist." && exit 1
    
    # Polishing
    polishHERVLoci.py \
        --resolve_file resolve.$FAM.json \
        $FAM/filtered.gtf | \
        sortgtf.py --chrom $CHROM > $FAM/polished.gtf
fi

### Step 8: Rename loci by genomic location ##############################################
nameHERV.py --cytoband ../cytoband.gtf $FAM/polished.gtf | \
    transferGTFAttr.py locus transcript_id | \
    transferGTFAttr.py locus gene_id > $FAM/$FAM.gtf

### Step 9: Visually inspect resulting annotation ########################################
if [[ $AUTO == false ]]; then
    # Compare conflicts before and after
    [[ $overlaps -ne 0 ]] && inspectHERVLoci.py --inspect overlap --snapshot --snapshot_dir $FAM/snapshots --compare_tracks $FAM/$FAM.gtf $FAM/tmp.igv.gtf

    # Snapshot all final annotations
    inspectHERVLoci.py --inspect all --snapshot --snapshot_dir $FAM/snapshots $FAM/$FAM.gtf
fi
