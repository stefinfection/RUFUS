# This script removes inherited variants called by rufus which lie on the same contig as
# a somatic variant. It accomplishes this task by performing a pileup and variant call
# in the control bam for each of the rufus-identified variants, then taking the complement
# of the control set in an intersection.

# parse command line arguments
REFERENCE_FILE=$1  # Reference file (.fa) - must be the same as RUFUS run that generated RUFUS_VCF
RUFUS_VCF=$2       # The 'FINAL' vcf generated by rufus, to be isec'd
OUT_VCF=$3
SRC_DIR=$4
ARG_LIST=("$@")
CONTROL_BAM_LIST=("${ARG_LIST[@]:4}") # Remaining args, all control bams

cd $SRC_DIR

# static vars
CONTROL_ALIGNED="temp_aligned.bam"
CONTROL_VCF="isec_control.vcf.gz"
NORMED_VCF="normed.${RUFUS_VCF}"
BWA=/opt/RUFUS/bin/externals/bwa/src/bwa_project/bwa
PILEUP_SCRIPT="/opt/RUFUS/post_process/single_pileup.sh"

# make intersection directory
ISEC_OUT_DIR="temp_isecs"
mkdir -p "${ISEC_OUT_DIR}"

#format final rufus vcf for intersections
vt normalize -n $RUFUS_VCF -r $REFERENCE_FILE | vt decompose_blocksub - | bgzip > $NORMED_VCF
bcftools index -t $NORMED_VCF

#for loop for each control file provided by user
for CONTROL in "${CONTROL_BAM_LIST[@]}"; do
    MADE_ALIGN_CONTROL=false
    
    #check to see if the provided bam file is aligned
    if [ "$(samtools view -H "$CONTROL" | grep -c '^@SQ')" -gt 0 ]; then
        CONTROL_BAM=$CONTROL
    else
        $BWA mem -t 40 $REFERENCE_FILE $CONTROL | samtools view -S -@ 12 -b - > $CONTROL_ALIGNED
        CONTROL_BAM=$CONTROL_ALIGNED
	    MADE_ALIGN_CONTROL=true
    fi
    
    #run pileup and call variants
	
	#bcftools query -f '%CHROM\t%POS0\t%POS\n' $NORMED_VCF > $NORMED_BED
	#PILEUP_VCF="pileup.vcf"
	bcftools query -f '%CHROM\n' $NORMED_VCF | sort | uniq > regions.out
	cat regions.out | parallel -j +0 $PILEUP_SCRIPT {} $CONTROL_BAM $REFERENCE
	
	# combine pileups	
	MERGED_PILEUP="merged_pileup.vcf.gz"
	bcftools concat -o $MERGED_PILEUP -Oz mpileup_*.vcf

    #bcftools mpileup -d 100 -r -f $REFERENCE_FILE -o $PILEUP_VCF $CONTROL_BAM
	echo "Finished pileups, starting call" >&2
	
	# call variants from merged pileup vcf
	bcftools call -cv -Oz -o $CONTROL_VCF $MERGED_PILEUP
	echo "Finished call, starting to index $CONTROL_VCF" >&2

	if [ -z "$CONTROL_VCF" ]; then
		"Could not find $CONTROL_VCF; exiting..." >&2
		exit
	fi
    bcftools index -t $CONTROL_VCF
    	
    #intersect the control vcf with formatted rufus vcf
	echo "Starting intersection..." >&2
    bcftools isec -Oz -w1 -n=1 -p $ISEC_OUT_DIR $NORMED_VCF $CONTROL_VCF    
 
    # save the new vcf as rufus final vcf
    OUTFILE="$ISEC_OUT_DIR/0000.vcf.gz"
	bcftools index -t $OUTFILE
    OUT_INDEX="$ISEC_OUT_DIR/0000.vcf.gz.tbi"

    cp "$OUTFILE" "${OUT_VCF}" 
    cp "$OUT_INDEX" "${OUT_VCF}.tbi" 

    # clean up aligned control file, if it exists
	if [ "$MADE_CONTROL_BAM" = true ]; then
    	rm $CONTROL_ALIGNED
	fi

	#TODO: delete this after testing instead of moving
	mv $ISEC_OUT_DIR rufus_supplementals/
	rm $NORMED_BED
done
