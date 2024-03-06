#assign necessary variables
REFERENCE_FILE=$0
RUFUS_VCF=$1
OUTFILE=$2
CONTROL_VCF=$3
CONTROL_ALIGNED=$4
control_bams=(list of control files given by user)

#format final rufus vcf for intersections
vt normalize -n -r $REFERENCE_FILE $RUFUS_VCF -| vt decompose_blocksub - | bgzip > $OUTFILE
bcftools index $OUTFILE

#for loop for each control file provided by user
for CONTROL in "${control_bams[@]}"; do
    #check to see if the provided bam file is aligned
    if [ "$(samtools view -H "$CONTROL" | grep -c '^@SQ')" -gt 0 ]; then
        CONTROL_BAM=$CONTROL
        #run pileup and call variants
        bcftools mpileup -Oz -d600 -f $REFERENCE_FILE -T $OUTFILE $CONTROL_BAM - | bcftools call -Oz -v > $CONTROL_VCF
        #intersect the control vcf with formatted rufus vcf
        bcftools isec -Oz -w1 -n=1 -p $OUTPUT_DIR $OUTFILE $CONTROL_VCF
        #save new rufus only vcf from intersection as the new outfile
        OUTFILE=$OUTPUT_DIR/0000.vcf.gz
    else
        bwa mem -t 40 $REFERENCE_FILE $CONTROL | samtools view -S -@ 12 -b - > $CONTROL_ALIGNED
        CONTROL_BAM=$CONTROL_ALIGNED
        bcftools mpileup -Oz -d600 -f $REFERENCE_FILE -T $OUTFILE $CONTROL_BAM - | bcftools call -Oz -v > $CONTROL_VCF
        bcftools isec -Oz -w1 -n=1 -p $OUTPUT_DIR $OUTFILE $CONTROL_VCF
        #save the new vcf as rufus final vcf
        OUTFILE=$OUTPUT_DIR/0000.vcf.gz
    fi
done
