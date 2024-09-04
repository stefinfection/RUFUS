echo "RUFUS post-process version C.0.1"
date
echo "RUFUS_$0 $@" >> .rufus.cmd

SUBJECT_STRING=$(basename $SUBJECT_FILE)
POST_PROCESS_DIR=/opt/RUFUS/post_process/

TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_STRING}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${SUBJECT_STRING}.combined.vcf.gz"

if [ "$WINDOW_SIZE" != "0" ]; then
	echo "Windowed run performed, trimming and combining region vcfs..."
	bash ${POST_PROCESS_DIR}trim_and_combine.sh $SUBJECT_STRING $CONTROL_STRING $WINDOW_SIZE
	mv $PREFILTERED_VCF rufus_supplementals/ # todo: am I doing this in trim and combine also?
fi

# Check for empty lines
echo "Checking vcf formatting..."
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_FINAL_VCF
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_PREFILTERED_VCF

# Sort
echo "Sorting..."
bcftools sort $TEMP_FINAL_VCF > "sorted.${TEMP_FINAL_VCF}"
bcftools sort $TEMP_PREFILTERED_VCF > "sorted.${TEMP_PREFILTERED_VCF}"

FINAL_VCF="RUFUS.Final.${subject_string}.combined.vcf.gz"
PREFILTERED_VCF="RUFUS.Prefiltered.${subject_string}.combined.vcf.gz"

# Remove coinheriteds
echo "Removing coinheriteds..."
SOURCE_DIR="."
bash ${POST_PROCESS_DIR}remove_coinheriteds.sh "$REFERENCE" "sorted.${TEMP_FINAL_VCF}" "$FINAL_VCF" "$SOURCE_DIR" "$CONTROL_STRING"

# Inject RUFUS command into header
echo "Composing final vcfs..."
bcftools view -h $COINHERITED_REMOVED_VCF | head -n -1 > $FINAL_VCF
cat .rufus.cmd >> $FINAL_VCF
bcftools view -h $COINHERITED_REMOVED_VCF | tail -n 1 >> $FINAL_VCF
bcftools view -H $COINHERITED_REMOVED_VCF >> $FINAL_VCF

bcftools view -h $TEMP_PREFILTERED_VCF | head -n -1 > $PREFILTERED_VCF
cat .rufus.cmd >> $PREFILTERED_VCF
bcftools view -h $TEMP_PREFILTERED_VCF | tail -n 1 >> $PREFILTERED_VCF
bcftools view -H $TEMP_PREFILTERED_VCF >> $PREFILTERED_VCF

# TODO: Separate SVs and SNV/Indels
#echo "Separating snvs/indels and SVs..."

# Cleanup
echo "Cleaning up intermediate post-processing files..."
#TODO: cleanup combined vcf
mv $PREFILTERED_VCF rufus_supplementals/
rm $TEMP_PREFILTERED_VCF
rm $TEMP_FINAL_VCF
rm "sorted.$TEMP_PREFILTERED_VCF"
rm "sorted.$TEMP_FINAL_VCF"

echo "Post-processing complete."

