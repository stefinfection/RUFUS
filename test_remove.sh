#!/bin/bash


if [ $_arg_simple_ouput_mode == "TRUE" ]; then
	rm -r Intermediates/
       	rm -r TempOverlap/
	rm mer_counts_merged.jf
#todo: left off removing temp files		
else
	control_files=(
		bam.generator
		bam.generator.Jelly.chr
		bam.generator.Jhash
	        bam.generator.Jhash.histo
		bam.generator.Jhash.histo.7.7.dist
		bam.generator.Jhash.histo.7.7.model
	        bam.generator.Jhash.histo.7.7.out
		bam.generator.Jhash.histo.7.7.prob 	
	)	

	# remove control files
	for control in "${_arg_controls[@]}"
	do
		for postfix in "${control_files[@]}"
		do
			if [ -e ${control}.${postfix} ];
				rm 
			fi
		done
	done
	
	# remove subject files
	subject_files=(
		"bam.generator"
		"bam.generator"       
		"bam.generator.V2.overlap.fastqd"
		"bam.generator.Jelly.chr"
		"bam.generator.V2.overlap.hashcount.fastq"
		"bam.generator.Jhash"
		"bam.generator.V2.overlap.hashcount.fastq.bam"
		"bam.generator.Jhash.histo"
		"bam.generator.V2.overlap.hashcount.fastq.bam.bai"
		"bam.generator.Jhash.histo.7.7.dist"
		"bam.generator.V2.overlap.hashcount.fastq.bam.coinherited.vcf.gz"
		"bam.generator.Jhash.histo.7.7.model"
		"bam.generator.V2.overlap.hashcount.fastq.bam.coinherited.vcf.gz.tbi"
		"bam.generator.Jhash.histo.7.7.out"
		"bam.generator.V2.overlap.hashcount.fastq.bam.vcf"
		"bam.generator.Jhash.histo.7.7.prob"
		"bam.generator.V2.overlap.hashcount.fastq.bam.vcf.bed"
		"bam.generator.Mutations.Mate1.fastq"    
		"bam.generator.filter.chr"
		"bam.generator.Mutations.Mate2.fastq"    
		"bam.generator.k25_c4.HashList"
		"bam.generator.Mutations.fastq.bam"      
		"bam.generator.temp"
		"bam.generator.Mutations.fastq.bam.bai"  
		"bam.generator.temp.mate1.fastq"
		"bam.generator.V2.overlap.fastq"        
		"bam.generator.temp.mate2.fastq"
	)
		for postfix in "${subject_files[@]}"
		do
			rm ${subject}.${postfix}
		done
fi

echo "done with everything"
exit 0
# ] <-- needed because of Argbash
