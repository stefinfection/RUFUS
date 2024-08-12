#!/bin/bash

_arg_subject=$1
_arg_ref=$2
_arg_threads=$3
_arg_kmersize=$4
_arg_min=$5
_arg_refhash=$6
_arg_saliva=$7
_arg_exome=$8
_MaxAlleleSize=$9
_arg_mosaic=${10}
_assemblySpeed=${11}
_parallel_jelly=${12}
_pairedEnd=${13}
_arg_region=${14}
_arg_filterK=${15}
_arg_ParLowK=${16}
_filterMinQ=${17}

#########__CREATE_ALL_GENERATOR_FILES_AND_VARIABLES__#############
ProbandFileName=$(basename "$_arg_subject")
ProbandExtension="${ProbandFileName##*.}"
#echo "proband extension is $ProbandExtension"

######## checking proband extension, FASTQ is not handled, need to add that, for the meantime generator dumping to SAM needs to be used #############
if [[ "$ProbandExtension" != "cram" ]] && [[ "$ProbandExtension" != "bam" ]] || [[ ! -e "$_arg_subject" ]] && [[ "$ProbandExtension" != "generator" ]]
then
    echo "The proband bam/generator file" "$_arg_subject" " was not provided or does not exist; killing run with non-zero exit status"
    kill -9 $$
elif [[ "$ProbandExtension" == "bam" ]]
then
#   echo "you provided the proband cram file" "$_arg_subject"
    ProbandGenerator="$ProbandFileName".generator
    echo "samtools view -F 3328 $_arg_subject $_arg_region" > "$ProbandGenerator"
elif [[ "$ProbandExtension" == "cram" ]]
then
#   echo "you provided the proband cram file" "$_arg_subject"
    ProbandGenerator="$ProbandFileName".generator
    if [ "$_arg_cramref" == "" ]
    then
         echo "ERROR cram reference not provided for cram input";
        kill -9 $$
     fi
    echo "samtools view -F 3328 -T $_arg_cramref $_arg_subject  $_arg_region" > "$ProbandGenerator"
elif [[ "$ProbandExtension" = "generator" ]]
then
#   echo "you provided the proband bam file" "$_arg_subject"
    ProbandGenerator="$ProbandFileName"
else
    echo "unknown error during generator generation, killing run with non-zero exit status"
fi

ParentGenerators=()
ParentJhash=()
ParentFileNames=""
space=" "

for parent in "${Parents[@]}"
do
    parentFileName=$(basename "$parent")
    ParentFileNames=$ParentFileNames$space$parent
#    echo "parent file name is" "$parentFileName"
    parentExtension="${parentFileName##*.}"
#    echo "parent file extension name is" "$parentExtension"

    if  [[ "$parentExtension" != "cram" ]] && [[ "$parentExtension" != "bam" ]]  && [[ "$parentExtension" != "generator" ]]
    then
	echo "The control bam/generator file" "$parent" " was not provided, or does not exist; killing run with non-zero exit status"
	kill -9 $$
    elif [[ "$parentExtension" == "bam" ]]
    then
	    parentGenerator="$parentFileName".generator
	    ParentGenerators+=("$parentGenerator")
	    echo "samtools view -F 3328 $parent  $_arg_region" > "$parentGenerator"
#	    echo "You provided the control bam file" "$parent"
    elif [[ "$parentExtension" == "cram" ]]
    then
            parentGenerator="$parentFileName".generator
            ParentGenerators+=("$parentGenerator")
	    if [ "$_arg_cramref" == "" ]
	    then
		echo "ERROR cram reference not provided for cram input";
		 kill -9 $$
	    fi
            echo "samtools view -F 3328 -T $_arg_cramref $parent  $_arg_region" > "$parentGenerator"
 #           echo "You provided the control cram file" "$parent"
    elif [[ "$parentExtension" = "generator" ]]
    then
	parentGenerator="$parentFileName"
        ParentGenerators+=("$parentGenerator")
#	echo "You provided the control bam file" "$parent"
    fi
done
#################################################################


################__COPY_ARG_BASH_VARIABLES_TO_SCRIPT_VARIABLES__##################

K=$_arg_kmersize
Threads=$_arg_threads
ref=$_arg_ref
#################################################################################

if [[ -z "$K" ]]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "kmer size ([-k|kmersize]) was not provided, killing run with non-zero exit status"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$

fi

if [[ -z "$Threads" ]]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "number of threads ([-t|--threads]) was not provided, killing run with non-zero exit status"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$
fi

if [[ -z "$ref" ]]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "reference file ([-r|--ref]) was not provided, killing run with non-zero exit status"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$

fi


###################__PRINT_VARIABLES_USED__######################################
#echo "~~~~~~~~~~~~ printing out parameter values used in script ~~~~~~~~~~~~~~~~"
#echo "value of ProbandGenerator $ProbandGenerator"
#echo "Value of ParentGenerators:"
#for parent  in "${ParentGenerators[@]}"
#do
#  echo " $parent"
#done
#echo "Value of K is: $K"
#echo "Value of Threads is: $Threads"
#echo "value of ref is: $ref"
#echo "value of min is: $_arg_min"
#echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#################################################################################


if [ -z "$_arg_refhash" ]
then
    echo "Did not provide refHash"
else
    echo "provided refHash of: " "$_arg_refhash"
fi

if ! [ -z "$_arg_min" ]
then
      echo "\$_arg_min is NOT empty"
      MutantMinCov=$_arg_min
fi
######################################


############__BUILD_JHASH_STRING__################
parentsString=""
parentsExcludeString=""
space=" "
jhash=".Jhash"

for parent in "${ParentGenerators[@]}"
do
  #echo "parent is  $parent "
  parentsString=$parentsString$space$parent$jhash
done
for exclude in "${_arg_exclude[@]}"
do
    parentsExcludeString=$parentsExcludeString$space$exclude
done

##################################################

## EXECUTABLE PATHS ##
RUFUSmodel=$RDIR/bin/ModelDist
RUFUSfilter=$RDIR/bin/RUFUS.Filter
RufAlu=$RDIR/bin/externals/rufalu/src/rufalu_project/src/aluDetect
RUFUSOverlap=$RDIR/scripts/Overlap.shorter.sh
RunJelly=$RDIR/scripts/RunJellyForRUFUS.sh
PullSampleHashes=$RDIR/scripts/CheckJellyHashList.sh
modifiedJelly=$RDIR/bin/externals/modified_jellyfish/src/modified_jellyfish_project/bin/jellyfish
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
RUFUSfilterFASTQ=$RDIR/bin/RUFUS.Filter
RUFUSfilterFASTQse=$RDIR/bin/RUFUS.Filter.single
fastp=$RDIR/bin/externals/fastp/src/fastp_project/fastp
samblaster=$RDIR/bin/externals/samblaster/src/samblaster_project/samblaster


####################__GENERATE_JHASH_FILES_FROM_JELLYFISH__#####################
if [ $_parallel_jelly == "yes" ]
then
	######## TODO instead of assuming 3 samples
	JThreads=$(( Threads / 3 ))
	if [ "$JThreads" -lt 3 ]
	then
	    JThreads=3
	fi
	#JThreads=$Threads

	for parent in "${ParentGenerators[@]}"
	do
	      bash $RunJelly $parent $K $(echo $JThreads -2 | bc) $_arg_ParLowK  &
	done

	bash $RunJelly $ProbandGenerator $K $(echo $JThreads -2 | bc) 2  &
	wait
else
        JThreads=$Threads
	if [ "$JThreads" -lt 3 ]
        then
            JThreads=3
        fi

        for parent in "${ParentGenerators[@]}"
        do
              bash $RunJelly $parent $K $(echo $JThreads -2 | bc) $_arg_ParLowK
        done

        # bash $RunJelly $ProbandGenerator $K  $Threads 2
         bash $RunJelly $ProbandGenerator $K $(echo $JThreads -2 | bc) 2
fi
##############################################################################


###########################_EMPTY_JHASH_CHECK##############################
########TODO just checking file size isn't a great idea, when jellyfish fails the fields arent zero size
for parent in "${ParentGenerators[@]}"
do
    ## Check Jhash files are not empty
     if [ ! -s "$parent".Jhash ]
     then
        echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
        echo "$parent.Jhash  is empty"
        echo "Killing run with exit status 1"
        echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        kill -9 $$
     fi
done

if [ ! -s "$ProbandGenerator".Jhash ]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "$ProbandGenerator.Jhash  is empty"
    echo "Killing run with exit status 1"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$
fi
##############################################################################



##################__GENERATE_JHASH_HISTOGRAMS__#################################
######TODO I can probably get rid of this if I just make model read either tab or space
perl -ni -e 's/ /\t/;print' "$ProbandGenerator".Jhash.histo
for parent in "${ParentGenerators[@]}"
do
  perl -ni -e 's/ /\t/;print' "$parent".Jhash.histo
done
##############################################################################



#######################__RUFUS_Model__############################################
#if [ $_arg_exome == "FALSE" ] #[	-z "$_arg_min" ]
if [ -z "$_arg_min" ]  && [ $_arg_exome == "FALSE" ]
then
	echo "exome not set, assuming data is whole genome, building model" #echo "min not provided, building model"
	if [ -e "$ProbandGenerator.Jhash.histo.7.7.model" ]
	then
	 	echo "skipping model"
	else
		echo "starting model"
		"$RUFUSmodel" "$ProbandGenerator".Jhash.histo $K 150 $Threads > "$ProbandGenerator".Jhash.histo.7.7.out
		for parent in "${ParentGenerators[@]}"
		do
			"$RUFUSmodel" "$parent".Jhash.histo $K 150 $Threads > "$parent".Jhash.histo.7.7.out &
		done
		echo "done with model"
	fi

	if [ -z "$_arg_min" ]
	then
		if [ -e "$ProbandGenerator".Jhash.histo.7.7.model ]
		then
			echo "$(grep Best\ Model "$ProbandGenerator".Jhash.histo.7.7.out)"
			MutantMinCov=$(head -2 "$ProbandGenerator".Jhash.histo.7.7.model | tail -1 )
			echo "INFO: mutant min coverage from generated model is $MutantMinCov"

			MutantSC=$(head -4 "$ProbandGenerator".Jhash.histo.7.7.model | tail -1 )
			echo "INFO: mutant SC coverage from generated model is $MutantSC"
			MaxHashDepth=$(echo "$MutantSC * 5" | bc)
			echo "INFO: MaxHashDepth = $MaxHashDepth"
		else
			echo "ERROR Model didnt run correctly, exiting"
			return -1
		fi
	else
		echo "min coverage provided of $_arg_min, setting min kmer to that"
		MutantMinCov="$_arg_min"
	fi

else
	if [	-z "$_arg_min" ]
	then
		echo "min coverage must be provided with an exome run"
		return -1;
	else
####TODO: check what im dond here
		echo "3" > "$ProbandGenerator".Jhash.histo.7.7.model;
		echo "$_arg_min" >> "$ProbandGenerator".Jhash.histo.7.7.model;
		echo "3.1392e+09" >> "$ProbandGenerator".Jhash.histo.7.7.model;
		echo "1000000" >> "$ProbandGenerator".Jhash.histo.7.7.model;
		echo "min was provided, min is $_arg_min"
		MutantMinCov="$_arg_min"
		#touch "$ProbandGenerator".Jhash.histo.7.7.model
	fi
fi
########################################################################################

if [ "$_arg_stop" = "jelly" ];
then
        echo "-StJ used, stopping run";
        exit 1;
fi
#######################################################################################

if [ -z $MutantMinCov ]; then
	echo "ERROR: No min coverage set, possible error in Model"
	exit 100
fi
if [ "$MutantMinCov" -lt "2" ]
then
	echo "ERROR, model couldn't pick a sensible lower cutoff, check your subject bam file"
        exit
fi
#################################__HASH_LIST_FILTER__#####################################

echo "########### Running Mutant Hash Identification ##############"

if [ -s "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList ]
then
    echo "skipping $ProbandGenerator.HashList pull "
else
    if [ -e "$ProbandGenerator".temp ]
    then
    	rm  "$ProbandGenerator".temp
    fi
    mkfifo "$ProbandGenerator".temp
    $modifiedJelly merge "$ProbandGenerator".Jhash $(echo $parentsString) $(echo $parentsExcludeString)  > "$ProbandGenerator".temp &
    bash $PullSampleHashes $ProbandGenerator.Jhash "$ProbandGenerator".temp $MutantMinCov $MaxHashDepth > "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList
    wait

fi

########################################################################################

if [ $(head  "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList | wc -l | awk '{print $1}') -eq "0" ]; then
	echo "ERROR: No mutant hashes identified, either the files are exactly the same of something went wrong in previous step"
	exit 100
fi
########################################################################################
if [ "$_arg_stop" = "hash" ];
then
        echo "-StH used, stopping run";
        exit 1;
fi
######################__RUFUS_FILTER__##################################################
echo "########### starting RUFUS filter ###########"

if [ $_pairedEnd == "true" ]
then
	if [ -e "$ProbandGenerator".Mutations.Mate1.fastq ]
	then
		echo "skipping filter"
	else
		if [ -z $_arg_fastqA ]
		then
		    if [ -e "$ProbandGenerator".temp.mate1.fastq ]; then
		    	rm  "$ProbandGenerator".temp.mate1.fastq
		    fi
		    if [ -e "$ProbandGenerator".temp.mate2.fastq ]; then
	                rm  "$ProbandGenerator".temp.mate2.fastq
	            fi
		    if [ -e "$ProbandGenerator".temp ]; then
			    rm "$ProbandGenerator".temp
	            fi
		    echo "running this one "
		    mkfifo "$ProbandGenerator".temp.mate1.fastq "$ProbandGenerator".temp.mate2.fastq
		    sleep 1
		      bash "$ProbandGenerator" | "$RDIR"/bin/PassThroughSamCheck.stranded "$ProbandGenerator".filter.chr  "$ProbandGenerator".temp >  "$ProbandGenerator".temp &
		       $RUFUSfilterFASTQ  "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$ProbandGenerator".temp.mate1.fastq "$ProbandGenerator".temp.mate2.fastq "$ProbandGenerator" "$K" $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)" &

		    wait
		else
			echo "Running RUFUS.filter from paired FASTQ files"
			FileName=$(basename $_arg_fastqA)
			Extension="${FileName##*.}"
			if [[ $Extension == 'gz' ]]
			then
				echo "Compressed fastq files found"
				$RUFUSfilterFASTQ "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList  <(zcat $_arg_fastqA) <(zcat $_arg_fastqB) "$ProbandGenerator" $K $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)"

			else
				echo "Uncompressed fastq files found"
				$RUFUSfilterFASTQ "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList  $_arg_fastqA $_arg_fastqB "$ProbandGenerator" $K $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)"
			fi
			wait
		fi
	fi

    if [ $(head "$ProbandGenerator".Mutations.Mate1.fastq | wc -l | awk '{print $1}') -eq "0" ]; then
		echo "ERROR: No mutant fastq reads idenfied.  Either the files are exactly the same of something went wrong in previous step"
		exit 100
	fi

	shortinsert="false"
	if [ -e "$ProbandGenerator".Mutations.fastq.bam ]
	then
		echo "skipping mapping mates"
	else
        # Sort fastq mates
        sortedMate1Fastq="$ProbandGenerator".sorted.Mutations.Mate1.fastq
        sortedMate2Fastq="$ProbandGenerator".sorted.Mutations.Mate2.fastq

        cat "$ProbandGenerator".Mutations.Mate1.fastq | paste - - - - | sort -k1 -S 8G | tr "\t" "\n" > $sortedMate1Fastq
        cat "$ProbandGenerator".Mutations.Mate2.fastq | paste - - - - | sort -k1 -S 8G | tr "\t" "\n" > $sortedMate2Fastq

		if [ $shortinsert = "false" ]
		then
			echo "skipping fastp fix"
	                $bwa mem -t $Threads $_arg_ref_bwa $sortedMate1Fastq $sortedMate2Fastq | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.bam
	                samtools index "$ProbandGenerator".Mutations.fastq.bam
		else
			echo "using fastp fix"
	        $fastp -i $sortedMate1Fastq -I $sortedMate2Fastq -m -o "$ProbandGenerator".Mutations.Mate1.fastq.fastp.fastq -O "$ProbandGenerator".Mutations.Mate2.fastq.fastp.fastq --merged_out "$ProbandGenerator".Mutations.Mate1.fastq.merged.fastq
			$bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.Mate1.fastq.fastp.fastq "$ProbandGenerator".Mutations.Mate2.fastq.fastp.fastq  | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.pared.bam
			$bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.Mate1.fastq.merged.fastq  | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.merged.bam
			samtools merge "$ProbandGenerator".Mutations.fastq.bam "$ProbandGenerator".Mutations.fastq.merged.bam "$ProbandGenerator".Mutations.fastq.pared.bam
			samtools index "$ProbandGenerator".Mutations.fastq.merged.bam
			samtools index "$ProbandGenerator".Mutations.fastq.pared.bam
			samtools index "$ProbandGenerator".Mutations.fastq.bam
		fi
	fi
else
#########put se pipe here
	if [ -e "$ProbandGenerator".Mutations.fastq ]
	then
		echo "skipping filter"
	else
		if [ -z $_arg_fastqA ]
		then

		    echo "running this one filer SE"
	            sleep 1
	            if [ -e "$ProbandGenerator".temp ]; then
	                            rm  "$ProbandGenerator".temp
	            fi
	                mkfifo "$ProbandGenerator".temp
	              bash "$ProbandGenerator" | "$RDIR"/bin/PassThroughSamCheck.stranded.se "$ProbandGenerator".filter.chr  "$ProbandGenerator".temp >  "$ProbandGenerator".temp &
	               $RUFUSfilterFASTQse  "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$ProbandGenerator".temp  "$ProbandGenerator" "$K" $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)" &
		    wait
		else
			echo "Running RUFUS.filter from single FASTQ files"
			echo "havent written this yet EXITing"
			exit
			#########WRITE THIS##########
			wait
		fi
	fi

	#if [ $(wc -l "$ProbandGenerator".Mutations.fastq | awk '{print $1}') -eq "0" ]; then
	if [ $(head "$ProbandGenerator".Mutations.fastq | wc -l  | awk '{print $1}') -eq "0" ]; then
		echo "ERROR: No mutant fastq reads idenfied.  Either the files are exactly the same of something went wrong in previous step"
		exit 100
	fi

	shortinsert="false"
	if [ -e "$ProbandGenerator".Mutations.fastq.bam ]
	then
		echo "skipping mapping mates"
	else
	                $bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.fastq | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.bam
	                samtools index "$ProbandGenerator".Mutations.fastq.bam

	fi

fi

########################################################################################
if [ $_arg_saliva == "TRUE" ]
then
	echo "saliva sample provided, only using aligned mutant contigs"
	if [ -e  "$ProbandGenerator".Mutations.fastq.FULL.bam ]
	then
		echo "skipping saliva filter"
	else

		mv "$ProbandGenerator".Mutations.fastq.bam "$ProbandGenerator".Mutations.fastq.FULL.bam
		samtools index "$ProbandGenerator".Mutations.fastq.FULL.bam
		rm "$ProbandGenerator".Mutations.fastq.bam.bai
		samtools view -F 12 -b "$ProbandGenerator".Mutations.fastq.FULL.bam > "$ProbandGenerator".Mutations.fastq.bam
		samtools index "$ProbandGenerator".Mutations.fastq.bam
	fi
fi



if [ $( samtools view "$ProbandGenerator".Mutations.fastq.bam | head | wc -l | awk '{print $1}') -eq "0" ]; then
        echo "ERROR: BWA failed on "$ProbandGenerator".Mutations.fastq.  Either the files are exactly the same of something went wrong in previous step"
        exit 100
fi
#################################################################################
if [ "$_arg_stop" = "filter" ];
then
        echo "-StF used, stopping run";
        exit 1;
fi
###################__RUFUS_OVERLAP__#############################################
if [ -e $ProbandGenerator.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz ]
then
    echo "########### Skipping overlap step ###########"
else
    echo "########### Starting RUFUS overlap ###########"
    echo " bash  $RUFUSOverlap "$_arg_ref" "$ProbandGenerator".Mutations.fastq $MutantMinCov $ProbandGenerator "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$K" "$Threads" "$_MaxAlleleSize" "$ProbandGenerator".Jhash "$parentsString" "$_arg_ref_bwa" "$_arg_refhash""
     bash  $RUFUSOverlap "$_arg_ref" "$ProbandGenerator".Mutations.fastq $MutantMinCov $ProbandGenerator "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$K" "$Threads" "$_MaxAlleleSize" "$_assemblySpeed" "$ProbandGenerator".Jhash "$parentsString" "$_arg_ref_bwa" "$_arg_refhash"
    #bash  $RUFUSOverlap "$_arg_ref" "$ProbandGenerator".Mutations.fastq 3 $ProbandGenerator "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$K" "$Threads" "$_MaxAlleleSize" "$_assemblySpeed" "$ProbandGenerator".Jhash "$parentsString" "$_arg_ref_bwa" "$_arg_refhash"
    echo "Done with RUFUS overlap"
fi
##############################################################################################


############################__RUFALU__#############################
#aluList=$RDIR/resources/primate_non-LTR_Retrotransposon.fasta
#fastaHackPath=$RDIR/bin/externals/fastahack/src/fastahack_project/bin/tools/fastahack
#jellyfishPath=$RDIR/src/externals/jellyfish-2.2.5/bin/jellyfish
#echo "$RufAlu $ProbandFileName $ProbandGenerator.V2.overlap.hashcount.fastq  $aluList $_arg_ref $jellyfishPath $(echo $ParentFileNames) "
#$RufAlu $_arg_subject $_arg_subject.generator.V2.overlap.hashcount.fastq  $aluList $_arg_ref $fastaHackPath $jellyfishPath  $(echo $ParentFileNames)
########################################################################


echo "cleaning up VCF"

PREFINAL_VCF="$ProbandGenerator.V2.overlap.hashcount.fastq.bam.coinherited.vcf"

grep ^# $ProbandGenerator.V2.overlap.hashcount.fastq.bam.vcf> ./Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf
grep -v  ^# $ProbandGenerator.V2.overlap.hashcount.fastq.bam.vcf | sort -k1,1V -k2,2n >> ./Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf
echo "arg_mosaic = $_arg_mosaic"
if [ "$_arg_mosaic" == "TRUE" ]
then
	echo "including mosaic";
	bash $RDIR/scripts/VilterAutosomeOnly ./Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf | perl $RDIR/scripts/ColapsDuplicateCalls.stream.pl > ./$PREFINAL_VCF
	#todo: guessing this is asynch because of stream in perl script title - which causes the next line to run before the file is created
	#todo: instead will incorporate 1mb mode, trim and combine, then filter inheriteds
else
	echo "excluding mosaic";
	bash $RDIR/scripts/VilterAutosomeOnly.withoutMosaic ./Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf | perl $RDIR/scripts/ColapsDuplicateCalls.stream.pl > ./$PREFINAL_VCF
fi