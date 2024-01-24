#!/bin/bash
File=$1
FinalCoverage=$2
NameStub=$3.V2
HashList=$4
HashSize=$5
Threads=$6
 
SampleJhash=$7
Parent1Jhash=$8
Parent2Jhash=$9
Parent3Jhash=$10
Parent3Jhash=$11 

echo " you gave
File=$1
FinalCoverage=$2
NameStub=$3.V2
HashList=$4
HashSize=$5
Threads=$6

SampleJhash=$7
Parent1Jhash=$8
Parent2Jhash=$9
Parent3Jhash=$10
Parent3Jhash=$11
"

mkdir ./TempOverlap/
echo "Overlaping $File"

RDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/RUFUS

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
gkno=$RDIR/bin/gkno_launcher/gkno
samtools=$RDIR/bin/gkno_launcher/tools/samtools/samtools
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
humanRef=$RDIR/bin/gkno_launcher/resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa
CheckHash=$RDIR/cloud/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam




#if [ -s $NameStub.overlap.hashcount.fastq ]
#then 
#	echo "Skipping Overlap"
#else

if [ -s ./TempOverlap/$NameStub.sam.fastqd ]  
then 
	echo "skipping sam assemble"
else
	$gkno bwa-se -ps human  -q $File -id $File -s $File -o $File.bam -p ILLUMINA
 	$OverlapSam <($samtools view $File.bam ) .95 50 5 ./TempOverlap/$NameStub.sam $NameStub 1 $Threads
fi
if [ -s ./TempOverlap/$NameStub.1.fastqd ]
then 	
	echo "skipping first overlap"
else
	time $OverlapHash ./TempOverlap/$NameStub.sam.fastqd .98 50 2 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out
fi 
if [ -s ./TempOverlap/$NameStub.2.fastqd ]
then 
	echo "skipping second overlap"
else
	time $OverlapHash ./TempOverlap/$NameStub.1.fastqd .98 50 2 FP 15 1 ./TempOverlap/$NameStub.2 1 $Threads #>>  $File.overlap.out
fi
$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.2.fastqd > ./TempOverlap/$NameStub.3.fastqd
if [ -s ./TempOverlap/$NameStub.4.fastqd ]
then 
	echo "skipping overlap 4"
else 
	time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .95 30 0 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads #>>  $File.overlap.out
fi
if [ -s ./TempOverlap/$NameStub.5.fastqd ]
then 
	echo "skipping ovelrap 5"
else
	time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .95 30 $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads #>>  $File.overlap.out
fi 

if [ -s ./$NameStub.overlap.hashcount.fastq ]
then 
	echo "skipping final overlap work"
else

	$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
	$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
	$AnnotateOverlap $HashList ./$NameStub.overlap.fastq $NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq 
fi


if [ -s ./$NameStub.overlap.hashcount.fastq.bam ]
then 
	echo "skipping contig alignment" 
else
	
	$gkno bwa-se -ps human  -q ./$NameStub.overlap.hashcount.fastq -id ./$NameStub.overlap.hashcount.fastq -s ./$NameStub.overlap.hashcount.fastq -o ./$NameStub.overlap.hashcount.fastq.bam -p ILLUMINA
fi 


$RDIR/bin/bedtools2/bin/fastaFromBed -bed <( $RDIR/bin/bedtools2/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $RDIR/bin/gkno_launcher/resources/homo_sapiens/current/human_reference_v37.fa -fo $NameStub.overlap.asembly.hash.fastq.ref.fastq 
$RDIR/bin/jellyfish/bin/jellyfish count -C -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq
$RDIR/bin/jellyfish/bin/jellyfish count -C -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq
$RDIR/bin/jellyfish/bin/jellyfish dump -c ./$NameStub.overlap.hashcount.fastq.Jhash > ./$NameStub.overlap.hashcount.fastq.Jhash.tab
$RDIR/bin/jellyfish/bin/jellyfish dump -c ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab

if [ -s $NameStub.overlap.asembly.hash.fastq.p1 ]
then
        echo "skipping hash lookup"
else
	echo "stargin hash lookup"
        bash $CheckHash $SampleJhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.sample
        bash $CheckHash $Parent1Jhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.p1
        bash $CheckHash $Parent2Jhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.p2
	bash $CheckHash $Parent3Jhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.p3
	echo "done with hash lookup"
fi


if [ -s $NameStub.overlap.asembly.hash.fastq.Ref.sample ]	
then 
	echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.sample
	bash $CheckHash $Parent1Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.p1
	bash $CheckHash $Parent2Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.p2
	bash $CheckHash $Parent3Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.p3
	 
fi 

wait

mkfifo check 
$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p1 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p2 -sR $NameStub.overlap.asembly.hash.fastq.Ref.sample -s $NameStub.overlap.asembly.hash.fastq.sample 

grep ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
grep -v  ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf | sort -k1,1 -k2,2n >> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
$RDIR/bin/gkno_launcher/tools/tabix/bgzip -f ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
$RDIR/bin/gkno_launcher/tools/tabix/tabix ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.gz
