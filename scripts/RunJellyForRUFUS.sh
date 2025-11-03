#!/bin/sh
set -e
GEN=$1
K=$2
T=$3
L=$4
HASH_SIZE=$5

RDIR=/opt/RUFUS/

JELLYFISH="$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish"

# If we're using a region-specific hash, adjust size accordingly (1MB hashes made w/ 1G)
if [ -e "$GEN.Jhash" ]
then
	echo "Skipping jelly, $GEN.Jhash alreads exists"
else
	echo "Running jellyfish for $GEN"
	if [ -e $GEN.Jhash.temp ]; then 
		rm $GEN.Jhash.temp
	fi
	mkfifo $GEN.Jhash.temp
	if [ -e $GEN.fq ]; then 
		rm $GEN.fq
	fi
	mkfifo $GEN.fq
	bash $GEN | $RDIR/bin/PassThroughSamCheck $GEN.Jelly.chr > $GEN.fq &

	# -C is canonical ("Count both strand, canonical representation")
	# -L is filtering out low frequency kmers ("Don't output k-mer with count < lower-count")
	# For subject, we keep kmers with 2+ counts
	# For controls, we keep kmers with 2+ counts OR the provided argument to rufus (_argParLowK)
	# These arguments combined, I interpret this as keeping kmers with a single read
	# --disk means hash will be written to disk if entire thing can't be held in memory
	# -s (intial hash size) is G + Gcek (genome size * coverage * error * kmer length) ~228G for 300x, 22.8G for 30x, etc
	# guessing this starting number is far too low and there's a lot of memory swapping happening here
	# good area of parallelization and possible merging after - will neeed to think through
	
	echo "about to count jellyfish for $GEN"
	echo "$JELLYFISH count --disk -m $K -L $L -s $HASH_SIZE -t $T -o $GEN.Jhash -C $GEN.fq"
	$JELLYFISH count --disk -m $K -L $L -s $HASH_SIZE -t $T -o $GEN.Jhash -C $GEN.fq
	rm $GEN.Jhash.temp
	rm $GEN.fq

	wait
fi

if [ ! -s  $GEN.Jhash.histo ]; then 
	$JELLYFISH histo -f -o $GEN.Jhash.histo $GEN.Jhash
fi
if [ $(awk '$2 > 0' $GEN.Jhash.histo | wc -l ) -eq "0" ]; then  
	echo "ERROR: jellyfish failed on the file $GEN"
	exit 100
fi

exit 0
