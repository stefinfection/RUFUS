#!/bin/sh

CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
RDIR=$CDIR/../

JellyFish=$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish
Jhash=$1 # List of kmers
HashList=$2 # File to count within
MinCov=$3
MaxCov=$4

# This searches the HashList for specific kmers listed in the Jhash argument, and counts them
timeout 1h $JellyFish query -s <(cat $HashList | awk '{print ">"$1"\n"$1}') $Jhash | awk -v var=$MinCov ' $2 >= var ' | awk -v var=$MaxCov ' $2 <= var '
#cat $HashList | awk '{print $1}' | $JellyFish query -i $Jhash | awk -v var=$MinCov ' $2 >= var ' | awk -v var=$MaxCov ' $2 <= var '
