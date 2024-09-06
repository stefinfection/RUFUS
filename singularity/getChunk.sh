#!/bin/bash

# Define arrays for chromosomes and lengths
chroms=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
lengths=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468 156040895 57227415)

# Define chunk size
chunkSize=1000000  # Adjust this value as needed

# Get chunk number from command line argument
chunkNum=$1

# Calculate chromosome and coordinates for the given chunk
chunkStart=0
chunkEnd=0
for ((i=0; i<${#chroms[@]}; i++)); do
  if (( chunkNum < lengths[i] / chunkSize )); then
    chunkStart=$((chunkNum * chunkSize))
    chunkEnd=$((chunkStart + chunkSize))
    echo "${chroms[i]}:${chunkStart}-${chunkEnd}"
    exit 0
  else
    chunkNum=$((chunkNum - lengths[i] / chunkSize))
  fi
done

echo "Chunk not found"