#!/bin/bash

HASH_TYPE="$1"
HASH_VERSION="$2"
FMTD_REGION="$3"

if [ ! -d "/mnt/rufus_resources" ]; then
    mkdir -p "/mnt/rufus_resources"
fi

chrom=$(echo "$FMTD_REG" | cut -d'_' -f1)

if [ "$REGION" == "wg" ]; then
    aws s3 sync "s3://rufus.marth.lab/public_access_data/${HASH_TYPE}_hashes/${HASH_VERSION}/wg_${HASH_TYPE}_${HASH_VERSION}.Jhash" "/mnt/rufus_resources/wg_${HASH_TYPE}_${HASH_VERSION}.Jhash"
else
    aws s3 sync "s3://rufus.marth.lab/public_access_data/${HASH_TYPE}_hashes/${HASH_VERSION}/${chrom}/${FMTD_REG}_${HASH_TYPE}_${HASH_VERSION}.Jhash" "/mnt/rufus_resources/${FMTD_REG}_${HASH_TYPE}_${HASH_VERSION}.Jhash"
fi