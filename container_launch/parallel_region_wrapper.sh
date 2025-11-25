#!/bin/bash
# A lightweight wrapper for feeding in to parallel execution point of call.
# Contains hash fetching, if optioned, and RUFUS execution.
# Only used for GNU parallel method, not SLURM.
# For now, only uses Docker for invocation.
# Run directly on host, external to container
# SJ Georges N

# Arguments
CONTAINER_ID="$1"
REGION="$2"

echo "Running RUFUS for $REGION..."

# Get hash args
HASH_CMD="/opt/RUFUS/get_hash_args.sh" # TODO: if this ends up being array, printf instead of simple assign
hash_args=$(docker exec "$CONTAINER_ID bash -c $HASH_CMD") # TODO: will this output as expected?

# RUFUS call
reg_agnostic_args=get_rufus_args
RUFUS_CMD="/opt/RUFUS/runRufus.sh \
    $reg_agnostic_args \
    $hash_args \
    -r $REGION"

# Convert region to proper format
fmtd_reg=$(echo "$REGION" | tr ':-' '_')

# Execute with Docker
docker exec "$CONTAINER_ID bash -c \
    $RUFUS_CMD \
    > /mnt/rufus_supplementals/logs/${fmtd_reg}.out \
    2> /mnt/rufus_supplementals/logs/${fmtd_reg}.err"

# TODO: do I need to do any file cleanup here? For downloaded hashses?