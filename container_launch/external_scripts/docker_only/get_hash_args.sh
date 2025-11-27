#!/bin/bash
# Returns single string of arguments of prebuilt hashes according to config yaml options.
# Run within Docker container.
# Called only for GNU Parallel execution for now.
# SJ Georges Nov2025

# Import cleaned env file
ENV_FILE="/mnt/rufus_temp/cleaned.env"
set -a
source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
set +a

source "/opt/RUFUS/container_launch/internal_hash_helpers.sh"

parallel_hash_args=""

# TODO: implement
# Get KG1 hash arg

# Get control hash arg

# Combine into single string (or array of strings to printf at PoC) and return
echo "$parallel_hash_args"