#!/bin/bash
# Returns single string of arguments for RUFUS invocation according to config yaml options.
# Run within Docker container.
# Called only for GNU Parallel execution for now.
# SJ Georges Nov2025

# Import cleaned env file
ENV_FILE="/mnt/rufus_temp/cleaned.env"
set -a
source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
set +a

source "$EXEC_HELPERS" # TODO: is this sourced correctly?

parallel_hash_args=""
parallel_hash_args=$(get_rufus_args) # TODO: if change this to array of strings update out here
echo "$parallel_hash_args"