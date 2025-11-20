#!/bin/bash

# A file of functions that will be called by launch scripts BEFORE the container starts running

# Checks all required variables
# Some variables are required for singularity, some are required for docker, some are required for both
# Returns container_type to be passed to other setup scripts
# TODO: add check for index here for subject + controls if they exist
# check_req_args


# Check for correct control setup - either prebuilt hashes or control array filled, or both
# Returns mount clause for container - which will differ for singularity vs docker
# set_up_controls


# Same as above, but for 1000G
# Returns mount clause for container which again will differ for sing vs docker
# set_up_kg1

# Checks if all BWA indexes are present in the same folder as the reference
# If one or more is missing, sets build_refs to true
# Returns mount clause for container and build_refs variable
#set_up_ref()

