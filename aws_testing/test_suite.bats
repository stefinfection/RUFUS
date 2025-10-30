#/usr/bin/env bats

# This is a BATS test suite for AWS-related functionalities in the RUFUS project.

@test "Make sure Singularity is installed" {
    run singularity --version
    [ "$status" -eq 0 ]
    [[ "$output" == Singularity* ]]
}

@test "Check RUFUS container exists" {
    RUFUS_CONTAINER_PATH="${RUFUS_CONTAINER_PATH:-./rufus.sif}"
    [ -f "$RUFUS_CONTAINER_PATH" ]
    [ "$status" -eq 0 ]
}

