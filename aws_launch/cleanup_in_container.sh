#!/bin/bash

dirs_to_cleanup=(
    "/mnt/rufus_supplementals/downlaoded_kg1_hashes"
    "/mnt/rufus_supplementals/downlaoded_control_hashes"
)

files_to_cleanup=(
    "/mnt/rufus_supplementals/process_region_worker.sh"
)

for cdir in "${dirs_to_cleanup[@]}"; do
    if [ -d "$cdir" ]; then
        rm -r $cdir
    fi
done

for file in "${files_to_cleanup[@]}"; do
    if [ -f "$file" ];
        rm $file
    fi
done