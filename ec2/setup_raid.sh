#!/bin/bash

mdadm --create /dev/md0 --level=0 --raid-devices=4 /dev/nvme[1-4]n1
mkfs.xfs /dev/md0
mount /dev/md0 /mnt
chown ubuntu /mnt
