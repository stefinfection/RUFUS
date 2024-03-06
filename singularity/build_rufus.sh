export SINGULARITY_TMPDIR=/home/ubuntu/singularity_tmp
export SINGULARITY_BIND=/home/ubuntu/rufus_data:/mnt
sudo -E singularity build --sandbox rufus.sif rufus.def
