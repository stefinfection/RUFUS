RUFUS Singularity Container
=====

K-mer based variant detection. 

Developed by Stephanie Georges, MS
Based on the thesis project of Andrew Farrell, PhD               

## Tool Overview

RUFUS is intended to run on a high performance computing cluster to take advantage of parallelism. At a high level, you'll need to download the pre-built singularity container (detailed below) and set up a batch script corresponding to your resource manager. If your HPC system uses SLURM, you can utilize the built-in helper functions to create your SBATCH scripts. 

RUFUS has two stages: a variant calling stage, and a post-processing stage. Separation of the stages is necessary because the calling stage may be run in a windowed fashion, requiring multiple parallel RUFUS jobs over all of the windows. The combination stage must wait to proceed until all calling jobs are complete.


## Obtaining the RUFUS Singularity image

#TODO
 **1) Download container**
```
git clone https://github.com/marthlab/RUFUS
```

**2) Set up resource management script**
The RUFUS container has two execution stages:
1) The calling stage, invoked by the following
```
singularity exec --bind {PATH_TO_LOCAL_DATA_DIR}:/mnt {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/runRufus.sh [OPTIONS]
```
With the following usage:

```
d
```

2) The post-processing stage, invoked by the following
```
singularity exec --bind {PATH_TO_LOCAL_DATA_DIR}:/mnt {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/post_process/post_process.sh [OPTIONS]
```

## RUFUS Requirements

## Providing a reference file.

After RUFUS has identified reads containing mutant kmers, the reads must be aligned to a reference fasta file.  Any fasta file can be used as a reference, as long as the fasta file has been indexed for BWA.  If a fasta has been indexed by bwa, there will be reference files with the following extensions: pac, .ann, .abm, .bwt, sa.  In order to prepare a reference fasta for bwa, simply type:

```
bwa index -a bwtsw reference.fa
samtools faidx reference.fa
```

This will produce the BWA index files, and the fasta file index respectively.  Make sure that the bwa index files and the fasta index file are in the same directory as reference.fa

## Ubuntu dependencies
In order for RUFUS to run on a fresh Ubuntu build, all of the following packages must be installed:

**General**
```
sudo apt-get update
sudo apt-get install python
sudo apt-get install cmake
sudo apt-get install wget
```

**GCC-4.9** (c/c++ compiler)
```
sudo apt-get install build-essential
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get install g++-4.9
```

**zlib** (file compression library)
```
sudo apt-get install zlib1g-dev
sudo apt-get install libbz2-dev
```

**bzlib** (bz2 file compression library)
```
sudo apt-get install libbz2-dev
sudo apt-get install liblzma-dev
```

**bc** (floating point precision library)
```
sudo apt-get install bc
```

**Curse** (terminal control library)
```
sudo apt-get install libncurses5-dev
```

## To Build Singularity Container on EC2 instance
c6i.8xl with Ubuntu
[Install Singularity](https://github.com/apptainer/singularity/blob/master/INSTALL.md)
```
mkdir singularity_temp
export SINGULARITY_TMPDIR=/home/ec2-user/singularity-temp
sudo -E singularity build rufus.sif rufus.def
```

=======
