RUFUS Singularity Container
=====

K-mer based variant detection. 

Developed by Stephanie Georges, MS\
Based on the thesis project of Andrew Farrell, PhD\
Marth Lab, University of Utah\
*Publication in preparation*

For questions and feature requests, please contact [stephanie.georges@genetics.utah.edu](stephanie.georges@genetics.utah.edu)

## RUFUS Overview

RUFUS is a reference-free, K-mer based variant detection algorithm, for short-read DNA sequence data. RUFUS is intended to run on a high performance computing (HPC) cluster with singularity installed to take advantage of parallelism. At a high level, you'll need to download the pre-built singularity container (detailed below) and set up a batch script corresponding to your resource manager. If your HPC system uses SLURM, you can utilize the provided helper functions to create your SBATCH scripts (detailed below).

RUFUS currently supports a single subject sample, multiple control samples, and only accepts GRCh38 as a reference genome. The subject sample may be in BAM, CRAM, FASTQ, or FASTQ pair format. The samples must be in BAM format (though may be unaligned). The reference genome must be in FASTA format, and must be indexed by BWA. If the BWA indexes are not detected in the same directory as the reference genome, RUFUS will create them.

RUFUS has two stages: a variant calling stage, and a post-processing stage. Separation of the stages is necessary because the calling stage may be run in a windowed fashion, requiring multiple parallel RUFUS jobs over all of the windows. The combination stage must wait to proceed until all calling jobs are complete. Algorithmic runtime increases roughly linearly with sample coverage. Generally with whole-genome mode, a 100x sample run will take 1 day. Windowed mode completes significantly faster.


## Running RUFUS

### Obtaining the RUFUS Singularity Image

 The pre-built RUFUS singularity container may be obtained from [Zenodo](https://zenodo.org/records/13694211). To download:
```
curl "https://zenodo.org/records/13694211/files/rufus.sif" -o rufus.sif
```

### Input Data

RUFUS requires the following data to run:
1) A subject sample in BAM format (this may be unaligned)
2) One or more control samples in BAM format (these may be unaligned)
3) A reference fasta file (this must be indexed by BWA) - for use in reporting the called variants. *It's recommended to provide the BWA indexes in the same data directory if you have them to save time creating them during the RUFUS run.*
   * To create the BWA indexes, run the following commands:
   ```
   bwa index -a bwtsw {REFERENCE.fa}
   samtools faidx {REFERENCE.fa}
   ```

**All of the above files must be located in a single directory, which will be mounted to the singularity container.**

### Output Data

RUFUS will, by default, output the following files *in the same bound directory containing the input data*:
1) A VCF file containing the called variants
2) A supplemental directory with:
    * A pre-filtered VCF file
    * A BAM file containing the raw reads containing the mutant kmers
    * A BAM file containing the assembled contigs from the raw reads containing the mutant kmers
    * A hash table containing the unique subject kmers and their counts

 
### The Two Stages of RUFUS

RUFUS has two execution stages:
1) The calling stage, invoked by the following
```
singularity exec --bind {PATH_TO_LOCAL_DATA_DIR}:/mnt {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/runRufus.sh [-s|--subject <arg>] [-r|--ref <arg>] [-t|--threads <arg>] [-k|--kmersize <arg>] [-m|--min <arg>] [-h|--help] [-c|<controls-1>] ... [-c|<controls-n>] ...OPTIONS
```
With the following usage:
```
    -s,--subject: bam/cram/fastq(or pair of fastq files)/generator file containing the subject of interest (no default, only one subject per run for now)"
    -c, --controls: bam/cram/fastq(or pair of fastq files)/generator file for the sequence data of the control sample (can be used multipe times)"
    -e,--exclude: Jhash file of kmers to exclude from mutation list, k must be  (no default, can be used multiple times)"
    -se,--single_end_reads: subject bam file is single end reads, not paired (default is to assume paired end data)"
    -r,--ref: file path to the desired reference file (no default)"
    -cr,--cramref: file path to the desired reference file to decompress input cram files (no default)"
    -t,--threads: number of threads to use (no default) (min 3)"
    -k,--kersize: size of k-mer to use (no default)"
    -m,--min: overwrites the minimum k-mer count to call variant (no default)"
    -i,--saliva: flag to indicate that the subject sample is a buccal swab and likely contains a significant fraction of contaminant DNA"
    -mx,--MaxAllele: Max size for insert/deletion events to put the entire alt sequence in. (default 1000)"
    -L,--Report_Low_Freq: Reprot Mosaic/Low Frequency/Somatic variants (default FALSE)"
    -h,--help: Print help"
```

2) The post-processing stage, invoked by the following
```
singularity exec --bind {PATH_TO_LOCAL_DATA_DIR}:/mnt {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/post_process/post_process.sh [-w window_size] [-r reference] [-subject] [-c control1,control2,control3...] [-d source_dir]
```
With the following usage:
```
    echo "Options:"
    echo " -w window_size   Required: The size of the window used in the RUFUS run"
    echo " -r reference Required: The reference used in the RUFUS run"
    echo " -c controls  Required: The control bam files used in the RUFUS run"
    echo " -s subject_file  Required: The name of the subject file: must be the same as that supplied to the RUFUS run"
    echo " -d source_dir    Required: The source directory where the RUFUS vcf(s) are located"
    echo " -h help  Print help message"
```


## Using the SLURM Helper Script & Executing the SLURM Batch Scripts

The SLURM helper script automatically creates the two SLURM batch scripts necessary to run RUFUS on a SLURM-managed HPC cluster, as well as a bash script to execute them. To use:
1) Execute the helper script (see full usage options below):
``` 
singularity exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/singularity/launch_container.sh [-s subject] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...OPTIONS
```

2) Then execute the generated bash script:
```
bash launch_rufus.sh
```

\
The full usage options for the helper script are as follows:
```
Required Arguments:
    -d data_directory The directory containing the subject, control, and reference files to be used in the run
    -s subject    The subject sample of interest; must be located in data_directory
    -c control(s) A single control or comma-delimited array of multiple controls; must be located in data_directory
    -b genome_build  The desired genome build; currently only supports GRCh38
    -r reference  The reference file matching the genome build; must be located in data_directory
    -a slurm_account  The account for the slurm job
    -p slurm_partition    The partition for the slurm job
    -l slurm_job_array_limit    The maximum amount of jobs slurm allows in an array
    
Optional Arguments:
    -r reference  If not provided, will automatically provide GRCh38 file
    -m kmer_depth_cutoff  The amount of kMers that must overlap the variant to be included in the final call set
    -w window_size    The size of the windows to run RUFUS on, in units of kilabases (KB); allowed range between 500-5000; defaults to single run of entire genome if not provided
    -e email  The email address to notify with slurm updates
    -q slurm_job_queue_lgimit    The maximum amount of jobs able to be ran at once; defaults to 20
    -t slurm_time_limit   The maximum amount of time to let the slurm job run; defaults to 7 days for full run, or one hour per window (DD-HH:MM:SS)
    -f path_to_rufus_container    If not provided, will look in current directory for rufus.sif
    -z rufus_threads  Number of threads provided to RUFUS; defaults to 36
    -h help   Print usage
```
\
*Notes on SLURM arguments*:\
This script utilizes SLURM arrays to batch RUFUS call runs and thus requires the SLURM job array limit to comply with user settings. To find your SLURM job array limit:
```
scontrol show config | grep "MaxArraySize"
```

To maximize parallelism, filling in the slurm job queue limit (-q) is recommended. You can find your limit by typing:
```
scontrol show config | grep "default_queue_depth"
```

#### Example Invocation of the helper script
```
singularity exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/launch_container.sh -d /home/my_data_dir/ -s subject.bam -c control_a.bam, control_b.bam -r GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -t "00:30:00" -m 5 -l 20 -z 36 -e "my_email@utah.edu" -f /home/my_container_path/
```
=======

