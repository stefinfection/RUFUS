RUFUS Singularity Container
=====

K-mer based variant detection. d1.1.8.

Developed by Stephanie Georges, MS\
Based on the thesis project of Andrew Farrell, PhD\
Marth Lab, University of Utah\
*Publication in preparation*

For questions and feature requests, please contact [stephanie.georges@genetics.utah.edu](stephanie.georges@genetics.utah.edu)

## RUFUS Overview

RUFUS is a reference-free, K-mer based variant detection algorithm, for short-read DNA sequence data. RUFUS is intended to run on a high performance computing (HPC) cluster with singularity installed. At a high level, you'll need to download the pre-built singularity container (detailed below) and set up a batch script corresponding to your resource manager. If your HPC system uses SLURM, you can utilize the provided helper functions to create your SBATCH scripts (detailed below).

RUFUS currently supports a single subject sample, multiple control samples, and only accepts GRCh38 as a reference genome. The samples must be in BAM format (though may be unaligned). The reference genome must be in FASTA format, and must be indexed by BWA. If the BWA indexes are not detected in the same directory as the reference genome, RUFUS will create them.

RUFUS has two stages: a variant calling stage, and a post-processing stage. Separation of the stages is necessary because the calling stage may be run in a windowed fashion, requiring multiple parallel RUFUS jobs over all of the windows. The combination stage must wait to proceed until all calling jobs are complete. Algorithmic runtime increases roughly linearly with sample coverage. Generally with whole-genome mode, a 100x sample run will take 1 day. Windowed mode completes significantly faster.


## Running RUFUS

### Obtaining the RUFUS Singularity Image

 The pre-built RUFUS singularity container may be obtained from [Zenodo](https://zenodo.org/records/18284901). To download:
```
curl "https://zenodo.org/records/18918769/files/rufus_d1.1.8.sif" -o rufus.sif
```

### Input Data

RUFUS requires the following data to run:
1) A subject sample in BAM/CRAM format (this may be unaligned if using whole genome mode)
2) One or more control samples in BAM/CRAM format (these may be unaligned)
3) A reference fasta file (this must be indexed by BWA) - for use in reporting the called variants. *It's recommended to provide the BWA indexes in the same directory as the reference to save time creating them during the RUFUS run.*\
\
To create the BWA indexes, run the following commands:
   ```
   bwa index -a bwtsw {REFERENCE.fa}
   samtools faidx {REFERENCE.fa}
   ```

All input files are specified by their full paths. The necessary host directories are automatically bind-mounted into the singularity container.

### Output Data

RUFUS will, by default, output the following files *in the current working directory*:
1) A VCF file containing the called variants
2) A supplemental directory with:
    * A pre-filtered VCF file
    * A BAM file containing the raw reads containing the mutant kmers
    * A BAM file containing the assembled contigs from the raw reads containing the mutant kmers
    * A hash table containing the unique subject kmers and their counts
3) In windowed (region) mode, a `region_status.log` file summarizing the outcome of each region (variants called, no variants found with reason, or error with exit code)

 
### The Two Stages of RUFUS

RUFUS has two execution stages:
1) The calling stage, invoked by the following
```
singularity exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/runRufus.sh [-s|--subject <arg>] [-r|--ref <arg>] [-t|--threads <arg>] [-k|--kmersize <arg>] [-m|--min <arg>] [-h|--help] [-c|<controls-1>] ... [-c|<controls-n>] ...OPTIONS
```
With the following usage:
```
Required Arguments:
    -s,--subject: single bam file (may be unaligned) containing the subject of interest
    -c,--controls: bam file (may be unaligned) for the sequence data of the control sample (can be used multiple times, e.g. -c control1 -c control2)
    -r,--ref: file path to the desired reference file
    -t,--threads: number of threads to use (min 3)

Optional Arguments:
    -k,--kmersize: length of k-mer to use (defaults to 25)
    -m,--min: overwrites the minimum k-mer depth count to call variant (defaults to 5)
    -e,--exclude: Jhash file of kmers to exclude from mutation list (can be used multiple times, e.g. -e Jhash1 -e Jhash2)
    -f,--refhash: Jhash file containing reference hashList
    -R,--region: genomic region to call variants on (e.g., chr1:1-1000000); used in windowed mode
    -h,--help: Print help
```

2) The post-processing stage, invoked by the following
```
singularity exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/post_process/post_process.sh -s <subject> -w <window_size>
```
With the following usage:
```
Required Arguments:
    -s subject_file  The name of the subject file: must be the same as that supplied to the RUFUS run
    -w window_size   The size of the window used in the RUFUS run (0 for whole genome mode)
Optional Arguments:
    -h help  Print help message
```

In windowed mode, the post-processing stage will print a region status summary showing how many regions called variants, how many had no variants (with breakdown by reason), and how many encountered errors. The final VCF header will include a `##RUFUS_runMode` line indicating region mode was used.


## Using the SLURM Helper Script & Executing the SLURM Batch Scripts

The SLURM helper script automatically creates the two SLURM batch scripts necessary to run RUFUS on a SLURM-managed HPC cluster, as well as a bash script to execute them. To use:
1) Execute the helper script (see full usage options below):
``` 
singularity exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh [-s subject] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...OPTIONS
```

2) Then execute the generated bash script:
```
bash launch_rufus.sh
```

The full usage options for the helper script are as follows:
```
Required Arguments:
    -s subject    Full path to the subject sample BAM/CRAM
    -c control(s) A single control or comma-delimited array of multiple controls (full paths)
    -b genome_build  The desired genome build; currently only supports GRCh38
    -r reference  Full path to the reference file matching the genome build
    -a slurm_account  The account for the slurm job
    -p slurm_partition    The partition for the slurm job
    -l slurm_job_array_limit    The maximum amount of jobs slurm allows in an array
    
Optional Arguments:
    -m kmer_depth_cutoff  The amount of kMers that must overlap the variant to be included in the final call set
    -w window_size    The size of the windows to run RUFUS on, in units of kilabases (KB); allowed range between 500-5000; defaults to single run of entire genome if not provided
    -f reference_hash   Jhash file containing reference kMer hash list
    -x exclude_hash     Single or comma-delimited list of Jhash file(s) containing kMers to exclude (static, same for all regions)
    -K kg1_hash_dir     Full path to directory of per-region KG1 Jhash files (files named *{region}*.Jhash)
    -G kg1_version      KG1 hash version to download from S3 (e.g., v3.0)
    -D ctrl_hash_dir    Full path to directory of per-region control Jhash files (files named *{region}*.Jhash)
    -V ctrl_version     Control hash version to download from S3 (e.g., v1.0)
    -M memory_per_call  How much memory to allot to the rufus calling stage job (e.g., 150G or 20G)
    -C cpus_per_call    How many cpus to allot to each rufus calling stage job
    -y path_to_rufus_container   If not provided, will look in current directory for rufus.sif
    -z rufus_threads  Number of threads provided to RUFUS; defaults to 36 for entire genome; 10 for 1MB windows
    -e email  The email address to notify with slurm updates
    -q slurm_job_queue_limit    The maximum amount of jobs able to be ran at once; defaults to 20
    -t slurm_time_limit   The maximum amount of time to let the slurm job run; defaults to 7 days for full run, or one hour per window (DD-HH:MM:SS)
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

#### Example Invocations of the helper script

Basic windowed mode:
```
singularity exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -c /home/controls/control_a.bam,/home/controls/control_b.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -t "00:30:00" -m 5 -l 20 -z 36 -e "my_email@utah.edu"
```

With local per-region hash directories:
```
singularity exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -c /home/controls/control_a.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -l 1000 -K /data/kg1_hashes/v3.0/ -D /data/ctrl_hashes/v1.0/
```

With S3-downloaded hashes (downloaded at setup time):
```
singularity exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -c /home/controls/control_a.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -l 1000 -G v3.0 -V v1.0
```

*Note*: `-K`/`-G` (KG1 hashes) and `-D`/`-V` (control hashes) are mutually exclusive per type. You may mix local and S3 across types (e.g., `-K /local/kg1/ -V v1.0`). Hash files in local directories must be named with the region string (e.g., `*chr1_1_1000000*.Jhash` for region `chr1:1-1000000`, or `*wg*.Jhash` for whole-genome mode).

=======
