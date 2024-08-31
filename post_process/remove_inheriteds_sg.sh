#!/bin/bash

OUT_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/HapMap/rufus_runs/yu/kmer_tests/whole_genome/k12_one_perc_test
REF=/scratch/ucgd/lustre/work/marth/shared/references/human/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
RUFUS_VCF=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/HapMap/rufus_runs/yu/kmer_tests/whole_genome/k12_one_perc_test/hd_af.exon_filtered.combined.vcf.gz
SAMPLE_NAME="NA12878-NA18517-1percent-lib1-umi-liquid_tumor"
CONTROL_BAM=/scratch/ucgd/lustre-work/marth/u0880188/smaht/uw/hapmap_data/NA12878.sorted.bam
INHERITED_SLURM=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/helper_scripts/remove_inheriteds_SG.slurm

sbatch $INHERITED_SLURM $REF $RUFUS_VCF $SAMPLE_NAME $OUT_DIR $CONTROL_BAM
