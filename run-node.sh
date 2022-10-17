#!/bin/bash -l

#SBATCH -A snic2022-5-42
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0:20:00
#SBATCH -J ont_barcoding
#SBATCH -C usage_mail
#SBATCH -M rackham
#SBATCH --mail-type=ALL
#SBATCH --output="logs/snakemake-%j.log"
#SBATCH --error="logs/snakemake-%j.log"

module load conda bioinfo-tools snakemake &&
snakemake -pr --jobs $SLURM_JOB_CPUS_PER_NODE\
    --use-envmodules\
    --use-conda\
    --conda-frontend conda\
    --shadow-prefix /scratch
