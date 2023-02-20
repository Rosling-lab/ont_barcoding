#!/bin/bash -l
#SBATCH -A naiss2023-5-37
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J ont_barcoding_2023
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
