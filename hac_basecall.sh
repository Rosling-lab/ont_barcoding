#!/bin/bash -l

#SBATCH -A snic2022-5-42
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -J ont_barcoding
#SBATCH -C usage_mail
#SBATCH -M snowy
#SBATCH --gpus=1
#SBATCH --gpus-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --output="logs/snakemake-%j.log"
#SBATCH --error="logs/snakemake-%j.log"

module load conda bioinfo-tools snakemake &&
snakemake -pr --jobs $SLURM_JOB_CPUS_PER_NODE\
    --use-envmodules\
    --use-conda\
    --conda-frontend conda\
    --shadow-prefix /scratch\
    all_hac
