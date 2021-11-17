# Snakemake workflow file
# This workflow demultiplexes custom "inner" barcodes from MinION plates.
# The intention is for "quick" checks during sequencing runs.

import subprocess
from glob import glob
import re
import snakemake

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = snakemake.utils.available_cpu_count()

def get_reads(wildcards):
    return glob(f"data/reads/barcode{wildcards.i}/*.fastq.gz")

def get_samples(wildcards):
    with open(checkpoints.sample_tags.get(i = wildcards.i).output.sample_names) as f:
        return([s.rstrip('\n') for s in f.readlines()])

def get_demuxes(wildcards):
    return(expand("data/demultiplex/barcode{i}/{sample}.fastq.gz",
        i = wildcards.i,
        sample = get_samples(wildcards)
    ))

def get_consensus(wildcards):
    checkpoints.demultiplex.get(i = wildcards.i)
    return(expand("data/consensus/barcode{i}/{sample}/final_clusters.tsv",
        i = wildcards.i,
        sample = get_samples(wildcards)
    ))

rule all:
    input:
        expand("data/consensus/barcode{i}",
            i = [n[1][0] for n in snakemake.utils.listfiles("samples/barcode{i}.xlsx")]),
        expand("data/demultiplex/barcode{i}/single/barcode{i}.summary",
            i = [n[1][0] for n in snakemake.utils.listfiles("tags/barcode{i}_single.fasta")])

# generate barcoding tag files for generic plates
rule plate_tags:
    output:
        tags_ITS1_fasta = "tags/ITS1_tags.fasta",
        tags_3NDf_LR5_fasta = "tags/3NDf_LR5_tags.fasta",
        tags_ITS1_table = "tags/ITS1_tags.tsv",
        tags_3NDf_LR5_table = "tags/3NDf_LR5_tags.tsv"
    input:
        tags_3NDf = "tags/3NDf_barcodes.fasta",
        tags_ITS1_LR5 = "tags/its1_lr5_barcodes.fasta",
        script = "scripts/tags.R"
    conda: "conda/tags.yaml"
    threads: 1
    script: "scripts/tags.R"

wildcard_constraints:
    i = "\d+"

# generate barcoding tag files for labeled samples
checkpoint sample_tags:
    output:
        sample_tags_fasta = "tags/barcode{i}.fasta",
        sample_tags_table = "tags/barcode{i}.tsv",
        sample_names = "samples/barcode{i}.txt"
    input:
        tags_3NDf = "tags/3NDf_barcodes.fasta",
        tags_ITS1_LR5 = "tags/its1_lr5_barcodes.fasta",
        tag_plate = "tags/3NDf-LR5_tagplate.xlsx",
        sample_plate = "samples/barcode{i}.xlsx",
        script = "scripts/tags.R"
    conda: "conda/tags.yaml"
    threads: 1
    script: "scripts/tags.R"

checkpoint demultiplex_cutadapt:
    output:
        summary = "data/demultiplex/barcode{i}/dual_cutadapt/barcode{i}.summary",
        unknowns = "data/demultiplex/barcode{i}/dual_cutadapt/unknown.fastq.gz"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}/dual_cutadapt/",
        samples = get_demuxes
    log: "logs/demux_barcode{i}_dual_cutadapt.log"
    threads: maxthreads
    conda: "conda/demultiplex.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        for s in {params.samples}; do
            echo "" | gzip -c - >"$s"
        done
        zcat {input.reads} |
        cutadapt -m 400\
                 -o - - |
        cutadapt -g file:{input.tags}\
                 --revcomp\
                 -o {params.outdir}/{{name}}.fastq.gz\
                 -j {threads}\
                 -\
                 >{log}
        grep -B2 "[Tt]rimmed: [^0]" {log} >{output.summary}
        """

checkpoint demultiplex_minibar:
    output:
        summary = "data/demultiplex/barcode{i}/dual_minibar/barcode{i}.summary"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.tsv",
        minibar = "scripts/minibar/minibar.py"
    params:
        outdir = "data/demultiplex/dual_minibar/barcode{i}",
        samples = get_demuxes,
        temp_fastq = "data/reads/barcode{i}/all.fastq.gz"
    log: "logs/demux_barcode{i}_dual_minibar.log"
    threads: maxthreads
    conda: "conda/minibar.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        for s in {params.samples}; do
            echo "" | gzip -c - >"$s"
        done
        trap "rm -f {params.temp_fastq}" EXIT
        cat {input.reads} >{params.temp_fastq} &
        python {input.minibar}\
               {input.tags}\
               {params.temp_fastq}\
               -e 2\
               -E 5\
               -l 120\
               -F\
               -P {params.outdir}/\
               -fh &>{log}
        cp {log} {output.summary}
        """

rule demultiplex_single:
    output:
        summary = "data/demultiplex/barcode{i}/single_cutadapt/barcode{i}.summary"
    input:
        reads = "data/demultiplex/barcode{i}/dual_cutadapt/unknown.fastq.gz",
        tags = "tags/barcode{i}_single.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}/single_cutadapt"
    log: "logs/demux_barcode{i}_single_cutadapt.log"
    threads: maxthreads
    conda: "conda/demultiplex.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        cutadapt -g file:{input.tags}\
                 --revcomp\
                 -o {params.outdir}/{{name}}.fastq.gz\
                 -j {threads}\
                 {input.reads}\
                 >{log}
        grep -B2 "[Tt]rimmed: [^0]" {log} >{output.summary}
        """

rule rDNA_consensus:
    output: "data/consensus/barcode{i}/{demux_algo}/{sample}/final_clusters.tsv"
    input: "data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/{demux_algo}/{sample}",
        unzipped = "data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq"
    log: "logs/consensus_barcode{i}_{demux_algo}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        gunzip -k {input}
        NGSpeciesID\
          --ont\
          --consensus\
          --racon\
          --racon_iter 3\
          --fastq {params.unzipped}\
          --outfolder {params.outdir}\
          --m 2800\
          --s 500\
          --sample_size 300\
          --t 1
          >{log}
        rm {params.unzipped}
        """

rule ITS_consensus:
    output: "data/consensus/barcode{i}/single/{sample}/final_clusters.tsv"
    input: "data/demultiplex/barcode{i}/single/{sample}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/{sample}",
        unzipped = "data/demultiplex/barcode{i}/single/{sample}.fastq"
    log: "logs/consensus_barcode{i}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        gunzip -k {input}
        NGSpeciesID\
          --ont\
          --consensus\
          --racon\
          --racon_iter 3\
          --fastq {params.unzipped}\
          --m 800\
          --s 200\
          --sample_size 300\
          --t 1
          >{log}
        rm {params.unzipped}
        """

rule all_consensus:
    output: directory("data/consensus/barcode{i}")
    input: get_consensus
