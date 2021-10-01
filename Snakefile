# Snakemake workflow file
# This workflow demultiplexes custom "inner" barcodes from MinION plates.
# The intention is for "quick" checks during sequencing runs.

import subprocess
from glob import glob
import re

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = int(subprocess.check_output("nproc").decode())

def get_reads(wildcards):
    return glob(f"data/reads/barcode{wildcards.i}/*.fastq.gz")

# generate barcoding tag files for generic plates
rule plate_tags:
    output:
        tags_ITS1 = "tags/ITS1_tags.fasta",
        tags_3NDf_LR5 = "tags/3NDf_LR5_tags.fasta"
    input:
        tags_3NDf = "tags/3NDf_barcodes.fasta",
        tags_ITS1_LR5 = "tags/its1_lr5_barcodes.fasta",
        script = "scripts/tags.R"
    conda: "conda/tags.yaml"
    threads: 1
    script: "{input.script}"

# generate barcoding tag files for labeled samples
rule sample_tags:
    output:
        sample_tags = "tags/barcode{i}.fasta"
    input:
        tags_3NDf = "tags/3NDf_barcodes.fasta",
        tags_ITS1_LR5 = "tags/its1_lr5_barcodes.fasta",
        tag_plate = "tags/3NDf-LR5_tagplate.xlsx",
        sample_plate = "samples/barcode{i}.xlsx",
        script = "scripts/tags.R"
    conda: "conda/tags.yaml"
    threads: 1
    script: "{input.script}"


rule demultiplex:
    output:
        summary = "data/demultiplex/barcode{i}/barcode{i}.summary",
        unknowns = "data/demultiplex/barcode{i}/unknown.fastq.gz"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}"
    log: "logs/demux_barcode{i}.log"
    threads: maxthreads
    shell:
        """
        mkdir -p {params.outdir}
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

rule demultiplex_single:
    output:
        summary = "data/demultiplex/barcode{i}/single/barcode{i}.summary"
    input:
        reads = "data/demultiplex/barcode{i}/unknown.fastq.gz",
        tags = "tags/barcode{i}_single.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}/single"
    log: "logs/demux_single_barcode{i}.log"
    threads: maxthreads
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
