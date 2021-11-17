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
    if wildcards.demux_algo == "single_cutadapt":
        with open(checkpoints.sample_tags.get(i = wildcards.i).output.sample_names_single) as f:
            return([s.rstrip('\n') for s in f.readlines()])
    else:
        with open(checkpoints.sample_tags.get(i = wildcards.i).output.sample_names) as f:
            return([s.rstrip('\n') for s in f.readlines()])

def get_demuxes(wildcards):
    return(expand("data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq.gz",
        i = wildcards.i,
        demux_algo = wildcards.demux_algo,
        sample = get_samples(wildcards)
    ))

def get_consensus(wildcards):
    if wildcards.demux_algo == "dual_minibar":
        print('Running get_consensus for dual_minibar.')
        my_checkpoint = checkpoints.demultiplex_minibar
    elif wildcards.demux_algo == "dual_cutadapt":
        print('Running get_consensus for dual_cutadapt.')
        my_checkpoint = checkpoints.demultiplex_cutadapt
    elif wildcards.demux_algo == "single_cutadapt":
        print('Running get_consensus for single_cutadapt.')
        my_checkpoint = checkpoints.demultiplex_single
    else:
        raise ValueError('unknown demux_algo:' + wildcards.demux_algo)
    my_checkpoint.get(i = wildcards.i, demux_algo = wildcards.demux_algo)
    samples = get_samples(wildcards)
    for s in samples:
        print(s)
    return(expand("data/consensus/barcode{i}/{demux_algo}/{sample}/consensus.fasta",
        i = wildcards.i,
        demux_algo = wildcards.demux_algo,
        sample = get_samples(wildcards)
    ))

rule all:
    input:
        expand("data/consensus/barcode{i}/{demux_algo}/barcode{i}.fasta",
            i = [n[1][0] for n in snakemake.utils.listfiles("samples/barcode{i}.xlsx")],
            demux_algo = ["dual_cutadapt", "dual_minibar", "single_cutadapt"],
        )

# generate barcoding tag files for generic plates
# rule plate_tags:
#     output:
#         tags_ITS1_fasta = "tags/ITS1_tags.fasta",
#         tags_3NDf_LR5_fasta = "tags/3NDf_LR5_tags.fasta",
#         tags_ITS1_table = "tags/ITS1_tags.tsv",
#         tags_3NDf_LR5_table = "tags/3NDf_LR5_tags.tsv"
#     input:
#         tags_3NDf = "tags/3NDf_barcodes.fasta",
#         tags_ITS1_LR5 = "tags/its1_lr5_barcodes.fasta",
#         script = "scripts/tags.R"
#     conda: "conda/tags.yaml"
#     threads: 1
#     script: "scripts/tags.R"

wildcard_constraints:
    i = "\d+"

# generate barcoding tag files for labeled samples
checkpoint sample_tags:
    output:
        sample_tags_fasta = "tags/barcode{i}.fasta",
        sample_tags_single_fasta = "tags/barcode{i}_single.fasta",
        sample_tags_table = "tags/barcode{i}.tsv",
        sample_names = "samples/barcode{i}.txt",
        sample_names_single = "samples/barcode{i}_single.txt"
    input:
        tags_3NDf = "tags/3NDf_barcodes.fasta",
        tags_ITS1_LR5 = "tags/its1_lr5_barcodes.fasta",
        sample_plate = "samples/barcode{i}.xlsx",
        script = "scripts/tags.R"
    conda: "conda/tags.yaml"
    threads: 1
    script: "scripts/tags.R"

checkpoint demultiplex_cutadapt:
    output:
        summary = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.summary",
        unknowns = "data/demultiplex/barcode{i}/{demux_algo}/unknown.fastq.gz"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}/{demux_algo}/",
        samples = get_demuxes
    wildcard_constraints:
        demux_algo = "dual_cutadapt"
    log: "logs/demux_barcode{i}_{demux_algo}.log"
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
        summary = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.summary"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.tsv",
        minibar = "scripts/minibar/minibar.py"
    params:
        outdir = "data/demultiplex/barcode{i}/{demux_algo}",
        samples = get_demuxes,
        temp_fastq = "data/reads/barcode{i}/all.fastq.gz"
    wildcard_constraints:
        demux_algo = "dual_minibar"
    log: "logs/demux_barcode{i}_{demux_algo}.log"
    threads: maxthreads
    conda: "conda/minibar.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        for s in {params.samples}; do
            echo "" | gzip -c - >"$s"
        done
        cat {input.reads} >{params.temp_fastq} &
        python {input.minibar}\
               {input.tags}\
               {params.temp_fastq}\
               -e 2\
               -E 5\
               -l 120\
               -F\
               -P {params.outdir}/\
               -fh\
               &>{log}
        gzip -f {params.outdir}/*.fastq
        cp {log} {output.summary}
        rm -f {params.temp_fastq}
        """

checkpoint demultiplex_single:
    output:
        summary = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.summary"
    input:
        reads = "data/demultiplex/barcode{i}/dual_cutadapt/unknown.fastq.gz",
        tags = "tags/barcode{i}_single.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}/{demux_algo}"
    wildcard_constraints:
        demux_algo = "single_cutadapt"
    log: "logs/demux_barcode{i}_{demux_algo}.log"
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
    output:
        cluster_map = "data/consensus/barcode{i}/{demux_algo}/{sample}/final_clusters.tsv",
        consensus = "data/consensus/barcode{i}/{demux_algo}/{sample}/consensus.fasta"
    input: "data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/{demux_algo}/{sample}",
        unzipped = "data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq"
    wildcard_constraints:
        demux_algo = "dual_(cutadapt|minibar)"
    log: "logs/consensus_barcode{i}_{demux_algo}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        gunzip -fk {input}
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
          --t 1\
          &>{log} || echo "Error in NGSpeciesID" >>{log}
        touch "{output.consensus}"
        for f in $(find {params.outdir} -path '*racon_cl_id_*/consensus.fasta') ; do
            echo "$f" | sed -r 's|.*/racon_cl_id(_[0-9]+)/consensus.fasta|>{wildcards.sample}\\1|' >>"{output.consensus}"
            tail -n+2 "$f" >>{output.consensus}
        done
        rm {params.unzipped}
        """

rule ITS_consensus:
    output:
        cluster_map = "data/consensus/barcode{i}/{demux_algo}/{sample}/final_clusters.tsv",
        consensus = "data/consensus/barcode{i}/{demux_algo}/{sample}/consensus.fasta"
    input: "data/demultiplex/barcode{i}/single_cutadapt/{sample}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/single_cutadapt/{sample}",
        unzipped = "data/demultiplex/barcode{i}/single_cutadapt/{sample}.fastq"
    wildcard_constraints:
        demux_algo = "single_(cutadapt|minibar)"
    log: "logs/consensus_barcode{i}_{demux_algo}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        gunzip -fk {input}
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
          &>{log} || echo "Error in NGSpeciesID" >>{log}
        touch "{output.consensus}"
        for f in $(find data/consensus/barcode01/dual_cutadapt/{wildcards.sample} -path '*racon_cl_id_*/consensus.fasta') ; do
            echo "$f" | sed -r 's|.*/racon_cl_id(_[0-9]+)/consensus.fasta|>{wildcards.sample}\\1|' >>"{output.consensus}"
            tail -n+2 "$f" >>{output.consensus}
        done
        rm {params.unzipped}
        """

rule all_consensus:
    output: "data/consensus/barcode{i}/{demux_algo}/barcode{i}.fasta"
    wildcard_constraints:
        demux_algo = "(single|dual)_(cutadapt|minibar)"
    input: get_consensus
    shell:
        """
        cat {input} >{output}
        """
