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
    print('Getting samples for ' + wildcards.demux_algo + ".")
    if wildcards.demux_algo == "single_cutadapt":
        with open(checkpoints.sample_tags.get(i = wildcards.i).output.sample_names_single) as f:
            return [s.rstrip('\n') for s in f.readlines()]
    else:
        with open(checkpoints.sample_tags.get(i = wildcards.i).output.sample_names) as f:
            return [s.rstrip('\n') for s in f.readlines()]

# def get_demuxes(wildcards):
#     return expand("data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq.gz",
#         i = wildcards.i,
#         demux_algo = wildcards.demux_algo,
#         sample = get_samples(wildcards)
#     )

def get_consensus(wildcards):
    # if wildcards.demux_algo == "dual_minibar":
    #     print('Running get_consensus for dual_minibar.')
    #     my_checkpoint = checkpoints.demultiplex_minibar
    # elif wildcards.demux_algo == "dual_cutadapt":
    #     print('Running get_consensus for dual_cutadapt.')
    #     my_checkpoint = checkpoints.demultiplex_cutadapt
    # elif wildcards.demux_algo == "single_cutadapt":
    #     print('Running get_consensus for single_cutadapt.')
    #     my_checkpoint = checkpoints.demultiplex_single
    # else:
    #     raise ValueError('unknown demux_algo:' + wildcards.demux_algo)
    # my_checkpoint.get(i = wildcards.i, demux_algo = wildcards.demux_algo)
    samples = get_samples(wildcards)
    print("Samples:")
    for s in samples:
        print(s)
    consensus_files = expand("data/consensus/barcode{i}/{demux_algo}/{sample}/consensus.fasta",
        i = wildcards.i,
        demux_algo = wildcards.demux_algo,
        sample = samples
    )
    print("Consensus files:")
    for f in consensus_files:
        print(f)
    return(consensus_files)

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

rule demultiplex_cutadapt:
    output:
        summary = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.summary",
        demux = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.fastq.gz",
        unknowns = "data/demultiplex/barcode{i}/{demux_algo}/unknown.fastq.gz"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.fasta"
    params:
        outdir = "data/demultiplex/barcode{i}/{demux_algo}/"
    wildcard_constraints:
        demux_algo = "dual_cutadapt"
    log: "logs/demux_barcode{i}_{demux_algo}.log"
    threads: maxthreads
    conda: "conda/demultiplex.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        zcat {input.reads} |
        cutadapt -m 400\\
                 -o - - |
        cutadapt -g file:{input.tags}\\
                 --revcomp\\
                 -o {output.demux}\\
                 --untrimmed-output {output.unknowns}\\
                 --rename='{{header}} sample:{{adapter_name}};'\\
                 -j {threads}\\
                 - |
                 tee {log}
        grep -B2 "[Tt]rimmed: [^0]" >{output.summary}
        """

rule demultiplex_minibar:
    output:
        summary = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.summary",
        demux = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.fastq.gz"
    input:
        reads = get_reads,
        tags = "tags/barcode{i}.tsv",
        minibar = "scripts/minibar/minibar.py"
    params:
        outdir = "data/demultiplex/barcode{i}/{demux_algo}",
        temp_fastq = "data/reads/barcode{i}/all.fastq.gz"
    wildcard_constraints:
        demux_algo = "dual_minibar"
    log: "logs/demux_barcode{i}_{demux_algo}.log"
    threads: 1
    conda: "conda/minibar.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        cat {input.reads} >{params.temp_fastq}
        python {input.minibar}\\
               {input.tags}\\
               {params.temp_fastq}\\
               -e 2\\
               -E 5\\
               -l 120\\
               -S\\
               -fh\\
               2>{log} |
        sed -r 's/^(@.+) [Hhx][+-]?\\([-, 0-9]+\\), ?[Hhx][+-]?\\([-, 0-9]+\\) (.+)$/\\1 sample:\\2;/' |
        gzip -c - >{output.demux}
        cp {log} {output.summary}
        rm -f {params.temp_fastq}
        """

rule demultiplex_single:
    output:
        summary = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.summary",
        demux = "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.fastq.gz"
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
        cutadapt -g file:{input.tags}\\
                 --revcomp\\
                 -o {output.demux}\\
                 --untrimmed-output /dev/null \\
                 --rename='{{header}} sample:{{adapter_name}};'\\
                 -j {threads}\\
                 {input.reads} |
                 tee {log} |
        grep -B2 "[Tt]rimmed: [^0]" >{output.summary}
        """

rule rDNA_consensus:
    output:
        cluster_map = "data/consensus/barcode{i}/{demux_algo}/{sample}/final_clusters.tsv",
        consensus = "data/consensus/barcode{i}/{demux_algo}/{sample}/consensus.fasta"
    input: "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/{demux_algo}/{sample}",
        filtered = "data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq",
        sample_label = "sample:{sample};"
    wildcard_constraints:
        demux_algo = "dual_(cutadapt|minibar)"
    log: "logs/consensus_barcode{i}_{demux_algo}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        vsearch --fastx_getseqs {input}\\
                --fastqout {params.filtered}\\
                --label '{params.sample_label}'\\
                --label_substr_match\\
                --notrunclabels
        NGSpeciesID\\
          --ont\\
          --consensus\\
          --racon\\
          --racon_iter 3\\
          --fastq {params.filtered}\\
          --outfolder {params.outdir}\\
          --m 2800\\
          --s 500\\
          --sample_size 1000\\
          --t 1\\
          &>{log} || echo "Error in NGSpeciesID" >>{log}
        touch "{output.consensus}"
        touch "{output.cluster_map}"
        for f in $(find {params.outdir} -path '*racon_cl_id_*/consensus.fasta') ; do
            echo "$f" | sed -r 's|.*/racon_cl_id(_[0-9]+)/consensus.fasta|>{wildcards.sample}\\1|' >>"{output.consensus}"
            tail -n+2 "$f" >>{output.consensus}
        done
        rm {params.filtered}
        """

rule ITS_consensus:
    output:
        cluster_map = "data/consensus/barcode{i}/{demux_algo}/{sample}/final_clusters.tsv",
        consensus = "data/consensus/barcode{i}/{demux_algo}/{sample}/consensus.fasta"
    input: "data/demultiplex/barcode{i}/{demux_algo}/barcode{i}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/{demux_algo}/{sample}",
        filtered = "data/demultiplex/barcode{i}/{demux_algo}/{sample}.fastq",
        sample_label = "sample:{sample};"
    wildcard_constraints:
        demux_algo = "single_(cutadapt|minibar)"
    log: "logs/consensus_barcode{i}_{demux_algo}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        vsearch --fastx_getseq {input}\\
                --fastqout {params.filtered}\\
                --label '{params.sample_label}'\\
                --label_substr_match\\
                --notrunclabels
        NGSpeciesID\\
          --ont\\
          --consensus\\
          --racon\\
          --racon_iter 3\\
          --fastq {params.filtered}\\
          --outfolder {params.outdir}\\
          --m 800\\
          --s 200\\
          --sample_size 1000\\
          --t 1
          &>{log} || echo "Error in NGSpeciesID" >>{log}
        touch "{output.consensus}"
        for f in $(find data/consensus/barcode01/{wildcards.demux_algo}/{wildcards.sample} -path '*racon_cl_id_*/consensus.fasta') ; do
            echo "$f" | sed -r 's|.*/racon_cl_id(_[0-9]+)/consensus.fasta|>{wildcards.sample}\\1|' >>"{output.consensus}"
            tail -n+2 "$f" >>{output.consensus}
        done
        rm {params.filtered}
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
