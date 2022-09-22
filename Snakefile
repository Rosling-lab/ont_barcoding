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

# input function to get all reads associated with a barcode
def get_reads(wildcards):
    return glob(f"data/reads/**/barcode{wildcards.i}/*_pass_barcode{wildcards.i}*.fastq.gz", recursive = True)

# get the tag files associated with a native barcode and locus (direct call)
def get_tags_(i, locus):
    checkpoints.loci_and_primers.get(i = i)
    with open(f"data/samples/barcode{i}/{locus}.primerlist") as f:
        primers = [s.rstrip('\n') for s in f.readlines()]
        return expand("tags/{primer}.fasta", primer = primers)

# get the tag files associated with a native barcode and locus (input function)
def get_tags(wildcards):
    return get_tags_(wildcards.i, wildcards.locus)

# get the loci associated with a native barcode (direct call)
def get_loci_(i):
    with open(checkpoints.loci_and_primers.get(i = i).output.loci) as f:
        loci = [s.rstrip('\n') for s in f.readlines()]
        return loci

# get the loci associated with a native barcode (input function)
def get_loci(wildcards):
    return get_loci_(i = wildcards.i)

# get the samples associated with a native barcode and locus (direct call)
def get_samples_(i, locus):
    checkpoints.loci_and_primers.get(i = i)
    with open(f"data/samples/barcode{i}/{locus}.samplelist") as f:
        return [s.rstrip('\n') for s in f.readlines()]

# get the samples associated with a native barcode and locus (input function)
def get_samples(wildcards):
    return get_samples_(wildcards.i, wildcards.locus)

# get the consensus files associated with a native barcode, locus, and demultiplexing algorithm (input function)
def get_consensus(wildcards):
    samples = get_samples(wildcards)
    consensus_files = expand("data/consensus/barcode{i}/{locus}/{demux_algo}/{sample}/consensus.fasta",
        i = wildcards.i,
        locus = wildcards.locus,
        demux_algo = wildcards.demux_algo,
        sample = samples
    )
    return(consensus_files)

# get all the outputs (input function)
def everything(wildcards):
    out = []
    for i in [x[1][0] for x in snakemake.utils.listfiles("samples/barcode{i}.xlsx")]:
        for locus in get_loci_(i):
            tags = get_tags_(i, locus)
            if len(tags) == 2:
                algos = ['cutadapt', 'minibar']
            elif len(tags) == 1:
                algos = ['cutadapt']
            else:
                raise ValueError("Too many tags: " + tags)
            for x in expand("output/barcode{i}_{locus}_{algo}.fasta", i = i, locus = locus, algo = algos):
                out.append(x)
    return out

rule all:
    input: everything

wildcard_constraints:
    i = "\d+"

# generate barcoding tag files for labeled samples
checkpoint loci_and_primers:
    output:
        loci = "data/samples/barcode{i}/locuslist"
    input:
        primers = "tags/primers.fasta",
        config = "samples/barcode{i}.xlsx",
        script = "scripts/loci_and_primers.R"
    envmodules: "R_packages/4.1.1"
    conda: "conda/tags.yaml"
    threads: 1
    script: "scripts/loci_and_primers.R"

rule tags:
    output:
        primers_fasta = "data/samples/barcode{i}/{locus}_primers.fasta",
        trim_fasta = "data/samples/barcode{i}/{locus}_trim.fasta",
        tags_fasta = "data/samples/barcode{i}/{locus}_tags.fasta",
        table = "data/samples/barcode{i}/{locus}.minibar",
        lengths = "data/samples/barcode{i}/{locus}.lengths"
    input:
        tags = get_tags,
        primers = "tags/primers.fasta",
        config = "samples/barcode{i}.xlsx",
        script = "scripts/tags.R"
    envmodules: "R_packages/4.1.1"
    conda: "conda/tags.yaml"
    threads: 1
    script: "scripts/tags.R"

def get_min_outer_length(wildcards):
    with open(f"data/samples/barcode{wildcards.i}/{wildcards.locus}.lengths") as f:
        return [s.rstrip('\n') for s in f.readlines()][0]

def get_max_outer_length(wildcards):
    with open(f"data/samples/barcode{wildcards.i}/{wildcards.locus}.lengths") as f:
        return [s.rstrip('\n') for s in f.readlines()][1]

def get_inner_length_center(wildcards):
    with open(f"data/samples/barcode{wildcards.i}/{wildcards.locus}.lengths") as f:
        return [s.rstrip('\n') for s in f.readlines()][2]

def get_inner_length_width(wildcards):
    with open(f"data/samples/barcode{wildcards.i}/{wildcards.locus}.lengths") as f:
        return [s.rstrip('\n') for s in f.readlines()][3]

rule demultiplex_cutadapt:
    output:
        summary = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/barcode{i}.summary",
        demux = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/barcode{i}.fastq.gz",
        unknowns = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/unknown.fastq.gz"
    input:
        reads = get_reads,
        primers = rules.tags.output.primers_fasta,
        tags = rules.tags.output.tags_fasta
    params:
        outdir = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/",
        min_length = get_min_outer_length,
        max_length = get_max_outer_length,
        temp1 = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/temp1.fastq.gz",
        temp2 = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/temp2.fastq.gz"
    wildcard_constraints:
        demux_algo = "cutadapt"
    log:
        filter = "logs/demux_barcode{i}_{locus}_{demux_algo}.log",
        demux = "logs/demux_barcode{i}_{locus}_{demux_algo}.log",
        orient = "logs/orient_barcode{i}_{locus}_{demux_algo}.log"
    threads: maxthreads
    conda: "conda/demultiplex.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        zcat {input.reads} |
        cutadapt --minimum-length {params.min_length}\\
                 --maximum-length {params.max_length}\\
                 -o -\\
                 -j {threads}\\
                 - 2>{log.filter} |
        cutadapt -g file:{input.primers}\\
                 --revcomp\\
                 -o -\\
                 --untrimmed-output {output.unknowns}\\
                 --action=retain\\
                 -j {threads}\\
                 - 2>{log.orient} |
        cutadapt -g file:{input.tags}\\
                 -o {output.demux}\\
                 -j {threads}\\
                 --untrimmed-output /dev/null\\
                 --action=retain\\
                 --rename='{{header}} sample:{{adapter_name}};'\\
                 - |
                 tee {log.demux} |
        grep -B2 "[Tt]rimmed: [^0]" >{output.summary}
        rm -f {params.temp2}
        """

rule demultiplex_minibar:
    output:
        summary = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/barcode{i}.summary",
        demux = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/barcode{i}.fastq.gz"
    input:
        reads = get_reads,
        tags = rules.tags.output.table,
        minibar = "scripts/minibar/minibar.py"
    params:
        outdir = "data/demultiplex/barcode{i}/{locus}/{demux_algo}",
        temp_fastq = "data/reads/barcode{i}/all.fastq.gz"
    wildcard_constraints:
        demux_algo = "minibar"
    log: "logs/demux_barcode{i}_{locus}_{demux_algo}.log"
    threads: 1
    conda: "conda/minibar.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        for f in {input.reads} ; do
            python {input.minibar}\\
                   {input.tags}\\
                   ${{f}}\\
                   -e 2\\
                   -E 5\\
                   -l 120\\
                   -T\\
                   -S\\
                   -fh
        done 2>{log} |
        sed -r 's/^(@.+) [Hhx][+-]?\\([-, 0-9]+\\), ?[Hhx][+-]?\\([-, 0-9]+\\) (.+)$/\\1 sample:\\2;/' |
        gzip -c - >{output.demux}
        cp {log} {output.summary}
        """

rule consensus:
    output:
        cluster_map = "data/consensus/barcode{i}/{locus}/{demux_algo}/{sample}/final_clusters.tsv",
        consensus = "data/consensus/barcode{i}/{locus}/{demux_algo}/{sample}/consensus.fasta"
    input: "data/demultiplex/barcode{i}/{locus}/{demux_algo}/barcode{i}.fastq.gz"
    params:
        outdir = "data/consensus/barcode{i}/{locus}/{demux_algo}/{sample}",
        filtered = "data/demultiplex/barcode{i}/{locus}/{demux_algo}/{sample}.fastq",
        sample_label = "sample:{sample};",
        center_length = get_inner_length_center,
        length_window = get_inner_length_width
    log: "logs/consensus_barcode{i}_{locus}_{demux_algo}/{sample}.log"
    conda: "conda/NGSpeciesID.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        vsearch --fastx_getseqs {input}\\
                --fastqout {params.filtered}\\
                --label '{params.sample_label}'\\
                --label_substr_match\\
                --notrunclabels &>{log}
        NGSpeciesID\\
          --ont\\
          --consensus\\
          --racon\\
          --racon_iter 3\\
          --fastq {params.filtered}\\
          --outfolder {params.outdir}\\
          --m {params.center_length}\\
          --s {params.length_window}\\
          --sample_size 1000\\
          --t 1\\
          &>{log} || echo "Error in NGSpeciesID" >>{log}
        touch "{output.consensus}"
        touch "{output.cluster_map}"
        for f in $(find {params.outdir} -path '*racon_cl_id_*/consensus.fasta') ; do
            sed -r 's/>consensus_cl_id(_[0-9]+)_total_supporting_reads_([0-9]+) .+/>{wildcards.sample}\\1_(\\2x_coverage)/' <"$f" >>"{output.consensus}"
        done
        rm {params.filtered}
        """

rule all_consensus:
    output: "output/barcode{i}_{locus}_{demux_algo}.fasta"
    wildcard_constraints:
        demux_algo = "(cutadapt|minibar)"
    input:
        consensus = get_consensus,
        primers = rules.tags.output.trim_fasta
    conda: "conda/demultiplex.yaml"
    threads: 1
    log: "logs/consensus_barcode{i}_{locus}_{demux_algo}.log"
    shell:
        """
        cat {input.consensus} |
        cutadapt -g file:{input.primers}\\
                 --revcomp\\
                 -o -\\
                 -O 10\\
                 --untrimmed-output /dev/null\\
                 -j {threads}\\
                 - 2>{log} |
        sed 's/ rc$//' >{output}
        """

rule itsx:
    input: "output/barcode{i}_rDNA_{demux_algo}.fasta"
    wildcard_constraints:
        demux_algo = "(cutadapt|minibar)"
    output:
        ITS = "output/barcode{i}_rDNA_{demux_algo}.ITS.fasta",
        LSU = "output/barcode{i}_rDNA_{demux_algo}.LSU.fasta",
        SSU = "output/barcode{i}_rDNA_{demux_algo}.SSU.fasta"
    params:
        prefix = "output/barcode{i}_rDNA_{demux_algo}",
        fullfile = "output/barcode{i}_rDNA_{demux_algo}.full.fasta"
    log: "logs/itsx{i}_{demux_algo}.log"
    conda: "conda/ITSx.yaml"
    threads: 4
    shell: """
    ITSx\\
        -i {input}\\
        -o {params.prefix}\\
        -t "."\\
        --cpu {threads}\\
        --complement F\\
        --save_regions LSU,SSU\\
        --graphical F\\
        --preserve T\\
        --positions F\\
        --summary F\\
        --not_found F\\
        >{log}
    mv {params.fullfile} {output.ITS}
    """

rule download_reunite_rdp_lsu:
    output: "references/rdp_train.LSU.sintax.fasta.gz"
    shell: """
    mkdir -p references
    wget -O {output} https://github.com/brendanf/reUnite/releases/download/v0.1/rdp_train.LSU.sintax.fasta.gz
    """

rule download_reunite_unite_its:
    output: "references/unite.ITS.sintax.fasta.gz"
    shell: """
    mkdir -p references
    wget -O {output} https://github.com/brendanf/reUnite/releases/download/v0.1/unite.ITS.sintax.fasta.gz
    """

def get_reference(wildcards):
    if wildcards.sublocus == "ITS":
        return "references/unite.ITS.sintax.fasta.gz"
    elif wildcards.sublocus == "LSU":
        return "references/rdp_train.LSU.sintax.fasta.gz"


rule sintax:
    output: touch("output/barcode{i}_{locus}_{demux_algo}_{sublocus}_sintax.tsv")
    wildcard_constraints:
        demux_algo = "(cutadapt|minibar)"
    input:
        consensus = "output/barcode{i}_{locus}_{demux_algo}.{sublocus}.fasta",
        reference = get_reference
    threads: maxthreads
    conda: "conda/vsearch.yaml"
    log: "logs/sintax_barcode{i}_{locus}_{demux_algo}_{sublocus}.log"
    shell: """
    vsearch\\
      --sintax {input.consensus}\\
      --db {input.reference}\\
      --sintax_cutoff 0.9\\
      --tabbedout {output}\\
      --threads {threads}\\
      &>{log}
    """

