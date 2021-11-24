## Bioinformatics pipeline for Minion barcoding

This Snakemake pipeline processed Minion reads from multiplexed barcoding experiments into high-quality consensus sequences. It has support for both single- and dual-indexing for multiplexing barcodes, as well as for libraries including multiple primer pairs (i.e., different genes).

## Overview

Sequences are first demultiplexed using [minibar](https://github.com/calacademy-research/minibar) for dual-indexed samples, and [cutadapt](https://cutadapt.readthedocs.io) for single-indexed samples.

Then, consensus sequences for each sample are generated using [NGspeciesID](https://github.com/ksahlin/NGSpeciesID).

## Configuration

Required configuration files are:

* `tags/primers.fasta` should include the sequences of all of the primer sequences in use. These should not include any sample-specific tagged/barcoded/indexed versions, only the basic primer sequences. These should be in fasta format. Primer names should include only the characters "`A-Za-z0-9_-`", and should not end in "`-tag`"- *Example:*
```
>ITS1
TCCGTAGGTGAACCTGC
>ITS4
TCCTCCGCTTATTGATATGC
```

* For each primer which is tagged/barcoded/indexed there should also be a file `tags/{primer}-tag.fasta` which includes all the tagged variants. The named should be in the format `{primer}_{ID}`. *Example:* (file `tags/ITS1-tag.fasta`)
```
>ITS1_1
GGTAGTCAGACGATGCGTCATTCCGTAGGTGAACCTGC
>ITS1_2
GGTAGCTATACATGACTCTGCTCCGTAGGTGAACCTGC
```

* For each sequencing library, a file `samples/barcode{n}.xlsx`.  This is an excel spreadsheet which defines which genes have been sequenced, which primers were used, and how the multiplexing tags/barcodes/indexes correspond to individual samples.  Different data are given in different sheets, many of which can be copied over between different sequencing runs using the same basic setup.

  * **Instructions** Gives instructions for using the spreadsheet.
  * **Loci** Defines the different loci (often genes) that have been sequenced.  For each locus, fill in a **Locus name**, the **Forward primer** and **Reverse primer**, as well as the expected **Min length** and **Max length** for amplicons, not including the length of the primers. Some experiments may have only one locus, but it should be defined here anyway. Primer names should match exactly the names given in `tags/primers.fasta` for non-tagged primers, or the name of a `tags/{primer}-tag.fasta` file (i.e., "`{primer}-tag`") for tagged primers. A locus is dual-indexed if both of its primers are tagged, and single-indexed if only one of its primers is tagged.  *Example:*
  
  | Locus name | Forward primer | Reverse primer | Min length | Max length |
  |------------|----------------|----------------|------------|------------|
  | ITS        | ITS1-tag       | ITS4           | 750        | 1000       |
  
  * It is assumed that both tagged primers and samples are stored/processed in plates, so that the primer (pair) which identifies a particular sample can be identified by finding the primer(s) which correspond to the same plate location as the sample.
  * A sheet named `{primer}-tag` should define the plate layout used for each indexed primer.
  The first row gives plate column number, and the first column gives plate row letters.  Each cell should give the ID of the primer tag/barcode/index for that well of the plate.  IDs should exactly match the names of the sequences in `tags/{primer}-tag.fasta`. If you do not use a whole plate worth of primers, leave some of the sheet blank. *Example:* (sheet `ITS1-tag`)
  
  |   | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
  |---|---|---|---|---|---|---|---|---|---|----|----|----|
  | A | 1 | 2 |
  | B |
  | C |
  | D |
  | E |
  | F |
  | G |
  | H |
  
  * A sheet named `{locus}` for each locus.  The layout is the same as in the tagged primer sheets, but now each cell gives a sample name.  Again, some of the cells can be blank, but every cell with a sample should have corresponding primers. *Example:* (sheet `ITS`)
  
  |   | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
  |---|---|---|---|---|---|---|---|---|---|----|----|----|
  | A | MySample | YourSample |
  | B |
  | C |
  | D |
  | E |
  | F |
  | G |
  | H |

## Input Files

Basecalled reads (`fastq.gz` files) should be placed in `data/reads/fastq_passed/barcode{n}`.  The pipeline can do several barcodes at once, as long as they all have their own reads and `samples/barcode{n}.xlsx` configuration file.

## Running

The pipeline is designed to be run on UPPMAX. To run the script on a single node,

```
sbatch run_node.sh
```

To run it on a local machine, you will need to have [`conda`](https://www.anaconda.org) and [`snakemake`](https://snakemake.readthedocs.io) installed.  `snakemake` will install the other dependencies via `conda`.  This will only work on OSX or Linux, sorry no Windows.

```
snakemake --use-conda -j{n}
```

where `{n}` is the desired number of cpus to use.

Depending on the version of `snakemake` and whether you have `mamba` installed, you may need to include `--conda-frontend conda`.

## Output

Consensus files are generated as `output/{locus}_barcode{n}.fasta`.  Sequence names are `{sample}_{m}`, where `m` is the cluster number; in some cases `NGSpeciesID determines that there is more than one sequence in the reads, and generates a consensus for each of them separately.  These may be due to contamination, sequencing a composite sample, or very high error rates.
