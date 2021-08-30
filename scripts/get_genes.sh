#!/usr/bin/env bash
# usage: get_genes.sh {basename} {regex} [options]
# inputs:
#    {basename}_individual.fasta - a fasta file of matching sequences fron GenBank, not including genomic contigs
#    {basename}_genomes.fasta - a gene features fasta file from GenBank for genomic contigs
#    {basename}_genomes.summary - text summary file of the same contig sequences from GenBank
#    {basename}_primers.fasta - file giving primer sequences to trim, in cutadapt format
#
# Extracts all gene features with a name matching {regex} (as for sed -r),
# and also full sequences with {regex} in the title, to {basename}_selected.fasta
# Then runs cutadapt to trim at the primer locations and writes the output to {basename}_trimmed.fasta
set -x

BASE=$1
shift
SUMMARY=${BASE}_genomes.summary
REPLACE=${BASE}.sed
# parse species names out of the contigs summary
paste - - - - <"$SUMMARY" |
sed -r 's/\s*[0-9]+\. (UNVERIFIED: |TSA: |[A-Za-z0-9-]+[0-9] |.*library |.*libraries |.*library during development )?(\[?[A-Z]\w+\]? (cf. |aff. )?([-a-z]+ |sp\.( [0-9]+ | ([A-Z]+|[A-Z][a-z]+) [A-Z]*[0-9]+a? | [A-Z0-9'"'"'][A-Za-z0-9_.'"'"'-]+ | )(\b|\())|Uncultured [A-Z][a-z]+).+linear (DNA|mRNA) *[\t ]+([A-Z]+_?[0-9]+).+/s@(lcl\\|)?(\9)\.[1-9]@\2\\2@/' >"$REPLACE"
GENOMES=${BASE}_genomes.fasta
INDIVIDUAL=${BASE}_individual.fasta
SELECTED=${BASE}_selected.fasta
CLUSTER=${BASE}_cluster.fasta
CLUSTERALIGN=${BASE}_cluster_aligned.fasta

TEMP=$(tempfile -s .fasta)
trap "rm -f $TEMP" EXIT

# select locus names which are one of the chosen genes
PATTERN=$1
shift
PRIMERS=${BASE}_primers.fasta
PAIRS=${BASE}_pairs.fasta
PRIMERALIGN=${BASE}_aligned_primers.fasta
vsearch --fastx_getseqs "$GENOMES" --notrunclabels --label "gene=" --label_substr_match --fastaout "$TEMP" 2>"${BASE}_select_gene.log"
sed -nr 's/>(.+) .*\[gene='"$PATTERN"'\].+/\1/ip' <"$GENOMES" |
tee "${BASE}_selected.txt" |
vsearch --fastx_getseqs "$TEMP" --fastaout $SELECTED --labels - 2>>"${BASE}_select_gene.log"

# modify titles to include species names
sed -i -r -f "$REPLACE" "$SELECTED"

# select sequences which have not already been chosen, and have the chosen gene name in the title
sed -nr 's/>([^ ]+) .*'"$PATTERN".*'/\1/ip' <"$INDIVIDUAL" |
vsearch --fastx_getseqs "$INDIVIDUAL" --notrunclabels --label_substr_match --fastaout - --labels - 2>"${BASE}_select_ind.log" |
# put taxonomy first
sed -r 's/>([A-Z]+_?[0-9]+\.?[0-9]) (UNVERIFIED: |TSA: |[A-Za-z0-9-]+[0-9] |.*library |.*libraries |.*library during development )?(\[?[A-Z]\w+\]? (cf. |aff. )?([-a-z]+ |sp\.( [0-9]+ | ([A-Z]+|[A-Z][a-z]+) [A-Z]*[0-9]+a? | [A-Z0-9'"'"'][A-Za-z0-9_.'"'"'-]+ | )(\b|\())|Uncultured [A-Z][a-z]+).+/>\3 \1/' >>"$SELECTED"

# cluster so that we don't have to look at lots of almost identical sequences
vsearch --cluster_fast "$SELECTED" --id 0.95 --notrunclabels --consout - 2>"${BASE}_cluster.log" |
sed 's/centroid=//' >"$CLUSTER"

# align
rm -f "${BASE}_align.log"
mafft --auto --reorder --thread 8 "$CLUSTER" >"$CLUSTERALIGN" 2>>"${BASE}_align.log"

mafft --auto --addfragments "$PRIMERS" --thread 8 "$CLUSTERALIGN" >"$PRIMERALIGN" 2>>"${BASE}_align.log"

# trim with the given primer pair(s)
rm -f "${BASE}_trim.log"
while IFS= read -r HEADER; do
  IFS= read -r SEQ
  PAIRBASE="${BASE}_trim_${HEADER#>}"
  cutadapt --revcomp -a "$SEQ" -O 15 -e 1 "$@" --action=retain --discard-untrimmed -o "${PAIRBASE}.fasta" "$CLUSTER" >>"${PAIRBASE}_trim.log"
  mafft --auto --reorder --thread 8 "${PAIRBASE}.fasta" >"${PAIRBASE}_aligned.fasta" 2>>"${PAIRBASE}_align.log"
done <"$PAIRS"

# just look for primers(s)
rm -f "${BASE}_primer.log"
while IFS= read -r HEADER; do
  IFS= read -r SEQ
  PRIMERBASE="${BASE}_primer_${HEADER#>}"
  echo ${HEADER} >"${PRIMERBASE}.fasta"
  echo ---${SEQ}--- >>"${PRIMERBASE}.fasta"
  cutadapt --revcomp -a "N{3}${SEQ}N{3}" -O 15 -e 0.2 "$@" --action=retain --discard-untrimmed -o - "$CLUSTER" |
  cutadapt --revcomp -g "N{3}${SEQ}N{3}" -O 15 -e 0.2 "$@" --action=retain --discard-untrimmed -o - >>"${PRIMERBASE}.fasta" -
done <"$PRIMERS"

LONG=${BASE}_long.fasta
LONGALIGN=${BASE}_long_align.fasta
LONGPRIMERALIGN=${BASE}_long_aligned_primers.fasta
# select just the long sequences
vsearch --fastx_filter "$CLUSTER" --fastq_minlen 1000 --fastaout "$LONG" 2>"${BASE}_long.log"
mafft --auto --reorder --thread 8 "$LONG" >"$LONGALIGN" 2>>"${BASE}_align.log"
mafft --auto --addfragments "$PRIMERS" --thread 8 "$LONGALIGN" >"$LONGPRIMERALIGN" 2>>"${BASE}_align.log"

