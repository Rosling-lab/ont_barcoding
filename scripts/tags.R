library(magrittr)

if (exists("snakemake")) {
  # input files
  # (some may be NULL)
  tags_3NDf <- snakemake@input$tags_3NDf
  tags_ITS1_LR5 <- snakemake@input$tags_ITS1_LR5
  tag_plate <- snakemake@input$tag_plate
  sample_plate <- snakemake@input$sample_plate

  # output files
  # (some may be NULL)
  tags_ITS1_fasta <- snakemake@output$tags_ITS1_fasta
  tags_ITS1_table <- snakemake@output$tags_ITS1_table
  tags_3NDf_LR5_fasta <- snakemake@output$tags_3NDf_LR5_fasta
  tags_3NDf_LR5_table <- snakemake@output$tags_3NDf_LR5_table
  sample_tags_fasta <- snakemake@output$sample_tags_fasta
  sample_tags_table <- snakemake@output$sample_tags_table
  sample_names <- snakemake@output$sample_names
} else {
  # defaults without Snakemake
  tags_3NDf <- "tags/3NDf_barcodes.fasta"
  tags_ITS1_LR5 <- "tags/its1_lr5_barcodes.fasta"
  tag_plate <- "tags/3NDf-LR5_tagplate.xlsx"
  sample_plate <- list.files("samples", "barcode[0-9]+\\.xlsx", full.names = TRUE)

  tags_ITS1_fasta <- "tags/ITS1_tags.fasta"
  tags_ITS1_table <- "tags/ITS1_tags.tsv"
  tags_3NDf_LR5_fasta <- "tags/3NDf_LR5_tags.fasta"
  tags_3NDf_LR5_table <- "tags/3NDf_LR5_tags.tsv"
  sample_tags_fasta <- file.path(
    "tags",
    sub("xlsx$", "fasta", basename(sample_plate))
  )
  sample_tags_table <- file.path(
    "tags",
    sub("xlsx$", "tsv", basename(sample_plate))
  )
  sample_names <- file.path(
    "samples",
    sub("xlsx$", "txt", basename(sample_plate))
  )
}

# read the files containing the barcodes
barcodes <- Biostrings::readDNAStringSet(tags_ITS1_LR5)
seq_ITS1 <- barcodes[startsWith(names(barcodes), "its1")]
seq_LR5 <- barcodes[startsWith(names(barcodes), "lr5")]
seq_3NDf <- Biostrings::readDNAStringSet(tags_3NDf)

# functions to remove common prefix and suffix
remove_lcsuffix <- function(seq) {
  lcs <- Biobase::lcSuffix(seq)
  sub(paste0(lcs, "$"), "", seq)
}

remove_lcprefix <- function(seq) {
  lcp <- Biobase::lcPrefix(seq)
  sub(paste0("^", lcp), "", seq)
}

# make a data frame for each tag set:
# name = name,
# seq = full sequence,
# revcomp = reverse complement of full sequence
# tag = only barcoding tag sequence
# primer = only primer sequence
enframe_tags <- function(tags) {
  revcomp <- as.character(Biostrings::reverseComplement(tags))
  tags <- as.character(tags)
  primer <- Biobase::lcSuffix(tags)
  tags |>
    tibble::enframe(value = "seq") |>
    dplyr::mutate(
      tag = seq |>
        remove_lcsuffix() |> #remove the primer
        remove_lcprefix(), #remove the "pad", if any
      primer = primer,
      revcomp = revcomp
    )
}

tagset_3NDf <- enframe_tags(seq_3NDf)
tagset_LR5 <- enframe_tags(seq_LR5)
tagset_ITS1 <- enframe_tags(seq_ITS1)

# find the minimum edit distance between tags in a tag set
mindist <- function(tags) {
  adist(tags$tag) |>
  purrr::keep(~.>0) |>
  min()
}

mindist(tagset_3NDf) # 5
mindist(tagset_LR5) # 7
mindist(tagset_ITS1) # 7

# make a set of tags in cutadapt fasta format for 3NDf--LR5

rDNA_tags_cutadapt <-
  tidyr::expand_grid(fwd = tagset_3NDf, rev = tagset_LR5) |>
  tidyr::unpack(c("fwd", "rev"), names_sep = "_") |>
  glue::glue_data(">{fwd_name}_{rev_name}\n{fwd_seq}...{rev_revcomp}")

# make a set of tags in minibar table format for 3NDf--LR5
rDNA_tags_minibar <-
  tidyr::expand_grid(fwd = tagset_3NDf, rev = tagset_LR5) |>
  tidyr::unpack(c("fwd", "rev"), names_sep = "_") |>
  dplyr::select(fwd_name, rev_name, fwd_tag, fwd_primer, rev_tag, rev_primer) |>
  tidyr::unite("name", c(fwd_name, rev_name))

# ITS(1) is already in a fine format for cutadapt,
# which can do single-end demultiplexing

# make a set of tags in minibar table format for ITS1-ITS4
# pretend that the ITS4 is a tag (with only one variant) with no primer

ITS_tags_minibar <-
  tibble::enframe(as.character(seq_ITS1)) |>
  dplyr::transmute(
    name = name,
    fwd_tag = remove_lcsuffix(value),
    fwd_primer = Biobase::lcSuffix(value),
    rev_tag = "GCATATCAATAAGCGGAGGA",
    rev_primer = ""
  )

# write out the tags; in this form they are not designed for a specific plate
# layout.

if (!is.null(tags_ITS1_fasta)) {
  Biostrings::writeXStringSet(seq_ITS1, tags_ITS1_fasta)
}
if (!is.null(tags_ITS1_table)) {
  write.table(ITS_tags_minibar, tags_ITS1_table, sep = "\t", row.names = FALSE,
              quote = FALSE)
}

if (!is.null(tags_3NDf_LR5_fasta)) {
  writeLines(rDNA_tags_cutadapt, tags_3NDf_LR5_fasta)
}

if (!is.null(tags_3NDf_LR5_table)) {
  write.table(rDNA_tags_minibar, tags_3NDf_LR5_table, sep = "\t", row.names = FALSE,
              quote = FALSE)
}

# if we have a plate, then convert tag (pairs) to wells
if (!is.null(tag_plate)) {
  platekey <- dplyr::left_join(
    readxl::read_xlsx(
      tag_plate,
      "3NDf",
      range = "B2:M9",
      col_names = as.character(1:12),
      col_types = "text"
    ) |>
      dplyr::mutate(row = LETTERS[1:8]) |>
      tidyr::pivot_longer(cols = 1:12, names_to = "col", values_to = "3NDf"),
    readxl::read_xlsx(
      tag_plate,
      "LR5",
      range = "B2:M9",
      col_names = as.character(1:12),
      col_types = "text"
    ) |>
      dplyr::mutate(row = LETTERS[1:8]) |>
      tidyr::pivot_longer(cols = 1:12, names_to = "col", values_to = "LR5"),
    c("row", "col")
  ) |>
    dplyr::mutate(tagname = glue::glue("3NDf_bc{`3NDf`}_lr5_{LR5}"))

# convert wells to samples for each plate
  for (i in seq_along(sample_plate)) {
    sample_data <- readxl::read_xlsx(
      sample_plate[i],
      range = "B2:M9",
      col_names = as.character(1:12),
      col_types = "text"
    ) |>
      dplyr::mutate(row = LETTERS[1:8]) |>
      tidyr::pivot_longer(
        cols = 1:12,
        names_to = "col",
        values_to = "sample"
      ) |>
      dplyr::left_join(platekey, by = c("row", "col"))
    sample_data %$%
      stringi::stri_replace_all_fixed(
        rDNA_tags_cutadapt,
        paste0(">", tagname),
        paste0(">", dplyr::coalesce(sample, tagname)),
        vectorize_all = FALSE
      ) |>
      writeLines(sample_tags_fasta[i])

    sample_data |>
      dplyr::left_join(rDNA_tags_minibar, by = c("tagname" = "name")) |>
      dplyr::transmute(
        name = dplyr::coalesce(sample, tagname),
        fwd_tag = fwd_tag,
        fwd_primer = fwd_primer,
        rev_tag = rev_tag,
        rev_primer = rev_primer
      ) |>
      write.table(sample_tags_table[i], sep = "\t", row.names = FALSE,
                  quote = FALSE)

    sample_data$sample |>
      purrr::discard(is.na) |>
      writeLines(sample_names[i])
  }
}
