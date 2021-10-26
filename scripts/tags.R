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
  tags_ITS1 <- snakemake@output$tags_ITS1
  tags_3NDf_LR5 <- snakemake@output$tags_3NDf_LR5
  sample_tags <- snakemake@output$sample_tags
  sample_names <- snakemake@output$sample_names
} else {
  # defaults without Snakemake
  tags_3NDf <- "tags/3NDf_barcodes.fasta"
  tags_ITS1_LR5 <- "tags/its1_lr5_barcodes.fasta"
  tag_plate <- "tags/3NDf-LR5_tagplate.xlsx"
  sample_plate <- list.files("samples", "barcode[0-9]+\\.xlsx", full.names = TRUE)

  tags_ITS1 <- "tags/ITS1_tags.fasta"
  tags_3NDf_LR5 <- "tags/3NDf_LR5_tags.fasta"
  sample_tags <- file.path(
    "tags",
    sub("xlsx$", "fasta", basename(sample_plate))
  )
  sample_names <- file.path(
    "samples",
    sub("xlsx$", "txt", basename(sample_plate))
  )
}

barcodes <- Biostrings::readDNAStringSet(tags_ITS1_LR5)
its1 <- barcodes[startsWith(names(barcodes), "its1")]
lr5 <- barcodes[startsWith(names(barcodes), "lr5")] |>
  Biostrings::reverseComplement()
threeNDf <- Biostrings::readDNAStringSet(tags_3NDf)

if (!is.null(tags_ITS1)) {
  Biostrings::writeXStringSet(its1, tags_ITS1)
}

rDNA_tags <-
    tidyr::expand_grid(
  fwd = tibble::enframe(as.character(threeNDf)),
  rev = tibble::enframe(as.character(lr5))
) |>
  tidyr::unpack(c("fwd", "rev"), names_sep = "_") |>
  glue::glue_data(">{fwd_name}_{rev_name}\n{fwd_value}...{rev_value}")

if (!is.null(tags_3NDf_LR5)) {
  writeLines(rDNA_tags, tags_3NDf_LR5)
}

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
        rDNA_tags,
        paste0(">", tagname),
        paste0(">", dplyr::coalesce(sample, tagname)),
        vectorize_all = FALSE
      ) |>
      writeLines(sample_tags[i])
    sample_data$sample |>
      purrr::discard(is.na) |>
      writeLines(sample_names[i])
  }
}
