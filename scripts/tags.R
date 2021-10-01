library(magrittr)

barcodes <- Biostrings::readDNAStringSet("tags/its1_lr5_barcodes.fasta")
its1 <- barcodes[startsWith(names(barcodes), "its1")]
lr5 <- barcodes[startsWith(names(barcodes), "lr5")] |>
  Biostrings::reverseComplement()
threeNDf <- Biostrings::readDNAStringSet("tags/3NDf_barcodes.fasta")

Biostrings::writeXStringSet(its1, "tags/ITS1_tags.fasta")

rDNA_tags <-
    tidyr::expand_grid(
  fwd = tibble::enframe(as.character(threeNDf)),
  rev = tibble::enframe(as.character(lr5))
) |>
  tidyr::unpack(c("fwd", "rev"), names_sep = "_") |>
  glue::glue_data(">{fwd_name}_{rev_name}\n{fwd_value}...{rev_value}")

platekey <- dplyr::left_join(
    readxl::read_xlsx("tags/3NDf-LR5_tagplate.xlsx", "3NDf",
                      col_names = as.character(1:12)) |>
        dplyr::mutate(row = LETTERS[1:8]) |>
        tidyr::pivot_longer(cols = 1:12, names_to = "col", values_to = "3NDf"),
    readxl::read_xlsx("tags/3NDf-LR5_tagplate.xlsx", "LR5",
                                      col_names = as.character(1:12)) |>
        dplyr::mutate(row = LETTERS[1:8]) |>
        tidyr::pivot_longer(cols = 1:12, names_to = "col", values_to = "LR5"),
    c("row", "col")
) |>
    dplyr::mutate(tagname = glue::glue("3NDf_bc{`3NDf`}_lr5_{LR5}"))

samplekeys <- list()
for (i in 1:12) {
    f <- sprintf("data/samples/barcode%02d.xlsx", i)
    if (file.exists(f)) {
        samplekeys[[i]] <-
            readxl::read_xlsx(f, range = "B2:M9", col_names = as.character(1:12),
                          col_types = "text") |>
            dplyr::mutate(row = LETTERS[1:8]) |>
            tidyr::pivot_longer(
                cols = 1:12,
                names_to = "col",
                values_to = "sample"
            ) |>
            dplyr::left_join(platekey, by = c("row", "col"))
        samplekeys[[i]] %$%
            stringi::stri_replace_all_fixed(
                rDNA_tags,
                paste0(">", tagname),
                paste0(">", dplyr::coalesce(sample, tagname)),
                vectorize_all = FALSE
            ) |>
            writeLines(sprintf("tags/barcode%02d.fasta", i))
    }
}
