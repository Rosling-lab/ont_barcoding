library(magrittr)

barcodes <- Biostrings::readDNAStringSet("tags/its1_lr5_barcodes.fasta")
its1 <- barcodes[startsWith(names(barcodes), "its1")]
lr5 <- barcodes[startsWith(names(barcodes), "lr5")] |>
  Biostrings::reverseComplement()
threeNDf <- Biostrings::readDNAStringSet("tags/3NDf_barcodes.fasta")

Biostrings::writeXStringSet(its1, "tags/ITS1_tags.fasta")

tidyr::expand_grid(
  fwd = tibble::enframe(as.character(threeNDf)),
  rev = tibble::enframe(as.character(lr5))
) |>
  tidyr::unpack(c("fwd", "rev"), names_sep = "_") |>
  glue::glue_data(">{fwd_name}_{rev_name}\n{fwd_value}...{rev_value}") |>
  writeLines("tags/3NDf_LR5_tags.fasta")
