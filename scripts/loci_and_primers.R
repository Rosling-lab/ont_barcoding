## Open xlsx configuration file, do a few consistency checks, and write simple
## lists for Snakemake
## Brendan Furneaux
## November 2021

library(magrittr)
if (exists("snakemake")) {
  configfile <- snakemake@input$config
  primerfile <- snakemake@input$primers
  locuslist <- snakemake@output$loci
} else {
  configfile <- "samples/mycoweek2022/barcode12.xlsx"
  primerfile <- "tags/primers.fasta"
  basename <- sub("\\.xlsx$", "", basename(configfile))
  expname <- basename(dirname(configfile))
  locuslist <- file.path("data", expname, "samples", basename, "locuslist")
}

primers <- Biostrings::readDNAStringSet(primerfile)

# check all the primer names are alphanumeric plus - and _,
# and do not end with "-tag"
primer_name_is_valid <-
  grepl("^[A-Za-z0-9_-]+$", names(primers)) & !endsWith(names(primers), "-tag")
if (!all(primer_name_is_valid)) {
  stop(
    "Invalid primer names: ",
    paste(names(primers)[!primer_name_is_valid], collapse = ", ")
  )
}

loci <- readxl::read_xlsx(configfile, sheet = "Loci")

# check we have sequences for all the primers
all_primers <- union(loci$`Forward primer`, loci$`Reverse primer`)
primers_untagged <- sub("-tag$", "", all_primers) %>% unique()
if (!all(primers_untagged %in% names(primers))) {
  stop(
    "Primer missing from ",
    primerfile,
    ": ",
    paste(setdiff(primers_untagged, names(primers)), collapse = ", ")
  )
}

# check we have sequences for all the tagged primers
primers_tagged <- all_primers[endsWith(all_primers, "-tag")]
missing_tag_file <- !file.exists(sprintf("tags/%s.fasta", primers_tagged))
if (any(missing_tag_file)) {
  stop("Primer tag files not found for primers: ", paste(primers_tagged[missing_tag_file], collapse = ", "))
}

# check we have a sheet defining samples for each locus
sheets <- readxl::excel_sheets(configfile)
missing_sheets <- !(loci$`Locus name` %in% sheets)
if (any(missing_sheets)) {
  stop(
    "Configuration file ",
    configfile,
    " is missing sheets for loci: ",
    paste(loci$`Locus name`[missing_sheets], collapse = ", ")
  )
}

missing_sheets <- !(primers_tagged %in% sheets)
if (any(missing_sheets)) {
  stop(
    "Configuration file ",
    configfile,
    " is missing sheets for tagged primers: ",
    paste(primers_tagged[missing_sheets], collapse = ", ")
  )
}

# write the names of the loci
outdir <- dirname(locuslist)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
writeLines(loci$`Locus name`, locuslist)

metadata_sheet <- grep("metadata", sheets, ignore.case = TRUE, value = TRUE)
if (length(metadata_sheet) > 1) {
  stop("Multiple sheets seem to contain metadata: ", paste(metadata_sheet, collapse = ", "))
} else if (length(metadata_sheet) == 1) {
  metadata <- readxl::read_xlsx(
    configfile,
    sheet = metadata_sheet
  ) %>%
    dplyr::rename(sample = 1) %>%
    tidyr::unite("metadata", -1) %>%
    tidyr::unite("long_sample", everything(), remove = FALSE) %>%
    dplyr::select(sample, long_sample)
} else {
  metadata <- tibble::tibble(sample = character(), long_sample = character())
}

# write the names of the primers and samples for each locus
# snakemake will use these to get dependencies for future steps
for (i in seq_len(nrow(loci))) {
  mytags <- intersect(loci[i, c("Forward primer", "Reverse primer")], primers_tagged)
  writeLines(
    mytags,
    file.path(outdir, sprintf("%s.primerlist", loci$`Locus name`[i]))
  )
  samplelist_file <- file.path(outdir, sprintf("%s.samplelist", loci$`Locus name`[i]))
  readxl::read_xlsx(
    configfile,
    sheet = loci$`Locus name`[i],
    range = "B2:M9",
    col_names = as.character(1:12),
    col_types = "text"
  ) %>%
    unlist() %>%
    na.omit() %>%
    tibble::tibble(sample = .) %>%
    dplyr::left_join(metadata, by = "sample") %>%
    dplyr::mutate(long_sample = dplyr::coalesce(long_sample, sample)) %>%
    readr::write_tsv(
      samplelist_file,
      col_names = FALSE,
      na = ""
    )
}
