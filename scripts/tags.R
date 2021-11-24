# generate tag files for demultiplexing for one library+locus combination

library(magrittr)
if (exists("snakemake")) {
  # input files
  tagfiles <- snakemake@input$tags
  primerfile <- snakemake@input$primers
  configfile <- snakemake@input$config

  # parameters
  barcode_id <- snakemake@wildcards$i
  locus <- snakemake@wildcards$locus

  # output files
  # (some may be NULL)
  primers_fasta <- snakemake@output$primers_fasta
  tags_fasta <- snakemake@output$tags_fasta
  tags_table <- snakemake@output$table
} else {
  # defaults without Snakemake
  # input files
  tagfiles <- c("tags/3NDf-tag.fasta", "tags/LR5-tag.fasta", "tags/ITS1-tag.fasta")
  primerfile <- "tags/primers.fasta"
  configfile <- "samples/barcode01.xlsx"

  # parameters
  basename <- sub("\\.xlsx$", "", basename(configfile))
  locus <- "rDNA"

  # output files
  primers_fasta <- file.path(
    "tags",
    basename,
    paste0(locus, "_primers.fasta")
  )
  tags_fasta <- file.path(
    "tags",
    basename,
    paste0(locus, "_tags.fasta")
  )
  tags_table <- file.path(
    "tags",
    basename,
    paste0(locus, ".minibar")
  )
}



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
enframe_tags <- function(tags, primer) {
  revcomp <- as.character(Biostrings::reverseComplement(tags))
  tags <- as.character(tags)
  stopifnot(all(endsWith(tags, primer)))
  primer <- Biobase::lcSuffix(tags)
  pad <- Biobase::lcPrefix(tags)
  tags %>%
    tibble::enframe(value = "seq") %>%
    dplyr::mutate(
      tag = seq %>%
        remove_lcsuffix() %>% #remove the primer
        remove_lcprefix(), #remove the "pad", if any
      primer = primer,
      pad = pad,
      revcomp = revcomp
    )
}

# find the minimum edit distance between tags in a tag set
mindist <- function(tags) {
  adist(tags) %>%
    purrr::keep(~.>0) %>%
    min()
}

load_platekey <- function(file, tagname, which) {
  primername <- sub("-tag$", "", tagname)
  readxl::read_xlsx(
    file,
    sheet = tagname,
    range = "B2:M9",
    col_names = as.character(1:12),
    col_types = "text"
  ) %>%
    dplyr::mutate(row = LETTERS[1:8]) %>%
    tidyr::pivot_longer(cols = 1:12, names_to = "col", values_to = "id") %>%
    dplyr::filter(!is.na(id)) %>%
    tibble::add_column(!!which := primername) %>%
    tidyr::unite(!!paste0(which, "_id"), !!which, id)
}

#### read input files ####
primers <- Biostrings::readDNAStringSet(primerfile)

names(tagfiles) <- sub("\\.fasta", "", basename(tagfiles))
tags <- lapply(tagfiles, Biostrings::readDNAStringSet)

loci <- readxl::read_xlsx(configfile, "Loci")
stopifnot(locus %in% loci$`Locus name`)
locus <- dplyr::filter(loci, `Locus name` == locus)
stopifnot(nrow(locus) == 1L)


fwd_is_tagged <- endsWith(locus$`Forward primer`, "-tag")
rev_is_tagged <- endsWith(locus$`Reverse primer`, "-tag")
stopifnot(fwd_is_tagged || rev_is_tagged)
if (fwd_is_tagged) {
  stopifnot(locus$`Forward primer` %in% names(tags))
  fwd_primer_name <- sub("-tag$", "", locus$`Forward primer`)
  fwd_primer_seq <- as.character(primers[fwd_primer_name])
  fwd_tags <- as.character(tags[[locus$`Forward primer`]])
  stopifnot(all(startsWith(names(fwd_tags), fwd_primer_name)))
  stopifnot(all(endsWith(fwd_tags, fwd_primer_seq)))
  fwd_primer <- Biobase::lcSuffix(fwd_tags)
  fwd_pad <- Biobase::lcPrefix(fwd_tags)
  fwd_tags <- fwd_tags %>%
    remove_lcprefix() %>%
    remove_lcsuffix()
  fwd_tag_nchar <- unique(nchar(fwd_tags))
  fwd_pad_spacer <- if (fwd_pad == "") "" else sprintf("N{%d}", nchar(fwd_pad))
  fwd_template <- sprintf(
    "%sN{%d}%s;o=%d;e=0.15",
    fwd_pad,
    fwd_tag_nchar,
    fwd_primer,
    nchar(fwd_primer) + fwd_tag_nchar
  ) %>%
    tibble::tibble(fwd_primer_name, fwd_template = .) %>%
    dplyr::mutate(fwd_primer_name = make.unique(fwd_primer_name, sep = "_"))
  fwd_tags <- sprintf(
    "X%s%s;o=%d;e=%.2f",
    fwd_pad,
    fwd_tags,
    nchar(fwd_tags),
    mindist(fwd_tags)/2/nchar(fwd_tags)
  ) %>%
    tibble::tibble(
      fwd_tagline = .,
      fwd_id = names(fwd_tags),
      fwd_tag = fwd_tags,
      fwd_primer = fwd_primer
    )
  fwd_platekey <- load_platekey(configfile, locus$`Forward primer`, "fwd")
} else {
  stopifnot(locus$`Forward primer` %in% names(primers))
  fwd_primer_name <- locus$`Forward primer`
  fwd_primer <- as.character(primers[fwd_primer_name])
  fwd_template <- sprintf("%s;o=%d;e=0.15", fwd_primer, nchar(fwd_primer))
  fwd_tags <- tibble::tibble(
    fwd_tagline = fwd_template,
    fwd_id = fwd_primer_name,
    fwd_primer,
    fwd_tag = ""
  )
  fwd_template <- tibble::tibble(fwd_primer_name, fwd_template = fwd_template)
  fwd_platekey <- tidyr::crossing(
    row = LETTERS[1:8],
    col = as.character(1:12),
    fwd_id = fwd_primer_name
  )
}
if (rev_is_tagged) {
  stopifnot(locus$`Reverse primer` %in% names(tags))
  rev_primer_name <- sub("-tag$", "", locus$`Reverse primer`)
  rev_primer_seq <- as.character(primers[rev_primer_name])
  rev_primer_rc <- as.character(Biostrings::reverseComplement(primers[rev_primer_name]))
  rev_tags <- as.character(tags[[locus$`Reverse primer`]])
  rev_tags_rc <- as.character(Biostrings::reverseComplement(tags[[locus$`Reverse primer`]]))
  stopifnot(all(startsWith(names(rev_tags), rev_primer_name)))
  rev_id <- sub(sprintf("^%s_", rev_primer_name), "", names(rev_tags))
  stopifnot(all(startsWith(rev_tags_rc, rev_primer_rc)))
  rev_primer_rc <- Biobase::lcPrefix(rev_tags_rc)
  rev_pad_rc <- Biobase::lcSuffix(rev_tags_rc)
  rev_primer <- Biobase::lcSuffix(rev_tags)
  rev_pad_spacer <- if (rev_pad == "") "" else sprintf("N{%d}", nchar(rev_pad))
  rev_tags <- rev_tags %>%
    remove_lcprefix() %>%
    remove_lcsuffix()
  rev_tags_rc <- rev_tags_rc %>%
    remove_lcprefix() %>%
    remove_lcsuffix()
  rev_tag_nchar <- unique(nchar(rev_tags_rc))
  rev_template <- sprintf(
    "%sN{%d}%s;o=%d;e=0.15",
    rev_primer_rc,
    rev_tag_nchar,
    rev_pad_rc,
    nchar(rev_primer_rc) + rev_tag_nchar
  ) %>%
    tibble::tibble(rev_primer_name, rev_template = .) %>%
    dplyr::mutate(rev_primer_name = make.unique(rev_primer_name, sep = "_"))
  rev_tags <- sprintf(
    "%s%sX;o=%d;e=%.2f",
    rev_tags_rc,
    rev_pad_spacer,
    nchar(rev_tags_rc),
    mindist(rev_tags_rc)/2/nchar(rev_tags_rc)
  ) %>%
    tibble::tibble(
      rev_tagline = .,
      rev_id = names(rev_tags),
      rev_tag = rev_tags,
      rev_tag_rc = rev_tags_rc,
      rev_primer = rev_primer,
      rev_primer_rc = rev_primer_rc
    )
  rev_platekey <- load_platekey(configfile, locus$`Reverse primer`, "rev")
} else {
  stopifnot(locus$`Reverse primer` %in% names(primers))
  rev_primer_name <- locus$`Reverse primer`
  rev_primer <- as.character(primers[rev_primer_name])
  rev_primer_rc <- as.character(Biostrings::reverseComplement(primers[rev_primer_name]))
  rev_template <- sprintf("%s;o=%d;e=0.15", rev_primer_rc, nchar(rev_primer_rc))
  rev_tags <- tibble::tibble(
    rev_tagline = rev_template,
    rev_id = rev_primer_name,
    rev_primer = rev_primer,
    rev_primer_rc = rev_primer_rc,
    rev_tag = "",
    rev_tag_rc = ""
  )
  rev_template <- tibble::tibble(rev_primer_name, rev_template = rev_template)
  rev_platekey <- tidyr::crossing(
    row = LETTERS[1:8],
    col = as.character(1:12),
    rev_id = rev_primer_name
  )
}

# write the template file
# this is for detecting the primers, orienting the reads, and trimming
# (but retain primers!)
if (!dir.exists(dirname(primers_fasta))) dir.create(dirname(primers_fasta), recursive = TRUE)
tidyr::crossing(fwd_template, rev_template) %$%
  sprintf(">%s_%s\n%s...%s", fwd_primer_name, rev_primer_name, fwd_template, rev_template) %>%
  writeLines(primers_fasta)

platekey <-
  dplyr::inner_join(fwd_platekey, rev_platekey, by = c("row", "col")) %>%
  dplyr::left_join(fwd_tags, by = "fwd_id") %>%
  dplyr::left_join(rev_tags, by = "rev_id")

sample_data <- readxl::read_xlsx(
  configfile,
  sheet = locus$`Locus name`,
  range = "B2:M9",
  col_names = as.character(1:12),
  col_types = "text"
) %>%
  dplyr::mutate(row = LETTERS[1:8]) %>%
  tidyr::pivot_longer(
    cols = 1:12,
    names_to = "col",
    values_to = "sample"
  ) %>%
  dplyr::filter(!is.na(sample)) %>%
  dplyr::full_join(platekey, by = c("row", "col")) %>%
  dplyr::mutate(
    sample = dplyr::coalesce(
      sample,
      paste(
        if(fwd_is_tagged) fwd_id else NULL,
        if (rev_is_tagged) rev_id else NULL
      ) %>%
        trimws() %>%
        chartr(old = " ", new = "_")
    )
  )

# write the tags file
# this is for actually demultiplexing
if (!dir.exists(dirname(tags_fasta))) dir.create(dirname(tags_fasta), recursive = TRUE)
sample_data %$%
  sprintf(">%s\n%s...%s", sample, fwd_tagline, rev_tagline) %>%
  writeLines(tags_fasta)

# write the minibar file
# this does both steps at once
if (!dir.exists(dirname(tags_table))) dir.create(dirname(tags_table), recursive = TRUE)
sample_data %>%
  dplyr::transmute(
    name = paste(fwd_primer_name, rev_primer_name, sep = "_"),
    fwd_tag = fwd_tag,
    fwd_primer = fwd_primer,
    rev_tag = rev_tag,
    rev_primer = rev_primer
  ) %>%
  write.table(tags_table, sep = "\t", row.names = FALSE, quote = FALSE)

# # make a set of tags in cutadapt fasta format for 3NDf--LR5
# rDNA_tags_cutadapt <-
#   tidyr::expand_grid(fwd = tagset_3NDf, rev = tagset_LR5) %>%
#   tidyr::unpack(c("fwd", "rev"), names_sep = "_") %>%
#   glue::glue_data(">{fwd_name}_{rev_name}\n{fwd_seq}...{rev_revcomp}")
#
# # make a set of tags in minibar table format for 3NDf--LR5
# rDNA_tags_minibar <-
#   tidyr::expand_grid(fwd = tagset_3NDf, rev = tagset_LR5) %>%
#   tidyr::unpack(c("fwd", "rev"), names_sep = "_") %>%
#   dplyr::select(fwd_name, rev_name, fwd_tag, fwd_primer, rev_tag, rev_primer) %>%
#   tidyr::unite("name", c(fwd_name, rev_name))
#
# # ITS(1) is already in a fine format for cutadapt,
# # which can do single-end demultiplexing
#
# # make a set of tags in minibar table format for ITS1-ITS4
# # pretend that half the ITS4 is a tag (with only one variant)
# # and half is the primer
#
# ITS_tags_minibar <-
#   tibble::enframe(as.character(seq_ITS1)) %>%
#   dplyr::transmute(
#     name = name,
#     fwd_tag = remove_lcsuffix(value),
#     fwd_primer = Biobase::lcSuffix(value),
#     rev_tag = "GCATATCAATAAGCGG",
#     rev_primer = "AGGA"
#   )
#
# # write out the tags; in this form they are not designed for a specific plate
# # layout.
#
# # if (!is.null(tags_ITS1_fasta)) {
# #   Biostrings::writeXStringSet(seq_ITS1, tags_ITS1_fasta)
# # }
# # if (!is.null(tags_ITS1_table)) {
# #   write.table(ITS_tags_minibar, tags_ITS1_table, sep = "\t", row.names = FALSE,
# #               quote = FALSE)
# # }
# #
# # if (!is.null(tags_3NDf_LR5_fasta)) {
# #   writeLines(rDNA_tags_cutadapt, tags_3NDf_LR5_fasta)
# # }
# #
# # if (!is.null(tags_3NDf_LR5_table)) {
# #   write.table(rDNA_tags_minibar, tags_3NDf_LR5_table, sep = "\t", row.names = FALSE,
# #               quote = FALSE)
# # }
#
# # if we have a plate, then convert tag (pairs) to wells
#
#
#
#
#   sample_data <- readxl::read_xlsx(
#     configfile,
#     sheet = locus$`Locus name`,
#     range = "B2:M9",
#     col_names = as.character(1:12),
#     col_types = "text"
#   ) %>%
#     dplyr::mutate(row = LETTERS[1:8]) %>%
#     tidyr::pivot_longer(
#       cols = 1:12,
#       names_to = "col",
#       values_to = "sample"
#     ) %>%
#     dplyr::filter(!is.na(sample)) %>%
#     dplyr::left_join(platekey, by = c("row", "col"))
#
#   sample_data %$%
#     stringi::stri_replace_all_fixed(
#       rDNA_tags_cutadapt,
#       paste0(">", tagname),
#       paste0(">", dplyr::coalesce(sample, tagname)),
#       vectorize_all = FALSE
#     ) %>%
#     writeLines(sample_tags_fasta[i])
#
#   sample_data %>%
#     dplyr::left_join(rDNA_tags_minibar, by = c("tagname" = "name")) %>%
#     dplyr::transmute(
#       name = dplyr::coalesce(sample, tagname),
#       fwd_tag = fwd_tag,
#       fwd_primer = fwd_primer,
#       rev_tag = rev_tag,
#       rev_primer = rev_primer
#     ) %>%
#     write.table(sample_tags_table[i], sep = "\t", row.names = FALSE,
#                 quote = FALSE)
#
#   sample_data$sample %>%
#     purrr::discard(is.na) %>%
#     writeLines(sample_names[i])
#
#   sample_data_single <-
#     readxl::read_xlsx(
#       sample_plate[i],
#       sheet = "Single",
#       col_names = TRUE,
#       col_types = "text"
#     )
#
#   set_names(
#     seq_ITS1,
#     stringi::stri_replace_all_regex(
#       names(seq_ITS1),
#       sprintf("^%s$", sample_data_single$Primer),
#       sample_data_single$Sample,
#       vectorize_all = FALSE
#     )
#   ) %>%
#     Biostrings::writeXStringSet(sample_tags_single_fasta[i])
#
#   writeLines(sample_data_single$Sample, sample_names_single[i])
# }
