library(rmelting)

max_replace <- function(sequence) {
    chartr("YRSWMKBDHVN", "CGGACGGGCGG", sequence)
}
max_comp <- function(sequence) {
    chartr("ACGTYRSWMKBDHVNI", "TGCAGCCTGCCCGCCC", sequence)
}
min_replace <- function(sequence) {
    chartr("YRSWMKBDHVN", "TAGAATTAAAA", sequence)
}
min_comp <- function(sequence) {
    chartr("ACGTYRSWMKBDHVNI", "TGCAATCTTAATTTTG", sequence)
}
mismatch_bases <- list(
    A = c("A", "C", "G"),
    C = c("A", "C", "T"),
    G = c("A", "G", "T"),
    T = c("C", "G", "T"),
    Y = c("C", "T"),
    R = c("A", "G"),
    S = c("A", "T"),
    W = c("C", "G"),
    K = c("G", "T"),
    M = c("A", "C"),
    B = "T",
    D = "G",
    H = "C",
    V = "A",
    N = character(),
    I = character()
)

oneoff_min <- function(sequence, n = NULL) {
    if (is.null(n)) {
        purrr::map_dfr(seq_len(nchar(sequence)), oneoff_min, sequence = sequence)
    } else {
        expand.grid(
            sequence = glue::glue("{min_replace(substr(sequence, 0, n-1))}{expand_bases[[substr(sequence, n, n)]]}{min_replace(substr(sequence, n+1, nchar(sequence)+1))}"),
            comp.sequence = glue::glue("{min_comp(substr(sequence, 0, n-1))}{mismatch_bases[[substr(sequence, n, n)]]}{min_comp(substr(sequence, n+1, nchar(sequence)+1))}"),
            stringsAsFactors = FALSE
        )
    }
}

expand_bases = list(
    A = "A",
    C = "C",
    G = "G",
    T = "T",
    Y = c("C", "T"),
    R = c("A", "G"),
    S = c("C", "G"),
    W = c("A", "T"),
    K = c("G", "T"),
    M = c("A", "C"),
    B = c("C", "G", "T"),
    D = c("A", "G", "T"),
    H = c("A", "C", "T"),
    V = c("A", "C", "G"),
    N = c("A", "C", "G", "T"),
    I = "I"
)
melt_range <- function(sequence, ...) {
    c(
        min = melting(min_replace(sequence), comp.sequence = min_comp(sequence), ...)$Results$`Melting temperature (C)`,
        max = melting(max_replace(sequence), comp.sequence = max_comp(sequence), ...)$Results$`Melting temperature (C)`
    )
}

enthalpy_range <- function(sequence, ...) {
    c(
        min = melting(min_replace(sequence), comp.sequence = min_comp(sequence), ...)$Results$`Enthalpy (J)`,
        max = melting(max_replace(sequence), comp.sequence = max_comp(sequence), ...)$Results$`Enthalpy (J)`
    )
}
entropy_range <- function(sequence, ...) {
    c(
        min = melting(min_replace(sequence), comp.sequence = min_comp(sequence), ...)$Results$`Entropy (J)`,
        max = melting(max_replace(sequence), comp.sequence = max_comp(sequence), ...)$Results$`Entropy (J)`
    )
}

thermo_range <- function(sequence, name = "", ...) {
    if (length(sequence) == 1) {
        stopifnot(length(name) == 1)
        rbind(
            min = unlist(melting(min_replace(sequence), comp.sequence = min_comp(sequence), ...)$Results),
            max = unlist(melting(max_replace(sequence), comp.sequence = max_comp(sequence), ...)$Results)
        ) |>
            as.data.frame() |>
            dplyr::select(H = `Enthalpy (J)`, S = `Entropy (J)`) |>
            tibble::rownames_to_column("optimum") |>
            tidyr::pivot_wider(names_from = "optimum", values_from = c("H", "S")) |>
            dplyr::transmute(
                Name = name,
                Sequence = sequence,
                dH = H_min,
                dHmax = H_max,
                dS = S_min,
                dSmax = S_max
            )
    } else {
        purrr::map2_dfr(sequence, name, thermo_range, ...)
    }
}

thermo_data <- function(sequ) {}

threeNDf <- "GGCAAGTCTGGTGCCAG"
LR5 <- "TCCTGAGGGAAACTTCG"
EF1_938f <- "GCYCCYGGHCAYCGTGAYTTYAT"
EF1_938f2 <- "GCYCCIGGICAYCGTGAYTTYAT"
EF1_2218r <- "ATGACACCRACRGCRACRGTYTG"
RPB1_aAf <- "GAGTGTCCGGGGCATTTYGG"
RPB1_aCr <- "ARAARTCBACHCGYTTBCCCAT"
gRPB1_A_for <- "GAKTGTCCKGGWCATTTTGG"
fRPB1_C_rev <- "CIGCDATITCRTTRTCCATRTA"
bRPB2_6F <- "TGGGGYATGGTITGYCCYGC"
bRPB2_7R <- "GAYTGRTTRTGRTCRGGGAAVGG"
fRPB2_5F <- "GAYGAYMGWGATCAYTTYGG"
fRPB2_5F2 <- "GAYGAYIGIGAYCAYTTYGG"
bRPB2_6R2 <- "GGRCAIACCATICCCCARTG"
bMCM7_F <- "TTYCARGARGTIAARATICARGARATGG"
bMCM7_R <- "TCCATYTTRTCRAAYTCRTC"
MCM7_709for <- "ACIMGIGTITCVGAYGTHAARCC"
MCM7_16r <- "GTYTGYTGYTCCATIACYTCRTG"
Bsens <- "ATYACICAYTCIYTIGGTGG"
Bsens2 <- "ATCACWCACTCICTIGGTGG"
Brev <- "CATGAAGAARTGIAGACGIGG"
Brev2 <- "AARAARTGIAGSCGIGGGAAIGG"

bMCM7_F_tagF <- "ACGACGTTGTAAAATTYCARGARGTIAARATICARGARATGG"
bMCM7_R_tagR <- "CATTAAGTTCCCATTATCCATYTTRTCRAAYTCRTC"
MCM7_709for_tagF <- "ACGACGTTGTAAAAACIMGIGTITCVGAYGTHAARCC"
MCM7_16r_tagR <- "CATTAAGTTCCCATTAGTYTGYTGYTCCATIACYTCRTG"
fRPB2_5F2_tagF <- "ACGACGTTGTAAAAGAYGAYIGIGAYCAYTTYGG"
bRPB2_7R_tagR <- "CATTAAGTTCCCATTAGAYTGRTTRTGRTCRGGGAAVGG"
gRPB1_A_for_tagF <- "ACGACGTTGTAAAAGAKTGTCCKGGWCATTTTGG"
fRPB1_C_rev_tagR <- "CATTAAGTTCCCATTACIGCDATITCRTTRTCCATRTA"
EF1_938f_tagF <- "ACGACGTTGTAAAAGCYCCYGGHCAYCGTGAYTTYAT"
EF1_938f2_tagF <- "ACGACGTTGTAAAAGCYCCIGGICAYCGTGAYTTYAT"
EF1_2218r_tagR <- "CATTAAGTTCCCATTAATGACACCRACRGCRACRGTYTG"
Bsens_tagF <- "ACGACGTTGTAAAAATYACICAYTCIYTIGGTGG"
Bsens2_tagF <- "ACGACGTTGTAAAAATCACWCACTCICTIGGTGG"
Brev_tagR <- "CATTAAGTTCCCATTACATGAAGAARTGIAGACGIGG"
Brev2_tagR <- "CATTAAGTTCCCATTAAARAARTGIAGSCGIGGGAAIGG"

min_comp(EF1_2218r)

melting(min_replace(EF1_2218r), "TACTGTGGTTGTCGTTGTCAAAC", nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
melt_range(RPB1_aCr, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
enthalpy_range(bRPB2_6R2, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
entropy_range(bRPB2_6R2, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
thermo_range(c(fRPB2_5F, bRPB2_6R2), c("fRPB2-5F", "bRPB2-6R2"),
             nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
purrr::pmap_dfr(oneoff_min(EF1_2218r), function(...) melting(...)$Results,
                nucleic.acid.conc = 4, hybridisation.type = "dnadna", Na.conc=1) |>
    dplyr::arrange(`Enthalpy (J)` / `Entropy (J)`) |>
    dplyr::slice(1)

melt_range(threeNDf, nucleic.acid.conc = 100e-9, hybridisation.type = "dnadna", Na.conc=1)
c(enthalpy_range(Brev2_tagR, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005),
entropy_range(Brev2_tagR, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005))[c(1,3,2,4)] |>
    paste(collapse="\t") |>
    chartr(old=".", new = ",") |>
    cat()
