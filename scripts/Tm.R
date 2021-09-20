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

threeNDf <- "GGCAAGTCTGGTGCCAG"
LR5 <- "TCCTGAGGGAAACTTCG"
EF1_938f <- "GCYCCYGGHCAYCGTGAYTTYAT"
EF1_938f2 <- "GCYCCIGGICAYCGTGAYTTYAT"
EF1_2218r <- "ATGACACCRACRGCRACRGTYTG"
RPB1_aAf <- "GAGTGTCCGGGGCATTTYGG"
RPB1_aCr <- "ARAARTCBACHCGYTTBCCCAT"
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
melt_range(threeNDf, nucleic.acid.conc = 100e-9, hybridisation.type = "dnadna", Na.conc=1)
c(enthalpy_range(EF1_938f2, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005),
entropy_range(EF1_938f2, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005))[c(1,3,2,4)] |>
    paste(collapse=" ") |>
    chartr(old=".", new = ",")
