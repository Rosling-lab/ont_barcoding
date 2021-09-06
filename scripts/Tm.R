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

EF1_938f <- "GCYCCYGGHCAYCGTGAYTTYAT"
EF1_2218r <- "ATGACACCRACRGCRACRGTYTG"
RPB1_aAf <- "GAGTGTCCGGGGCATTTYGG"
RPB1_aCr <- "ARAARTCBACHCGYTTBCCCAT"
bRPB2_6F <- "TGGGGYATGGTITGYCCYGC"
bRPB2_7R <- "GAYTGRTTRTGRTCRGGGAAVGG"
fRPB2_5F <- "GAYGAYMGWGATCAYTTYGG"
bRPB2_6R2 <- "GGRCAIACCATICCCCARTG"

melt_range(RPB1_aCr, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
enthalpy_range(bRPB2_6R2, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
entropy_range(bRPB2_6R2, nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
