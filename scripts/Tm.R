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

melt_range("ACGACGTTGTAAAAACIMGIGTITCVGAIGTHAARCC", nucleic.acid.conc = 500e-9, hybridisation.type = "dnadna", Mg.conc=0.005)
