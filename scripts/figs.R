library(magrittr)
library(ggplot2)

theme_set(theme_bw())

# draw the rDNA, primers, and amplicons
draw_amplicon_plot <- function(amplicons_file, ylim = c(-25, 40), arrowlen = 0.1) {
  readr::read_csv(amplicons_file) %>%
    dplyr::mutate(xmid = (x + xend) / 2,
                  label = tidyr::replace_na(label, "")) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend, color = color,
               linetype = linetype, size = size)) +
    geom_segment(data = ~dplyr::filter(., !arrow, color == "black")) +
    geom_segment(data = ~dplyr::filter(., !arrow, color != "black")) +
    geom_segment(
      data = ~dplyr::filter(., arrow),
      arrow = arrow(length = unit(arrowlen, "inches")),
      linejoin = "mitre"
    ) +
    geom_text(aes(x = xmid, y = y_label, label = label, color = labelcolor,
                  size = size_label)) +
    scale_linetype_identity() +
    scale_color_identity() +
    scale_size_identity() +
    ylim(ylim) +
    coord_equal(ratio = 20) +
    theme(axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()
    )
}

draw_amplicon_plot("figures/rDNA_amplicons.csv")
ggsave("figures/rDNA_amplicons.png", height = 3, width = 8)
draw_amplicon_plot("figures/barcoding_strategy.csv", ylim = c(-25, 7), arrowlen = 0.3)
ggsave("figures/barcoding_strategy.png", height = 4, width = 8)

col2hex <- function(color)
{
  clr <- col2rgb(color)
  sprintf("#%02X%02X%02X", clr[1],clr[2],clr[3])
}
