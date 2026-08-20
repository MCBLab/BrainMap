#!/usr/bin/env Rscript
# Placeholder figure PDFs so main.tex compiles end to end while the real
# renders are pending. Each placeholder carries the caption slug and the
# command that will replace it, so a printed draft is self-documenting.
#
#   Rscript figures/make_placeholders.R
#
# Replace with the real thing via figures/make_figures.R (see IMPLEMENTATION.md).

specs <- list(
  list(file = "fig2_interface.pdf",  w = 9, h = 6.5,
       title = "Figure 2 -- VEGABrain web interface",
       body  = c("Screenshots: the three query modes (gene list / single gene / ontology),",
                 "the micro-macro region toggle and the reference-map modal.",
                 "Source: screen capture of https://mcblab.github.io/BrainMap/")),
  list(file = "fig3_gene_sox10.pdf", w = 9, h = 7,
       title = "Figure 3 -- Single-gene expression map (SOX10)",
       body  = c("DK lateral / DK medial / aseg sagittal / aseg coronal,",
                 "faceted over the five developmental windows.",
                 "Source: make_figures.R -> fig_gene('SOX10')")),
  list(file = "fig4_ontology.pdf",   w = 9, h = 7,
       title = "Figure 4 -- Ontology composite ssGSEA map",
       body  = c("GOBP_FOREBRAIN_GENERATION_OF_NEURONS, micro and macro scales.",
                 "Headline figure for the per-region enrichment contribution.",
                 "Source: make_figures.R -> fig_ontology()")),
  list(file = "fig5_genelist.pdf",   w = 9, h = 7,
       title = "Figure 5 -- User-submitted gene-list map",
       body  = c("On-demand ssGSEA for GAD1, GAD2, SLC32A1, GABRB3 (GABAergic set),",
                 "computed per request and mapped across development.",
                 "Source: make_figures.R -> fig_genelist()"))
)

for (s in specs) {
  pdf(s$file, width = s$w, height = s$h)
  op <- par(mar = c(0, 0, 0, 0))
  plot.new()
  rect(0.02, 0.02, 0.98, 0.98, border = "grey40", lty = 2, lwd = 2)
  text(0.5, 0.70, s$title, cex = 1.5, font = 2, col = "grey20")
  text(0.5, 0.58, "PLACEHOLDER", cex = 2.2, font = 2, col = "firebrick")
  text(0.5, seq(0.44, by = -0.06, length.out = length(s$body)),
       s$body, cex = 1.0, col = "grey30")
  par(op)
  dev.off()
  message("wrote ", s$file)
}
