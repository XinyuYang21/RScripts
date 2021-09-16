extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}

# Extract wanted legend
#wantedLegend <- extractLegend(template)
# ####### Conver igraph to ggraph

