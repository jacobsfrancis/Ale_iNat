# ------------------------------------------------------------
# Theoretical bipartite networks varying in nestedness
# Plants on top, pollinators on bottom
# Exports individual panels for PowerPoint/poster assembly
# ------------------------------------------------------------

library(bipartite)

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

plant_col <- "#003366"       # plants on top, blue
pollinator_col <- "#CC0000"  # pollinators on bottom, red

# Use slightly darker grey + alpha because transparency lightens on white
link_col <- "#999999"
link_alpha <- 0.45

# ------------------------------------------------------------
# Networks
# ------------------------------------------------------------
# Important:
# In bipartite::plotweb():
#   rows    = lower trophic level  = bottom
#   columns = higher trophic level = top
#
# For this figure:
#   rows    = pollinators = bottom = red
#   columns = plants      = top    = blue
# ------------------------------------------------------------

low_nested <- matrix(
  c(
    1, 1, 0, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 1
  ),
  nrow = 8,
  byrow = TRUE
)

medium_nested <- matrix(
  c(
    1, 1, 1, 1, 1, 0, 0, 0,
    1, 1, 1, 1, 0, 1, 0, 0,
    1, 1, 1, 0, 1, 0, 1, 0,
    1, 1, 0, 1, 0, 1, 0, 1,
    1, 0, 1, 0, 1, 0, 0, 0,
    0, 1, 0, 1, 0, 1, 0, 0,
    1, 0, 0, 0, 1, 0, 1, 0,
    0, 1, 0, 0, 0, 1, 0, 1
  ),
  nrow = 8,
  byrow = TRUE
)

high_nested <- matrix(
  c(
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 0,
    1, 1, 1, 1, 1, 1, 0, 0,
    1, 1, 1, 1, 1, 0, 0, 0,
    1, 1, 1, 1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0
  ),
  nrow = 8,
  byrow = TRUE
)

network_list <- list(
  low_nestedness = low_nested,
  medium_nestedness = medium_nested,
  high_nestedness = high_nested
)

# Add row/column names
# These are hidden in the plot, but useful for checking and analysis
for (i in seq_along(network_list)) {
  rownames(network_list[[i]]) <- paste0("Pollinator_", 1:8)
  colnames(network_list[[i]]) <- paste0("Plant_", 1:8)
}

# ------------------------------------------------------------
# Optional: check nestedness values
# ------------------------------------------------------------

nodf_results <- sapply(network_list, function(web) {
  bipartite::networklevel(web, index = "NODF")
})

nestedness_results <- sapply(network_list, function(web) {
  bipartite::networklevel(web, index = "nestedness")
})

print(nodf_results)
print(nestedness_results)

# ------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------

plot_nested_web <- function(web, title = NULL, label_cex = 2) {
  
  bipartite::plotweb(
    web,
    sorting = "normal",
    empty = FALSE,
    
    lower_labels = FALSE,
    higher_labels = FALSE,
    
    lower_color = pollinator_col,
    lower_border = "same",
    
    higher_color = plant_col,
    higher_border = "same",
    
    link_color = link_col,
    link_border = "same",
    link_alpha = link_alpha,
    
    spacing = 0.25,
    box_size = 0.08,
    curved_links = FALSE,
    
    # increased margin room for larger labels
    mar = c(4.2, 1, 4.2, 1)
  )
  
  if (!is.null(title)) {
    title(main = title, line = 1.3, cex.main = 1.1)
  }
  
  mtext(
    "Plants",
    side = 3,
    line = 0.7,
    cex = label_cex,
    font = 2,
    col = plant_col
  )
  
  mtext(
    "Pollinators",
    side = 1,
    line = 1.2,
    cex = label_cex,
    font = 2,
    col = pollinator_col
  )
}

# ------------------------------------------------------------
# Preview all three together in R
# ------------------------------------------------------------

par(mfrow = c(1, 3))

plot_nested_web(network_list$low_nestedness)
plot_nested_web(network_list$medium_nestedness)
plot_nested_web(network_list$high_nestedness)

# ------------------------------------------------------------
# Export individual panels
# ------------------------------------------------------------

dir.create("nestedness_network_exports", showWarnings = FALSE)

# PDF export
for (nm in names(network_list)) {
  
  pdf(
    file = file.path("nestedness_network_exports", paste0(nm, ".pdf")),
    width = 4,
    height = 3.3,
    useDingbats = FALSE
  )
  
  plot_nested_web(network_list[[nm]])
  
  dev.off()
}

