#' Generate plots of results
#' @param results Output from summarise_results
#' @return List of ggplot objects
plot_result <- function(results) {
  require(ggplot2)
  require(ggtree)
  
  plots <- list()
  
  # AICc Plot
  plots$aicc_plot <- ggplot(results$single_models$AICc, aes(x = Model, y = AICc, fill = Model)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = round(AICc, 1)), vjust = -0.5) +
    labs(title = "Model Comparison", y = "AICc") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Tree Plot
  if (!is.null(results$data$regime)) {
    regimes <- na.omit(unique(results$data$regime))
    colors <- setNames(rainbow(length(regimes)), regimes)
    
    plots$tree_plot <- ggtree(results$surface$named_tree) %<+% results$data +
      geom_tree(aes(color = regime), size = 0.7) +
      scale_color_manual(values = colors, name = "Regime") +
      theme(legend.position = "right") +
      labs(title = "Phylogenetic Regimes")
  }
  
  return(plots)
}