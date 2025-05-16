#' Generate plots of results
#' @param results Output from summarise_results
#' @return List of ggplot objects
plot_result <- function(results) {
  require(ggplot2)
  require(ggtree)
  
  plots <- list()
  
  # AIC plot
  plots$aic_plot <- ggplot(results$single_models$AIC, aes(x = Model, y = AIC, fill = Model)) +
    geom_col() +
    geom_text(aes(label = round(AIC, 1)), vjust = -0.5) +
    labs(title = "AIC Comparison") +
    theme_minimal()
  
  # Only create tree plot if <= 50 tips
  if (length(results$data$community) <= 50) {
    n_regimes <- length(unique(results$regimes$regime))
    pal <- if (n_regimes == 1) {
      setNames("#1f77b4", results$regimes$regime[1])
    } else {
      setNames(rainbow(n_regimes), results$regimes$regime)
    }
    
    tip_data <- data.frame(
      label = results$data$community,
      regime = results$data$regime,
      stringsAsFactors = FALSE
    )
    
    theta_labels <- paste0(
      results$regimes$regime, 
      ": θ=", 
      round(if ("theta_original" %in% colnames(results$regimes)) {
        results$regimes$theta_original
      } else {
        results$regimes$theta
      }, 3)
    )
    
    plots$tree_plot <- ggtree(results$surface$named_tree) %<+% tip_data +
      geom_tree(aes(color = regime), size = 0.8) +
      scale_color_manual(
        values = pal,
        labels = theta_labels,
        name = "Regime (θ)"
      ) +
      theme_tree2() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 12, face = "bold")
      )
  }
  
  return(plots)
}