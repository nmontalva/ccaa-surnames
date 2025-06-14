plot_result <- function(results) {
  require(ggplot2)
  require(ggtree)
  require(viridis)  # For better color scales
  
  plots <- list()
  
  # 1. AIC Comparison Plot (always created)
  plots$aic_plot <- ggplot(results$single_models$AIC, aes(x = Model, y = AIC, fill = Model)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_text(aes(label = round(AIC, 1)), vjust = -0.5, size = 3.5) +
    labs(title = paste("Model Comparison for", results$prepped$original_var),
         subtitle = "Akaike Information Criterion",
         y = "AIC", x = "") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10))
  
  # 2. Tree Plot (now without tip number restriction)
  if (!is.null(results$surface$named_tree)) {
    # Create a color palette that works well for multiple regimes
    n_regimes <- length(unique(results$regimes$regime))
    regime_colors <- viridis(n_regimes, option = "D")
    names(regime_colors) <- results$regimes$regime
    
    # Prepare tip labels with regime information
    tip_data <- results$data
    tip_data$theta_label <- sapply(1:nrow(tip_data), function(i) {
      r <- tip_data$regime[i]
      theta <- results$regimes$theta_original[results$regimes$regime == r]
      sprintf("%s (θ=%.2f)", r, theta)
    })
    
    # Create the tree plot with optimized settings for large trees
    plots$tree_plot <- ggtree(results$surface$named_tree, size = 0.4, layout = "rectangular") %<+% tip_data +
      geom_tree(aes(color = regime), size = 0.5) +
      scale_color_manual(values = regime_colors,
                        name = "Regime",
                        labels = unique(tip_data$theta_label)) +
      theme_tree2() +
      theme(legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            plot.title = element_text(size = 11)) +
      ggtitle(paste("Phylogenetic Regimes for", results$prepped$original_var)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    # Add tip points only for trees under 150 tips
    if (length(results$data$community) <= 150) {
      plots$tree_plot <- plots$tree_plot + 
        geom_tippoint(aes(color = regime), size = 1.5, alpha = 0.7)
    }
  }
  
  # 3. Regime Distribution Plot
  if (!is.null(results$data$regime)) {
    reg_summary <- as.data.frame(table(results$data$regime))
    colnames(reg_summary) <- c("regime", "count")
    reg_summary$percent <- round(100 * reg_summary$count / sum(reg_summary$count), 1)
    
    plots$regime_dist <- ggplot(reg_summary, aes(x = regime, y = count, fill = regime)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = paste0(count, "\n(", percent, "%)")), 
                vjust = -0.3, size = 3) +
      scale_fill_manual(values = regime_colors) +
      labs(title = paste("Regime Distribution for", results$prepped$original_var),
           x = "Regime", y = "Count") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none",
            plot.title = element_text(size = 12))
  }
  
  return(plots)
}