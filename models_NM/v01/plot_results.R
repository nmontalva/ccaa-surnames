plot_result <- function(results) {
  require(ggplot2)
  require(ggtree)
  #require(paletteer)  # For better color scales
  require(RColorBrewer)
  plots <- list()

  color_pal <- c("#E69F00","#56B4E9","#009E73","#D55E00","#0072B2" ) 
  # 1. AICc Comparison Plot (always created)
  n_models <- if (!is.null(nrow(results$single_models$AICc))) nrow(results$single_models$AICc) else length(results$single_models$AICc)
  
  plots$aicc_plot <- ggplot(results$single_models$AICc, aes(x = Model, y = AICc, fill = Model)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_text(aes(label = round(AICc, 1)), vjust = -0.5, size = 3.5) +
    #scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.2, end = 0.8)+
    scale_fill_manual(values = color_pal[seq_len(n_models)])+
    labs(title = paste("Model Comparison for", results$prepped$original_var),
         subtitle = "Corrected Akaike Information Criterion (AICc)",
         y = "AICc", x = "") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10))
  if (!is.null(results$surface$named_tree)) {
    tip_data <- results$data %>%
      dplyr::filter(!is.na(regime)) %>% 
      left_join(results$regimes, by = "regime") %>%
      mutate(
        theta_fmt = ifelse(!is.na(theta_original), sprintf("%.2f", theta_original), "NA"),
        alpha_fmt = ifelse(!is.na(alpha), sprintf("%.2f", alpha), "NA"),
        label = sprintf("%s (θ=%s, α=%s)", regime, theta_fmt, alpha_fmt)
      )
    n_tips <- length(results$data$community)
    line_size <- ifelse(n_tips > 100, 0.5, 0.8)  # thicker lines
  # 2. Tree Plot (now without tip number restriction)
  if (!is.null(results$surface$named_tree)) {
    # Create a color palette that works well for multiple regimes
    regimes_unique <- unique(results$regimes$regime)
    regimes_unique <- regimes_unique[!is.na(regimes_unique)]  # Quitar NA
    n_regimes <- length(unique(results$regimes$regime))
    #regime_colors <- viridis(n_regimes, option = "D",begin = 0.1, end = 0.9)
    regime_colors <- setNames(color_pal[1:n_regimes], unique(results$regimes$regime))
    names(regime_colors) <- results$regimes$regime
    
    # Prepare tip labels with regime information
    tip_data <- results$data %>% dplyr::filter(!is.na(regime)) %>%
      left_join(results$regimes, by = "regime")
    tip_data$theta_label <- sapply(1:nrow(tip_data), function(i) {
      r <- tip_data$regime[i]
      theta <- results$regimes$theta_original[results$regimes$regime == r]
      sprintf("%s (θ=%.2f)", r, theta)
    })
    
    # Create the tree plot with optimized settings for large trees
    plots$tree_plot <- ggtree(results$surface$named_tree, size = line_size, layout = "rectangular") %<+% tip_data +
      geom_tree(aes(color = regime), size = 0.5) +
      scale_color_manual(values = regime_colors,
                        name = "Regime",
                        labels = unique(tip_data$label)) +
      theme_tree2() +
      theme(legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 11)) +
      ggtitle(paste("Phylogenetic Regimes for", results$prepped$original_var)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    # Add tip points only for trees under 150 tips
    if (length(results$data$community) <= 150) {
      plots$tree_plot <- plots$tree_plot + 
        geom_tippoint(aes(color = regime), size = 1.5, alpha = 0.7)
    }
  }
    if (!interactive()) {
      plots$tree_plot <- plots$tree_plot + 
        theme(plot.background = element_rect(fill = "white", color = NA))
    }
  } # Better export visualization
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
