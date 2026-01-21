source("models_NM/v01/prep_traits.R")
source("models_NM/v01/fit_single_regime.R")
source("models_NM/v01/run_surface.R")
source("models_NM/v01/summarise_results.R")
source("models_NM/v01/compare_models.R")
source("models_NM/v01/plot_results.R")
source("models_NM/v01/evolutionary_analysis.R")

# Run analysis (now with large tree)
results <- evolutionary_analysis(
  data = GM_df,
  variables = c("G", "M"),
  tree = y_total,  # Your large tree
  steps = 1:6
)


# View the tree plot (will work for 200+ tips)
print(results$G$plots$tree_plot)

# To customize further:
results$G$plots$tree_plot +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1")

# To customize further:
results$M$plots$tree_plot +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1")

# Access results for each variable