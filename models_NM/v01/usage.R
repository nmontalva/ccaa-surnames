# Run analysis (now with large tree)
results <- evolutionary_analysis(
  data = GM_df,
  variables = c("G", "M"),
  tree = y_total,  # Your large tree
  steps = 1:5
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
print(results$G$summary$regimes)
print(results$M$summary$regimes)

# View plots for a specific variable
print(results$G$plots$tree_plot)
print(results$M$plots$tree_plot)
