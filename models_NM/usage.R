# Test with your data
results <- evolutionary_analysis(
  data = GM_df,
  variables = c("G"),
  tree = y_total,
  steps = 1:5,
  verbose = TRUE
)

# 2. Check results
print(results$G$plots$tree_plot)
print(results$G$plots$aic_plot)

# 3. View all parameters
results$G$summary$regimes
results$G$summary$parameters