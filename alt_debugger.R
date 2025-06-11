# For quick testing with small tree
test_results <- run_evolutionary_analysis(
  data = GM_df,
  variable = "G",
  tree = consensus_tree,
  verbose = TRUE
)

# Check results
print(test_results$regimes)
print(test_results$plots$tree_plot)

# For full analysis (run overnight)
#full_results <- run_evolutionary_analysis(
#  data = GM_df,
#  variable = "G",
#  tree = y_total,
#  verbose = TRUE
#)

# Save results
#saveRDS(full_results, "surface_results_G.rds")