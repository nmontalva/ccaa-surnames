source("models_NM/v01/prep_traits.R")
source("models_NM/v01/fit_single_regime.R")
source("models_NM/v01/run_surface.R")
source("models_NM/v01/summarise_results.R")
source("models_NM/v01/compare_models.R")
source("models_NM/v01/plot_results.R")
source("models_NM/v01/evolutionary_analysis.R")
# Test with your data
results <- evolutionary_analysis(
  data = GM_df,
  variables = c("G","M"),
  tree = y_total,
  steps = 1:6,
  verbose = TRUE
)
save(results, file="models_NM/sandbox/model_results.RData")

# 2. Check results
print(results$G$plots$tree_plot)
print(results$G$plots$aicc_plot)

# 3. View all parameters
results$G$summary$regimes
results$G$summary$parameters
