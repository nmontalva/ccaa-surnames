## Incremental Testing debugger

# 1. First test just data preparation:

result <- evolutionary_analysis(
  data = GM_df,
  variable = "G",
  tree = consensus_tree,
  steps = 1
)
str(result$prepped)

# 2. Add single regime models:

result <- evolutionary_analysis(
  data = GM_df,
  variable = "G",
  tree = consensus_tree,
  steps = 1:2
)
print(result$single_regime$AIC)

# 3. Add SURFACE analysis:

result <- evolutionary_analysis(
data = GM_df,
variable = "G", 
tree = consensus_tree,
steps = 1:3
)

# 4. Add summary

result <- evolutionary_analysis(
  data = GM_df,
  variable = "G",
  tree = consensus_tree,
  steps = 1:4
)
print(result$summary$regimes)

# 5. Complete analysis with plots:

result <- evolutionary_analysis(
  data = GM_df,
  variable = "G",
  tree = consensus_tree,
  steps = 1:5
)
print(result$plots$tree_plot)