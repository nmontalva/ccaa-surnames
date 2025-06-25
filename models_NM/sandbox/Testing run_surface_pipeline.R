# debugger.R
#— Iteration 3: test prep_traits(), fit_single_regime(), run_surface()

library(geiger)
library(surface)

# 1) Source your analysis file
#source("surface_analysis.R")

# assume you have GM_df, consensus_tree in your environment

# 2) Test prep_traits
cat("\n--- prep_traits() ---\n")
dfG  <- prep_traits(GM_df, "G")
print(head(dfG[, c("community","G","G_adj","G_logit")]))

# 3) Test fit_single_regime
cat("\n--- fit_single_regime() ---\n")
vecG <- setNames(dfG$G_logit, dfG$community)
aicsG <- fit_single_regime(vecG, consensus_tree)
print(aicsG)

# 4) Test run_surface
cat("\n--- run_surface() ---\n")
surfG <- run_surface(vecG, consensus_tree)
cat(" forward steps:", length(surfG$fwd), "\n")
cat(" backward steps:", length(surfG$bwd), "\n")
