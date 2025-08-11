# Selected results

## Idea:
# 1. After running evolutionary analysis.R
# 2. Pick the resulting object (i.e. "results")
# 3. Use it to build tailored outputs.

## A table to compare models # LISTO
# Compare BM, single OU, multi OU
# For both G and M
# Include relevant numbers such as:
# - AIC
# - What else?
# Comparison G
Publish::publish(results$G$comparison)
# Comparison M
Publish::publish(results$M$comparison)

## Nicer trees
# Better colours
# I don't like the colours
# Titles lack the pointer to the trait name.
# Remove "NA" from legend
# Considering to add alpha besides theta
# Eliminate the numeric scale at the bottom 
# For some reason, the smaller display sinked to RStudio visual device looks better thant the sinked output. No idea why.
# View plots for a specific variable
print(results$G$plots$tree_plot)
print(results$M$plots$tree_plot)

## Regimes
# Summaries already have thetha and logit(theta).
# I think we may also need alpha and other numbers. Maybe checkout what papers normally report.
# Curiously, it seems like ´$parameters´ is empty (return NAs)
# Read surface documentation to figure out what are does values stored at n_regimes: k, kprime, deltak, etc.
# Regimes dist (plot) is nice, but that should just go on a table (?)
print(results$G$summary$regimes)
print(results$M$summary$regimes)
print(as.table(as.matrix(dist(plot))))

# What's next?
# I need to list actual things to get and sort them.
# Resulting number I need
# Tables
# Plots (includin tree(s))
