#TESTING STEPS OUTSIDE A SCRIPT

# Step 1. Single-fit models
# 1. Logit transform for M
n       <- nrow(GM_df)
eps     <- 0.5 / n
GM_df$M_adj   <- (GM_df$M * (n - 1) + 0.5) / n
GM_df$M_logit <- log(GM_df$M_adj / (1 - GM_df$M_adj))

# 2. Named vector
M_logit_vec <- setNames(GM_df$M_logit, GM_df$community)

# 3. Fit BM and OU
library(geiger)
bm.M <- fitContinuous(y_total, M_logit_vec, model="BM")
ou.M <- fitContinuous(y_total, M_logit_vec, model="OU")

# 4. Print AICs
cat("BM AIC:", bm.M$opt$aic, "\n")
cat("OU AIC:", ou.M$opt$aic, "\n")

# Step 2: SURFACE forward for M_logit

library(surface)

# a) Ensure unique node names
named_tree <- nameNodes(y_total)

# b) Prepare data.frame of M_logit
datM <- data.frame(M_logit = GM_df$M_logit,
                   row.names = GM_df$community)

# c) Convert to OUCH objects
convM  <- convertTreeData(named_tree, datM)
otreeM <- convM[[1]]
odataM <- convM[[2]]

# d) Forward phase (verbose)
fwd_M <- surfaceForward(
  otreeM, odataM,
  aic_threshold = 0,
  verbose       = FALSE
)

fwd_M
surfaceSummary(fwd_M)

# Step 3: SURFACE backward for M_logit

# Assuming you already have fwd_M and odataM/otreeM from the forward phase:
bwd_M <- surfaceBackward(
  otreeM, 
  odataM, 
  fwd_M[[ length(fwd_M) ]], 
  aic_threshold = 0, 
  verbose       = TRUE
)

# Then get its summary:
ss_M <- surfaceSummary(bwd_M)
print(ss_M)

# Step 4. Summary
# 1. Parameter table (logit scale)
param_tab_M <- data.frame(
  alpha    = ss_M$alpha["M_logit"],
  sigma2   = ss_M$sigma_squared["M_logit"]
)
param_tab_M$stat_var <- param_tab_M$sigma2 / (2 * param_tab_M$alpha)
rownames(param_tab_M) <- "M_logit"

cat("\nModel parameters (logit scale):\n")
print(param_tab_M)

# 2. Regime optima table (logit + original M scale)
theta_logit_M <- ss_M$theta[,1]
theta_tab_M <- data.frame(
  regime      = names(theta_logit_M),
  theta_logit = theta_logit_M,
  theta_M     = round(plogis(theta_logit_M), 4)
)
cat("\nRegime optima:\n")
print(theta_tab_M)

# Step 5: Tip assignments for M

# Extract the hansentree
hmod_M   <- bwd_M[[length(bwd_M)]]$fit$M_logit

# regimes factor (length = #nodes), tips are the first 170 entries
regs_fac <- hmod_M@regimes$regs

# Pull out tip regimes
tip_regs_M <- as.character(regs_fac[1:length(GM_df$community)])

# Name them by community
names(tip_regs_M) <- GM_df$community

# Build a data.frame
tip_df_M <- data.frame(
  community = names(tip_regs_M),
  regime    = tip_regs_M,
  stringsAsFactors = FALSE
)

# Show counts and percentages
library(dplyr)
reg_summary_M <- tip_df_M %>%
  count(regime) %>%
  mutate(pct = round(100 * n / sum(n), 1))

print(reg_summary_M)

# Step 6: Bar chart of communities per regime for M

# Make a manual palette matching theta_tab_M
pal_M <- setNames(rainbow(nrow(theta_tab_M)), theta_tab_M$regime)

# Plot
library(ggplot2)
ggplot(reg_summary_M, aes(x = regime, y = n, fill = regime)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = pal_M) +
  labs(
    title = "Number of Communities per Regime (M)",
    x     = "Regime",
    y     = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 14),
    plot.title   = element_text(size = 16, face = "bold")
  )

surfaceTreePlot(named_tree,reg_summary_M)

# Step 7: Tree plot with regime colours and θ legend for M

# 1. Build palette with θ labels
theta_lab_M <- paste0(
  theta_tab_M$regime, 
  ": θ=", 
  theta_tab_M$theta_M
)
pal_M <- setNames(rainbow(nrow(theta_tab_M)), theta_lab_M)

# 2. Create a named vector of tip regimes (with θ‐label names)
tip_regs_lab <- factor(
  paste0(tip_df_M$regime, ": θ=", 
         theta_tab_M$theta_M[match(tip_df_M$regime, theta_tab_M$regime)]),
  levels = theta_lab_M
)

## Alternative tree No 1 ##
# 3. Plot
library(phytools)
plotTree(y_total, fsize=0.6, lwd=1,
         main = "SURFACE Regimes on Tree (M)")
tiplabels(pch=21, bg = pal_M[as.character(tip_regs_lab)], cex=0.5)

# 4. Legend
legend(
  "topright", 
  legend = theta_lab_M, 
  pt.bg  = pal_M, 
  pch    = 21, 
  pt.cex = 1, 
  bty    = "n",
  title  = "Regime (θ)"
)

# The tree above: the branches aren't coloured and the names of the communities are to cluttered

## Alternative tree No 2 ##
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)

# Prepare a data frame mapping tips → regimes (use the theta-label palette if you like)
tip_info <- data.frame(
  label  = tip_df_M$community,
  regime = tip_df_M$regime,
  stringsAsFactors = FALSE
)

# Build a ggtree plot with coloured tip points, no tip labels
p <- ggtree(y_total, size=0.3) %<+% tip_info +
  geom_tippoint(aes(color = regime), size = 1.5) +
  scale_color_manual(
    values = setNames(
      rainbow(nrow(theta_tab_M)), 
      theta_tab_M$regime
    ),
    name = "Regime"
  ) +
  theme_tree2() +
  ggtitle("Regimes on Tree (M)") +
  theme(
    legend.position = "right",
    plot.title      = element_text(size=14, face="bold")
  )

print(p)
# does not show colours *on the branches* And does not have theta values any longer.
# What does OK: Omit the messy tip text, Color each tip’s point by its regime, Show a clean legend (with regime names) on the right.

## Alternative tree No 3 ##
library(phytools)

# 1. Build a named vector of tip regimes (no θ in names)
tip_regs_simple <- tip_df_M$regime
names(tip_regs_simple) <- tip_df_M$community

# 2. Build a palette keyed by regime only
pal_simple <- setNames(
  rainbow(nrow(theta_tab_M)), 
  theta_tab_M$regime
)

# 3. Color the branches
bst <- colorBranches(
  y_total, 
  tipRegs = tip_regs_simple, 
  col = pal_simple
)

# 4. Plot the tree (branches colored)
plot(bst, 
     show.tip.label = FALSE, 
     main = "Branches by SURFACE Regime (M)")

# 5. Add tip points (optional, to mark tips)
tiplabels(pch = 21, 
          bg  = pal_simple[tip_regs_simple],
          cex = 0.5)

# 6. Legend with θ values
legend(
  "topright",
  legend = paste0(theta_tab_M$regime, ": θ=", theta_tab_M$theta_M),
  fill   = pal_simple,
  border = NA,
  bty    = "n",
  title  = "Regime (θ)"
)

#The above tree does not draw and returns and error:
# Error en colorBranches(y_total, tipRegs = tip_regs_simple, col = pal_simple): 
#   no se pudo encontrar la función "colorBranches"

## Alternative tree No 4 ##
# One‐step tree plotting with surfaceTreePlot

# 1. Define your color vector keyed by regime
pal_simple <- setNames(
  rainbow(nrow(theta_tab_M)), 
  theta_tab_M$regime
)

# 2. Plot with no tip‐label text
surfaceTreePlot(
  named_tree,
  bwd_M[[ length(bwd_M) ]],
  cols        = pal_simple,
  cex         = 0.00001,       # tip labels hidden
  labelshifts = FALSE    # no shift labels
)

# 3. Legend showing θ
legend(
  "topright",
  legend = paste0(theta_tab_M$regime, ": θ=", theta_tab_M$theta_M),
  fill   = pal_simple,
  border = NA,
  bty    = "n",
  title  = "Regime (θ)"
)

# Almost there!! But cex = 0 is invalid *and* the legend overlaps the plot.

## Alternative tree No 5 ##
# 1. Expand right margin to make space for legend
op <- par(mar = c(4, 1, 4, 10))  # bottom, left, top, right

# 2. Tiny tip labels (effectively hidden) and no shift labels
surfaceTreePlot(
  named_tree,
  bwd_M[[length(bwd_M)]],
  cols        = pal_simple,
  cex         = 0.01,    # almost zero
  labelshifts = FALSE
)

# 3. Legend positioned in the outer right margin
legend(
  x       = 1.02,       # 100% + a bit on the plotting region (in npc)
  y       = 1,          # top
  legend  = paste0(theta_tab_M$regime, ": θ=", theta_tab_M$theta_M),
  fill    = pal_simple,
  border  = NA,
  bty     = "n",
  xpd     = TRUE,       # allow drawing in the margin
  title   = "Regime (θ)"
)

# Restore original par
par(op)

#Now the lenged is not showing. And the tiny labels are still there.

## Alternative tree No 6 ##
library(ggtree)
library(ggplot2)

# 1. Create a data frame with tip & internal node regimes.
#    We’ll assign each internal node the majority regime of its descendants.
tip_df_M$node <- 1:length(tip_df_M$community)  # temporary placeholder

# Use ggtree’s groupClade helper to propagate regimes upward.
p <- ggtree(named_tree) %<+% tip_df_M +
  aes(color = regime, subset = (isTip)) +
  geom_tippoint(size=1.5) +
  geom_tree(aes(color = regime), size=0.5) +
  scale_color_manual(
    values = setNames(pal_simple, names(pal_simple)),
    name   = "Regime (θ)",
    labels = paste0(theta_tab_M$regime, ": θ=", theta_tab_M$theta_M)
  ) +
  theme_tree2() +
  theme(
    legend.position = c(1.05, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(5, 100, 5, 5)  # make room on the right
  ) +
  ggtitle("SURFACE Regimes on Tree (M)")

print(p)

## Alternative tree No 7 ##

# 1. Prepare tip‐to‐regime data
tip_info <- data.frame(
  label  = tip_df_M$community,
  regime = tip_df_M$regime,
  stringsAsFactors = FALSE
)

# 2. Build rounded θ labels
theta_lab <- setNames(
  paste0(theta_tab_M$regime, ": θ=", format(round(theta_tab_M$theta_M, 3), nsmall=3)),
  theta_tab_M$regime
)

# 3. Palette keyed by regime
pal <- setNames(rainbow(length(theta_lab)), names(theta_lab))

# 4. Plot: only colored branches (no points), fatter lines
p <- ggtree(named_tree, size = 0.8) %<+% tip_info +
  geom_tree(aes(color = regime), size = 0.8) +
  scale_color_manual(
    values = pal,
    labels = theta_lab,
    name   = "Regime (θ)"
  ) +
  theme_tree2() +
  theme(
    plot.title            = element_text(size=14, face="bold"),
    legend.position.inside = c(1.02, 0.5),  # push legend into right margin
    legend.justification = c(0, 0.5),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_text(size=12),
    legend.text          = element_text(size=10)
  ) +
  ggtitle("SURFACE Regimes on Tree (M)")

print(p)

#Drop NA from legend

# Remove any NA entries
keep <- !is.na(names(pal))
pal_clean      <- pal[keep]
theta_lab_clean <- theta_lab[keep]

# After building pal_clean and theta_lab_clean as before:

p <- ggtree(named_tree, size = 0.8) %<+% tip_info +
  geom_tree(aes(color = regime), size = 0.8) +
  scale_color_manual(
    values        = pal_clean,
    labels        = theta_lab_clean,
    name          = "Regime (θ)",
    na.translate  = FALSE    # <-- drop NA from the legend
  ) +
  theme_tree2() +
  theme(
    plot.title             = element_text(size=14, face="bold"),
    legend.position.inside = c(1.02, 0.5),
    legend.justification   = c(0, 0.5),
    legend.background      = element_blank(),
    legend.key             = element_blank(),
    legend.title           = element_text(size=12),
    legend.text            = element_text(size=10)
  ) +
  ggtitle("SURFACE Regimes on Tree (M)")

print(p)