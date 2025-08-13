## ================================================================
## Inferred trees — summary + tanglegram + entanglement + tests
## Drop-in: source() this file after your trees/distances exist.
## Requires: ape, dendextend, vegan, treestats (optional)
## ================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dendextend)
  library(vegan)
  suppressWarnings(try(library(treestats), silent = TRUE))
  suppressWarnings(try(library(conflicted), silent = TRUE))
})

## ---------- global settings ----------

## ---------- prefer dendextend where names clash ----------
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dendextend::ladderize)
  conflicted::conflicts_prefer(dendextend::untangle)
  conflicted::conflicts_prefer(dendextend::prune)
  conflicted::conflicts_prefer(dendextend::rotate)
}

## ---------- helpers ----------
tree_height   <- function(tr) max(node.depth.edgelength(tr))
mean_bl       <- function(tr) mean(tr$edge.length)

colless_norm  <- function(tr) {
  if (requireNamespace("treestats", quietly = TRUE)) {
    treestats::colless(tr, normalization = "pda")
  } else {
    suppressWarnings(as.numeric(colless.phylo(tr))) # fallback (may differ)
  }
}

sackin_norm   <- function(tr) {
  if (requireNamespace("treestats", quietly = TRUE)) {
    treestats::sackin(tr, normalization = "pda")
  } else {
    NA_real_  # ape::sackin is unnormalized; better to skip than mix scales
  }
}

# Align two square matrices by common rows/cols (names or dimnames)
align_mats <- function(A, B) {
  if (inherits(A, "dist")) A <- as.matrix(A)
  if (inherits(B, "dist")) B <- as.matrix(B)
  ra <- rownames(A); rb <- rownames(B)
  common <- intersect(ra, rb)
  A2 <- A[common, common, drop = FALSE]
  B2 <- B[common, common, drop = FALSE]
  list(A = A2, B = B2, labs = common)
}

# Cophenetic vs given distance (matrix or dist), aligned by tree tip labels
coph_corr <- function(tr, D) {
  C <- cophenetic(tr)
  if (inherits(D, "dist")) D <- as.matrix(D)
  labs <- intersect(rownames(C), rownames(D))
  cor(C[labs, labs], D[labs, labs])
}

# Mantel wrapper that prints r, p, n (after alignment)
mantel_quick <- function(A, B, perms = iter) {
  AB <- align_mats(A, B)
  m  <- mantel(as.dist(AB$A), as.dist(AB$B), permutations = perms, method = "pearson")
  list(r = as.numeric(m$statistic),
       p = as.numeric(m$signif),
       n = nrow(AB$A) - 2L,
       perms = perms)
}

# ---------- Baker’s gamma with permutation p (only replace THIS function) ----------

bg_gamma_with_p <- function(t1, t2, perms = iter, force_ultra = FALSE, seed = 9) {
  if (force_ultra) {
    if (!is.ultrametric(t1)) t1 <- chronoMPL(t1)
    if (!is.ultrametric(t2)) t2 <- chronoMPL(t2)
  }
  # keep common tips and ladderize
  keep <- intersect(t1$tip.label, t2$tip.label)
  t1k  <- ladderize(reorder(keep.tip(t1, keep), "postorder"))
  t2k  <- ladderize(reorder(keep.tip(t2, keep), "postorder"))

  # to dendrograms
  d1 <- as.dendrogram(t1k)
  d2 <- as.dendrogram(t2k)

  # Align label sets and rotate d2 to d1's order (cosmetic)
  labs <- intersect(labels(d1), labels(d2))
  d1   <- prune(d1, setdiff(labels(d1), labs))
  d2   <- prune(d2, setdiff(labels(d2), labs))
  d2   <- rotate(d2, order = labels(d1))

  # observed gamma
  g_obs <- cor_bakers_gamma(d1, d2)

  # permutation: RANDOMIZE LABELS of d2 (not just rotate order)
  set.seed(seed)
  g_perm <- numeric(perms)
  base_labels <- labels(d2)
  for (i in seq_len(perms)) {
    d2p <- d2
    labels(d2p) <- sample(base_labels, length(base_labels), replace = FALSE)
    g_perm[i] <- cor_bakers_gamma(d1, d2p)
  }

  # two-sided p
  p_two <- (1 + sum(abs(g_perm) >= abs(g_obs))) / (perms + 1)
  list(gamma = unname(g_obs), p = p_two, perms = perms, tips = length(labs))
}
## ---------- 1) BASIC TREE STATS ----------
cat("\n# === TREE STATS ===\n")
trees <- list(Ta = Ta, Ts = Ts, Tg = Tg, Tc = Tc)

tips <- sapply(trees, Ntip)
cat("Tips per tree:\n"); print(tips)

ultra <- sapply(trees, is.ultrametric)
cat("\nUltrametric?\n"); print(ultra)

heights <- sapply(trees, tree_height)
cat("\nTree heights:\n"); print(heights)

mbl <- sapply(trees, mean_bl)
cat("\nMean branch lengths:\n"); print(mbl)

cat("\nColless (normalized) — Ta & Tc:\n")
print(sapply(list(Ta = Ta, Tc = Tc), colless_norm))

cat("\nSackin (normalized) — Ta & Tc:\n")
print(sapply(list(Ta = Ta, Tc = Tc), sackin_norm))

## ---------- 2) COPHENETIC vs SOURCE DISTANCES ----------
cat("\n# === COPHENETIC vs SOURCE DISTANCES ===\n")
cat("Ta vs Da (all surnames):     ", coph_corr(Ta, Da), "\n")
cat("Ts vs Ds (subset surnames):  ", coph_corr(Ts, Ds), "\n")
cat("Tg vs Dg_dps (genetics):     ", coph_corr(Tg, Dg_dps), "\n")

## ---------- 3) MANTEL TESTS ----------
cat("\n# === MANTEL TESTS (Pearson; iter perms) ===\n")
cat("All tips: Da vs Da_geo\n");     print(mantel_quick(Da,     Da_geo, perms = iter))
cat("\nSubset:  Ds vs D_geo\n");      print(mantel_quick(Ds,     D_geo,  perms = iter))
cat("\nSubset:  Dg_dps vs D_geo\n");  print(mantel_quick(Dg_dps, D_geo,  perms = iter))
cat("\nSubset:  Ds vs Dg_dps\n");     print(mantel_quick(Ds,     Dg_dps, perms = iter))

## ---------- 4) TANGLEGRAM & ENTANGLEMENT (Ts vs Tg) ----------
cat("\n# === TANGLEGRAM & ENTANGLEMENT (Ts vs Tg) ===\n")
keep <- intersect(Ts$tip.label, Tg$tip.label)
Ts2  <- ladderize(reorder(keep.tip(Ts, keep), "postorder"))
Tg2  <- ladderize(reorder(keep.tip(Tg, keep), "postorder"))

# numeric labels for clean crossing lines + legend printed below
Ts_num <- Ts2; Tg_num <- Tg2
Ts_num$tip.label <- seq_along(Ts2$tip.label)
Tg_num$tip.label <- seq_along(Tg2$tip.label)

dTs <- as.dendrogram(Ts_num)
dTg <- as.dendrogram(Tg_num)

dl_step1 <- dendlist(dTs, dTg) %>% untangle(method = "step1side")
dl_step2 <- dendlist(dTs, dTg) %>% untangle(method = "step2side")

ent1 <- entanglement(dl_step1)
ent2 <- entanglement(dl_step2)

if (ent2 < ent1) {
  best_dl <- dl_step2; best_ent <- ent2; best_method <- "step2side"
} else {
  best_dl <- dl_step1; best_ent <- ent1; best_method <- "step1side"
}

cat(sprintf("Entanglement (method=%s): %.6f  (lower is better)\n", best_method, best_ent))

cat("\nLegend (number -> community):\n")
lab_map <- paste(seq_along(Ts2$tip.label), Ts2$tip.label, sep = "    ")
cat(paste(lab_map, collapse = "\n"), "\n")

# If you want to plot the tanglegram, uncomment:
# tanglegram(
#   best_dl,
#   common_subtrees_color_lines = FALSE,
#   highlight_distinct_edges    = FALSE,
#   lab.cex = 0.7, margin_inner = 4,
#   color_lines = "grey40", lwd = 1
# )

## ---------- 5) BAKER’S GAMMA (gamma + permutation p) ----------
cat("\n# === BAKER'S GAMMA (gamma + permutation p) ===\n")
keep_all <- Reduce(intersect, list(Ts$tip.label, Tg$tip.label, Tc$tip.label))
Ts3 <- keep.tip(Ts, keep_all)
Tg3 <- keep.tip(Tg, keep_all)
Tc3 <- keep.tip(Tc, keep_all)

bg_Ts_Tg <- bg_gamma_with_p(Ts3, Tg3, perms = iter, force_ultra = FALSE, seed = 9)
bg_Ts_Tc <- bg_gamma_with_p(Ts3, Tc3, perms = iter, force_ultra = FALSE, seed = 9)
bg_Tg_Tc <- bg_gamma_with_p(Tg3, Tc3, perms = iter, force_ultra = FALSE, seed = 9)

set.seed(9)
bg_raw  <- list(
  Ts_Tg = bg_gamma_with_p(Ts3, Tg3, perms = iter, force_ultra = FALSE),
  Ts_Tc = bg_gamma_with_p(Ts3, Tc3, perms = iter, force_ultra = FALSE),
  Tg_Tc = bg_gamma_with_p(Tg3, Tc3, perms = iter, force_ultra = FALSE)
)

set.seed(9)
bg_ultra <- list(
  Ts_Tg = bg_gamma_with_p(Ts3, Tg3, perms = iter, force_ultra = TRUE),
  Ts_Tc = bg_gamma_with_p(Ts3, Tc3, perms = iter, force_ultra = TRUE),
  Tg_Tc = bg_gamma_with_p(Tg3, Tc3, perms = iter, force_ultra = TRUE)
)

print(bg_raw); print(bg_ultra)

cat(sprintf("Ts vs Tg:  gamma = %.3f, p = %.4f, tips = %d, perms = %d\n",
            bg_Ts_Tg$gamma, bg_Ts_Tg$p, bg_Ts_Tg$tips, bg_Ts_Tg$perms))
cat(sprintf("Ts vs Tc:  gamma = %.3f, p = %.4f, tips = %d, perms = %d\n",
            bg_Ts_Tc$gamma, bg_Ts_Tc$p, bg_Ts_Tc$tips, bg_Ts_Tc$perms))
cat(sprintf("Tg vs Tc:  gamma = %.3f, p = %.4f, tips = %d, perms = %d\n",
            bg_Tg_Tc$gamma, bg_Tg_Tc$p, bg_Tg_Tc$tips, bg_Tg_Tc$perms))

## ================================================================
