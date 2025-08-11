# Objective: Describe the resulting trees and correlations between them

# Names

# Ta  <-  surnames tree (all tips, UPGMA, ultrametric)
# Ts  <-  surnames tree (16 tips subset, UPGMA)
# Tg  <-  genetic tree (16 tips, UPGMA on Dps)  # final choice per Methods
# Tc  <-  consensus tree (15/16 tips — verify)
# D_s_all <- 
# D_s_16 <-   # surname distances (1 - Hedrick’s H)
# D_ps <-            # genetic distances (pick your chosen, e.g., D_ps)
# D_geo_all
# D_geo_16 <-   # geography distances (km)

# Tips
sapply(list(Ta=Ta, Ts=Ts, Tg=Tg, Tc=Tc), Ntip)
Reduce(intersect, list(Ts$tip.label, Tg$tip.label, Tc$tip.label))  # common set
setdiff(Ts$tip.label, Tg$tip.label)  # just to be sure
 
# Shape and ultrametricity
sapply(list(Ta=Ta, Ts=Ts, Tg=Tg, Tc=Tc), is.ultrametric)

# Heights (tree depth) & mean branch length
tree_height <- function(tr) max(node.depth.edgelength(tr))
mean_bl     <- function(tr) mean(tr$edge.length)

sapply(list(Ta=Ta, Ts=Ts, Tg=Tg, Tc=Tc), tree_height)
sapply(list(Ta=Ta, Ts=Ts, Tg=Tg, Tc=Tc), mean_bl)

# Imbalance

colless_norm <- function(tr) colless.phylo(tr, norm=TRUE)
sackin_norm  <- function(tr) sackin.phylo(tr, norm=TRUE)

sapply(list(Ta=Ta, Tc=Tc), colless_norm)
sapply(list(Ta=Ta, Tc=Tc), sackin_norm)

# Cophenetic correlation (tree vs the distance it was built from)

 coph_corr <- function(tr, D) {
  labs <- tr$tip.label
  C <- cophenetic(tr)[labs, labs]
  cor(as.numeric(C[lower.tri(C)]), as.numeric(D[labs, labs][lower.tri(D[labs, labs])]))
}

# Ta vs surname distances (all)
coph_corr(Ta, D_s_all)

# Ts vs surname distances (subset)
coph_corr(Ts, D_s_16)

# Tg vs chosen genetic distance (subset; e.g., D_ps)
coph_corr(Tg, D_ps)

# Tc: just compute the cophenetic correlation (no native “source” D); you can still report this number
cor(as.numeric(cophenetic(Tc)[lower.tri(cophenetic(Tc))]),
    as.numeric(cophenetic(Tc)[lower.tri(cophenetic(Tc))]))  # will be 1 by definition; skip

# Mantel tests (matrix ↔ matrix; report r, p, n, perms)
mantel_quick <- function(A,B,perms=iter){
  m <- mantel(as.dist(A), as.dist(B), permutations=perms, method="pearson")
  list(r=m$statistic, p=m$signif, n=attr(as.dist(A),"Size"), perms=perms)
}

# All
mantel_quick(D_s_all, D_geo_all)

# 16‑tip subset
mantel_quick(D_s_16, D_geo_16)
mantel_quick(D_ps,    D_geo_16)
mantel_quick(D_s_16,  D_ps)

# (Optional) Tc cophenetic vs Geo
C_tc <- cophenetic(Tc)
mantel_quick(C_tc, D_geo_16[rownames(C_tc), colnames(C_tc)])

# Baker’s Gamma (tree ↔ tree; 16‑tip shared set only)
# Make sure tip sets match exactly
keep <- Reduce(intersect, list(Ts$tip.label, Tg$tip.label, Tc$tip.label))
Ts2 <- keep.tip(Ts, keep); Tg2 <- keep.tip(Tg, keep); Tc2 <- keep.tip(Tc, keep)

bg <- function(t1,t2, nperm=iter){
  g  <- cor_bakers_gamma(as.hclust(t1), as.hclust(t2))
  gp <- cor_bakers_gamma(as.hclust(t1), as.hclust(t2), nperm=nperm)
  list(gamma=g, p=gp$p.value)
}

bg(Ts2, Tg2)  # report γ and p
bg(Ts2, Tc2)
bg(Tg2, Tc2)

# (Optional) Tanglegram entanglement (you have ~0.12; recompute to confirm)
t1 <- as.dendrogram(as.hclust(Ts2)); t2 <- as.dendrogram(as.hclust(Tg2))
entanglement(dendlist(t1, t2))
