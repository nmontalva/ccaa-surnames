library(tidyverse)
library(ade4) # for mantel.randtest

# compare surnames distance matrix with S, A, R or G matrices
compare_dist <- function(commoners,
                         trait,
                         dist_method="euclid",
                         nrepet=99999) {
  tc <- traits(commoners)
  trait_dist <- dist(tc[[trait]], dist_method)
  hedkin_dist <- surname_distance_matrix(commoners)
  # Ric: I don't understand why one would use quasieuclid or not
  mantel.randtest(hedkin_dist, # quasieuclid(hedkin_dist),
                  trait_dist,
                  nrepet=nrepet)
}


k_clusters <- function(commoners,
                       k,
                       hclust_method=hclust_default_method) {
  sur_clust <- surname_clustering(commoners, hclust_method)
  kgroups <- cutree(sur_clust, k=k)
  tibble(community=names(kgroups),
         branch=as.character(kgroups))
}

k_histogram <- function(commoners,
                        k,
                        save_as=NULL,
                        hclust_method=hclust_default_method) {
  tc <- traits(commoners) %>%
    left_join(k_clusters(commoners, k, hclust_method),
              by="community")
  if (!is.null(save_as))
    pdf(save_as, width=5, height=4)
  p <- ggplot(data=tc) +
    geom_histogram(aes(x=G, y=..density.., fill=branch),
                   bins=10, alpha=1/k,
                   position="identity", boundary=0)
  if (!is.null(save_as)) {
    print(p)
    dev.off()
  } else p
}

# anova of a trait by clusters
clust_pval <- function(commoners,
                       ks,
                       trait="G",
                       # hclust_method=hclust_default_method) {
                       hclust_method="ward.D2") {
  tc <- traits(commoners)
  ks %>% map_dfr(function(k) {
    tck <- left_join(tc, k_clusters(commoners, k, hclust_method),
                     by="community")
    clusters <- tck %>% group_by(branch) %>%
      summarise(size=n())
    # TODO: make code parametric on trait
    # gini_mean=mean(G, na.rm=TRUE),
    # gini_sd=sd(G, na.rm=TRUE))
    tibble(
      num_clusters=k,
      # plot=ggboxplot(tck, x="branch", y=trait, xlab="Cluster", ylab=trait),
      # anova=summary(aov(tck[[trait]] ~ tck$branch)),
      pval=kruskal.test(tck[[trait]] ~ as.factor(tck$branch))$p.value,
      n=mean(clusters$size),
      sd=sd(clusters$size))
  })
}

# TODO: write function similar to clust_pval that computes p-value
# for a pair of groups (t-test?), not all groups at the same time
# then draw the dendrogram with the p-values at the joints.
