library(tidyverse)
library(ade4) # for mantel.randtest
# library(ggpubr) # for ggboxplot

# compare surnames distance matrix with S, A, R or G matrices
compare_dist <- function(commoners,
                         trait,
                         dist_method="euclid",
                         nrepet=99999) {
  ts <- traits(commoners)
  trait_dist <- dist(ts[[trait]], dist_method)
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
                        nbins=10,
                        save_as=NULL,
                        hclust_method=hclust_default_method) {
  tc <- traits(commoners) %>%
    left_join(k_clusters(commoners, k, hclust_method),
              by="community")
  if (!is.null(save_as))
    pdf(save_as, width=5, height=4)
  p <- ggplot(data=tc) +
    geom_histogram(aes(x=G, y=..density.., fill=branch),
                   bins=nbins, alpha=1/k,
                   position="identity", boundary=0)
  if (!is.null(save_as)) {
    print(p)
    dev.off()
  } else p
}

# anova of a trait by clusters
clusters_pval <- function(commoners,
                          trait,
                          ks=seq(2,169), # vector with the number of groups in which the tree is split
                          hclust_method=hclust_default_method) {
                          # hclust_method="ward.D2") {
  ts <- traits(commoners)
  ks %>% map_dfr(function(k) {
    tks <- left_join(ts, k_clusters(commoners, k, hclust_method),
                     by="community")
    clusters <- tks %>% group_by(branch) %>%
      summarise(size=n())
    # TODO: make code parametric on trait
    # gini_mean=mean(G, na.rm=TRUE),
    # gini_sd=sd(G, na.rm=TRUE))
    tibble(
      num_clusters=k,
      # plot=ggboxplot(tck, x="branch", y=trait, xlab="Cluster", ylab=trait),
      # anova=summary(aov(tks[[trait]] ~ tks$branch)),
      pval_aov=summary(aov(tks[[trait]] ~ tks$branch))[[1]][["Pr(>F)"]][[1]],
      pval_kw=kruskal.test(tks[[trait]] ~ as.factor(tks$branch))$p.value,
      n=mean(clusters$size),
      sd=sd(clusters$size))
  })
}

plot_pvalues_kw <- function(pvalues) {
  ggplot(pvalues, aes(num_clusters, pval_kw)) +
    geom_line() +
    scale_x_continuous("Number of clusters") +
    scale_y_continuous("Kruskal-Wallis p-value")
}

plot_pvalues_anova <- function(pvalues) {
  ggplot(pvalues, aes(num_clusters, pval_aov)) +
    geom_line() +
    scale_x_continuous("Number of clusters") +
    scale_y_continuous("ANOVA p-value")
}

colourbar_pvalues <- function(pvalues) {
  ggplot(pvalues, aes(1, 1, fill=pval_aov)) +
    geom_bar(stat="identity") +
    scale_x_discrete("", expand=expand_scale()) +
    scale_y_continuous("Number of clusters", expand=expand_scale()) +
    scale_fill_viridis_c(option = "plasma", direction=-1) +
    guides(fill = guide_colourbar(title="p-value",
                                  barwidth=0.5, barheight=10)) +
    coord_flip()
}

# TODO: write function similar to clust_pval that computes p-value
# for a pair of groups (t-test?), not all groups at the same time
# then draw the dendrogram with the p-values at the joints.
