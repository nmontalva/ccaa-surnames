library(tidyverse)

k_clades <- function(commoners,
                     k,
                     hclust_method=hclust_default_method) {
  sur_clust <- surname_clustering(commoners, hclust_method)
  kgroups <- cutree(sur_clust, k=k)
  tibble(community=names(kgroups),
         branch=kgroups)
}

k_histogram <- function(commoners,
                        k,
                        nbins=10,
                        save_as=NULL,
                        hclust_method=hclust_default_method) {
  tc <- traits(commoners) %>%
    left_join(k_clades(commoners, k, hclust_method),
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

# anova of a trait by clades
# ks is a vector with the number of groups in which the tree is split
clades_pval <- function(commoners,
                        trait,
                        ks=seq(2,169),
                        hclust_method=hclust_default_method) {
                        # hclust_method="ward.D2") {
  ts <- traits(commoners)
  ks %>% map_dfr(function(k) {
    tks <- left_join(ts, k_clades(commoners, k, hclust_method),
                     by="community")
    # clades <- tks %>% group_by(branch) %>%
    #   summarise(size=n(),
    #             trait_mean=mean(!!sym(trait), na.rm=TRUE),
    #             trait_sd=sd(!!sym(trait), na.rm=TRUE))
    # size_mean <- mean(clades$size)
    # size_sd <- sd(clades$size)
    anova <- summary(aov(tks[[trait]] ~ tks$branch))[[1]]
    pval_aov <- if ("Pr(>F)" %in% colnames(anova))
      anova[["Pr(>F)"]][[1]]
    else 0
    # kw <- kruskal.test(tks[[trait]] ~ as.factor(tks$branch))
    kw <- kruskal.test(tks[[trait]] ~ tks$branch)
    tibble(
      num_clades=k,
      pval_aov=pval_aov,
      pval_kw=kw$p.value)
  })
}

plot_pvalues_kw <- function(pvalues) {
  ggplot(pvalues, aes(num_clades, pval_kw)) +
    geom_line() +
    scale_x_continuous("Number of clades") +
    scale_y_continuous("Kruskal-Wallis p-value")
}

plot_pvalues_anova <- function(pvalues) {
  ggplot(pvalues, aes(num_clades, pval_aov)) +
    geom_line() +
    scale_x_continuous("Number of clades") +
    scale_y_continuous("ANOVA p-value")
}

colourbar_pvalues <- function(pvalues) {
  ggplot(pvalues, aes(1, 1, fill=pval_aov)) +
    geom_bar(stat="identity") +
    scale_x_discrete("", expand=expand_scale()) +
    scale_y_continuous("Number of clades", expand=expand_scale()) +
    scale_fill_viridis_c(option = "plasma", direction=-1) +
    guides(fill = guide_colourbar(title = "p-value")) +
    coord_flip()
}

# TODO: write function similar to clust_pval that computes p-value
# for a pair of groups (t-test?), not all groups at the same time
# then draw the dendrogram with the p-values at the joints.
