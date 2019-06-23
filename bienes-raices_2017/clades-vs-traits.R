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
    else stop('p-value not found in anova of k=', k)
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

butlast <- function(v) head(v, n=-1)

# pvalues could be computed from commoners using clades_pval
# but it takes too long so we defer the computation to the user
pvalue_tiles <- function(commoners,
                         pvalues,
                         pvalue_method = "kw",
                         hclust_method = hclust_default_method) {
  # fill_var <- if (pvalue_method == "anova") sym("pval_aov")
  # else if (pvalue_method == "kw") sym("pval_kw")
  # else stop("parameter 'pvalue_method' isn't recognised")
  pvalue <- if (pvalue_method == "anova") pvalues$pval_aov
  else if (pvalue_method == "kw") pvalues$pval_kw
  else stop("parameter 'pvalue_method' isn't recognised")
  sur_clust <- surname_clustering(commoners, hclust_method)
  # compute the position and width of the tiles
  inv_height <- rev(1 - sur_clust$height)
  dx <- inv_height - c(0, butlast(inv_height))
  xmin <- inv_height - dx/2
  xmax <- c(butlast(inv_height) + dx[-1]/2, 1)
  stopifnot(xmax == c(xmin[-1], 1))
  width <- xmax - xmin
  # there is nothing between 0 and xmin[1] because
  # there is no p-value for the case when there is only one clade
  x <- xmin + width / 2
  pv <- tibble(pvalue = pvalue, x = x, width = width)
  ggplot(pv) +
    geom_tile(aes(x = x, y = 0, width = width, fill = pvalue)) +
    scale_x_continuous("", breaks = NULL, expand = expand_scale()) +
    scale_y_continuous("", breaks = NULL, expand = expand_scale()) +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    guides(fill = guide_colourbar(title = "p-value"))
}

# TODO: write function similar to clust_pval that computes p-value
# for a pair of groups (t-test?), not all groups at the same time
# then draw the dendrogram with the p-values at the joints.
