library(tidyverse)
library(Biodem)

surname_distance_matrix <- function(commoners,
                                    group_by_col="community") {
  # cross tabulate
  surnames_freq <- table(commoners$surname_father,
                         commoners[[group_by_col]])
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # hedrick returns values of similarity
  # transform them into values of dissimilarity (distance)
  as.dist(1-hedkin)
}

hclust_default_method <- "complete"
surname_clustering <- function(commoners,
                               hclust_method=hclust_default_method,
                               group_by_col="community") {
  commoners %>%
    surname_distance_matrix(group_by_col) %>%
    # hierarchical clustering of the distance matrix
    hclust(method=hclust_method)
}
