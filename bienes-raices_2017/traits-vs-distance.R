library(tidyverse)
library(ade4) # for mantel.randtest

# compare surnames distance matrix with S, A, R or G delta matrices
mantel <- function(commoners,
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

