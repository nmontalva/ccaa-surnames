library(tidyverse)
library(ade4) # for mantel.randtest
library(latex2exp)

dist_vs_delta <- function(commoners,
                          trait,
                          dist_method="euclid") {
  ts <- traits(commoners)
  trait_deltas <- dist(ts[[trait]], dist_method)
  hedkin_dist <- surname_distance_matrix(commoners)
  # assume trait_dist and hedkin_dist are indexed equally
  stopifnot(all(ts$community == labels(hedkin_dist)))
  dd <- tibble(surname_dist = hedkin_dist,
               trait_dist = trait_deltas)
  m <- lm(trait_dist ~ surname_dist, dd) # y ~ x
  # r2 <- paste("r^2 = ", format(summary(m)$r.squared, digits = 3))
  r2 <- substitute(italic(r)^2~"="~r2, list(
    r2 = format(summary(m)$r.squared, digits = 3)))
  ggplot(dd, aes(x = surname_dist,
                 y = trait_dist)) +
    geom_point(size = 0.1, alpha = 0.25) +
    geom_smooth(method = 'lm', se = FALSE, size = 0.5) +
    annotate("text", x = 0.103, y = 0.55, label = deparse(r2),
             size = 3, colour = "#3366FF", parse = TRUE) +
    # geom_text(x = 0.15, y = 0.58, label = r2, parse = TRUE) +
    scale_x_continuous("Surname distance") +
    # scale_y_continuous(paste("Δ", trait))
    # scale_y_continuous(paste("Delta ", trait))
    scale_y_continuous(TeX(paste("$\\Delta \\,", trait, "$")))
}

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

