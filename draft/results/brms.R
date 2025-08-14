library(ape)
library(brms)
# Ensure GM_df has a column 'community' matching tree tip labels
stopifnot("community" %in% names(GM_df))

# Make sure tree tips and data align
Ta <- keep.tip(Ta, intersect(Ta$tip.label, GM_df$community))
Tc <- keep.tip(Tc, intersect(Tc$tip.label, GM_df$community))

d_Ta <- subset(GM_df, community %in% Ta$tip.label)
d_Tc <- subset(GM_df, community %in% Tc$tip.label)

# Order data to match the correlation matrices row order
d_Ta <- d_Ta[match(Ta$tip.label, d_Ta$community), ]
d_Tc <- d_Tc[match(Tc$tip.label, d_Tc$community), ]

# Phylogenetic correlation matrices (DON'T call these 'A'; use K_*)
K_Ta <- ape::vcv(Ta, corr = TRUE)
K_Tc <- ape::vcv(Tc, corr = TRUE)
stopifnot(identical(rownames(K_Ta), d_Ta$community))
stopifnot(identical(rownames(K_Tc), d_Tc$community))

# Scale predictors to help sampling (optional but recommended)
scale_std <- function(x) as.numeric(scale(x))
d_Ta$M_z <- scale_std(d_Ta$M)
d_Ta$A_z <- scale_std(d_Ta$A)
d_Ta$S_z <- scale_std(d_Ta$S)

d_Tc$M_z <- scale_std(d_Tc$M)
d_Tc$A_z <- scale_std(d_Tc$A)
d_Tc$S_z <- scale_std(d_Tc$S)

# Family handles exact 0 and 1: zero_one_inflated_beta()
# Formula: mean part; also model phi (precision) simply with intercept
bform_G <- bf(G ~ M_z + A_z + S_z + (1 | gr(community, cov = K_Ta)),
              phi ~ 1,
              zoi ~ 1, coi ~ 1)   # model zero/one inflation as intercepts

priors_G <- c(
  prior(normal(0, 1), class = "b"),            # slopes
  prior(normal(0, 1.5), class = "Intercept"),  # mean intercept
  prior(exponential(1), class = "sd"),         # phylo RE sd
  prior(exponential(1), class = "phi"),        # precision
  prior(beta(1, 9), class = "zi"),             # weakly favor low zero mass
  prior(beta(1, 9), class = "coi")             # weakly favor low one mass
)

fit_G_Ta <- brm(
  formula = bform_G,
  data = d_Ta,
  family = zero_one_inflated_beta(),
  data2 = list(K_Ta = K_Ta),
  prior = priors_G,
  chains = 4, iter = 4000, cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 123
)

bform_G_Tc <- bf(G ~ M_z + A_z + S_z + (1 | gr(community, cov = K_Tc)),
                 phi ~ 1, zoi ~ 1, coi ~ 1)

fit_G_Tc <- brm(
  bform_G_Tc, data = d_Tc, family = zero_one_inflated_beta(),
  data2 = list(K_Tc = K_Tc),
  prior = priors_G,
  chains = 4, iter = 4000, cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 123
)

library(bayesplot)
summary(fit_G_Ta)
bayes_R2(fit_G_Ta)        # Report as Bayes R^2
pp_check(fit_G_Ta)        # Posterior predictive fit
posterior_interval(fit_G_Ta, prob = 0.95)  # 95% CrIs

# A Bayesian zero–one‑inflated beta regression with a phylogenetic random effect was fitted for G as a function of M, A, and S. The phylogenetic covariance was derived from Ta (Tc in the subset analysis). The effect of M was positive (estimate = XXXX, 95% CrI XXXX–XXXX); effects for A and S were XXXX. Model fit was adequate by posterior predictive checks; Bayes R^2 = XXXX.”