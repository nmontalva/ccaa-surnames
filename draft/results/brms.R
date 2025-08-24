
library(ape)
library(brms)
library(bayesplot)

## ---- Align data with tree tips ----
stopifnot("community" %in% names(result_traits))

Ta <- keep.tip(Ta, intersect(Ta$tip.label, result_traits$community))
Tc <- keep.tip(Tc, intersect(Tc$tip.label, result_traits$community))

d_Ta <- subset(result_traits, community %in% Ta$tip.label)
d_Tc <- subset(result_traits, community %in% Tc$tip.label)

# Order rows to match covariance matrices
d_Ta <- d_Ta[match(Ta$tip.label, d_Ta$community), ]
d_Tc <- d_Tc[match(Tc$tip.label, d_Tc$community), ]

## ---- Phylogenetic correlation matrices (use K_*; don't clash with A=agnatic) ----
K_Ta <- ape::vcv(Ta, corr = TRUE)
K_Tc <- ape::vcv(Tc, corr = TRUE)
stopifnot(identical(rownames(K_Ta), d_Ta$community))
stopifnot(identical(rownames(K_Tc), d_Tc$community))

## ---- Scale predictors (helps sampling; response stays in [0,1]) ----
scale_std <- function(x) as.numeric(scale(x))
for (nm in c("M","A","S")) {
  d_Ta[[paste0(nm,"_z")]] <- scale_std(d_Ta[[nm]])
  d_Tc[[paste0(nm,"_z")]] <- scale_std(d_Tc[[nm]])
}

## ---- ziB family handles exact 0 and 1; model zi/coi as intercepts ----
bform_G_Ta <- bf(
  G ~ M_z + A_z + S_z + (1 | gr(community, cov = K_Ta)),
  phi ~ 1,     # precision (log link)
  zi ~ 1     # zero-or-one inflation (logit)
)

# Example weakly-informative choices:

logit <- function(p) qlogis(p)
priors_G <- c(
  prior(normal(0, 1.5), class = "Intercept"),              # mu intercept
  prior(normal(0, 1),   class = "b"),                      # mu slopes
  prior(exponential(2), class = "sd"),                     # phylo RE sd (regularizing)
  prior(normal(log(5), 1),        class = "Intercept", dpar = "phi"),
  # Option A: generic weak prior near 10% zeros on logit scale:
  prior(normal(logit(0.10), 1.5), class = "Intercept", dpar = "zi")
)

# Option B (data-guided zi prior): replace the last line with this per dataset:
# prior(normal(logit((sum(d_Ta$G==0)+0.5)/(nrow(d_Ta)+1)), 1.0), class="Intercept", dpar="zi")
# prior(normal(logit((sum(d_Tc$G==0)+0.5)/(nrow(d_Tc)+1)), 1.0), class="Intercept", dpar="zi")

#Families:
fam_Ta <- zero_inflated_beta()
fam_Tc <- zero_inflated_beta()   # no ones → no 'coi' term


ctrl <- list(adapt_delta = 0.999, max_treedepth = 15)

fit_G_Ta <- brm(
  bform_G_Ta, data = d_Ta, family = fam_Ta,
  data2 = list(K_Ta = K_Ta), prior = priors_G,
  chains = 4, iter = 6000, warmup = 2000, cores = 4,
  control = ctrl, init_r = 0.1, seed = 123
)

## ---- Same model on Tc ----
bform_G_Tc <- bf(
  G ~ M_z + A_z + S_z + (1 | gr(community, cov = K_Tc)),
  phi ~ 1, zi ~ 1
)


fit_G_Tc <- brm(
  bform_G_Tc, data = d_Tc, family = fam_Tc,
  data2 = list(K_Tc = K_Tc), prior = priors_G,
  chains = 4, iter = 8000, warmup = 3000, cores = 4,
  control = ctrl, init_r = 0.1, seed = 124
)

#Cheks
get_prior(bform_G_Ta, d_Ta, fam_Ta, data2 = list(K_Ta = K_Ta))
summary(fit_G_Ta); bayes_R2(fit_G_Ta); pp_check(fit_G_Ta)

get_prior(bform_G_Tc, d_Tc, fam_Tc, data2 = list(K_Tc = K_Tc))
summary(fit_G_Tc); bayes_R2(fit_G_Tc); pp_check(fit_G_Tc)

## ---- Diagnostics and reporting (same for both fits) ----
summary(fit_G_Ta)
bayes_R2(fit_G_Ta)                               # Bayes R^2
posterior_interval(fit_G_Ta, prob = 0.95)        # 95% CrIs (all params)
pp_check(fit_G_Ta)                               # posterior predictive checks
conditional_effects(fit_G_Ta)                    # marginal effects (mu)

summary(fit_G_Tc)
# bayes_R2(fit_G_Tc)
posterior_interval(fit_G_Tc, prob = 0.95)
pp_check(fit_G_Tc)
conditional_effects(fit_G_Tc)

#To remeber for methods:
#We accounted for phylogenetic non-independence among communities by including a random intercept with a covariance matrix derived from the phylogenetic tree. This approach models the correlation in residuals expected under shared ancestry, analogous to phylogenetic generalized least squares (PGLS), but in a hierarchical Bayesian framework. By using the full phylogenetic covariance as a random effect, we capture the entire phylogenetic signal without reducing the tree to a set of fixed predictors, thus avoiding arbitrary dimensionality reduction and preserving interpretability of variance components.