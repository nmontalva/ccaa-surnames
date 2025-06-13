################################################################################
################################################################################
######### Project 111160402: Cultural phylogenetics and coevolution of wealth###
###### inheritance and land tenure norms in agropastoralist communities ########
################################################################################
################################################################################

# OBJETIVO 5
# Fit models of trait evolution (e.g., Brownian motion, Ornstein-Uhlenbeck)
# to G and M. Compare AIC for model selection and perform simulations.

# Cargar las librerías
library(dplyr)
library(geiger)
library(phytools)
library(picante)
library(surface)
library(ape)
library(phangorn)

# Trait data
G_values <- GM_logit %>%
  dplyr::select(G_logit)
M_values <- GM_logit %>%
  dplyr::select(M_logit)

# Randomise tips (with their respective traits, G and M) keeping the underlying tree topology.
# Configuración de los datos
trait_data_G <- setNames(G_values$G_logit, row.names(G_values))
trait_data_M <- setNames(M_values$M_logit, row.names(M_values))
tree_tips <- consensus_tree$tip.label
total_tips <- y_total$tip.label
# Randomize tip labels
## Randomize trees
howmanytrees(length(tree_tips)) #número de árboles filogenéticos posibles

# Etiquetas aleatorizadas
#consenso
set.seed(150)
random.trees <- lapply(1:1000, function(x) {
  tree_random <- consensus_tree  # Copia el árbol original
  tree_random$tip.label <- sample(tree_tips)  # Reordena las etiquetas
  tree_random
})
class(random.trees) <- "multiPhylo"

#total
set.seed(150)
random.total <- lapply(1:1000, function(x) {
  tree_random <- y_total  # Copia el árbol original
  tree_random$tip.label <- sample(total_tips)  # Reordena las etiquetas
  tree_random
})
class(random.total) <- "multiPhylo"

# Compare phylogenetic signal of real tree vs randomised tree.
# Asegurar correspondencia entre rasgos y puntas
phy_tree_s <- consensus_tree
phy_tree <- y_total
trait_data_G_s <- trait_data_G[names(trait_data_G) %in% tree_tips]
trait_data_M_s <- trait_data_M[names(trait_data_M) %in% tree_tips]

##Random.trees
##Todas las comunidades
# Señal filogenética observada (Pagel lambda)
observed_G_lambda <- phylosig(phy_tree, trait_data_G, method = "lambda", test = TRUE, nsim = 9999)
observed_M_lambda <- phylosig(phy_tree, trait_data_M, method = "lambda", test = TRUE, nsim = 9999)
# Señal filogenética observada (Blomberg K)
observed_G_K <- phylosig(phy_tree, trait_data_G, method = "K", test = TRUE, nsim = 9999)
observed_M_K <- phylosig(phy_tree, trait_data_M, method = "K", test = TRUE, nsim = 9999)

#G
random_lambdas_G <- sapply(random.total, function(tree) {
  phylosig(tree, trait_data_G, method = "lambda", test = FALSE)
})
random_K_G <- sapply(random.total, function(tree) {
  phylosig(tree, trait_data_G, method = "K", test = FALSE)
})
lambda_values_G <- as.numeric(random_lambdas_G["lambda",])
p_value_G_lambda <- mean(lambda_values_G>= observed_G_lambda$lambda)
cat("P-value random λ para G:", p_value_G_lambda)
cat("P-value de λ observado =", round(observed_G_lambda$P, 3))

#K_values_G <- as.numeric(random_K_G ["K",])
p_value_G_K <- mean(random_K_G  >= observed_G_K$K)
cat("P-value random K para G:", p_value_G_K)
cat("p-valor K observado =", round(observed_G_K$P, 3))

hist(lambda_values_G, main = "Distribución nula de λ para G", xlab = "λ")
abline(v = observed_G_lambda$lambda, col = "red", lwd = 2)
legend("topright", legend = paste("λ observado =", round(observed_G_lambda$lambda, 3)), 
       col = "red", lwd = 2)

hist(random_K_G, main = "Distribución nula de Blomerg K para G", xlab = "K")
abline(v = observed_G_K$K, col = "red", lwd = 2)
legend("topright", legend = paste("K observado =", round(observed_G_K$K, 3)), 
       col = "red", lwd = 2)

#M
random_lambdas_M <- sapply(random.total, function(tree) {
  phylosig(tree, trait_data_M, method = "lambda", test = TRUE)
})
random_K_M <- sapply(random.total, function(tree) {
  phylosig(tree, trait_data_M, method = "K", test = TRUE)
})
lambda_values_M <- as.numeric(random_lambdas_M["lambda",])
p_value_M_lambda <- mean(lambda_values_M  >= observed_M_lambda$lambda)
cat("P-value random λ para M:", p_value_M_lambda)
cat("P-value de λ observado para M =", round(observed_M_lambda$P, 3))

K_values_M <- as.numeric(random_K_M["K",])
p_value_M_K <- mean(K_values_M>= observed_M_K$K)
cat("P-value random K para M:", p_value_M_K)
cat("P-value de K observado para M=", round(observed_M_lambda$P, 3))

hist(lambda_values_M, main = "Distribución nula de λ para M", xlab = "λ")
abline(v = observed_M_lambda$lambda, col = "red", lwd = 2)
legend("topright", legend = paste("λ observado =", round(observed_M_lambda$lambda, 3)), 
       col = "red", lwd = 2)

hist(K_values_M, main = "Distribución nula de Blomerg K para M", xlab = "K")
abline(v = observed_M_K$K, col = "red", lwd = 2)
legend("topright", legend = paste("K observado =", round(observed_M_K$K, 3)), 
       col = "red", lwd = 2)

##Árbol de consenso
# Señal filogenética observada (Pagel lambda)
G_physignal_s_pagel <- phylosig(phy_tree_s, trait_data_G_s, method = "lambda",test = T)
M_physignal_s_pagel <- phylosig(phy_tree_s, trait_data_M_s, method = "lambda",test = T)
# Señal filogenética observada (Blomberg K)
G_physignal_s_K <- phylosig(phy_tree_s, trait_data_G_s, method = "K",test = T)
M_physignal_s_K <- phylosig(phy_tree_s, trait_data_M_s, method = "K",test = T) 

#G Pagel
random_signals_G_sample <- sapply(random.trees, function(tree) {
  phylosig(tree, trait_data_G_s, method = "lambda")
})
lambda_values_G_sample <- as.numeric(random_signals_G_sample["lambda",])
p_value_random.trees_G_sample_Pagel  <- mean(lambda_values_G_sample >= G_physignal_s_pagel$lambda)
cat("P-valor G usando random.trees sample Pagel's lambda:", p_value_random.trees_G_sample_Pagel, "\n")
cat("P-value de λ observado para M =", round(G_physignal_s_pagel$P, 3))

#G Blomberg K
random_signals_G_sample_K <- sapply(random.trees, function(tree) {
  phylosig(tree, trait_data_G_s, method = "K")
})
K_values_G_sample <- as.numeric(random_signals_G_sample_K)
p_value_random.trees_G_sample  <- mean(K_values_G_sample >= G_physignal_s_K$K)
cat("P-valor G usando random.trees sample Blomberg K:", p_value_random.trees_G_sample, "\n")
cat("P-value de λ observado para M =", round(G_physignal_s_K$P, 3))

#M Paglel lambda
random_signals_M <- sapply(random.trees, function(tree) {
  phylosig(tree, trait_data_M_s, method = "lambda")$lambda
})
p_value_random.trees.s_M <- mean(random_signals_M >= M_physignal_s_pagel$lambda)
cat("P-valor M usando random.trees.s:", p_value_random.trees.s_M, "\n")
cat("P-valor M usando random.trees.s:", M_physignal_s_pagel$P, "\n")

#M Blomberg K
random_signals_M_sample_K <- sapply(random.trees, function(tree) {
  phylosig(tree, trait_data_M_s, method = "K")
})
p_value_random.trees_M_sample  <- mean(random_signals_M_sample_K >= M_physignal_s_K$K)
cat("P-valor M usando random.trees sample Blomberg K:", p_value_random.trees_M_sample, "\n")
cat("P-valor M usando random.trees sample Blomberg K:", M_physignal_s_K$K, "\n")

