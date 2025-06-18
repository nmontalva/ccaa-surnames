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

## ---------------------------
## PREPARACION DE DATOS
## ---------------------------

# Trait data
G_values <- GM_logit %>%
  dplyr::select(G_logit)
M_values <- GM_logit %>%
  dplyr::select(M_logit)

# Randomise tips (with their respective traits, G and M) keeping the underlying tree topology.
# Configuración de los datos
<<<<<<< HEAD
#total
trait_data_G <- setNames(G_values$G_logit, row.names(G_values))
trait_data_M <- setNames(M_values$M_logit, row.names(M_values))
total_tips <- y_total$tip.label
#consenso
tree_tips <- consensus_tree$tip.label
=======
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
>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01
trait_data_G_s <- trait_data_G[names(trait_data_G) %in% tree_tips]
trait_data_M_s <- trait_data_M[names(trait_data_M) %in% tree_tips]

## ---------------------------
## GENERACION DE ARBOLES RANDOMIZADOS
## ---------------------------
#howmanytrees(length(tree_tips)) #numero de arboles filogeneticos posibles
# Generar arboles aleatorizados
generate_random_trees <- function(phy, n=1000) {
  random_trees <- replicate(n, {
    rt <- phy
    rt$tip.label <- sample(phy$tip.label)  # Permutación de tips
    return(rt)
  }, simplify = FALSE)
  class(random_trees) <- "multiPhylo"
  return(random_trees)
}

set.seed(150)  # Para reproducibilidad
random_trees_consensus <- generate_random_trees(consensus_tree, 1000) #consenso
random_trees_complete <- generate_random_trees(y_total, 1000) #total

<<<<<<< HEAD
## ---------------------------
## CALCULO DE SEÑAL FILOGENETICA
## ---------------------------
#Funcion para calculo de señal filogenetica
phylo_signal_analysis <- function(phy, trait, tree_set, trait_name="Trait") {
  
  # Señal observada
  obs_lambda <- phylosig(phy, trait, method="lambda", test=TRUE, nsim = 9999)
  obs_K <- phylosig(phy, trait, method="K", test=TRUE, nsim = 9999)
  
  # Señal en árboles randomizados (solo cálculo, sin test para mayor velocidad)
  rand_lambda <- sapply(tree_set, function(x) phylosig(x, trait, method="lambda",test=TRUE, nsim = 9999))
  rand_K <- sapply(tree_set, function(x) phylosig(x, trait, method="K",test=TRUE, nsim = 9999))
                
                return(list(obs_lambda=obs_lambda, obs_K=obs_K, 
                            rand_lambda=rand_lambda, rand_K=rand_K))
                }
## ---------------------------
## ANALISIS PARA CADA TRAIT Y ARBOL
## ---------------------------
# Análisis para el árbol completo
results_G_complete <- phylo_signal_analysis(y_total, trait_data_G, random_trees_complete, "G (árbol completo)")
results_M_complete <- phylo_signal_analysis(y_total, trait_data_M, random_trees_complete, "M (árbol completo)")
# Análisis para el árbol de consenso
results_G_consensus <- phylo_signal_analysis(consensus_tree, trait_data_G_s, random_trees_consensus, "G (consenso)")
results_M_consensus <- phylo_signal_analysis(consensus_tree, trait_data_M_s, random_trees_consensus, "M (consenso)")



# Comparaciones
comparar_phylo_signal <- function(results, nombre = "", archivo = NULL) {
  if (!is.null(archivo)) {
    png(archivo, width = 1000, height = 1000, res = 150)  # Espacio amplio
  }
  layout(matrix(1:4, 2, 2, byrow = TRUE))  # Panel 2x2 ordenado por filas
  par(mar = c(5, 4, 4, 2) + 0.1)           # Márgenes cómodos
  ## LAMBDA
  obs_lambda <- results$obs_lambda$lambda
  p_lambda <- results$obs_lambda$P
  sim_lambda <- as.numeric(results$rand_lambda["lambda", ])
  sim_p_lambda <- as.numeric(results$rand_lambda["P", ])
  
  # Histograma de lambda simulado
  hist(sim_lambda, breaks = 30, col = "lightblue",
       main = paste("Lambda simulado\n", nombre),
       xlab = "Lambda simulado")
  abline(v = obs_lambda, col = "red", lwd = 2)
  legend("topright",
         legend = paste("Lambda observado =", round(obs_lambda, 3)),
         col = "red", lwd = 2)
  
  # Histograma de p-valores simulados lambda
  hist(sim_p_lambda, breaks = 30, col = "pink",
       main = paste("p-valores simulados (Lambda)\n", nombre),
       xlab = "p-valor")
  abline(v = p_lambda, col = "red", lwd = 2)
  legend("topright",
         legend = paste("p observado =", round(p_lambda, 4)),
         col = "red", lwd = 2)
  
  ## K
  obs_K <- results$obs_K$K
  p_K <- results$obs_K$P
  sim_K <- as.numeric(results$rand_K["K", ])
  sim_p_K <- as.numeric(results$rand_K["P", ])
  
  # Histograma de K simulado
  hist(sim_K, breaks = 30, col = "lightblue",
       main = paste("K simulado\n", nombre),
       xlab = "K simulado")
  abline(v = obs_K, col = "red", lwd = 2)
  legend("topright",
         legend = paste("K observado =", round(obs_K, 3)),
         col = "red", lwd = 2)
  
  # Histograma de p-valores simulados K
  hist(sim_p_K, breaks = 30, col = "pink",
       main = paste("p-valores simulados (K)\n", nombre),
       xlab = "p-valor")
  abline(v = p_K, col = "red", lwd = 2)
  legend("topright",
         legend = paste("p observado =", round(p_K, 4)),
         col = "red", lwd = 2)
  if (!is.null(archivo)) {
    dev.off()
  }
  }
comparar_phylo_signal(results_G_complete, "G (Árbol completo)",archivo = "Figures/phylo_signal_G_total.png")
comparar_phylo_signal(results_M_complete, "M (Árbol completo)",archivo = "Figures/phylo_signal_M_total.png")
comparar_phylo_signal(results_G_consensus, "G (Árbol consenso)",archivo = "Figures/phylo_signal_G_consenso.png")
comparar_phylo_signal(results_M_consensus, "M (Árbol consenso)",archivo = "Figures/phylo_signal_M_consenso.png")
=======
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

>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01
