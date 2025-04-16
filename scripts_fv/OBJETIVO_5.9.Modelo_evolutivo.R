################################################################################
################################################################################
######### Project 111160402: Cultural phylogenetics and coevolution of wealth###
###### inheritance and land tenure norms in agropastoralist communities #########
################################################################################
################################################################################

# OBJETIVO 5
# Fit models of trait evolution (e.g., Brownian motion, Ornstein-Uhlenbeck)
# to G and M. Compare AIC for model selection and perform simulations.

# Cargar las librerías
library(dplyr)
library(geiger)
library(phytools)
library(surface)
library(ape)
library(phangorn)

# Trait data
G_values <- GM_df %>%
  dplyr::select(community, G)
M_values <- GM_df %>%
  dplyr::select(community, M)

########################## SAMPLED COMMUNITIES
###############################
########################## G #################################
# Filtrar los datos para que coincidan con las etiquetas del árbol
tree_tips <- consensus_tree$tip.label
data_tips <- G_values$community
setdiff(data_tips, tree_tips)  # Nombres en los datos que no están en el árbol
G_df_filtered_s <- G_values[rownames(G_values) %in% tree_tips, ]
G_df_filtered_s <- dplyr::select(G_df_filtered_s,G)
G_df_filtered_s

# Fit BM model for G
fit_bm_G_s <- fitContinuous(consensus_tree, G_df_filtered_s, model = "BM")
print(fit_bm_G_s)

# Fit OU model for G
fit_ou_G_s <- fitContinuous(consensus_tree, G_df_filtered_s, model = "OU")
print(fit_ou_G_s)

# Comparar AIC
AIC(fit_bm_G_s) ## ESTE?
AIC(fit_ou_G_s)
AIC(fit_bm_G_s,fit_ou_G_s)
cat("Delta AIC (BM - OU):", AIC(fit_bm_G_s) - AIC(fit_ou_G_s), "\n")

# Fase forward-backward con surface
consensus_tree_named <- nameNodes(consensus_tree)
olist <- convertTreeData(consensus_tree_named, G_df_filtered_s)
otree <- olist[[1]]
data <- olist[[2]]
forward_model_G_s <- surfaceForward(otree, data)
summary(forward_model_G_s)
k <- length(forward_model_G_s)
backward_model_G_s <- surfaceBackward(otree, data, starting_model = forward_model_G_s[[k]])
surfaceSummary(fwd = forward_model_G_s, bwd = backward_model_G_s)

# Graficar resultados
dataframe_to_plot <- G_df_filtered_s
names(dataframe_to_plot) <- G_df_filtered_s$community
surfaceTreePlot(consensus_tree, backward_model_G_s[[1]], labelshifts = T)
surfaceTraitPlot(G_df_filtered_s, backward_model_G_s[[1]])

########################## M #################################
# Filtrar los datos para que coincidan con las etiquetas del árbol
tree_tips <- consensus_tree$tip.label
data_tips <- M_values$community
setdiff(data_tips, tree_tips)  # Nombres en los datos que no están en el árbol
M_df_filtered_s <- M_values[rownames(M_values) %in% tree_tips, ]
M_df_filtered_s <- dplyr::select(M_df_filtered_sample,M)

# Fit BM model for M
fit_bm_M_s <- fitContinuous(consensus_tree, M_df_filtered_s, model = "BM")
print(fit_bm_M_s)

# Fit OU model for M
fit_ou_M_s <- fitContinuous(consensus_tree, M_df_filtered_s, model = "OU")
print(fit_ou_M_s)

# Comparar AIC
AIC(fit_bm_M_s) ## Este??
AIC(fit_ou_M_s)
AIC(fit_bm_M_s,fit_ou_M_s)
cat("Delta AIC (BM - OU):", AIC(fit_bm_G_s) - AIC(fit_ou_G_s), "\n")


# Fase forward-backward con surface
olist <- convertTreeData(consensus_tree_named, M_df_filtered_s)
otree <- olist[[1]]
data <- olist[[2]]
forward_model_M_s <- surfaceForward(otree, data)
summary(forward_model_M_s)
k <- length(forward_model_M_s)
backward_model_M_s <- surfaceBackward(otree, data, starting_model = forward_model_M_s[[k]])
surfaceSummary(fwd = forward_model_M_s, bwd = backward_model_M_s)

# Graficar resultados
dataframe_to_plot <- M_df_filtered_s
names(dataframe_to_plot) <- M_df_filtered_s$community
surfaceTreePlot(consensus_tree, backward_model_M_s[[1]], labelshifts = TRUE)
surfaceTraitPlot(dataframe_to_plot, backward_model_M_s[[1]])


########################## SIMULACIONES ARBOL ##################################
#Los modelos nulos simulados evalúan si los patrones de agrupamiento observados en los datos son estadísticamente significativos al compararlos con una distribución nula generada a partir de datos aleatorios o simulados. Este enfoque está relacionado en cierta medida con la prueba de Mantel (ambos comparan los patrones observados con las expectativas nulas), pero es más flexible y se adapta a hipótesis específicas.

library(picante)
library(phytools)
library(ape)

#11. Randomise tips (with their respective traits, G and M) keeping the underlying tree topology.
# Configuración de los datos
trait_data_G <- setNames(G_values$G, G_values$community)
trait_data_M <- setNames(M_values$M, G_values$community)
tree_tips <- consensus_tree$tip.label
total_tips <- y_total$tip.label
# Randomize tip labels
## Randomize trees
howmanytrees(length(tree_tips)) #número de árboles filogenéticos posibles

# Topología y ramas aleatorias
#random.trees <- rmtree(1000, n = length(tree$tip.label))

# Topología aleatoria
#random.trees <- lapply(1:1000, function(x) rtree(n = length(tree$tip.label)))
#class(random.trees) <- "multiPhylo"

# Etiquetas aleatorizadas
#consenso
random.trees <- lapply(1:1000, function(x) {
  tree_random <- consensus_tree  # Copia el árbol original
  tree_random$tip.label <- sample(tree_tips)  # Reordena las etiquetas
  tree_random
})
class(random.trees) <- "multiPhylo"
#total
random.total <- lapply(1:1000, function(x) {
  tree_random <- y_total  # Copia el árbol original
  tree_random$tip.label <- sample(total_tips)  # Reordena las etiquetas
  tree_random
})
class(random.total) <- "multiPhylo"


# Calcular la métrica observada (puede cambiarse por otra, aquí usamos la varianza)
observed_metric_G <- var(trait_data_G)
observed_metric_M <- var(trait_data_M)

# Calcular métricas nulas usando los árboles aleatorios
null_metrics_G <- sapply(random.trees, function(tree) {
  var(trait_data_G)  # Métrica basada en el rasgo G
})

null_metrics_M <- sapply(random.trees, function(tree) {
  var(trait_data_M)  # Métrica basada en el rasgo M
})

# Calcular p-valores
p_value_G <- mean(null_metrics_G >= observed_metric_G)
p_value_M <- mean(null_metrics_M >= observed_metric_M)

# Imprimir resultados
cat("P-valor para G (usando random trees):", p_value_G, "\n")
cat("P-valor para M (usando random trees):", p_value_M, "\n")


#12. Compare phylogenetic signal of real tree vs randomised tree.
# Asegurar correspondencia entre rasgos y puntas
trait_data_G <- setNames(G_values$G, G_values$community)
trait_data_M <- setNames(M_values$M, M_values$community)
phy_tree_s <- consensus_tree
phy_tree <- y_total

# Señal filogenética observada
G_physignal <- phylosig(phy_tree, trait_data_G, method = "lambda",test = T)
M_physignal <- phylosig(phy_tree, trait_data_M, method = "lambda",test = T)
G_physignal_s <- phylosig(phy_tree_s, trait_data_G, method = "lambda",test = T)
M_physignal_s <- phylosig(phy_tree_s, trait_data_M, method = "lambda",test = T)

##Random.trees
#Todas las comunidades
#random.total <- lapply(1:1000, function(x) rtopology(length(total_tips), tip.label = total_tips))
#class(random.total)<-"multiPhylo"
#Árbol de consenso
#random.trees <- lapply(1:1000, function(x) rtopology(length(tree_tips), tip.label = tree_tips))
#class(random.trees)<-"multiPhylo"

###Señal filogenética
##Todas las comunidades
#G
random_signals_G <- sapply(random.total, function(tree) {
  phylosig(tree, trait_data_G, method = "lambda")$lambda
})
p_value_random_trees_G <- mean(random_signals_G >= G_physignal$lambda)
cat("P-valor G usando random trees:", p_value_random_trees_G, "\n")

null_signals_G_bm <- apply(simulated_traits_G, 2, function(trait) {
  phylosig(phy_tree, trait, method = "lambda")$lambda
})
p_value_G_bm <- mean(null_signals_G_bm >= G_physignal$lambda)
cat("P-valor G (Brownian motion):", p_value_G_bm, "\n")

hist(random_G, main = "Distribución nula para G", xlab = "Lambda")
abline(v = G_physignal, col = "red", lwd = 2)

#M
random_signals_M <- sapply(random.total, function(tree) {
  phylosig(tree, trait_data_M, method = "lambda")$lambda
})
p_value_random_trees_M <- mean(random_signals_M >= M_physignal$lambda)
cat("P-valor M usando random trees:", p_value_random_trees_M, "\n")

null_signals_M_bm <- apply(simulated_traits_M, 2, function(trait) {
  phylosig(phy_tree, trait, method = "lambda")$lambda
})
p_value_M_bm <- mean(null_signals_M_bm >= M_physignal)
cat("P-valor M (Brownian motion):", p_value_M_bm, "\n")

hist(random_M, main = "Distribución nula para M", xlab = "Lambda")
abline(v = M_physignal, col = "red", lwd = 2)

##Árbol de consenso
#G
random_signals_G_sample <- sapply(random.trees, function(tree) {
  phylosig(tree, trait_data_G, method = "lambda")$lambda
})
p_value_random.trees_G_sample  <- mean(random_signals_G >= G_physignal$lambda)
cat("P-valor G usando random.trees.s:", p_value_random.trees_G_sample, "\n")

null_signals_G_bm_sample <- apply(simulated_traits_G, 2, function(trait) {
  phylosig(phy_tree_s, trait, method = "lambda")$lambda
})
p_value_G_bm_sample  <- mean(null_signals_G_bm_sample  >= G_physignal$lambda)
cat("P-valor G (Brownian motion):", p_value_G_bm_sample , "\n")

hist(random_G, main = "Distribución nula para G", xlab = "Lambda")
abline(v = G_physignal, col = "red", lwd = 2)

#M
random_signals_M <- sapply(random.trees, function(tree) {
  phylosig(tree, trait_data_M, method = "lambda")$lambda
})
p_value_random.trees.s_M <- mean(random_signals_M >= M_physignal$lambda)
cat("P-valor M usando random.trees.s:", p_value_random.trees.s_M, "\n")

null_signals_M_bm <- apply(simulated_traits_M, 2, function(trait) {
  phylosig(phy_tree_s, trait, method = "lambda")$lambda
})
p_value_M_bm <- mean(null_signals_M_bm >= M_physignal$lambda)
cat("P-valor M (Brownian motion):", p_value_M_bm, "\n")

hist(random_M, main = "Distribución nula para M", xlab = "Lambda")
abline(v = M_physignal, col = "red", lwd = 2)