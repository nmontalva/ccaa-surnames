################################################################################
################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth####
######inheritance and land tenure norms in agropastoralist communities##########
################################################################################
################################################################################

#### OBJETIVO 5 ####
### Part 1 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###
##IMPORTANTE: CORRER LOS SCRIPTS DE LOS OBJETIVOS 2 y 4 ##

#### CARGAR E INSTALLAR LIBRERIAS ####
library(tidyverse)
library(ape)
library(car)
library(caper)
library(conflicted)
library(dplyr)
library(geiger)
library(nlme)
library(phytools)
library(Publish)
library(BiocManager) #esto se demora mucho en compilar cosas. No sé si se puede instalar en modo binario ejecutable
library("treeio") #Ver "packages.R" para instalar este paqute
library(ggplot2)
library(gridExtra)
library(viridisLite) #colores
library(rr2)
library(grDevices)
library(png)

# Resolver conflictos de funciones
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("treedata", "geiger")
conflict_prefer("as.phylo", "phylogram")

#### CREACION DE ARBOLES PARA CADA TRAIT ####
#Ocuparemos (1) arbol de apellidos(y_total), (2) todas las tablas del objetivo 2, y (3) consenso STR/Apellido(consensus_3)

###Explorar los datos
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
selected_communities <- unique(STR$pop)

##y_total
plotTree(y_total,type="phylogram",ftype="off")

##Traits
select_variable <- function(result_traits, selected_communities = NULL, variable) {
  df <- as.data.frame(result_traits)
  row.names(df) <- result_traits$community
  if (!is.null(selected_communities)) {
    df <- df[df$community %in% selected_communities, ]
  }
  df <- df[, c("community", variable), drop = FALSE]
  row.names(df) <- df$community
  return(df)
}

# S
S_trait <- select_variable(result_traits, selected_communities, "S")
S_trait2 <- select_variable(result_traits, NULL, "S")
S_trait$community <- NULL
S_trait2$community <- NULL
S_logit <- select_variable(result_traits, selected_communities, "S_logit")
S_logit2 <- select_variable(result_traits, NULL, "S_logit")
S_logit$community <- NULL
S_logit2$community <- NULL

#N_trait2 <- select_variable(result_traits, selected_communities, "N")

# G
G_trait <- select_variable(result_traits, selected_communities, "G")
G_trait2 <- select_variable(result_traits, NULL, "G")
G_trait$community <- NULL
G_trait2$community <- NULL
G_logit <- select_variable(result_traits, selected_communities, "G_logit")
G_logit2 <- select_variable(result_traits, NULL, "G_logit")
G_logit$community <- NULL
G_logit2$community <- NULL
# A
A_trait <- select_variable(result_traits, selected_communities, "A")
A_trait2 <- select_variable(result_traits, NULL, "A")
A_trait$community <- NULL
A_trait2$community <- NULL
A_logit <- select_variable(result_traits, selected_communities, "A_logit")
A_logit2 <- select_variable(result_traits, NULL, "A_logit")
A_logit$community <- NULL
A_logit2$community <- NULL
# M
M_trait <- select_variable(result_traits, selected_communities, "M")
M_trait2 <- select_variable(result_traits, NULL, "M")
M_trait$community <- NULL
M_trait2$community <- NULL
M_logit <- select_variable(result_traits, selected_communities, "M_logit")
M_logit2 <- select_variable(result_traits, NULL, "M_logit")
M_logit$community <- NULL
M_logit2$community <- NULL
#Cluster
#C_trait <- select_variable(result_traits, selected_communities, "cluster")
#C_trait2 <- select_variable(result_traits, NULL, "cluster")
#C_trait$community <- NULL
#C_trait2$community <- NULL

#Exponenciar logit
logit_to_prob <- function(logit) {
  prob <- exp(logit) / (1 + exp(logit))
  return(prob)
}
# Transformar los valores logit a probabilidades (exponencial)
S_prob <- logit_to_prob(S_logit)
S_prob2 <- logit_to_prob(S_logit2)
G_prob <- logit_to_prob(G_logit)
G_prob2 <- logit_to_prob(G_logit2)
A_prob <- logit_to_prob(A_logit)
A_prob2 <- logit_to_prob(A_logit2)
M_prob <- logit_to_prob(M_logit)
M_prob2 <- logit_to_prob(M_logit2)

#Final_table
ft <- result_traits
ft <- ft %>% column_to_rownames(var = "community")
ft2 <- result_traits %>% filter(community %in% selected_communities)
ft2 <- ft2 %>% column_to_rownames(var = "community")

##Consensus##
consensus_tree <-as.phylo(consensus_tree)
plotTree(consensus_tree,type="phylogram", ftype="i",lwd=1)

###Estimar estados ancestrales
# S
sv1 <- as.matrix(S_logit)[,1] # con logit
sv2 <- as.matrix(S_logit2)[,1]
sv1p <- as.matrix(S_prob)[,1] # Con exponencial
sv2p <- as.matrix(S_prob2)[,1]

# G
gv1 <- as.matrix(G_logit)[,1]
gv2 <- as.matrix(G_logit2)[,1]
gv1p <- as.matrix(G_prob)[,1]
gv2p <- as.matrix(G_prob2)[,1]
# A
av1 <- as.matrix(A_logit)[,1]
av2 <- as.matrix(A_logit2)[,1]
av1p <- as.matrix(A_prob)[,1]
av2p <- as.matrix(A_prob2)[,1]
# M
mv1 <- as.matrix(M_logit)[,1]
mv2 <- as.matrix(M_logit2)[,1]
mv1p <- as.matrix(M_prob)[,1]
mv2p <- as.matrix(M_prob2)[,1]

# Funcion para estimar estados ancestrales y crear el arbol
estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = FALSE)
  obj$lims <- c(-6, 6)
  obj <- setMap(obj, viridisLite::viridis(n=100))
  plot(obj,lwd = 10, type = "phylogram", legend = 0.7 * max(nodeHeights(tree)),outline=FALSE, ftype = "i",fsize = c(0.7, 1.5), leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}

# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
#svg("outputs/Figures/S_muestra.svg")
#svc <- estimacion_estados_ancestrales(consensus_tree, sv1, "S")
#dev.off()

svg("outputs/Figures/G_muestra.svg")
gvc <- estimacion_estados_ancestrales(consensus_tree, gv1, "G")
dev.off()

#svg("outputs/Figures/A_muestra.svg")
#avc <- estimacion_estados_ancestrales(consensus_tree, av1, "A")
#dev.off()

svg("outputs/Figures/M_muestra.svg")
mvc <- estimacion_estados_ancestrales(consensus_tree, mv1, "M")
dev.off()

save(gvc, mvc, file = "outputs/Figures/ancestral_states_sampled.RData")

estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = FALSE)
  obj$lims <- c(-6,6)
  obj <- setMap(obj, viridisLite::viridis(n=100))
  plot(obj,lwd=4, type = "phylogram", legend = 0.5 * max(nodeHeights(tree)),outline=FALSE, ftype = "off",fsize = c(0, 1.5), leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}
# Estimar estados ancestrales para los ?rboles de apellidos totales
#svg("outputs/Figures/S_total.svg")
#svt <- estimacion_estados_ancestrales(y_total, sv2, "S")
#dev.off()
svg("outputs/Figures/G_total.svg")
gvt <- estimacion_estados_ancestrales(y_total, gv2, "G")
dev.off()
#svg("outputs/Figures/A_total.svg")
#avt <- estimacion_estados_ancestrales(y_total, av2, "A")
#dev.off()
svg("outputs/Figures/M_total.svg")
mvt <- estimacion_estados_ancestrales(y_total, mv2, "M")
dev.off()

save(gvt,mvt, file = "outputs/Figures/ancestral_states_all.RData")

###Comparaciones
##Comunidades muestreadas
plot_comparison <- function(obj1, obj2, xlabel) {
  par(mfrow = c(1, 2), mar = c(5, 5, 5, 5) + 1, oma = c(1, 1, 1, 1))
  # First plot
  obj1$lims <- c(-6,6)
  plot(obj1, lwd = 10, ftype = "off", legend = 0.7 * max(nodeHeights(consensus_tree)), outline = FALSE, fsize = c(0,1.5), leg.txt = xlabel[1])
  
  # Second plot
  obj2$lims <- c(-6,6)
  plot(obj2, lwd = 10, outline = FALSE, direction = "leftwards", ftype = "off", legend = 0.7 * max(nodeHeights(consensus_tree)), fsize = c(0, 1.5), leg.txt = xlabel[2])
}

# Comparaciones
#svg("outputs/Figures/S_A_muestra.svg",width = 12, height = 12)
#plot_comparison(svc$obj, avc$obj, c("S", "A"))
#dev.off()
#svg("outputs/Figures/S_G_muestra.svg",width = 12, height = 12)
#plot_comparison(svc$obj, gvc$obj, c("S", "G"))
#dev.off()
#svg("outputs/Figures/S_M_muestra.svg",width = 12, height = 12)
#plot_comparison(svc$obj, mvc$obj, c("S", "M"))
#dev.off()
#svg("outputs/Figures/A_G_muestra.svg",width = 12, height = 12)
#plot_comparison(avc$obj, gvc$obj, c("G", "A"))
#dev.off()
#svg("outputs/Figures/A_M_muestra.svg",width = 12, height = 12)
#plot_comparison(avc$obj, mvc$obj, c("A", "M"))
#dev.off()
svg("outputs/Figures/G_M_muestra.svg",width = 14, height = 14)
plot_comparison(gvc$obj, mvc$obj, c("G", "M"))
dev.off()

##Comunidades totales
plot_comparison2 <- function(obj1, obj2, xlabel) {
  par(pin=c(15, 15))
  par(mai = c(2, 10, 2, 10))
  layout(matrix(1:3, 1, 3), widths=c(0.48, 0.01, 0.48))
  obj1$lims <- c(-6,6)
  plot(obj1, lwd=7, ftype="off", outline=FALSE,legend=1*max(nodeHeights(y_total)),fsize=c(0,3),leg.txt=xlabel[1])
  plot.new()
  plot.window(xlim=c(1,2), ylim=get("last_plot.phylo", envir=.PlotPhyloEnv)$y.lim)
  obj2$lims <- c(-6,6)
  plot(obj2, lwd=7,outline=FALSE, direction="leftwards", legend=1*max(nodeHeights(y_total)),
       ftype="off",fsize=c(0,3),leg.txt=xlabel[2])
}

# Comparaciones
dev.off()
#svg("outputs/Figures/S_G_total.svg",width = 20, height = 20)
#plot_comparison2(svt$obj, gvt$obj, c("S", "G"))
#dev.off()
#svg("outputs/Figures/S_A_total.svg",width = 20, height = 20)
#plot_comparison2(svt$obj, avt$obj, c("S", "A"))
#dev.off()
#svg("outputs/Figures/S_M_total.svg",width = 20, height = 20)
#plot_comparison2(svt$obj, mvt$obj, c("S", "M"))
#dev.off()
#svg("outputs/Figures/A_G_total.svg",width = 20, height = 20)
#plot_comparison2(avt$obj,gvt$obj, c("G", "A"))
#dev.off()
#svg("outputs/Figures/A_M_total.svg",width = 20, height = 20)
#plot_comparison2(avt$obj, mvt$obj, c("A", "M"))
#dev.off()
svg("outputs/Figures/G_M_total.svg",width = 14, height = 14)
plot_comparison2(gvt$obj, mvt$obj, c("G", "M"))
dev.off()

### Phylosignal
safe_phylosignal <- function(tree, trait_vector, trait_name = "Trait") {
  
  # 1. Asegurar que el trait tenga nombres
  if (is.null(names(trait_vector))) {
    stop(paste("El vector", trait_name, "no tiene nombres de especies asignados."))
  }
  
  # 2. Intersección estricta: Quedarse solo con lo que existe en ambos lados
  common_taxa <- intersect(tree$tip.label, names(trait_vector))
  
  if (length(common_taxa) < 3) {
    stop(paste("Menos de 3 especies coinciden para", trait_name, "- No se puede calcular señal."))
  }
  
  # 3. ALINEACIÓN FORZOSA (La clave para evitar el error de lambda)
  # Recortar el árbol
  tree_clean <- keep.tip(tree, common_taxa)
  # Reordenar el trait para que coincida EXACTAMENTE con el orden del árbol recortado
  trait_clean <- trait_vector[tree_clean$tip.label] 
  
  # Asegurar numérico
  trait_clean <- as.numeric(trait_clean)
  names(trait_clean) <- tree_clean$tip.label
  
  # 4. Cálculo de Lambda (Pagel)
  # test=TRUE hace un Likelihood Ratio Test (LRT) comparando el modelo con Lambda libre vs Lambda=0
  res_lambda <- phylosig(tree_clean, trait_clean, method = "lambda", test = TRUE)
  
  # 5. Cálculo de K (Blomberg)
  # nsim=1000 es suficiente usualmente, sube a 9999 si necesitas p-valores muy precisos
  res_K <- phylosig(tree_clean, trait_clean, method = "K", test = TRUE, nsim = 1000)
  
  return(data.frame(
    Trait = trait_name,
    N_Taxa = length(common_taxa),
    Lambda = round(res_lambda$lambda, 5),
    p_Lambda = round(res_lambda$P, 5), # p-valor del LRT
    K = round(res_K$K, 5),
    p_K = round(res_K$P, 5), # p-valor de la permutación
    stringsAsFactors = FALSE
  ))
}

# ==============================================================================
# ARBOL DE CONSENSO
# ==============================================================================
list_muestra <- list(S = sv1, G = gv1, A = av1, M = mv1)

# Aplicar función
results_muestra <- do.call(rbind, lapply(names(list_muestra), function(n) {
  safe_phylosignal(consensus_tree, list_muestra[[n]], trait_name = n)
}))

print(results_muestra)

# Guardar Gráfico
svg("outputs/Figures/Phylosignal_muestra_check.svg", width = 7, height = 3)
grid.table(results_muestra, rows = NULL) # rows=NULL quita los números de fila
dev.off()

# ==============================================================================
# ÁRBOL COMPLETO
# ==============================================================================
list_total <- list(S = sv2, G = gv2, A = av2, M = mv2)

results_total <- do.call(rbind, lapply(names(list_total), function(n) {
  safe_phylosignal(y_total, list_total[[n]], trait_name = n)
}))

print(results_total)

# Guardar Gráfico
svg("outputs/Figures/Phylosignal_total_check.svg", width = 7, height = 3)
grid.table(results_total, rows = NULL)
dev.off()

###Regresion PGLS
##Phylo-regression for sampled data
sampled_matched <- geiger::treedata(consensus_tree, ft2, sort=F, warnings=TRUE)
spc <- sampled_matched$phy$tip.label
V<-corPagel(1,phy=sampled_matched$phy,form=~spc, fixed=FALSE)
C <- vcv.phylo(phy = sampled_matched$phy)
#Assign traits
obj <- ft2
x <- write.tree(consensus_tree)
tree_consensus.tree <- read.tree(text=x)
G <- setNames(obj[,"G_logit"],rownames(obj))
M <- setNames(obj[,"M_logit"],rownames(obj))
S <- setNames(obj[,"S_logit"],rownames(obj))
A <- setNames(obj[,"A_logit"],rownames(obj))
N <- setNames(obj[,"N"],rownames(obj))

# Verifica que N, M, A, S, y G sean vectores numéricos
N <- as.vector(N)
M <- as.vector(M)
A <- as.vector(A)
S <- as.vector(S)
G <- as.vector(G)

# Crea un dataframe con estos vectores
vector.data <- data.frame(N = N, M = M, A = A, S = S, G = G)
get_model_metrics <- function(model, model_name) {
  require(rr2)
  
  lambda_val <- if(inherits(model$modelStruct$corStruct, "corPagel")) {
    model$modelStruct$corStruct[[1]]
  } else {
    1 # Para corBrownian
  }
  
  metrics <- list(
    model_name = model_name,
    formula = deparse(formula(model)),
    coefficients = summary(model)$tTable,
    anova = anova(model),
    logLik = as.numeric(logLik(model)),
    AICc = AIC(model),
    AICc = AIC(model) + (2 * length(coef(model)) * (length(coef(model)) + 1)) / 
      (nobs(model) - length(coef(model)) - 1),
    BIC = BIC(model),
    R2 = R2_lik(model),
    lambda = lambda_val
  )
  
  return(metrics)
}

# Modelo N (ahora capturamos los resultados)
bm_glsN <- gls(N~M+A+S+G, correlation = V, data=vector.data)
metrics_N <- get_model_metrics(bm_glsN, "N ~ M + A + S + G")

# Modelo G
bm_glsG <- gls(G~M+A+S, correlation = V, data=vector.data)
metrics_G <- get_model_metrics(bm_glsG, "G ~ M + A + S")
print(metrics_G)
publish(bm_glsG)
# Modelo M
bm_glsM <- gls(M~G+A+S, correlation = V, data=vector.data)
metrics_M <- get_model_metrics(bm_glsM, "M ~ G + A + S")

##Phylo-regression for all communities
sampled_matched <- geiger::treedata(y_total, ft, sort=FALSE, warnings=TRUE)
spc <- sampled_matched$phy$tip.label
ft <- ft[spc, ]
V<-corPagel(1,phy=sampled_matched$phy,form = ~spc,fixed=FALSE)
C <- vcv.phylo(phy = sampled_matched$phy)

#Now we run the analysis:
#Assign traits
trait_names <- c("G_logit", "M_logit", "S_logit", "A_logit", "N")
traits <- lapply(trait_names, function(x) setNames(ft[,x], rownames(ft)))
names(traits) <- substr(trait_names, 1, 1)  # G, M, S, A, N
list2env(traits, envir = .GlobalEnv)

# Crear dataframe (asegurando orden)
vector.data <- data.frame(N, M, A, S, G, row.names = spc)

# Modelo N (ahora capturamos los resultados)
bm_glsNt <- gls(N~M+A+S+G, correlation = V, data=vector.data)
metrics_Nt <- get_model_metrics(bm_glsNt, "N ~ M + A + S + G")

# Modelo G
bm_glsGt <- gls(G~M+A+S, correlation = V, data=vector.data)
metrics_Gt <- get_model_metrics(bm_glsGt, "G ~ M + A + S")
print(metrics_Gt)
publish(bm_glsGt)
# Modelo M
bm_glsMt <- gls(M~G+A+S, correlation = V, data=vector.data)
metrics_Mt <- get_model_metrics(bm_glsM, "M ~ G + A + S")

##Generate PICs and test while conditioning on phylogeny
#Prepare the tree
consensus_tree_rooted <- midpoint.root(consensus_tree) 
# Resuelve las politomías
consensus_tree_dicotomous <- multi2di(consensus_tree_rooted) 
#Sampled communities
S_pic1<-pic(x = ft2$S_logit, phy = consensus_tree_dicotomous)
A_pic1<-pic(x = ft2$A_logit, phy = consensus_tree_dicotomous)
G_pic1<-pic(x=ft2$G_logit, phy = consensus_tree_dicotomous)
M_pic1<-pic(x=ft2$M_logit, phy = consensus_tree_dicotomous)

data1 <- data.frame(S_pic1,A_pic1,G_pic1,M_pic1)

#Total communities
S_pic<-pic(x = ft$S_logit, phy = y_total)
A_pic<-pic(x = ft$A_logit, phy = y_total)
G_pic<-pic(x=ft$G_logit, phy = y_total)
M_pic<-pic(x=ft$M_logit, phy = y_total)

data <- data.frame(S_pic,A_pic,G_pic,M_pic)