##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 5 ####
### Part 1 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###
##IMPORTANTE: CORRER LOS SCRIPTS DE LOS OBJETIVOS 2 y 4 ##

#### ESPACIO DE TRABAJO ####
#setwd("C:/Users/Kibif/Desktop/Proyecto desigualdad agropastores/Directorio_proyecto")

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
#library(treeio) #No lo puedo instalar
library(ggplot2)
library(gridExtra)
library(geomorph) #No lo puedo instalar
#install.packages("geomorph") 
#library(devtools)
#install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)

# Resolver conflictos de funciones
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("treedata", "geiger")
conflict_prefer("as.phylo", "phylogram")
#### CREACI�N DE �RBOLES PARA CADA TRAIT ####
#Ocuparemos (1) �rbol de apellidos(y_total), (2) todas las tablas del objetivo 2, y (3) consenso STR/Apellido(consensus_3)

###Explorar los datos
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
selected_communities <- unique(STR$pop)

##y_total
plotTree(y_total,type="phylogram",ftype="off")

##Traits
select_variable <- function(result, selected_communities = NULL, variable) {
  df <- as.data.frame(result)
  row.names(df) <- result$community
  if (!is.null(selected_communities)) {
    df <- df[df$community %in% selected_communities, ]
  }
  df <- df[, c("community", variable), drop = FALSE]
  row.names(df) <- df$community
  return(df)
}


# S
S_trait <- select_variable(result, selected_communities, "S")
S_trait2 <- select_variable(result, NULL, "S")
S_trait$community <- NULL
S_trait2$community <- NULL
N_trait2 <- select_variable(result, selected_communities, "N")

# R
R_trait <- select_variable(result, selected_communities, "R")
R_trait2 <- select_variable(result, NULL, "R")
R_trait$community <- NULL
R_trait2$community <- NULL
# G
G_trait <- select_variable(result, selected_communities, "G")
G_trait2 <- select_variable(result, NULL, "G")
G_trait$community <- NULL
G_trait2$community <- NULL
# A
A_trait <- select_variable(result, selected_communities, "A")
A_trait2 <- select_variable(result, NULL, "A")
A_trait$community <- NULL
A_trait2$community <- NULL
# M
M_trait <- select_variable(result, selected_communities, "M")
M_trait2 <- select_variable(result, NULL, "M")
M_trait$community <- NULL
M_trait2$community <- NULL


#Final_table
ft <- result
ft <- ft %>% column_to_rownames(var = "community")
ft2 <- result %>% filter(community %in% selected_communities)
ft2 <- ft2 %>% column_to_rownames(var = "community")

##Consensus##
consensus_tree <-consensus.edges(mphy_cR, method=c("least.squares"), rooted = TRUE)
consensus_tree <- as.phylo(consensus_tree)
plotTree(consensus_tree,type="phylogram", ftype="i",lwd=1)

###Estimar estados ancestrales
# S
sv1 <- as.matrix(S_trait)[,1]
sv2 <- as.matrix(S_trait2)[,1]

# R
rv1 <- as.matrix(R_trait)[,1]
rv2 <- as.matrix(R_trait2)[,1]

# G
gv1 <- as.matrix(G_trait)[,1]
gv2 <- as.matrix(G_trait2)[,1]

# A
av1 <- as.matrix(A_trait)[,1]
av2 <- as.matrix(A_trait2)[,1]

# M
mv1 <- as.matrix(M_trait)[,1]
mv2 <- as.matrix(M_trait2)[,1]

# Funci?n para estimar estados ancestrales y crear el ?rbol
estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram", legend = 0.7 * max(nodeHeights(tree)), ftype = "i", leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}

# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
svc <- estimacion_estados_ancestrales(consensus_tree, sv1, "S")
rvc <- estimacion_estados_ancestrales(consensus_tree, rv1, "R")
gvc <- estimacion_estados_ancestrales(consensus_tree, gv1, "G")
avc <- estimacion_estados_ancestrales(consensus_tree, av1, "A")
mvc <- estimacion_estados_ancestrales(consensus_tree, mv1, "M")

#Escribir estados ancestrales en el árbol de consenso
writeAncestors(consensus_tree, Anc=svc$anc, file="Figures/S_nodos_muestra.nex", digits=6, format=c("nexus"))
writeAncestors(consensus_tree, Anc=rvc$anc, file="Figures/R_nodos_muestra.nex", digits=6, format=c("nexus"))
writeAncestors(consensus_tree, Anc=gvc$anc, file="Figures/G_nodos_muestra.nex", digits=6, format=c("nexus"))
writeAncestors(consensus_tree, Anc=avc$anc, file="Figures/A_nodos_muestra.nex", digits=6, format=c("nexus"))
writeAncestors(consensus_tree, Anc=mvc$anc, file="Figures/M_nodos_muestra.nex", digits=6, format=c("nexus"))

estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram", legend = 0.5 * max(nodeHeights(tree)), ftype = "off", leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}
# Estimar estados ancestrales para los ?rboles de apellidos totales
svt <- estimacion_estados_ancestrales(y_total, sv2, "S")
rvt <- estimacion_estados_ancestrales(y_total, rv2, "R")
gvt <- estimacion_estados_ancestrales(y_total, gv2, "G")
avt <- estimacion_estados_ancestrales(y_total, av2, "A")
mvt <- estimacion_estados_ancestrales(y_total, mv2, "M")


###Comparaciones
##Comunidades muestreadas
plot_comparison <- function(obj1, obj2, xlabel) {
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 4) + 1, oma = c(1, 1, 1, 1))
  # First plot
  plot(obj1, lwd = 5, ftype = "off", legend = 0.7 * max(nodeHeights(consensus_tree)), outline = TRUE, fsize = c(0.3), leg.txt = xlabel[1])
  
  # Second plot
  plot(obj2, lwd = 5, outline = TRUE, direction = "leftwards", ftype = "off", legend = 0.7 * max(nodeHeights(consensus_tree)), fsize = c(0.3), leg.txt = xlabel[2])
}

# Comparaciones
png("Figures/S_R_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, rvc$obj, c("S", "R"))
dev.off()
png("Figures/S_G_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, gvc$obj, c("S", "G"))
dev.off()
png("Figures/S_A_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, avc$obj, c("S", "A"))
dev.off()
png("Figures/S_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, mvc$obj, c("S", "M"))
dev.off()
png("Figures/R_G_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(rvc$obj, gvc$obj, c("R", "G"))
dev.off()
png("Figures/R_A_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(rvc$obj, avc$obj, c("R", "A"))
dev.off()
png("Figures/R_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(rvc$obj, mvc$obj, c("R", "M"))
dev.off()
png("Figures/G_A_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(gvc$obj, avc$obj, c("G", "A"))
dev.off()
png("Figures/G_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(gvc$obj, mvc$obj, c("G", "M"))
dev.off()
png("Figures/A_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(avc$obj, mvc$obj, c("A", "M"))
dev.off()


##Comunidades totales
plot_comparison2 <- function(obj1, obj2, xlabel) {
  par(pin=c(10, 10))
  par(mai = c(1.5, 4, 1.5, 4))
  layout(matrix(1:3, 1, 3), widths=c(0.48, 0.04, 0.48)) 
  plot(obj1, lwd=6, ftype="off", outline=TRUE,legend=0.4*max(nodeHeights(y_total)),fsize=c(0.3,1.2),leg.txt=xlabel[1])
  plot.new()
  plot.window(xlim=c(1,2), ylim=get("last_plot.phylo", envir=.PlotPhyloEnv)$y.lim)
  plot(obj2, lwd=6,outline=TRUE, direction="leftwards", legend=0.4*max(nodeHeights(y_total)),ftype="off",fsize=c(0.3,1.2),leg.txt=xlabel[2])
}


# Comparaciones
dev.off()
png("Figures/S_R_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, rvt$obj, c("S", "R"))
dev.off()
png("Figures/S_G_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, gvt$obj, c("S", "G"))
dev.off()
png("Figures/S_A_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, avt$obj, c("S", "A"))
dev.off()
png("Figures/S_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, mvt$obj, c("S", "M"))
dev.off()
png("Figures/R_G_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(rvt$obj, gvt$obj, c("R", "G"))
dev.off()
png("Figures/R_A_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(rvt$obj, avt$obj, c("R", "A"))
dev.off()
png("Figures/R_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(rvt$obj, mvt$obj, c("R", "M"))
dev.off()
png("Figures/G_A_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(gvt$obj, avt$obj, c("G", "A"))
dev.off()
png("Figures/G_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(gvt$obj, mvt$obj, c("G", "M"))
dev.off()
png("Figures/A_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(avt$obj, mvt$obj, c("A", "M"))
dev.off()

###Phylogenetic signal
##Muestra
calculate_physignal_plot <- function(trait_vector, file_name) {
  class(trait_vector) <- "numeric"
  K <- physignal(A = trait_vector, phy = consensus_tree, iter = 1000, print.progress = FALSE)
  png(file_name, width = 800, height = 600)
  plot(K)
  dev.off()
  return(K)
}
sv1
# Vectores de datos
trait_vectors <- list(sv1, rv1, gv1, av1, mv1)
file_names <- c("Figures/Phylosignal_S_muestra2.png", "Figures/Phylosignal_R_muestra2.png",
                "Figures/Phylosignal_G_muestra2.png", "Figures/Phylosignal_A_muestra2.png",
                "Figures/Phylosignal_M_muestra2.png")

# Calcular la se?al f?sica y generar los gr?ficos en un bucle
K_values <- list()
for (i in seq_along(trait_vectors)) {
  K_values[[i]] <- calculate_physignal_plot(trait_vectors[[i]], file_names[i])
}

# Calcular los valores de se?al f?sica y p-value
Phylogenetic_signal <- sapply(K_values, function(K) K$phy.signal)
P_value <- sapply(K_values, function(K) K$pvalue)

# Combinar resultados en un marco de datos
Phy_sig <- data.frame(Phylogenetic_signal, P_value)
row.names(Phy_sig) <- c("S", "R", "G", "A", "M")

# Redondear los resultados
Phy_sig <- round(Phy_sig, digits = 4)
publish(Phy_sig)
# Generar la tabla y guardarla como imagen
png("Figures/Phylosignal_muestra2.png", width = 300, height = 200)
grid.table(Phy_sig)
dev.off()

##Total
calculate_physignal_plot <- function(trait_vector, file_name) {
  class(trait_vector) <- "numeric"
  K <- physignal(A = trait_vector, phy = y_total, iter = 1000, print.progress = FALSE)
  png(file_name, width = 800, height = 600)
  plot(K)
  dev.off()
  return(K)
}
# Vectores de datos
trait_vectors2 <- list(sv2, rv2, gv2, av2, mv2)
file_names <- c("Figures/Phylosignal_S_total.png", "Figures/Phylosignal_R_total.png",
                "Figures/Phylosignal_G_total.png", "Figures/Phylosignal_A_total.png",
                "Figures/Phylosignal_M_total.png")

# Calcular la se?al f?sica y generar los gr?ficos en un bucle
K_values <- list()
for (i in seq_along(trait_vectors2)) {
  K_values[[i]] <- calculate_physignal_plot(trait_vectors2[[i]], file_names[i])
}

# Calcular los valores de se?al f?sica y p-value
Phylogenetic_signal2 <- sapply(K_values, function(K) K$phy.signal)
P_value2 <- sapply(K_values, function(K) K$pvalue)

# Combinar resultados en un marco de datos
Phy_sig2 <- data.frame(Phylogenetic_signal2, P_value2)
row.names(Phy_sig2) <- c("S", "R", "G", "A", "M")

# Redondear los resultados
Phy_sig2 <- round(Phy_sig2, digits = 4)
publish(Phy_sig2)
# Generar la tabla y guardarla como imagen
png("Figures/Phylosignal_total.png", width = 300, height = 200)
grid.table(Phy_sig2)
dev.off()


###Regresi?n PGLS
##Phylo-regression for sampled data
sampled_matched <- treedata(consensus_tree, ft2, sort=F, warnings=TRUE)
spc <- sampled_matched$phy$tip.label
V<-corBrownian(phy=sampled_matched$phy,form = ~spc)
C <- vcv.phylo(phy = sampled_matched$phy)
#Assign traits
obj <- ft2
x <- write.tree(consensus_tree)
tree_consensus.tree <- read.tree(text=x)
G <- setNames(obj[,"G"],rownames(obj))
M <- setNames(obj[,"M"],rownames(obj))
S <- setNames(obj[,"S"],rownames(obj))
A <- setNames(obj[,"A"],rownames(obj))
R <- setNames(obj[,"R"],rownames(obj))
N <- setNames(obj[,"N"],rownames(obj))

# Verifica que N, M, A, S, y G sean vectores numéricos
N <- as.vector(N)
M <- as.vector(M)
A <- as.vector(A)
S <- as.vector(S)
G <- as.vector(G)

# Crea un dataframe con estos vectores
vector.data <- data.frame(N = N, M = M, A = A, S = S, G = G)

#model with N ( efecto del N en el indice)
bm_glsN<-gls(N~M+A+S+G,correlation = V, data= vector.data)
summary(bm_glsN) 
anova(bm_glsN)

#model with G 
bm_glsG<-gls(G~M+A+S,correlation = V, data=data.frame(N, M, A, S, G))
summary(bm_glsG) 
anova(bm_glsG)

#model with M 
bm_glsM<-gls(M~G+A+S,correlation = V, data=data.frame(N, M, A, S, G))
summary(bm_glsM) 
anova(bm_glsM)

publish(bm_glsN)
publish(bm_glsG)
publish(bm_glsM)

##Phylo-regression for all communities
sampled_matched <- treedata(y_total, ft, sort=FALSE, warnings=TRUE)
spc <- sampled_matched$phy$tip.label
V<-corBrownian(phy=sampled_matched$phy,form = ~spc)
C <- vcv.phylo(phy = sampled_matched$phy)

#Now we run the analysis:
#Assign traits
obj <- ft
x <- write.tree(y_total)
tree_consensus.tree <- read.tree(text=x)
G <- setNames(obj[,"G"],rownames(obj))
M <- setNames(obj[,"M"],rownames(obj))
S <- setNames(obj[,"S"],rownames(obj))
A <- setNames(obj[,"A"],rownames(obj))
R <- setNames(obj[,"R"],rownames(obj))
N <- setNames(obj[,"N"],rownames(obj))

#model with N ( efecto del N en el ?ndice)
bm_glsN<-gls(N~M+A+S+G,correlation = V, data=data.frame(N, M, A, S, G))
summary(bm_glsN) 
anova(bm_glsN)

#model with G 
bm_glsG<-gls(G~M+A+S,correlation = V, data=data.frame(N, M, A, S, G))
summary(bm_glsG) 
anova(bm_glsG)

#model with M 
bm_glsM<-gls(M~G+A+S,correlation = V, data=data.frame(N, M, A, S, G))
summary(bm_glsM) 
anova(bm_glsM)

publish(bm_glsN)
publish(bm_glsG)
publish(bm_glsM)


##Generate PICs and test while conditioning on phylogeny
#Sampled communities
S_pic1<-pic(x = ft2$S, phy = consensus_tree)
R_pic1<-pic(x= ft2$R, phy = consensus_tree)
A_pic1<-pic(x = ft2$A, phy = consensus_tree)
G_pic1<-pic(x=ft2$G, phy = consensus_tree)
M_pic1<-pic(x=ft2$M, phy = consensus_tree)

calc_r <- function(x, y) {
  r <- cor(x, y)  # Calcular el coeficiente de correlaci?n de Pearson
  return(r)
}
data1 <- data.frame(S_pic1,A_pic1,R_pic1,G_pic1,M_pic1)

#Total communities
S_pic<-pic(x = ft$S, phy = y_total)
R_pic<-pic(x= ft$R, phy = y_total)
A_pic<-pic(x = ft$A, phy = y_total)
G_pic<-pic(x=ft$G, phy = y_total)
M_pic<-pic(x=ft$M, phy = y_total)

calc_r <- function(x, y) {
  r <- cor(x, y)  # Calcular el coeficiente de correlaci?n de Pearson
  return(r)
}
data <- data.frame(S_pic,A_pic,R_pic,G_pic,M_pic)
