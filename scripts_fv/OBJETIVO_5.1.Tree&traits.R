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
library(geomorph) #No lo puedo instalar

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
#S_trait <- select_variable(result_traits, selected_communities, "S")
#S_trait2 <- select_variable(result_traits, NULL, "S")
#S_trait$community <- NULL
#S_trait2$community <- NULL
S_logit <- select_variable(result_traits, selected_communities, "S_logit")
S_logit2 <- select_variable(result_traits, NULL, "S_logit")
S_logit$community <- NULL
S_logit2$community <- NULL

#N_trait2 <- select_variable(result_traits, selected_communities, "N")

# G
#G_trait <- select_variable(result_traits, selected_communities, "G")
#G_trait2 <- select_variable(result_traits, NULL, "G")
#G_trait$community <- NULL
#G_trait2$community <- NULL
G_logit <- select_variable(result_traits, selected_communities, "G_logit")
G_logit2 <- select_variable(result_traits, NULL, "G_logit")
G_logit$community <- NULL
G_logit2$community <- NULL
# A
#A_trait <- select_variable(result_traits, selected_communities, "A")
#A_trait2 <- select_variable(result_traits, NULL, "A")
#A_trait$community <- NULL
#A_trait2$community <- NULL
A_logit <- select_variable(result_traits, selected_communities, "A_logit")
A_logit2 <- select_variable(result_traits, NULL, "A_logit")
A_logit$community <- NULL
A_logit2$community <- NULL
# M
#M_trait <- select_variable(result_traits, selected_communities, "M")
#M_trait2 <- select_variable(result_traits, NULL, "M")
#M_trait$community <- NULL
#M_trait2$community <- NULL
M_logit <- select_variable(result_traits, selected_communities, "M_logit")
M_logit2 <- select_variable(result_traits, NULL, "M_logit")
M_logit$community <- NULL
M_logit2$community <- NULL
#Cluster
#C_trait <- select_variable(result_traits, selected_communities, "cluster")
#C_trait2 <- select_variable(result_traits, NULL, "cluster")
#C_trait$community <- NULL
#C_trait2$community <- NULL

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
sv1 <- as.matrix(S_logit)[,1]
sv2 <- as.matrix(S_logit2)[,1]

# G
gv1 <- as.matrix(G_logit)[,1]
gv2 <- as.matrix(G_logit2)[,1]

# A
av1 <- as.matrix(A_logit)[,1]
av2 <- as.matrix(A_logit2)[,1]

# M
mv1 <- as.matrix(M_logit)[,1]
mv2 <- as.matrix(M_logit2)[,1]

# Cluster
#cv1 <- as.matrix(C_trait)[,1]
#cv2 <- as.matrix(C_trait2)[,1]

# Funci?n para estimar estados ancestrales y crear el ?rbol
estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram", legend = 0.7 * max(nodeHeights(tree)), ftype = "i", leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}

# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("outputs/Figures/S_muestra.png")
svc <- estimacion_estados_ancestrales(consensus_tree, sv1, "S")
dev.off()

#png("outputs/Figures/R_muestra.png")
#rvc <- estimacion_estados_ancestrales(consensus_tree, rv1, "R")
#dev.off()

png("outputs/Figures/G_muestra.png")
gvc <- estimacion_estados_ancestrales(consensus_tree, gv1, "G")
dev.off()

png("outputs/Figures/A_muestra.png")
avc <- estimacion_estados_ancestrales(consensus_tree, av1, "A")
dev.off()

png("outputs/Figures/M_muestra.png")
mvc <- estimacion_estados_ancestrales(consensus_tree, mv1, "M")
dev.off()

#png("outputs/Figures/C_muestra.png")
#cvc<- estimacion_estados_ancestrales(consensus_tree,cv1,"C")
#dev.off()

# COMO LO QUE VIENE DEPENDE DE LO ANTERIOR, NO LO SEGUÍ PROBANDO##

writeAncestors(consensus_tree, Anc=svc$anc, file="outputs/Figures/S_nodos_muestra.phy", digits=6, format=c("phylip"))
#writeAncestors(consensus_tree, Anc=rvc$anc, file="outputs/Figures/R_nodos_muestra.phy", digits=6, format=c("phylip"))
writeAncestors(consensus_tree, Anc=gvc$anc, file="outputs/Figures/G_nodos_muestra.phy", digits=6, format=c("phylip"))
writeAncestors(consensus_tree, Anc=avc$anc, file="outputs/Figures/A_nodos_muestra.phy", digits=6, format=c("phylip"))
writeAncestors(consensus_tree, Anc=mvc$anc, file="outputs/Figures/M_nodos_muestra.phy", digits=6, format=c("phylip"))

#Cluster discrete visualization
#trees<-make.simmap(consensus_tree,cv1,nsim=100)
#obj<-summary(trees,plot=FALSE)
#cols<-setNames(palette()[1:400],mapped.states(trees))
#cols<-cols[1:4]
#plot(obj,cols,type="phylogram",fsize=0.8,cex=c(0.5,0.3))
#add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
#                  y=0.9*par()$usr[4],prompt=FALSE,fsize=0.9)

estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram", legend = 0.5 * max(nodeHeights(tree)), ftype = "off", leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}
# Estimar estados ancestrales para los ?rboles de apellidos totales
png("outputs/Figures/S_total.png")
svt <- estimacion_estados_ancestrales(y_total, sv2, "S")
dev.off()
#png("outputs/Figures/R_total.png")
#rvt <- estimacion_estados_ancestrales(y_total, rv2, "R")
#dev.off()
png("outputs/Figures/G_total.png")
gvt <- estimacion_estados_ancestrales(y_total, gv2, "G")
dev.off()
png("outputs/Figures/A_total.png")
avt <- estimacion_estados_ancestrales(y_total, av2, "A")
dev.off()
png("outputs/Figures/M_total.png")
mvt <- estimacion_estados_ancestrales(y_total, mv2, "M")
dev.off()
#png("outputs/Figures/C_total.png")
#cvt <- estimacion_estados_ancestrales(y_total, cv2, "C")
#dev.off()

#Cluster discrete visualization
#trees<-make.simmap(y_total,cv2,nsim=100)
#obj<-summary(trees,plot=FALSE)
#cols<-setNames(palette()[1:400],mapped.states(trees))
#cols<-cols[1:4]
#plot(obj,cols,type="phylogram",fsize=0.8,cex=c(0.5,0.3))
#add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
#                  y=0.9*par()$usr[4],prompt=FALSE,fsize=0.9)

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
#png("outputs/Figures/S_R_muestra.png",width = 1042, height = 534, res = 300)
#plot_comparison(svc$obj, rvc$obj, c("S", "R"))
#dev.off()
png("outputs/Figures/S_A_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, avc$obj, c("S", "A"))
dev.off()
png("outputs/Figures/S_G_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, gvc$obj, c("S", "G"))
dev.off()
png("outputs/Figures/S_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(svc$obj, mvc$obj, c("S", "M"))
dev.off()
#png("outputs/Figures/R_A_muestra.png",width = 1042, height = 534, res = 300)
#plot_comparison(rvc$obj, avc$obj, c("R", "A"))
#dev.off()
#png("outputs/Figures/R_G_muestra.png",width = 1042, height = 534, res = 300)
#plot_comparison(rvc$obj, gvc$obj, c("R", "G"))
#dev.off()
#png("outputs/Figures/R_M_muestra.png",width = 1042, height = 534, res = 300)
#plot_comparison(rvc$obj, mvc$obj, c("R", "M"))
#dev.off()
png("outputs/Figures/A_G_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(avc$obj, gvc$obj, c("G", "A"))
dev.off()
png("outputs/Figures/A_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(avc$obj, mvc$obj, c("A", "M"))
dev.off()
png("outputs/Figures/G_M_muestra.png",width = 1042, height = 534, res = 300)
plot_comparison(gvc$obj, mvc$obj, c("G", "M"))
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
#png("outputs/Figures/S_R_total.png",width = 2084, height = 1068, res = 300)
#plot_comparison2(svt$obj, rvt$obj, c("S", "R"))
#dev.off()
png("outputs/Figures/S_G_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, gvt$obj, c("S", "G"))
dev.off()
png("outputs/Figures/S_A_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, avt$obj, c("S", "A"))
dev.off()
png("outputs/Figures/S_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(svt$obj, mvt$obj, c("S", "M"))
dev.off()
#png("outputs/Figures/R_G_total.png",width = 2084, height = 1068, res = 300)
#plot_comparison2(rvt$obj, gvt$obj, c("R", "G"))
#dev.off()
#png("outputs/Figures/R_A_total.png",width = 2084, height = 1068, res = 300)
#plot_comparison2(rvt$obj, avt$obj, c("R", "A"))
#dev.off()
#png("outputs/Figures/R_M_total.png",width = 2084, height = 1068, res = 300)
#plot_comparison2(rvt$obj, mvt$obj, c("R", "M"))
#dev.off()
png("outputs/Figures/A_G_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(avt$obj,gvt$obj, c("G", "A"))
dev.off()
png("outputs/Figures/A_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(avt$obj, mvt$obj, c("A", "M"))
dev.off()
png("outputs/Figures/G_M_total.png",width = 2084, height = 1068, res = 300)
plot_comparison2(gvt$obj, mvt$obj, c("G", "M"))
dev.off()

###Phylogenetic signal
### Cambiar por physignal.z? 
z<-phylosig(consensus_tree,sv1,method = "lambda",test=T)
z

##Muestra
analyze_phylosignal <- function(tree, trait) {
  # 1. Verificar correspondencia entre árbol y datos
  if(!all(names(trait) %in% tree$tip.label)) {
    stop("Los nombres en el trait no coinciden con las puntas del árbol")
  }
  
  # 2. Asegurar que el trait es numérico
  trait_num <- setNames(as.numeric(trait), names(trait))
  
  # 3. Calcular lambda de Pagel con más simulaciones
  res <- phylosig(tree, trait_num, method = "lambda", test = TRUE, nsim = 9999)
  
  # 4. Calcular K de Blomberg como alternativa
  res_k <- phylosig(tree, trait_num, method = "K", test = TRUE, nsim = 9999)
  
  list(Pagel_lambda = res, Blomberg_K = res_k)
}

# Aplicar a todos los traits
Phy_sigx <- lapply(list(S = sv1, G = gv1, A = av1, M = mv1), 
                  function(t) analyze_phylosignal(consensus_tree, t))  
# Tabla resumen mejorada
Phy_sig <- data.frame(
  Lambda = format(round(sapply(Phy_sigx, function(x) x$Pagel_lambda$lambda), 5), nsmall = 5),
  pvalue_lambda = format(round(sapply(Phy_sigx, function(x) x$Pagel_lambda$P), 5), nsmall = 5),
  Blomberg_K = format(round(sapply(Phy_sigx, function(x) x$Blomberg_K$K), 5), nsmall = 5),
  pvalue_K = format(round(sapply(Phy_sigx, function(x) x$Blomberg_K$P), 5), nsmall = 5)
)
publish(Phy_sig)
png("outputs/Figures/Phylosignal_muestra2.png", width = 400, height = 200)
grid.table(Phy_sig)
dev.off()

##Total
calculate_phylosignal <- function(tree, trait_vector, trait_name = "") {
  # Verificar y ajustar nombres
  trait_numeric <- as.numeric(trait_vector)
  names(trait_numeric) <- names(trait_vector)
  
  # Filtrar para coincidencia exacta entre árbol y datos
  matched_data <- match_tree_and_traits(tree, trait_numeric)
  
  # Calcular métricas
  lambda_res <- phylosig(matched_data$tree, matched_data$trait, 
                         method = "lambda", test = TRUE, nsim = 9999)
  K_res <- phylosig(matched_data$tree, matched_data$trait,
                    method = "K", test = TRUE, nsim = 9999)
  
  list(
    trait = trait_name,
    lambda = lambda_res$lambda,
    lambda_p = lambda_res$P,
    K = K_res$K,
    K_p = K_res$P
  )
}

# Función auxiliar para emparejar árbol y traits
match_tree_and_traits <- function(tree, trait) {
  common_taxa <- intersect(tree$tip.label, names(trait))
  if(length(common_taxa) == 0) stop("No hay coincidencia entre árbol y traits")
  
  list(
    tree = keep.tip(tree, common_taxa),
    trait = trait[common_taxa]
  )
}

# Análisis para todos los traits totales
trait_list_total <- list(
  "S" = sv2,
  "G" = gv2,
  "A" = av2,
  "M" = mv2
)

Phy_sig_t <- lapply(names(trait_list_total), function(trait_name) {
  calculate_phylosignal(y_total, trait_list_total[[trait_name]], trait_name)
})

# Crear tabla de resultados con redondeo a 5 decimales
Phy_sig_total <- do.call(rbind, lapply(Phy_sig_t , function(x) {
  data.frame(
    Trait = x$trait,
    Lambda = round(x$lambda, 5),
    pvalue_Lambda = round(x$lambda_p, 5),
    Blomberg_K = round(x$K, 5),
    pvalue_K = round(x$K_p, 5),
    stringsAsFactors = FALSE
  )
}))
publish(Phy_sig_total)
png("outputs/Figures/Phylosignal_total.png", width = 400, height = 200)
grid.table(Phy_sig_total)
dev.off()

###Regresion PGLS
##Phylo-regression for sampled data
sampled_matched <- treedata(consensus_tree, ft2, sort=F, warnings=TRUE)
spc <- sampled_matched$phy$tip.label
V<-corBrownian(phy=sampled_matched$phy,form = ~spc)
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
G <- setNames(obj[,"G_logit"],rownames(obj))
M <- setNames(obj[,"M_logit"],rownames(obj))
S <- setNames(obj[,"S_logit"],rownames(obj))
A <- setNames(obj[,"A_logit"],rownames(obj))
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
#Prepare the tree
consensus_tree_rooted <- midpoint.root(consensus_tree) 
# Resuelve las politomías
consensus_tree_dicotomous <- multi2di(consensus_tree_rooted) 
#Sampled communities
S_pic1<-pic(x = ft2$S_logit, phy = consensus_tree_dicotomous)
A_pic1<-pic(x = ft2$A_logit, phy = consensus_tree_dicotomous)
G_pic1<-pic(x=ft2$G_logit, phy = consensus_tree_dicotomous)
M_pic1<-pic(x=ft2$M_logit, phy = consensus_tree_dicotomous)

calc_r <- function(x, y) {
  r <- cor.test(x, y)  # Calcular el coeficiente de correlaci?n de Pearson
  r_squared <- r$estimate
  p_value<-r$p.value
  return(list(r_squared = r_squared, p_value = p_value))
}
data1 <- data.frame(S_pic1,A_pic1,G_pic1,M_pic1)

#Total communities
S_pic<-pic(x = ft$S_logit, phy = y_total)
A_pic<-pic(x = ft$A_logit, phy = y_total)
G_pic<-pic(x=ft$G_logit, phy = y_total)
M_pic<-pic(x=ft$M_logit, phy = y_total)

calc_r <- function(x, y) {
  r <- cor.test(x, y)  # Calcular el coeficiente de correlaci?n de Pearson
  r_squared <- r$estimate
  p_value<-r$p.value
  return(list(r_squared = r_squared, p_value = p_value))
}
data <- data.frame(S_pic,A_pic,G_pic,M_pic)

