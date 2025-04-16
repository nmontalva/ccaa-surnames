##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### HIPOTESIS 1 ####
#### OBJETIVO 5 ####
### Part 7 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

###There was a change in wealth inheritance cultural norms, transiting from agnatic primogeniture
###(Spanish majorat) to non partible multigeniture originating the common tenure norm of property of right.
###Predictions of this hypothesis are high values of S and A, and low values for R at the basal internal
###nodes of the tree if agnatic primogeniture is the ancestral state. If that is followed by an early stage of
###non-partible multigeniture, R and S should be both low at the following internal nodes of the tree.

#### ESPACIO DE TRABAJO ####

### Estados basales sampled communities###

crear_grafico <- function(vc,v1, label,filename) {
  png(filename,width = 1400, height = 1000, res = 200)
  # Ajustar margen y tama?o de texto para evitar colapso
  par(mar = c(1, 1, 4, 1) + 0.1)
  
  # Ajustar la separacion entre nodos y otros parametros
  plot(vc$obj, type = "phylogram",offset = 3, legend = 0.7 * max(nodeHeights(vc$obj$tree)), ftype = "reg", leg.txt = label, no.margin = F)
  title(main = paste(label, "ancestral tree"), line = 2)
  nodelabels(text = round(vc$anc$ace, 4), cex = 0.5, bg = "lightblue")
  tip_values <- v1[consensus_tree$tip.label] 
  tiplabels(text = round(tip_values, 4), cex = 0.6, bg = "lightpink", offset = 0.06)
  dev.off()
}
# Llamar a la funci?n para cada comunidad muestreada
crear_grafico(svc,sv1, "S", "Figures/S_ancestral_tree_muestra.png")

#=======
#TODO: REVISAR. error in evaluating the argument 'x' in selecting a method for function 'plot': object 'svc' not found
# SUENA COMO A QUE SVC ES UN METODO QUE PROVIENE DE UN PAQUETE NO INSTALADO (?)
#=======

crear_grafico(gvc,gv1,"G", "Figures/G_ancestral_tree_muestra.png")
crear_grafico(avc,av1, "A", "Figures/A_ancestral_tree_muestra.png")
crear_grafico(rvc,rv1, "R", "Figures/R_ancestral_tree_muestra.png")
crear_grafico(mvc,mv1, "M", "Figures/M_ancestral_tree_muestra.png")


### Estados basales total communities###
crear_grafico2 <- function(vt,v2, label,filename) {
  pdf(filename,
      width=200,
      height=250)
  # Ajustar margen y tama?o de texto para evitar colapso
  par(mar = c(4, 4, 2, 2) + 0.1)
  plot(vt$obj, type = "phylogram",offset = 20, legend = 0.7 * max(nodeHeights(vt$obj$tree)), ftype = "reg",cex=8, leg.txt = label, no.margin = F)
  title(main = paste(label, "ancestral tree"), line = 2)
  nodelabels(text = round(vt$anc$ace, 4), cex = 5, bg = "lightblue")
  tip_values <- v2[y_total$tip.label] 
  tiplabels(text = round(tip_values, 4), cex = 5, bg = "lightpink", offset = 0.005)
  # Cerrar el dispositivo gr?fico PDF
  dev.off()
}

# Llamar a la funcion para el total de comunidades
crear_grafico2(svt, sv2, "S", "Figures/S_ancestral_tree_total.pdf")
crear_grafico2(gvt, gv2, "G", "Figures/G_ancestral_tree_total.pdf")
crear_grafico2(avt, av2, "A", "Figures/A_ancestral_tree_total.pdf")
crear_grafico2(rvt, rv2, "R", "Figures/R_ancestral_tree_total.pdf")
crear_grafico2(mvt, mv2, "M", "Figures/M_ancestral_tree_total.pdf")

### CREACION DE POBLACION EN LA BASE (FOSIL) ###
## Hombre,una sola persona, con un solo derecho, un apellido.

##COMUNIDADES MUESTREADAS (consensus_tree)

incorporacion_fosil <- function(fosil,valor,or_tree,label,filename, valor_raiz = 0)
  {
  #incorporar fosil
  pm<-setNames(c(1000,rep(fosil,or_tree$Nnode)),
               c("sig2",1:or_tree$Nnode+length(or_tree$tip.label)))
  #incorporacion en la raiz
  nn<-as.character((length(or_tree$tip.label)+1))
  pm[nn]<- valor_raiz
  # Varianza previa
  pv <- setNames(c(1000^2, rep(1000, length(pm) - 1)), names(pm))
  pv[as.character(nn)] <- 1e-100
  # Ejecutar MCMC
  mcmc <- anc.Bayes(or_tree, valor, ngen = 100000,
                    control = list(pr.mean = pm, pr.var = pv,
                                   a = pm[nn],
                                   y = pm[as.character(2:or_tree$Nnode + length(or_tree$tip.label))]))
  mcmc_tree <- mcmc$tree
  # Obtener valores ancestrales estimados
  w <- as.data.frame(mcmc$mcmc)
  exclude_cols <- c("gen", "sig2", "logLik")
  existing_cols <- intersect(exclude_cols, colnames(w))
  w <- colMeans(w[, !(colnames(w) %in% existing_cols)])
  # Agregar valor tip.labels
  trait_name <- paste0(label, "_trait2")
  tip.community <- paste0 (mcmc_tree$tip.label)
  name <- get(trait_name)[tip.community,]
  
  # Guardar el grafico como una imagen PNG

  png(filename,width = 1400, height = 1000, res = 200)
  # Ajustar la separacion entre nodos y otros parametros
  par(mar = c(1, 1, 4, 1) + 0.1)
  sorted_trait_vector <- valor[sort(or_tree$tip.label)]
  obj <- contMap(mcmc_tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram",offset = 3, legend = 0.7 * max(nodeHeights(obj$tree)), ftype = "reg", leg.txt = label, no.margin = F)
  title(main = paste(label, "ancestral tree"), line = 2)
  nodelabels(text = round(w, 4), cex = 0.5, bg = "lightblue")
  tip_values <- valor[or_tree$tip.label] 
  tiplabels(text = round(tip_values, 4), cex = 0.6, bg = "lightpink", offset = 0.06)
  #plot(mcmc_tree, show.tip.label = TRUE, cex = 0.6, edge.width = 2, label.offset = 0.15, direction = "rightwards", mar = c(6, 6, 3, 3) + 0.1)
  #title(main = paste(label, "fossil ancestral tree"), line = 2)
  #nodelabels(text = round(w, 4), cex = 0.5, bg = "lightblue")
  #tiplabels(text = round(name, 4), cex = 0.6, bg = "lightpink", offset = 0.005)
  dev.off()
}

incorporacion_fosil(0,sv1,consensus_tree,"S","Figures/S_fosil_muestra.png", valor_raiz = 0) #S
incorporacion_fosil(1,rv1,consensus_tree,"R", "Figures/R_fosil_muestra.png",valor_raiz = 1) #R
incorporacion_fosil(1,av1,consensus_tree,"A", "Figures/A_fosil_muestra.png",valor_raiz = 1) #A
incorporacion_fosil(0,gv1,consensus_tree, "G", "Figures/G_fosil_muestra.png",valor_raiz = 0) #G
incorporacion_fosil(0,mv1,consensus_tree, "M", "Figures/M_fosil_muestra.png",valor_raiz = 0) #M


##COMUNIDADES TOTALES (y_total)
incorporacion_fosil2 <- function(fosil,valor,tree, label, filename)
{
  # Verificar longitudes de ramas
  branch_lengths <- tree$edge.length
  # Ajustar ramas con longitud cero
  min_branch_length <- 1e-6  # Define un valor m?nimo peque?o
  branch_lengths[branch_lengths == 0] <- min_branch_length
  tree$edge.length <- branch_lengths
  #incorporar fosil
  pm<-setNames(c(1000,rep(fosil,tree$Nnode)),
               c("sig2",1:tree$Nnode+length(tree$tip.label)))
  #incorporaci?n en la raiz
  nn<-as.character((length(tree$tip.label)+1))
  pm[nn]<-0
  # Varianza previa
  pv <- setNames(c(1000^2, rep(1000, length(pm) - 1)), names(pm))
  pv[as.character(nn)] <- 1e-100
  # Ejecutar MCMC
  mcmc <- anc.Bayes(tree, valor, ngen = 100000,
                    control = list(pr.mean = pm, pr.var = pv,
                                   a = pm[nn],
                                   y = pm[as.character(2:tree$Nnode + length(tree$tip.label))]))
  mcmc_tree <- mcmc$tree
  # Obtener valores ancestrales estimados
  w <- as.data.frame(mcmc$mcmc)
  exclude_cols <- c("gen", "sig2", "logLik")
  existing_cols <- intersect(exclude_cols, colnames(w))
  w <- colMeans(w[, !(colnames(w) %in% existing_cols)])
  # Agregar valor tip.labels
  trait_name <- paste0(label, "_trait2")
  tip.community <- paste0 (mcmc_tree$tip.label)
  name <- get(trait_name)[tip.community,]
  
  # Guardar el gr?fico como una imagen PNG
  pdf(filename,
      width=200,
      height=250)
  # Ajustar la separaci?n entre nodos y otros par?metros
  par(mar = c(4, 4, 2, 2) + 0.1)
  plot(mcmc_tree, show.tip.label =TRUE, cex = 2, edge.width = 2, label.offset = 0.001, direction = "rightwards", mar = c(6, 6, 3, 3) + 0.1)
  title(main = paste(label, "fossil ancestral tree.",label,"in root =",fosil), line = -15, cex.main = 12)
  nodelabels(text = round(w, 4), cex = 5, bg = "lightblue")
  tiplabels(text = round(name, 4), cex = 5, bg = "lightpink", offset = 0.02)
  dev.off()
}

incorporacion_fosil2(0,sv2,y_total,"S", "Figures/S_f?sil_total.pdf") #S
incorporacion_fosil2(1,rv2,y_total,"R", "Figures/R_f?sil_total.pdf") #R
incorporacion_fosil2(1,av2,y_total,"A", "Figures/A_f?sil_total.pdf") #A
incorporacion_fosil2(1,gv2,y_total,"G", "Figures/G_f?sil_total.pdf") #G
incorporacion_fosil2(0,mv2,y_total,"M", "Figures/M_f?sil_total.pdf") #M
