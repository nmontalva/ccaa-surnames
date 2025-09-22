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

### Estados basales sampled communities###
load("outputs/ancestral_states_sampled.RData")
load("outputs/ancestral_states_all.RData")
crear_grafico <- function(vc,v1, label,filename) {
  svg(filename,width = 2800, height = 1600, res = 300)
  # Ajustar margen y tama?o de texto para evitar colapso
  par(mar = c(1, 10, 10, 1) + 0.1)
  # Ajustar la separacion entre nodos y otros parametros
  plot(vc$obj,lwd=10, type = "phylogram",offset = 3.2, legend = 0.7 * max(nodeHeights(vc$obj$tree)),outline=FALSE,fsize=0.8, ftype = "reg", leg.txt = label, no.margin = F)
  title(main = paste(label, "ancestral tree"), line = 2)
  nodelabels(text = round(vc$anc$ace, 4), cex = 0.8, bg = "white")
  tip_values <- v1[consensus_tree$tip.label] 
  tiplabels(text = round(tip_values, 4), cex = 0.8, bg = "white", offset = 0.03)
  dev.off()
}
# Llamar a la funcion para cada comunidad muestreada
crear_grafico(svcp,sv1p, "S", "outputs/Figures/S_ancestral_tree_muestra.svg")
crear_grafico(gvcp,gv1p,"G", "outputs/Figures/G_ancestral_tree_muestra.svg")
crear_grafico(avcp,av1p, "A", "outputs/Figures/A_ancestral_tree_muestra.svg")
crear_grafico(mvcp,mv1p, "M", "outputs/Figures/M_ancestral_tree_muestra.svg")

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
crear_grafico2(svtp, sv2p, "S", "outputs/Figures/S_ancestral_tree_total.pdf")
crear_grafico2(gvtp, gv2p, "G", "outputs/Figures/G_ancestral_tree_total.pdf")
crear_grafico2(avtp, av2p, "A", "outputs/Figures/A_ancestral_tree_total.pdf")
crear_grafico2(mvtp, mv2p, "M", "outputs/Figures/M_ancestral_tree_total.pdf")

################################################################################
### CREACION DE POBLACION EN LA BASE (FOSIL) ###
################################################################################
## Hombre,una sola persona, con un solo derecho, un apellido.

##COMUNIDADES MUESTREADAS (consensus_tree)
incorporacion_fosil <- function(fosil, valor, or_tree, label, filename, valor_raiz = 0) {
  # Ajustar ramas de longitud 0
  or_tree$edge.length[or_tree$edge.length == 0] <- 1e-8
  
  # Preparar priors
  pm <- setNames(c(1000, rep(fosil, or_tree$Nnode)),
                 c("sig2", 1:or_tree$Nnode + length(or_tree$tip.label)))
  nn <- as.character(length(or_tree$tip.label) + 1)  # nodo raíz
  pm[nn] <- valor_raiz
  
  pv <- setNames(c(1000^2, rep(1000, length(pm) - 1)), names(pm))
  pv[nn] <- 1e-100
  
  # Ejecutar MCMC
  mcmc <- anc.Bayes(
    or_tree, valor, ngen = 100000,
    control = list(
      pr.mean = pm,
      pr.var = pv,
      a = pm[nn],
      y = pm[as.character(2:or_tree$Nnode + length(or_tree$tip.label))]
    )
  )
  
  mcmc_tree <- mcmc$tree
  
  # Extraer medias posteriores
  w <- as.data.frame(mcmc$mcmc)
  exclude_cols <- c("gen", "sig2", "logLik")
  w <- colMeans(w[, !(colnames(w) %in% intersect(exclude_cols, colnames(w)))])
  
  # Valores de las puntas
  tip_trait <- setNames(valor[mcmc_tree$tip.label], mcmc_tree$tip.label)
  
  # Nodo raíz
  nt <- length(mcmc_tree$tip.label)
  w_nodes <- w[as.character((nt + 1):(nt + mcmc_tree$Nnode))]
  w_nodes[as.character(nt + 1)] <- valor_raiz  # fijar raíz
  
  # Crear contMap sin plot
  obj <- contMap(mcmc_tree, tip_trait, plot = FALSE)
  obj$anc <- w_nodes
  obj$lims <- range(c(tip_trait, w_nodes, valor_raiz), na.rm = TRUE)
  obj <- setMap(obj, viridisLite::viridis(n = 30))
  
  # Guardar svg
  svg(filename, width = 2600, height = 1600, res = 200)

  par(mar = c(1, 1, 4, 1) + 0.1)
  plot(obj, lwd = 10, type = "phylogram", outline = FALSE, offset = 3.2,
       legend = 0.7 * max(nodeHeights(obj$tree)), ftype = "reg", leg.txt = label,
       no.margin = FALSE)
  
  # Dibujar explícitamente la rama que parte de la raíz con color exacto
  # Paleta viridis
  pal <- viridisLite::viridis(n = 100)
  # Normalizar valor de raíz al rango
  root_norm <- (valor_raiz - min(obj$lims)) / (max(obj$lims) - min(obj$lims))
  root_color <- pal[ceiling(root_norm * 99) + 1]
  root_edge <- which(obj$tree$edge[,1] == nt + 1)[1]
  #root_color <- obj$cols[root_edge]
  edge_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  segments(
    x0 = edge_coords$xx[obj$tree$edge[root_edge,1]],
    y0 = edge_coords$yy[obj$tree$edge[root_edge,1]],
    x1 = edge_coords$xx[obj$tree$edge[root_edge,2]],
    y1 = edge_coords$yy[obj$tree$edge[root_edge,2]],
    col = root_color, lwd = 10
  )
  
  # Etiquetas de nodos y puntas
  nodelabels(text = round(w_nodes, 4), cex = 1, bg = "white")
  tiplabels(text = round(tip_trait, 4), cex = 1, bg = "white", offset = 0.02)
  title(main = paste(label, "fossil ancestral tree. Root =", valor_raiz), line = 2)
  dev.off()
}
consensus_tree$edge.length <- consensus_tree$edge.length + 1e-8 # Se le agrega una distancia mínima
incorporacion_fosil(0,sv1p,consensus_tree,"S","outputs/Figures/S_fosil_muestra.svg", valor_raiz = 0) #S
incorporacion_fosil(1,av1p,consensus_tree,"A", "outputs/Figures/A_fosil_muestra.svg",valor_raiz = 1) #A
incorporacion_fosil(0,gv1p,consensus_tree, "G", "outputs/Figures/G_fosil_muestra.svg",valor_raiz = 0) #G
incorporacion_fosil(0,mv1p,consensus_tree, "M", "outputs/Figures/M_fosil_muestra.svg",valor_raiz = 0) #M

##COMUNIDADES TOTALES (y_total)
incorporacion_fosil2 <- function(fosil, valor, tree, label, filename, valor_raiz = 0) {
  # Verificar y corregir ramas de longitud 0
  tree$edge.length[tree$edge.length == 0] <- 1e-6
  
  # Preparar priors
  pm <- setNames(c(1000, rep(fosil, tree$Nnode)),
                 c("sig2", 1:tree$Nnode + length(tree$tip.label)))
  
  # Nodo de la raíz
  nn <- as.character(length(tree$tip.label) + 1)
  pm[nn] <- valor_raiz
  
  # Varianzas previas (muy baja para fijar la raíz)
  pv <- setNames(c(1000^2, rep(1000, length(pm) - 1)), names(pm))
  pv[nn] <- 1e-100
  
  # Ejecutar MCMC para reconstrucción ancestral
  mcmc <- anc.Bayes(
    tree, valor, ngen = 100000,
    control = list(
      pr.mean = pm,
      pr.var = pv,
      a = pm[nn],
      y = pm[as.character(2:tree$Nnode + length(tree$tip.label))]
    )
  )
  
  mcmc_tree <- mcmc$tree
  
  # Extraer medias posteriores de estados ancestrales
  w <- as.data.frame(mcmc$mcmc)
  exclude_cols <- c("gen", "sig2", "logLik")
  existing_cols <- intersect(exclude_cols, colnames(w))
  w <- colMeans(w[, !(colnames(w) %in% existing_cols)])
  
  # Vector nombrado de rasgos para tips
  tip_trait <- setNames(valor[mcmc_tree$tip.label], mcmc_tree$tip.label)
  
  # Identificadores de nodos internos
  nt <- length(mcmc_tree$tip.label)
  nnodes <- mcmc_tree$Nnode
  node_ids <- as.character((nt + 1):(nt + nnodes))
  
  # Filtrar valores ancestrales solo para nodos internos
  w_nodes <- w[node_ids]
  
  # Fijar explícitamente el valor de la raíz
  w_nodes[as.character(nt + 1)] <- valor_raiz
  
  # Crear contMap SIN plotear
  obj <- contMap(mcmc_tree, tip_trait, plot = FALSE)
  
  # Sobrescribir valores ancestrales por los estimados de MCMC
  obj$anc <- w_nodes
  
  # Ajustar límites del gradiente para abarcar puntas + nodos
  obj$lims <- range(c(tip_trait, w_nodes), na.rm = TRUE)
  
  # Elegir paleta de colores
  obj <- setMap(obj, viridisLite::viridis(n = 30))
  
  # Guardar en PDF
  pdf(filename, width = 12, height = 10)
  par(mar = c(4, 4, 4, 4) + 0.1)
  
  # Graficar el árbol recoloreado
  plot(
    obj,
    type = "phylogram",
    lwd = 6,
    offset = 3,
    legend = 0.7 * max(nodeHeights(obj$tree)),
    ftype = "reg",
    leg.txt = label,
    no.margin = FALSE
  )
  
  # Título del gráfico
  title(main = paste(label, "fossil ancestral tree. Root =", valor_raiz), line = 2)
  
  # Etiquetas para nodos y tips
  nodelabels(text = round(w_nodes, 4), cex = 1.2, bg = "lightblue")
  tiplabels(text = round(tip_trait, 4), cex = 1.2, bg = "lightpink", offset = 0.02)
  
  dev.off()
}

incorporacion_fosil2(0,sv2p,y_total,"S", "outputs/Figures/S_fosil_total.pdf") #S: 1 apellido
incorporacion_fosil2(1,av2p,y_total,"A", "outputs/Figures/A_fosil_total.pdf") #A: 1 hombre
incorporacion_fosil2(1,gv2p,y_total,"G", "outputs/Figures/G_fosil_total.pdf") #G: 
incorporacion_fosil2(0,mv2p,y_total,"M", "outputs/Figures/M_fosil_total.pdf") #M
