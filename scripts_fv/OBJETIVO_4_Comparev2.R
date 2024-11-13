################################################################################
################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth#### ##########inheritance and land tenure norms in agropastoralist communities######
################################################################################
################################################################################

#### OBJETIVO 4 ################################################################
### To compare trees builds with different data and their deviations from the consensus tree #################################################################

##########IMPORTANTE: CORRER LOS SCRIPTS DE LOS OBJETIVOS 1 Y 3#################

###LIBRERIAS
library(ape)
library(dendextend)
library(geiger)
library(phytools)

################ TANGLEGRAMA:COMPARACION APELLIDOS/STR #########################
# Etiquetas están en el mismo orden
hy <- reorder(hy, "postorder")
plotTree(hy)
phyDSW <- reorder(phyDSW,"postorder")
plotTree(phyDSW)
# Escalonar las ramas en el mismo sentido (izquierda o derecha)
hy <- ape::ladderize(hy)
phyDSW <- ape::ladderize(phyDSW)

## Tanglegrama
# Hacer dos árboles en dónde las poblaciones sean números
hy_num <- hy
hy_num$tip.label <- seq_along(hy$tip.label)
phyDSW_num <- phyDSW
phyDSW_num$tip.label <- seq_along(phyDSW$tip.label)
# Crear correspondencia entre número y comunidad
comunidades <- hy$tip.label 
# Crear leyenda como un vector de caracteres
leyenda <- paste(seq_along(comunidades), "  ", comunidades)
#Entanglement
e <- dendlist(as.dendrogram(hy_num), as.dendrogram(phyDSW_num)) %>%
  dendextend::untangle(method = "step1side") %>% # Find the best alignment layout
  entanglement()                     # Alignment?quality
#PNG
png("Figures/Tanglegram.png",width = 800, height = 400, units = "px", pointsize = 12,bg = "white")
dendlist(as.dendrogram(hy_num),as.dendrogram(phyDSW_num))%>%
  dendextend::untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(highlight_distinct_edges = FALSE, # Turn-off dashed lines,
             main = paste("entanglement =", round(e,4)),
             main_left = "Surnames",
             main_right = "STR",
             sort = F,
             edge.lwd = T,
             color_lines = T,
             intersecting = T,
             rank_branches = T,
             common_subtrees_color_branches = T
  ) 
# Agregar leyenda en el gráfico
legend(x = 5.8, y = 7.5,legend = leyenda, cex = 0.3, pt.cex = 0.5,inset=c(0.05,0.05),text.width = 0.8,box.lty = 2)

dev.off() 

###################### ARBOL DE CONSENSO #######################################
# Combina los árboles en una lista de clase "multiPhylo"
#combined_trees1 <- as.multiPhylo(hy,phyDSW)
combined_trees2 <- as.multiPhylo(phyDSW,hy)
#combined_trees3 <- as.multiPhylo(hy2,phyDSW2)
#combined_trees4 <- as.multiPhylo(phyDSW2,hy2)

# Genera el árbol de consenso
#consensus_tree1 <- consensus(combined_trees1,p=0.5,check.labels = TRUE, rooted =F)
#consensus_tree1 <- reorder(consensus_tree1, "postorder")
consensus_tree2 <- consensus(combined_trees2,p=0.5,check.labels = TRUE, rooted =F)
consensus_tree2 <- reorder(consensus_tree2, "postorder")
#consensus_tree3 <- consensus(combined_trees3,p=0.5,check.labels = TRUE, rooted =F)
#consensus_tree3 <- reorder(consensus_tree3, "postorder")
#consensus_tree4 <- consensus(combined_trees4,p=0.5,check.labels = TRUE, rooted =F)
#consensus_tree4 <- reorder(consensus_tree4, "postorder")
#plotTree(consensus_tree1)
plotTree(consensus_tree2)
#plotTree(consensus_tree3)
#plotTree(consensus_tree4)

#par(mfrow=c(1,2))
#plot(consensus_tree1, main="Apellidos primero, s/PUCLARO")
#plot(consensus_tree2, main="DSW primero, s/PUCLARO")
#par(mfrow=c(1,1))
#dev.off()
#par(mfrow=c(1,2))
#plot(consensus_tree3, main="Apellidos primero, c/PUCLARO")
#plot(consensus_tree4, main="DSW primero, c/PUCLARO")
#par(mfrow=c(1,1))
#dev.off()


png("Figures/Consensus_tree.png",width = 800, height = 800, units = "px", pointsize = 12,bg = "white")
plot.phylo(consensus_tree2)
dev.off()

#Preparacion de árbol con largo de ramas (Todavía está por verse)
largo.ramas <- ((hy$edge.length) + (phyDSW$edge.length))/2
consensus_tree2$edge.length <- largo.ramas
micalibracion <- makeChronosCalib(consensus_tree2)
mytimetree <- chronos(consensus_tree2, lambda = 1, model = "discrete", 
                      calibration = micalibracion, 
                      control = chronos.control(nb.rate.cat = 1))
plot.phylo(mytimetree)
consensus_tree<-mytimetree

###################### DENDROPLOT CONSENSO & TRAITS ############################
## Generar dendroplot con el ?rbol de consenso y los traits anotados
traits <- function(comuneros, group_by_cols = c("community","commune")) {
  # Asegurarse de que group_by_cols es un vector
  if (!is.vector(group_by_cols)) {
    group_by_cols <- as.vector(group_by_cols)
  }
  
  # Calcular los �ndices
  result <- comuneros %>%
    group_by(across(all_of(group_by_cols))) %>%
    summarise(
      N = n(),
      S = n_distinct(surname_father) / N,
      R = mean(rights, na.rm = TRUE),
      G = gini(shares),
      A = mean(sex == "M", na.rm = TRUE),
      M = sum(rights < 1, na.rm = TRUE) / N,
    )
  
  return(result)
}

result <-traits(comuneros) 
consensus_tree <- as.dendrogram(consensus_tree)
plot(consensus_tree)
comuneros$commune[comuneros$commune == "VICUÃ‘A"] <-"VICUÑA"
select_comuneros <- comuneros %>% dplyr::filter(comuneros$community %in% selected_communities)

dendroplot <- function(consensus_tree, save_as = NULL, group_by_col = "community") {
  library(ggdendro)
  library(ggplot2)
  
  # Comprobar si el objeto es un dendrograma
  if (!inherits(consensus_tree, "dendrogram")) {
    stop("consensus_tree debe ser un objeto de tipo dendrogram")
  }
  
  # Generar los datos del dendrograma
  hcd <- dendro_data(consensus_tree, type = "rectangle")
  # Obtener las etiquetas como caracteres
  hcd$labels$label <- as.character(hcd$labels$label)
  
  # Definir el container
  container <- if (group_by_col == "community") "commune"
  else if (group_by_col == "commune") "province"
  else if (group_by_col == "province") "region"
  else NULL
  
  if (is.null(container)) {
    stop("Container is NULL")
  }
  
  # Calcular traits
  tc <- traits(select_comuneros, c(group_by_col, container))
  
  # Funci?n auxiliar para obtener vectores
  vector_of <- function(target_col) {
    if (!target_col %in% colnames(tc)) {
      stop(paste("Column", target_col, "is not found in traits data"))
    }
    v <- tc[[target_col]]
    if (is.null(v)) {
      stop(paste("Column", target_col, "is NULL in traits data"))
    }
    names(v) <- tc[[group_by_col]]
    v
  }
  
  location <- if (!is.null(container)) vector_of(container) else NULL
  N <- vector_of("N")
  S <- vector_of("S")
  R <- vector_of("R")
  G <- vector_of("G")
  A <- vector_of("A")
  M <- vector_of("M")
  
  # Coordenadas ?tiles
  lastrow <- nrow(hcd$labels)
  x0 <- hcd$labels$x[[lastrow]]
  y0 <- hcd$labels$y[[lastrow]]
  x1 <- x0 + 1 + 0.5 * lastrow / 170
  ydiff <- if (is.null(container)) 0.4 else 0
  
  # Verificar y crear el directorio si no existe
  if (!is.null(save_as)) {
    dir_path <- dirname(save_as)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    png(filename = save_as, width = 10 + 4 * lastrow / 170, height = 1 + 40 * lastrow / 170, units = "in", res = 300)
  }
  
  size <- function(xs) {
    xs * 1.3 / max(xs) + (1.7 + lastrow / 170)
  }
  #Ajustar tama?o de las ramas
  scale_factor <- 0.2  # Ajusta este valor para cambiar el tama?o de las ramas
  hcd$segments$y <- hcd$segments$y * scale_factor
  hcd$segments$yend <- hcd$segments$yend * scale_factor
  
  # Graficar
  p <- ggplot() +
    geom_segment(data = segment(hcd), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = label, hjust = 0), nudge_y = 0.01, nudge_x = 0.05, size = 2) +
    annotate("text", x = x1, y = y0 - 0.215, label = str_to_title(group_by_col), fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 1.84 + ydiff, label = "#", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.06 + ydiff, label = "S", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.26 + ydiff, label = "R", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.46 + ydiff, label = "G", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.66 + ydiff, label = "A", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.86 + ydiff, label = "M", fontface = "bold", size = 4) +
    scale_size_identity() +
    coord_flip() +
    scale_y_reverse(expand = c(0.2, 0)) +
    theme_dendro()
  
  # Agregar condicionalmente el container
  if (!is.null(container)) {
    p <- p +
      annotate("text", x = x1, y = y0 - 1.265, label = str_to_title(container), fontface = "bold", size = 4) +
      geom_text(data = label(hcd), aes(x = x, y = y, hjust = 0, label = location[label], colour = location[label]), nudge_y = 1.1, size = 3, show.legend = FALSE)
  }
  
  # Agregar el resto de las anotaciones
  p <- p +
    geom_text(data = label(hcd), aes(x = x, y = y, label = N[label]), nudge_y = 1.8 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(S[label], format = "f", digits = 3)), nudge_y = 2.0 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(R[label], format = "f", digits = 3)), nudge_y = 2.2 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(G[label], format = "f", digits = 3)), nudge_y = 2.4 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(A[label], format = "f", digits = 3)), nudge_y = 2.6 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(M[label], format = "f", digits = 3)), nudge_y = 2.8 - ydiff, hjust = 0, size = 3)
  
  # Guardar el gr?fico
  if (!is.null(save_as)) {
    print(p)
    dev.off()
  } else {
    p
  }
}

# Llamar a la funci?n consensus_dendrogram
library(conflicted)
conflict_prefer("theme_dendro","ggdendro")
conflict_prefer("label", "ggdendro")

consensus_dendrogram <- function(select_comuneros, save_as=NULL,
                               hclust_method="complete",
                               group_by_col="community") {
  hc_total <- surname_clustering(select_comuneros, hclust_method, group_by_col)
  dendroplot(consensus_tree, save_as, group_by_col)
}

consensus_dendrogram(select_comuneros, save_as = "Figures/dendrograma_consenso.png")
