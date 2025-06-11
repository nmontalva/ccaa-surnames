################################################################################
################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth#### 
##########inheritance and land tenure norms in agropastoralist communities######
################################################################################
################################################################################
#### OBJETIVO 4 ################################################################
### To compare trees builds with different data and their deviations from the consensus tree #################################################################

##########IMPORTANTE: CORRER LOS SCRIPTS DE LOS OBJETIVOS 1 Y 3 #################

###LIBRERIAS
library(ape)
library(dendextend)
library(geiger)
library(phangorn)
library(phytools)
library(stringr)

################ TANGLEGRAMA:COMPARACION APELLIDOS/STR #########################
# Etiquetas están en el mismo orden
hy <- reorder(hy, "postorder")
plotTree(hy)
phyDPS <- reorder(phyDPS,"postorder")
plotTree(phyDPS)
# Escalonar las ramas en el mismo sentido (izquierda o derecha)
hy <- ape::ladderize(hy)
phyDPS <- ape::ladderize(phyDPS)

## Tanglegrama
# Hacer dos árboles en dónde las poblaciones sean números
hy_num <- hy
hy_num$tip.label <- seq_along(hy$tip.label)
phyDPS_num <- phyDPS
phyDPS_num$tip.label <- seq_along(phyDPS$tip.label)
# Crear correspondencia entre número y comunidad
comunidades <- hy$tip.label 
# Crear leyenda como un vector de caracteres
leyenda <- paste(seq_along(comunidades), "  ", comunidades)
#Entanglement
e <- dendlist(as.dendrogram(hy_num), as.dendrogram(phyDPS_num)) %>%
  dendextend::untangle(method = "step1side") %>% # Find the best alignment layout
  entanglement()                     # Alignment

dendlist(as.dendrogram(hy_num),as.dendrogram(phyDPS_num))%>% dendextend::untangle(method = "step1side") %>% entanglement # lower is better
dendlist(as.dendrogram(hy_num),as.dendrogram(phyDPS_num))%>% dendextend::untangle(method = "step1side") %>% 
  tanglegram(common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE)
x <- dendlist(as.dendrogram(hy_num),as.dendrogram(phyDPS_num))%>% dendextend::untangle(method = "step2side") 

#PNG DPS
png("Figures/Tanglegram_DPS.png",width = 800, height = 400, units = "px", pointsize = 12,bg = "white")
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)),common_subtrees_color_branches = T, highlight_distinct_edges  = FALSE,
  main_left = "Surnames",
  main_right = "STR",
  sort = F,
  edge.lwd = T,
  color_lines = T,
  intersecting = T,
  rank_branches = T,
)
legend(x = 5.8, y = 7.5,legend = leyenda, cex = 0.3, pt.cex = 0.5,inset=c(0.05,0.05),text.width = 0.8,box.lty = 2)
dev.off()

################ BAKER GAMMA:COMPARACION APELLIDOS/STR #########################
# Encontrar p-valor
set.seed(10000)
the_cor <- cor_bakers_gamma(hy,hy)
the_cor2 <- cor_bakers_gamma(as.dendrogram(phyDPS), as.dendrogram(hy))
R <- 1000
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- hd
for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = F)
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(hd, dend_mixed)
}
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft", legend = c("cor", "cor2"), fill = c(2,1))
sum(the_cor2 < cor_bakers_gamma_results)/ R #p-valor = 0.013
the_cor2 #Baker's gamma correlation coeff = 0.5660764 (Va de -1 a 1, 0 significa que NO son estadisticamente similares)
set.seed(NULL)

###################### ARBOL DE CONSENSO #######################################
# Combina los árboles en una lista de clase "multiPhylo"
combined_trees1 <- as.multiPhylo(hy,phyDPS)
combined_trees2 <- as.multiPhylo(phyDPS,hy)


# Genera el árbol de consenso
consensus_tree1 <- consensus(combined_trees1,p=0.5,check.labels = TRUE, rooted =F)
consensus_tree1 <- reorder(consensus_tree1, "postorder")
consensus_tree2 <- consensus(combined_trees2,p=0.5,check.labels = TRUE, rooted =F)
consensus_tree2 <- reorder(consensus_tree2, "postorder")

plotTree(consensus_tree1)
plotTree(consensus_tree2)

par(mfrow=c(1,2))
plot(consensus_tree1, main="Apellidos primero")
plot(consensus_tree2, main="DPS primero")
par(mfrow=c(1,1))
dev.off()

png("Figures/Consensus_tree.png",width = 800, height = 800, units = "px", pointsize = 12,bg = "white")
plot.phylo(consensus_tree2)
dev.off()

#Preparacion de árbol con largo de ramas (Todavía está por verse)
is.binary(consensus_tree2)
consensus_tree2$root.edge <- 0
# Promediar longitudes de ambos árboles (si topologías son compatibles)
consensus_tree2$edge.length <- (hy$edge.length + phyDPS$edge.length) / 2
consensus_tree2$edge.length[consensus_tree2$edge.length == 0] <- 1e-6  # Reemplazar ceros para evitar la división por 0
consensus_tree2 <- nnls.tree(cophenetic(consensus_tree2), consensus_tree2, rooted = TRUE)
consensus_tree2 <- multi2di(consensus_tree2)
plot.phylo(consensus_tree2)
is.rooted(consensus_tree2)
is.ultrametric(consensus_tree2)
consensus_tree<-(consensus_tree2)

###################### DENDROPLOT CONSENSO & TRAITS ############################
<<<<<<< HEAD
## Generar dendroplot con el ?rbol de consenso y los traits anotados
result_dendro2
=======
## Generar dendroplot con el árbol de consenso y los traits anotados
#TODO
# result_dendro2 no existe (!)
# Se genera en 2.2.
# Pero 2.2. no se puede correr si se quieren correr los del objetivo 5.
# Por eso, tenemos que solucionar esto

>>>>>>> master
consensus_tree <- as.dendrogram(consensus_tree2)

plot(consensus_tree)
#comuneros$commune <-gsub("Ñ", "N", comuneros$commune)
select_comuneros <- comuneros %>% dplyr::filter(comuneros$community %in% selected_communities)

dendroplot <- function(consensus_tree,result_dendro2, save_as = NULL, group_by_col = "community") {
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
  tc <- result_dendro2
  
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
  A <- vector_of("A")
  G <- vector_of("G")
  M <- vector_of("M")
#  cluster <- vector_of("cluster")
  
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
  scale_factor <- 0.5  # Ajusta este valor para cambiar el tama?o de las ramas
  hcd$segments$y <- hcd$segments$y * scale_factor
  hcd$segments$yend <- hcd$segments$yend * scale_factor
  
  # Graficar
  p <- ggplot() +
    geom_segment(data = hcd$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = hcd$labels, aes(x = x, y = y, label = label, hjust = 0), nudge_y = 0.01, nudge_x = 0.05, size = 2) +
    annotate("text", x = x1, y = y0 - 0.215, label = str_to_title(group_by_col), fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 1.84 + ydiff, label = "#", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.06 + ydiff, label = "S", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.26 + ydiff, label = "A", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.46 + ydiff, label = "G", fontface = "bold", size = 4) +
    annotate("text", x = x1, y = y0 - 2.66 + ydiff, label = "M", fontface = "bold", size = 4) +
 #   annotate("text", x = x1, y = y0 - 2.90 + ydiff, label = "Cluster", fontface = "bold", size = 4) +
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
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(A[label], format = "f", digits = 3)), nudge_y = 2.2 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(G[label], format = "f", digits = 3)), nudge_y = 2.4 - ydiff, hjust = 0, size = 3) +
    geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(M[label], format = "f", digits = 3)), nudge_y = 2.6 - ydiff, hjust = 0, size = 3)
  #  geom_text(data = label(hcd), aes(x = x, y = y, label = formatC(cluster[label], format = "f", digits = 0)), nudge_y = 2.85 - ydiff, hjust = 0, size = 3)
  
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

consensus_dendrogram <- function(select_comuneros, save_as=NULL,group_by_col="community") {
  raised_tree <- raise.dendrogram(as.dendrogram(consensus_tree), max(consensus_tree$edge.length))
  dendroplot(raised_tree, save_as, group_by_col)
}

consensus_dendrogram(select_comuneros, save_as = "Figures/dendrograma_consenso_DPS.png")

consensus_tree<-(consensus_tree2)

