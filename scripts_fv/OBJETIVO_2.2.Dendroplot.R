##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 2 ####
### To estimate the traits of surnames' diversity, concentration of commoners' rights and inheritance's agnatic bias for each community based on the distributions of surnames within communities ###

## Cargar paquetes y librerias ##
library(dendextend)
library(ggplot2)
library(gridExtra)
library(stringr)

##Conflict preference
conflict_prefer("label","ggdendro")
conflict_prefer("theme_dendro","ggdendro")
## ADVERTENCIA: NO CORRER LO SIGUIENTE SI QUIERE CORRER LOS SCRIPTS DEL OBJETIVO_5 ##

# Mostrar y escribir la tabla final 
result_dendro<-result
row.names(result_dendro) <-result_dendro$community
result_dendro <- result_dendro%>%select(-community)
write.table(result_dendro, file='Figures/Tabla_indices.txt', sep = '\t', row.names = T, quote = FALSE)

# Escribir la tabla de las comunidades muestreadas
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
selected_communities <- unique(STR$pop)
result_dendro2 <- result_dendro %>% dplyr::filter(row.names(result_dendro) %in% selected_communities)
write.table(result_dendro, file='Figures/Tabla_indices.txt', sep = '\t', row.names = T, quote = FALSE)

# Imprimir im?gen en png y en pdf
png("Figures/Tabla_indices.png", height=4000, width=1500)
p<-tableGrob(result_dendro)
grid.arrange(p)
dev.off()

pdf("Figures/Tabla_indices.pdf", height=60, width=20)
grid.table(result_dendro)
dev.off()


## Traits y comunidades pdf anotados en un arbol
dendroplot <- function(hc, save_as=NULL,
                       group_by_col="community") {
  # generate dendrogram from hclust data
  hcd <- dendro_data(hc, type="rectangle")
  # get rid of those factors
  hcd$labels$label <- as.character(hcd$labels$label)
  # traits
  container <- if (group_by_col == "community") "commune"
  else if (group_by_col == "commune") "province"
  else if (group_by_col == "province") "region"
  # Verificar si container es NULL
  if (is.null(container)) {
    stop("Container is NULL")
  }
  #tc <- traits(comuneros, c(group_by_col, container))
  tc <- result_dendro
  # vectors to obtain commune, S, R, G, A, cluster of communities
  vector_of <- function(target_col) {
    if (!target_col %in% colnames(tc)) {
      stop(paste("Column", target_col, "is not found in traits data"))
    }
    v <- tc[[target_col]]
    if (is.null(v)) {
      stop(paste("Column", target_col, "is NULL in traits data"))
    }
    print(v)  # Agregado para depurar
    names(v) <- tc[[group_by_col]]
    v
  }
  location <- if (!is.null(container))
    vector_of(container)
  else NULL
  N <- vector_of("N")
  S <- vector_of("S")
  A <- vector_of("A")
  G <- vector_of("G")
  M <- vector_of("M")
  #cluster <- vector_of("cluster")
  # useful coordinates
  lastrow <- nrow(hcd$labels)
  x0 <- hcd$labels$x[[lastrow]]
  y0 <- hcd$labels$y[[lastrow]]
  x1 <- x0 + 1 + 0.5 * lastrow / 170
  ydiff <- if (is.null(container)) 0.4 else 0
  # output to pdf
  if (!is.null(save_as)) {
    pdf(save_as,
        width=8 + 4 * lastrow / 170,
        height=1 + 40 * lastrow / 170)
  }
  size <- function(xs) {
    xs * 1.3 / max(xs) + (1.7 + lastrow / 170)
  }
  # plot
  p <- ggplot() +
    geom_segment(data=segment(hcd),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    # because of the coord_flip at the end, x and y are flipped
    geom_text(data=ggdendro::label(hcd),
              aes(x=x, y=y, label=label, hjust=0),
              nudge_y=0.01,
              size=3) +
    annotate("text", x=x1, y=y0-0.215, # y=y0-0.30,
             label=str_to_title(group_by_col),
             fontface="bold", size=4) +
    (if (!is.null(container))
      annotate("text", x=x1, y=y0-1.265, # y=y0-1.30,
               label=str_to_title(container),
               fontface="bold", size=4)) +
    annotate("text", x=x1, y=y0-1.64+ydiff,
             label="#", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-1.86+ydiff,
             label="S", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-2.06+ydiff,
             label="A", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-2.26+ydiff,
             label="G", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-2.46+ydiff,
             label="M", fontface="bold", size=4) +
    #annotate("text", x=x1, y=y0-2.66+ydiff,
    #         label="Cluster", fontface="bold", size=4) +
    (if (!is.null(container))
      geom_text(data=ggdendro::label(hcd),
                aes(x=x, y=y, hjust=0,
                    label=location[label],
                    colour=location[label]),
                nudge_y=1.1,
                size=3,
                show.legend=FALSE)) +
    geom_text(data=ggdendro::label(hcd),
              aes(x=x, y=y,
                  label=N[label],
                  size=size(N[label])),
              nudge_y=1.6-ydiff,
              hjust=0) +
    geom_text(data=ggdendro::label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(S[label], digits=3),
                  size=size(S[label])),
              nudge_y=1.8-ydiff,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(A[label], digits=3),
                  size=size(A[label])),
              nudge_y=2.0-ydiff,
              hjust=0) +
    geom_text(data=ggdendro::label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(G[label], digits=3),
                  size=size(G[label])),
              nudge_y=2.2-ydiff,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(M[label], digits=3),
                  size=size(M[label])),
              nudge_y=2.4-ydiff,
              hjust=0)  +
    #geom_text(data=label(hcd),
    #          aes(x=x, y=y,
    #              label=prettyNum(cluster[label], digits=1),
    #              size=size(cluster[label])),
    #          nudge_y=2.6-ydiff,
    #          hjust=0)  +
    scale_size_identity() +
    coord_flip() +
    scale_y_reverse(expand=c(0, 0.6)) +
    theme_dendro()
  # save plot
  if (!is.null(save_as)) {
    # ggsave(save_as, width=12, height=60, units="cm")
    print(p)
    dev.off()
  } else p
}

surname_dendrogram <- function(comuneros, save_as=NULL,
                               hclust_method=hclust_default_method,
                               group_by_col="community") {
  hc_total <- surname_clustering(comuneros, hclust_method, group_by_col)
  dendroplot(hc_total, save_as, group_by_col)
}

# Llamar a la funci?n surname_dendrogram
surname_dendrogram(comuneros, save_as = "Figures/dendrograma_total.pdf")
