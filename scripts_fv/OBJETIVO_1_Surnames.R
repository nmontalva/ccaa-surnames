##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### ESPACIO DE TRABAJO ####
#getwd()
#setwd("C:/Users/Kibif/Desktop/Proyecto desigualdad agropastores/Directorio_proyecto")

#### OBJETIVO 1 ####
###  To build a phylogenetic tree showing relationships between communities based on the distributions of surnames within and between communities. ###

## Cargar DATOS ##
comuneros <- read.csv("scripts_fv/Datos/commoners.csv")
comuneros$community <- gsub(" ", "_", comuneros$community)
comuneros$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , comuneros$community)] <- "RINCONADA_DE_PUNITAQUI"

## Cargar paquetes y librerias ##
library(Biodem)
library(reldist)
library(dplyr)
library(ape)
library(ggdendro)
library(ggpubr)
library(phylogram)
library(phytools)

##Metodo por cluster jerarquico de una matriz pairwise. Hedrick's standarized kinship coef.##

##1. TODAS LAS COMUNIDADES
surname_distance_matrix <- function(comuneros,
                                    group_by_col="community") {
  # cross tabulate
  surnames_freq <- table(comuneros$surname_father,
                         comuneros[[group_by_col]])
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # hedrick returns values of similarity
  # transform them into values of dissimilarity (distance)
  as.dist(1-hedkin) #IMPORTANT: This is likely wrong. We figured out a better way.
  #as.dist(Biodem::Fst(hedkin, sum(colSums(surnames_freq) )))
}
surname_matrix <- surname_distance_matrix(comuneros)

hclust_default_method <- "average"
surname_clustering <- function(comuneros,
                               hclust_method=hclust_default_method,
                               group_by_col="community") {
  comuneros %>%
    surname_distance_matrix(group_by_col) %>%
    # hierarchical clustering of the distance matrix
    hclust(method=hclust_method)
}
hc_total <- surname_clustering(comuneros)
plot(hc_total)
y_total <- as.phylo(hc_total) #Phylo format

#Personalizar ?rboles

plot(hc_total) #Dendrograma
ddata <- dendro_data(hc_total, type = "rectangle")
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0))+
  theme_dendro()
p

##2. COOMUNIDADES MUESTREADAS
# Seleccionar los datos que se necesitan
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
selected_communities <- unique(STR$pop)

# Funci?n de creaci?n de matriz y ?rbol
surnames <- comuneros %>% filter(community %in% selected_communities)
surnames$community <- factor(surnames$community, levels = selected_communities)
surname_distance_muestra <- function(surnames,
                                    group_by_col="community") {
  # cross tabulate
  surnames_freq <- table(surnames$surname_father,
                         surnames[[group_by_col]])
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # hedrick returns values of similarity
  # transform them into values of dissimilarity (distance)
  as.dist(1-hedkin) #IMPORTANT: This is likely wrong. We figured out a better way.
  #as.dist(Biodem::Fst(hedkin, sum(colSums(surnames_freq) )))
}
surname_matrix_muestra <- surname_distance_muestra(surnames)

hclust_default_method <- "average"
surname_clustering <- function(surnames,
                               hclust_method=hclust_default_method,
                               group_by_col="community") {
  surnames %>%
    surname_distance_muestra(group_by_col) %>%
    # hierarchical clustering of the distance matrix
    hclust(method=hclust_method)
}
hc <- surname_clustering(surnames)
plot(hc)
hd <- as.dendrogram(hc)
hy <- as.phylo(hc) #Phylo format
plotTree(hy)


#### Guardar dendrogramas
write.dendrogram(as.dendrogram(hc_total), file = "Figures/Apellidos.phy", edges = FALSE) #Guardar dendrograma archivo phy
write.nexus(y_total, file = "Figures/Apellidos.nex", translate = TRUE) #Guardar archivo nexus desde phy
write.dendrogram(hd,file = "Figures/Apellidos_muestra.phy", edges = FALSE) #Guardar dendrograma archivo phy
write.nexus(hy, file = "Figures/Apellidos_muestra.nex", translate = TRUE) #Guardar archivo nexus desde phy

