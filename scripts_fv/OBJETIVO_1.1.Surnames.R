##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

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
  as.dist(1-hedkin)
  #as.dist(1-hedkin) #IMPORTANT: This is likely wrong. We figured out a better way.
  #as.dist(Biodem::Fst(hedkin, sum(colSums(surnames_freq) )))
}
#TODO: REVISAR. Arriba hay una nota (muuuy antigua) que sugiere que el método es incorrecto.
# Deberíamos revisarlo

# RESOLUCION: NO CAMBIARLO AHORA, SI NO HASTA ARREGLAR TODOS LOS OTROS SCRIPTS
#=======


## REVISIÓN: Me di cuenta que (1-hedkin)  no parece tener sentido, al menos teórico. ¿Qué se intenta hacer con este paso? ¿Para qué se está invirtiendo la escala?
# Alternativas: 
# (1)Usar simplemente as.dist(hedkin). ¿Por qué no se optó por esto?
# (2)Usar sqrt(1 - hedkin).
# (3)Si lo que se busca es normalizar hedkin_normalizado <- hedkin / max(hedkin) y luego distancia <- as.dist(hedkin_normalizado)
# (4) Si se quiere si o si invertir la escala: distancia <- as.dist(sqrt(1 - (hedkin / max(hedkin))) 
#(5) as.dist(Biodem::Fst(hedkin, sum(colSums(surnames_freq) )))
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
y_total <- ape::as.phylo(hc_total) #Phylo format

#Personalizar ?rboles

plot(hc_total) #Dendrograma
ddata <- dendro_data(hc_total, type = "rectangle")
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0))+
  ggdendro::theme_dendro()
p

##2. COOMUNIDADES MUESTREADAS
# Seleccionar los datos que se necesitan
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
selected_communities <- unique(STR$pop)

STR2 <- STR
STR2$pop[STR2$pop %in% c("GUALLIGUAICA", "LA_POLVADA")] <- "PUCLARO"
selected_communities2<- unique(STR2$pop)
comuneros2 <- comuneros
comuneros2$community[comuneros2$community %in% c("GUALLIGUAICA", "LA_POLVADA")] <- "PUCLARO"

# Funcion de creacion de matriz y arbol (SIN PUCLARO)
surnames <- comuneros %>% dplyr::filter(community %in% selected_communities)
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

#=======
#TODO: REVISAR. Es la misma nota que hay arriba respecto al método de conversión de matriz de similitud a matriz de distnacia.
# Deberíamos revisarlo
#=======

surname_matrix_muestra <- surname_distance_muestra(surnames)
hc <-hclust(surname_matrix_muestra,method = "average")
hd <- as.dendrogram(hc)
hy <- ape::as.phylo(hc) #Phylo format
is.rooted(hy)
plotTree(hy)

#### Guardar dendrogramas
write.dendrogram(as.dendrogram(hc_total), file = "Figures/Apellidos.phy", edges = FALSE) #Guardar dendrograma archivo phy
write.nexus(y_total, file = "Figures/Apellidos.nex", translate = TRUE) #Guardar archivo nexus desde phy
write.dendrogram(hd,file = "Figures/Apellidos_muestra.phy", edges = FALSE) #Guardar dendrograma archivo phy
write.nexus(hy, file = "Figures/Apellidos_muestra.nex", translate = TRUE) #Guardar archivo nexus desde phy

