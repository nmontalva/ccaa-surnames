##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

## Cargar paquetes y librerias
library(ape)
library(conflicted)
library(graph4lg)
library(phylogram)
library(phytools)

##Conflict preference
conflict_prefer("as.phylo","ape")

#Dps
DPS<-mat_pw_dps(STR_genind)
DPS<-as.dist(DPS)
#Arbol con UPGMA
dend.DPS<-as.dendrogram(hclust(DPS,method = "complete"))
phyDPS <- as.phylo(dend.DPS)
is.rooted(phyDPS)
plot.phylo(phyDPS)

write.dendrogram(dend.DPS, file = "Figures/treeDPS.phy", edges = FALSE)
write.nexus(phyDPS, file = "Figures/treeDPS.nex", translate = TRUE)


### Revisar cluster con este árbol
C_trait <- select_variable(result_traits, selected_communities, "cluster")
C_trait2 <- select_variable(result_traits, NULL, "cluster")
C_trait$community <- NULL
C_trait2$community <- NULL
cv1 <- as.matrix(C_trait)[,1]
cv2 <- as.matrix(C_trait2)[,1]

estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram", legend = 0.7 * max(nodeHeights(tree)), ftype = "i", leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}
png("Figures/C_STR_muestra.png")
cvDPS<- estimacion_estados_ancestrales(phyDPS,cv1,"C")
dev.off()
