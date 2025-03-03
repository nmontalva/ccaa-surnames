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