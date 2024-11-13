##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

## Cargar paquetes y librerias
library(ape)
library(conflicted)
library(phangorn)
library(phylogram)
library(phytools)


##UPGMA preference
conflict_prefer("upgma", "phangorn")
conflict_prefer("as.phylo","ape")

#Dsw con populations
DSW <-read.table("Archives/Dsw.txt", header = FALSE, sep = "\t", dec = ".")
DSW <-as.data.frame(DSW)
rownames(DSW) <- as.factor(unique(STR$pop))
DSW<-DSW[,-1]
colnames(DSW) <- as.factor(unique(STR$pop))
DSW <- as.dist(DSW)

#Arbol con UPGMA
phyDSW <- upgma(DSW)

#phyDSW<-root(phyDSW, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyDSW)
plot.phylo(phyDSW)
dend.DSW<-as.dendrogram(upgma(DSW))

#Arbol con Neighbor Joining
njDSW<-nj(DSW)
is.rooted(njDSW)
njDSW$edge.length[njDSW$edge.length<0]<-0
plot.phylo(njDSW, use.edge.length = F)

#Arbol con hclust
hcDSW<-hclust(DSW)
plot(hcDSW)

write.dendrogram(dend.DSW, file = "Figures/treeDSW.phy", edges = FALSE)
write.nexus(phyDSW, file = "Figures/treeDSW.nex", translate = TRUE)

