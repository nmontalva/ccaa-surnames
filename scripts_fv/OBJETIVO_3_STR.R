##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### ESPACIO DE TRABAJO ####
#setwd("C:/Users/Kibif/Desktop/Proyecto desigualdad agropastores/Directorio_proyecto")

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

## Cargar paquetes y librerias
#library(assignPOP)
library(tidyverse)
library(ade4)
library(adegenet)
library(ape)
library(caper)
library(conflicted)
library(data.table)
library(dendextend)
library(geiger)
library(genepop)
library(graph4lg)
library(hierfstat)
library(pegas)
library(poppr)
library(phangorn)
library(phylogram)
library(phytools)
library(polysat)
library(mmod)

##UPGMA preference
conflict_prefer("upgma", "phangorn")
conflict_prefer("as.phylo","ape")

## Cargar bases de datos ##
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
colnames(STR) <- gsub("\\.1$", " ", colnames(STR))
write.table(STR, file = "scripts_fv/Datos/STR1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
STR [, 3:ncol(STR)] <- lapply(STR[, 3:ncol(STR)], function(x) sprintf("%02d", x)) #Creo que ese "%02d" es también un error de codificación de caracteres
locis <- colnames(STR[,3:32])
locis_unicos <- unique(trimws(locis))


## Generar archivos e inputs para crear matrices de distancia ##

## Funci?n para generar el contenido Genepop
generateGenepopContent <- function(df, filename, title = "Genepop file") {
  # Extraer los nombres de las columnas de los loci
  loci_names <- locis_unicos
  # Crear el contenido inicial con los nombres de los loci en una sola l?nea
  content <- paste0(title,"\n",paste(loci_names, collapse = ", "), "\n")
  # Iterar sobre cada poblaci?n
  unique_pops <- unique(df$pop)
  for (pop in unique_pops) {
    # A?adir la etiqueta 'Pop' antes de cada poblaci?n
    content <- paste0(content, "pop\n")
    # Filtrar los individuos de la poblaci?n actual
    pop_df <- df[df$pop == pop,]
    # Crear las l?neas correspondientes a cada individuo
    for (i in 1:nrow(pop_df)) {
      ind_name <- pop_df$ind[i]
      ind_data <- pop_df[i, 3:ncol(df)]
    #Concatenar los pares de alelos para cada locus sin tabulaciones, solo espacios
      loci_data <- apply(matrix(ind_data, ncol = 2, byrow = TRUE), 1, paste, collapse = "")
    #A?adir el nombre del individuo y sus datos al contenido
      content <- paste0(content, ind_name,", ", paste(loci_data, collapse = " "), "\n")
    }
  }
  #Escribir el contenido en un archivo .txt
  writeLines(content, con = filename)
}

# Generar el contenido
genepop_content <- generateGenepopContent(STR,"Archives/STR_GENEPOP.txt",title= "no comment")
genepop_content <- generateGenepopContent(STR,"Archives/STR_GENEPOP.gen",title= "no comment")

# Crear archivo genind
STR_genind <- read.genepop(file = "Archives/STR_GENEPOP.gen", ncode = 2L,quiet = F)
levels(STR_genind@pop) <-selected_communities

#Crear archivo Hierfstat 
STR_hierfstat <- genind2hierfstat(STR_genind,pop=STR_genind@pop)

#Crear archivo loci
STR_loci <- as.loci(STR_genind,ploidy=2,checkSNP = F)
STR_loci <- na.omit(STR_loci)

#Crear archivo genambig
STR_unique<-STR[,3:32]
mygen <- new("genambig", samples=STR$ind,
             loci=colnames(STR_loci[,-1]))
mygen
Description(mygen) <- "STR para Rst"
Usatnts(mygen)<-STR_genind@loc.n.all
PopNames(mygen)<-unique(STR$pop)
PopInfo(mygen) <- as.integer(factor(STR$pop, levels = unique(STR$pop)))
Ploidies(mygen) <-2
genotype_lists <- list()
for (locus_index in 1:length(locis_unicos)) {
  locus <- locis_unicos[locus_index]
  locus_list <- list()
  for (i in 1:nrow(STR_unique)) {
  # Extraer los valores de las dos columnas correspondientes al locus actual
  alelo1 <- as.numeric(STR_unique[i, (2 * locus_index - 1)])
  alelo2 <- as.numeric(STR_unique[i, (2 * locus_index)])
    
    # A?adir los alelos como un par a la lista del locus
    locus_list[[i]] <- c(alelo1, alelo2)
  }
  
  # A?adir la lista del locus a la lista principal usando el nombre del locus
  genotype_lists[[locus]] <- locus_list
}
for (locus in names(genotype_lists)) {
  genotype_lists[[locus]] <- lapply(genotype_lists[[locus]], function(x) {
    if (all(x == 0)) {
      return(c(-9))
    } else {
      return(x)
    }
  })
}
#Automatizar la asignaci?n de los genotipos a mygen
for (locus in locis_unicos) {
  Genotypes(mygen, loci = locus) <- genotype_lists[[locus]]
}
summary(mygen)
Missing(mygen) <- -9

#Correr s?lo si se quiere excluir datos con NA #
samToUse <- Samples(mygen)
exclude <- c((x<-find.missing.gen(mygen))$Sample)
mygen2 <- deleteSamples(mygen, exclude)
summary(mygen2)

#Crear archivo de frecuencuas
rrr <-deSilvaFreq(mygen,self=0,tol=10)# Porque es necesario hacer esto? tolerancia es un n?mero aleatorio
class(rrr) # Es necesario calcular las frecuencias

simple_STR<-simpleFreq(mygen) # Archivo de frecuencias simple
#class(simple_freq)

#Archivo de frecuencias funci?n homemade
allele_counts <- tab(STR_genind, freq=FALSE)
STR_freq <- matrix(0, nrow=length(levels(STR_genind@pop)), ncol=ncol(allele_counts))
rownames(STR_freq) <- levels(STR_genind@pop)
colnames(STR_freq) <- colnames(allele_counts)
genomes <- numeric(length(levels(STR_genind@pop)))
names(genomes) <- levels(STR_genind@pop)
## Calcular frecuencias al?licas para cada poblaci?n
for (i in seq_along(levels(STR_genind@pop))) {
  pop <- levels(STR_genind@pop)[i]
  pop_indices <- which(STR_genind@pop == pop)
  pop_counts <- colSums(allele_counts[pop_indices, , drop=FALSE], na.rm=TRUE)
  STR_freq[pop, ] <- pop_counts / sum(pop_counts, na.rm=TRUE)
}
STR_freq <- as.data.frame(STR_freq)
STR_freq$Genomes <- rrr$Genomes
class(STR_freq)

##Matrices de distancia ##
##Goldstein et al. (1995) Gst
#Gst por pairwise Gst-Hedrick
GST <- pairwise_Gst_Hedrick(STR_genind, linearized = F)
is.dist(GST)
GST <- as.matrix(GST)
GST<-ifelse (GST < 0, 0,GST)
GST<-as.dist(GST)
dend.gst <- as.dendrogram(upgma(GST))
phyGST <- as.phylo.dendrogram(dend.gst)
plot(dend.gst)
plot.phylo(phyGST)

##NOTA: GST.txt, Nei.txt, ASD.txt y otos fueron generados en Genepop

#Gst calculado por populations
GST2 <-read.table("Archives/GST.txt", header = FALSE, sep = "\t", dec = ".")
GST2 <-as.data.frame(GST2)
rownames(GST2) <- as.factor(selected_communities)
GST2<-GST2[,-1]
colnames(GST2) <- as.factor(selected_communities)
GST2 <- as.dist(GST2)

#Gst por calcPopDiff con archivos de frecuencia distintos
GST3<-calcPopDiff(simple_STR,metric = "Gst") #Da un resultado similar a GST
GST4<-calcPopDiff(STR_freq,metric = "Gst") #Da un resultado demasiado dispar con los otros

#Arbol con Neighbor-Joining
GSTsolo <- read.table("Archives/GST.txt", header = FALSE, sep = "\t", dec = ".")
GSTsolo <-as.data.frame(GSTsolo)
rownames(GSTsolo) <- as.factor(selected_communities)
GSTsolo<-GSTsolo[,-1]
colnames(GSTsolo) <- as.factor(selected_communities)
GSTsolo <- as.dist(GSTsolo)

njGST<-bionj(GST)

njGST$edge.length[njGST$edge.length<0]<-0
plot.phylo(njGST, use.edge.length = F)



#ArbodesignTree()#Arbol con hclust
conflicts_prefer(ape::as.phylo)
hcGST<-hclust(GST)
plotTree(as.phylo(hcGST))

## Nei Ds
Nei <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Ds"), labels=T)
colnames(Nei) <- rownames(Nei) <- levels(STR_hierfstat$pop)
Nei <- as.dist(Nei)

#Nei con Populations
Nei2 <-read.table("Archives/Nei.txt", header = FALSE, sep = "\t", dec = ".")
Nei2 <-as.data.frame(Nei2)
rownames(Nei2) <- as.factor(selected_communities)
Nei2<-Nei2[,-1]
colnames(Nei2) <- as.factor(selected_communities)
Nei2 <- as.dist(Nei2)

#Arbol con UPGMA
dend.ds <- as.dendrogram(upgma(Nei))
phyNei <- as.phylo.dendrogram(dend.ds)
plot(dend.ds)

#Arbol con Neighbor-Joining
njNei <- nj(Nei)
njNei$edge.length[njNei$edge.length<0]<-0
plot.phylo(njNei, use.edge.length = F)

#Arbol con hclust
hcNei<-hclust(Nei)
plot(hcNei)

# Cavalli Sforza
cs <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Dch"), labels=T)
colnames(cs) <- rownames(cs) <- levels(STR_hierfstat$pop)
cs <- as.dist(cs)
#Arbol con UPGMA
dend.dch <- as.dendrogram(upgma(cs))
phyCS <- as.phylo.dendrogram(dend.dch)
plot(dend.dch)
plot.phylo(phyCS)
#Arbol con Neighbor joining
njCS<-nj(cs)
njCS$edge.length[njCS$edge.length<0]<-0
plot.phylo(njCS, use.edge.length = F)
#Arbol con hclust
hccs<-as.phylo(hclust(cs))
dend.cs<- as.dendrogram(hccs)
plot(hccs)

#RST 1
RST <- calcPopDiff(simple_STR,metric = "Rst", object = mygen)
diag(RST) <- 0
RST <- as.matrix (RST)
RST<-ifelse (RST < 0, 0,RST) 
RST <- as.dist(RST)
#Arbol con Upgma
phyRST <- upgma(RST)
dend.rst <- as.dendrogram(phyRST)
plotTree(phyRST)
phyRST$tip.label<-gsub(" ","_",phyRST$tip.label)
is.rooted(phyRST)
#Arbol con Neighbor joining
njRST<-nj(RST)
njRST$edge.length[njRST$edge.length<0]<-0
plot.phylo(njRST, use.edge.length = F)
#Arbol con hclust
hcrst<-hclust(RST)
plotTree(as.phylo(hcrst))

#RST2 con Arlequin
RST2 <- read.table("Archives/Arlequin_RST.txt", header = T, fill = T, dec = ".")
RST2 <-as.matrix(RST2)
RST2<-RST2[,-1]
RST2<-RST2[-16,]
rownames(RST2) <- as.factor(selected_communities)
colnames(RST2) <- as.factor(selected_communities)
RST2<-ifelse (RST2 < 0, 0,RST2) 
RST2 <- as.dist(RST2)
#Arbol upgma
phyRST2 <- upgma(RST2)
plotTree(phyRST2)
#Arbol con Neighbor joining
njRST2<-nj(RST2)
njRST2$edge.length[njRST2$edge.length<0]<-0
plot.phylo(njRST2, use.edge.length = F)
#Arbol con hclust
hcrst2<-hclust(RST2)
plotTree(as.phylo(hcrst2))

#(Rst con otro nombre) con populations
ASD <-read.table("Archives/ASD.txt", header = FALSE, sep = "\t", dec = ".")
ASD <-as.data.frame(ASD)
rownames(ASD) <- as.factor(selected_communities)
ASD<-ASD[,-1]
colnames(ASD) <- as.factor(selected_communities)
ASD <- as.dist(ASD)
#Arbol con UPGMA
phyASD <- upgma(ASD)
plot.phylo(phyASD)
dend.asd<-as.dendrogram(upgma(ASD))
#Arbol con Neighbor joining
njASD<-nj(ASD)
njASD$edge.length[njASD$edge.length<0]<-0
plot.phylo(njASD, use.edge.length = F)
#Arbol con hclust
hcASD<-hclust(ASD)
plot(hcASD)

#Dsw con populations
DSW <-read.table("Archives/Dsw.txt", header = FALSE, sep = "\t", dec = ".")
DSW <-as.data.frame(DSW)
rownames(DSW) <- as.factor(selected_communities)
DSW<-DSW[,-1]
colnames(DSW) <- as.factor(selected_communities)
DSW <- as.dist(DSW)
#Arbol con UPGMA
phyDSW <- upgma(DSW)
plot.phylo(phyDSW)
dend.DSW<-as.dendrogram(upgma(DSW))
#Arbol con Neighbor Joining
njDSW<-nj(DSW)
njDSW$edge.length[njDSW$edge.length<0]<-0
plot.phylo(njDSW, use.edge.length = F)
#Arbol con hclust
hcDSW<-hclust(DSW)
plot(hcDSW)

#Delta-mu al cuadrado con Populations
Dmu2 <-read.table("Archives/dm2.txt", header = FALSE, sep = "\t", dec = ".")
Dmu2 <-as.data.frame(Dmu2)
rownames(Dmu2) <- as.factor(selected_communities)
Dmu2<-Dmu2[,-1]
colnames(Dmu2) <- as.factor(selected_communities)
Dmu2 <- as.matrix (Dmu2)
Dmu2<-ifelse (Dmu2 < 0, 0,Dmu2) 
Dmu2 <- as.dist(Dmu2)
#Arbol con UPGMA
phyDmu2 <- upgma(Dmu2)
plot.phylo(phyDmu2)
dend.Dmu2<-as.dendrogram(upgma(Dmu2))
#Arbol con Neighbor Joining
njDmu2<-nj(Dmu2)
njDmu2$edge.length[njDmu2$edge.length<0]<-0
plot.phylo(njDmu2, use.edge.length = F)
#Arbol con hclust
hcDmu2<-hclust(Dmu2)
plot(hcDmu2)

### Escribir �rboles
write.dendrogram(dend.gst, file = "Figures/treeGST.phy", edges = FALSE)
write.nexus(phyGST, file = "Figures/treeGST.nex", translate = TRUE)

write.dendrogram(dend.ds, file = "Figures/treeNei.phy", edges = FALSE)
write.nexus(phyNei, file = "Figures/treeNei.nex", translate = TRUE)

write.dendrogram(dend.dch, file = "Figures/treeCS.phy", edges = FALSE)
write.nexus(phyCS, file = "Figures/treeCS.nex", translate = TRUE)

write.dendrogram(dend.rst, file = "Figures/treeRST.phy", edges = FALSE)
write.nexus(phyRST, file = "Figures/treeRST.nex", translate = TRUE)

#write.dendrogram(dend.smm, file = "Figures/treeSMM.phy", edges = FALSE) #Objeto no encontrado
#write.nexus(phySMM, file = "Figures/treeSMM.nex", translate = TRUE) #Objeto no encontrado

###VOY A USAR EL RST POR EL MOMENTO ###
lambda_transformed_tree <- rescale(phyRST,model="lambda",lambda=0.5)
plotTree(lambda_transformed_tree)

