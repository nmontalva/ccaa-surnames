##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

## Cargar paquetes y librerias
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
STR2 <- STR
STR2$pop[STR$pop %in% c("GUALLIGUAICA", "LA_POLVADA")] <- "PUCLARO"
colnames(STR) <- gsub("\\.1$", " ", colnames(STR))
colnames(STR2) <- gsub("\\.1$", " ", colnames(STR2))
write.table(STR, file = "scripts_fv/Datos/STR1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(STR2, file = "scripts_fv/Datos/STR2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
STR [, 3:ncol(STR)] <- lapply(STR[, 3:ncol(STR)], function(x) sprintf("%02d", x)) #Creo que ese "%02d" es también un error de codificación de caracteres
STR2 [, 3:ncol(STR2)] <- lapply(STR2[, 3:ncol(STR2)], function(x) sprintf("%02d", x)) #Creo que ese "%02d" es también un error de codificación de caracteres
locis <- colnames(STR[,3:32])
locis_unicos <- unique(trimws(locis))


## Generar archivos e inputs para crear matrices de distancia ##

## Funcion para generar el contenido Genepop
generateGenepopContent <- function(df, filename, title = "Genepop file") {
  # Extraer los nombres de las columnas de los loci
  loci_names <- locis_unicos
  # Crear el contenido inicial con los nombres de los loci en una sola l?nea
  content <- paste0(title,"\n",paste(loci_names, collapse = ", "), "\n")
  # Iterar sobre cada poblacion
  unique_pops <- unique(df$pop)
  for (pop in unique_pops) {
    # A?adir la etiqueta 'Pop' antes de cada poblaci?n
    content <- paste0(content, "pop\n")
    # Filtrar los individuos de la poblacion actual
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
genepop_content <- generateGenepopContent(STR2,"Archives/STR2_GENEPOP.txt",title= "no comment")
genepop_content <- generateGenepopContent(STR2,"Archives/STR2_GENEPOP.gen",title= "no comment")
# Crear archivo genind
STR_genind <- read.genepop(file = "Archives/STR_GENEPOP.gen", ncode = 2L,quiet = F)
levels(STR_genind@pop) <-unique(STR$pop)
levels(STR_genind@pop)
STR2_genind <- read.genepop(file = "Archives/STR2_GENEPOP.gen", ncode = 2L,quiet = F)
levels(STR2_genind@pop) <-unique(STR2$pop)
levels(STR2_genind@pop)

#Crear archivo Hierfstat 
STR_hierfstat <- genind2hierfstat(STR_genind,pop=STR_genind@pop)
STR2_hierfstat <- genind2hierfstat(STR2_genind,pop=STR2_genind@pop)
#Crear archivo loci
STR_loci <- as.loci(STR_genind,ploidy=2,checkSNP = F)
STR_loci <- na.omit(STR_loci)
STR2_loci <- as.loci(STR2_genind,ploidy=2,checkSNP = F)
STR2_loci <- na.omit(STR2_loci)
#Crear archivo genambig
###STR###
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
#Automatizar la asignacion de los genotipos a mygen
for (locus in locis_unicos) {
  Genotypes(mygen, loci = locus) <- genotype_lists[[locus]]
}
summary(mygen)
Missing(mygen) <- -9
###STR2###
STR2_unique<-STR2[,3:32]
mygen2 <- new("genambig", samples=STR2$ind,
             loci=colnames(STR2_loci[,-1]))
mygen2
Description(mygen2) <- "STR para Rst"
Usatnts(mygen2)<-STR2_genind@loc.n.all
PopNames(mygen2)<-unique(STR2$pop)
PopInfo(mygen2) <- as.integer(factor(STR2$pop, levels = unique(STR2$pop)))
Ploidies(mygen2) <-2
genotype_lists <- list()
for (locus_index in 1:length(locis_unicos)) {
  locus <- locis_unicos[locus_index]
  locus_list <- list()
  for (i in 1:nrow(STR2_unique)) {
    # Extraer los valores de las dos columnas correspondientes al locus actual
    alelo1 <- as.numeric(STR2_unique[i, (2 * locus_index - 1)])
    alelo2 <- as.numeric(STR2_unique[i, (2 * locus_index)])
    
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
#Automatizar la asignacion de los genotipos a mygen
for (locus in locis_unicos) {
  Genotypes(mygen2, loci = locus) <- genotype_lists[[locus]]
}
summary(mygen2)
Missing(mygen2) <- -9

#Correr solo si se quiere excluir datos con NA #
samToUse <- Samples(mygen)
exclude <- c((x<-find.missing.gen(mygen))$Sample)
mygen_exclude <- deleteSamples(mygen, exclude)
summary(mygen_exclude)

#Crear archivo de frecuencuas
rrr <-deSilvaFreq(mygen,self=0,tol=10)# Porque es necesario hacer esto? tolerancia es un numero aleatorio
class(rrr) # Es necesario calcular las frecuencias

simple_STR<-simpleFreq(mygen) # Archivo de frecuencias simple
simple_STR2<-simpleFreq(mygen2) # Archivo de frecuencias simple

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
plot.phylo(phyGST)

GST2 <- pairwise_Gst_Hedrick(STR2_genind, linearized = F)
is.dist(GST2)
GST2 <- as.matrix(GST2)
GST2<-ifelse (GST2 < 0, 0,GST2)
GST2<-as.dist(GST2)
dend.gst <- as.dendrogram(upgma(GST2))
phyGST2 <- as.phylo.dendrogram(dend.gst)
#phyGST<-root(phyGST, outgroup = "PUCLARO", resolve.root = TRUE)
plot.phylo(phyGST)
##NOTA: ASD.txt, dm2.txt y Dsw.txt fueron generados en Populations

#Gst por calcPopDiff con archivos de frecuencia distintos
GST3<-calcPopDiff(simple_STR,metric = "Gst") #Da un resultado similar a GST
GST4<-calcPopDiff(STR_freq,metric = "Gst") #Da un resultado demasiado dispar con los otros

#Arbol con Neighbor-Joining
njGST<-bionj(GST)
#njGST<-root(njGST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njGST)
njGST$edge.length[njGST$edge.length<0]<-0
plot.phylo(njGST, use.edge.length = F)

#Arbol con hclust
conflicts_prefer(ape::as.phylo)
hcGST<-hclust(GST)
plotTree(as.phylo(hcGST))

## Nei Ds
Nei <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Ds"), labels=T)
colnames(Nei) <- rownames(Nei) <- levels(STR_hierfstat$pop)
Nei <- as.dist(Nei)
#Arbol con UPGMA
dend.ds <- as.dendrogram(upgma(Nei))
phyNei <- as.phylo.dendrogram(dend.ds)
#phyNei<-root(phyNei, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyNei)
plotTree(phyNei)

#Arbol con Neighbor-Joining
njNei <- nj(Nei)
njNei<-root(njNei, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njNei)
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
#phyCS<-root(phyCS, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyCS)
plotTree(phyCS)
#Arbol con Neighbor joining
njCS<-nj(cs)
njCS<-root(njCS, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njCS)
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
#phyRST<-root(phyRST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyRST)
plotTree(phyRST)
phyRST$tip.label<-gsub(" ","_",phyRST$tip.label)
#Arbol con Neighbor joining
njRST<-nj(RST)
njRST<-root(njRST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njRST)
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
#phyRST2<-root(phyRST2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyRST2)
dend.rst2<-as.dendrogram(phyRST2)
plotTree(phyRST2)
#Arbol con Neighbor joining
njRST2<-nj(RST2)
njRST2<-root(njRST2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njRST2)
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
#phyASD<-root(phyASD, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyASD)
plot.phylo(phyASD)
dend.asd<-as.dendrogram(upgma(ASD))
#Arbol con Neighbor joining
njASD<-nj(ASD)
njASD<-root(njASD, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njASD)
njASD$edge.length[njASD$edge.length<0]<-0
plot.phylo(njASD, use.edge.length = F)
#Arbol con hclust
hcASD<-hclust(ASD)
plot(hcASD)

#Dsw con populations
DSW <-read.table("Archives/Dsw.txt", header = FALSE, sep = "\t", dec = ".")
DSW <-as.data.frame(DSW)
rownames(DSW) <- as.factor(unique(STR$pop))
DSW<-DSW[,-1]
colnames(DSW) <- as.factor(unique(STR$pop))
DSW <- as.dist(DSW)

DSW2 <-read.table("Archives/Dsw2.txt", header = FALSE, sep = "\t", dec = ".")
DSW2 <-as.data.frame(DSW2)
DSW2
rownames(DSW2) <- as.factor(unique(STR2$pop))
DSW2<-DSW2[,-1]
colnames(DSW2) <- as.factor(unique(STR2$pop))
DSW2 <- as.dist(DSW2)

#Arbol con UPGMA
phyDSW <- upgma(DSW)
phyDSW2 <- upgma(DSW2)
#phyDSW<-root(phyDSW, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyDSW)
plotTree(phyDSW)
plotTree(phyDSW2)
dend.DSW<-as.dendrogram(upgma(DSW))
#Arbol con Neighbor Joining
njDSW<-nj(DSW)
#njDSW<-root(njDSW, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njDSW)
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
#phyDmu2<-root(phyDmu2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyDmu2)
plot.phylo(phyDmu2)
dend.Dmu2<-as.dendrogram(upgma(Dmu2))
#Arbol con Neighbor Joining
njDmu2<-nj(Dmu2)
njDmu2<-root(njDmu2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njDmu2)
njDmu2$edge.length[njDmu2$edge.length<0]<-0
plot.phylo(njDmu2, use.edge.length = F)
#Arbol con hclust
hcDmu2<-hclust(Dmu2)
plot(hcDmu2)

#FST
FST <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Fst"), labels=T)
colnames(FST) <- rownames(FST) <- levels(STR_hierfstat$pop)
FST <- as.dist(FST)
#Arbol con UPGMA
dend.Fst <- as.dendrogram(upgma(FST))
phyFST <- as.phylo.dendrogram(dend.Fst)
#phyFST<-root(phyFST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyFST)
plot(dend.Fst)
plot.phylo(phyFST)
#Arbol con Neighbor joining
njFST<-nj(FST)
#njFST<-root(njFST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njFST)
njFST$edge.length[njFST$edge.length<0]<-0
plot.phylo(njFST, use.edge.length = F)
#Arbol con hclust
hcFST<-as.phylo(hclust(FST))
dend.FST<- as.dendrogram(hcFST)
plot(hcFST)


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

