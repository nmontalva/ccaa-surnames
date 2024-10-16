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
plot.phylo(phyGST2)
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

njGST2<-bionj(GST2)
#njGST<-root(njGST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njGST2)
njGST2$edge.length[njGST2$edge.length<0]<-0
plot.phylo(njGST2, use.edge.length = F)

#Arbol con hclust
hcGST<-hclust(GST)
plotTree(as.phylo(hcGST))
hcGST2<-hclust(GST2)
plot.phylo(as.phylo(hcGST2))

## Nei Ds
Nei <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Ds"), labels=T)
colnames(Nei) <- rownames(Nei) <- levels(STR_hierfstat$pop)
Nei <- as.dist(Nei)
#Arbol con UPGMA
dend.ds <- as.dendrogram(upgma(Nei))
phyNei <- as.phylo.dendrogram(dend.ds)
#phyNei<-root(phyNei, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyNei)
plot.phylo(phyNei)

Nei2 <- as.matrix(genet.dist(STR2_hierfstat, diploid=T, method = "Ds"), labels=T)
colnames(Nei2) <- rownames(Nei2) <- levels(STR2_hierfstat$pop)
Nei2 <- as.dist(Nei2)
#Arbol con UPGMA
dend.ds2 <- as.dendrogram(upgma(Nei2))
phyNei2 <- as.phylo.dendrogram(dend.ds2)
#phyNei<-root(phyNei, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyNei2)
plot.phylo(phyNei2)

#Arbol con Neighbor-Joining
njNei <- nj(Nei)
njNei<-root(njNei, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njNei)
njNei$edge.length[njNei$edge.length<0]<-0
plot.phylo(njNei, use.edge.length = F)
njNei2 <- nj(Nei2)
#njNei2<-root(njNei2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njNei2)
njNei2$edge.length[njNei2$edge.length<0]<-0
plot.phylo(njNei2, use.edge.length = F)
#Arbol con hclust
hcNei<-hclust(Nei)
plot(hcNei)
hcNei2<-hclust(Nei2)
plot(hcNei2)

# Cavalli Sforza
cs <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Dch"), labels=T)
colnames(cs) <- rownames(cs) <- levels(STR_hierfstat$pop)
cs <- as.dist(cs)
cs2 <- as.matrix(genet.dist(STR2_hierfstat, diploid=T, method = "Dch"), labels=T)
colnames(cs2) <- rownames(cs2) <- levels(STR2_hierfstat$pop)
cs2 <- as.dist(cs2)
#Arbol con UPGMA
dend.dch <- as.dendrogram(upgma(cs))
phyCS <- as.phylo.dendrogram(dend.dch)
plot(dend.dch)
dend.dch2 <- as.dendrogram(upgma(cs2))
phyCS2 <- as.phylo.dendrogram(dend.dch2)
plot(dend.dch2)
#phyCS<-root(phyCS, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyCS)
plot.phylo(phyCS)
is.rooted(phyCS2)
plot.phylo(phyCS2)
#Arbol con Neighbor joining
njCS<-nj(cs)
#njCS<-root(njCS, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njCS)
njCS$edge.length[njCS$edge.length<0]<-0
plot.phylo(njCS, use.edge.length = F)
njCS2<-nj(cs2)
#njCS2<-root(njCS2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njCS2)
njCS2$edge.length[njCS2$edge.length<0]<-0
plot.phylo(njCS2, use.edge.length = F)
#Arbol con hclust
hccs<-as.phylo(hclust(cs))
dend.cs<- as.dendrogram(hccs)
plot(hccs)
hccs2<-as.phylo(hclust(cs2))
dend.cs2<- as.dendrogram(hccs2)
plot(hccs2)

#RST
RST <- calcPopDiff(simple_STR,metric = "Rst", object = mygen)
diag(RST) <- 0
RST <- as.matrix (RST)
RST<-ifelse (RST < 0, 0,RST) 
RST <- as.dist(RST)
RST2 <- calcPopDiff(simple_STR2,metric = "Rst", object = mygen2)
diag(RST2) <- 0
RST2 <- as.matrix (RST2)
RST2<-ifelse (RST2 < 0, 0,RST2) 
RST2 <- as.dist(RST2)
#Arbol con Upgma
phyRST <- upgma(RST)
dend.rst <- as.dendrogram(phyRST)
is.rooted(phyRST)
plot.phylo(phyRST)
phyRST2$tip.label<-gsub(" ","_",phyRST2$tip.label)

phyRST2 <- upgma(RST2)
dend.rst <- as.dendrogram(phyRST2)
#phyRST<-root(phyRST, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyRST2)
plot.phylo(phyRST2)
phyRST2$tip.label<-gsub(" ","_",phyRST2$tip.label)

#Arbol con Neighbor joining
njRST<-nj(RST)
is.rooted(njRST)
njRST$edge.length[njRST$edge.length<0]<-0
plot.phylo(njRST, use.edge.length = F)

njRST2<-nj(RST2)
#njRST2<-root(njRST2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njRST2)
njRST2$edge.length[njRST2$edge.length<0]<-0
plot.phylo(njRST2, use.edge.length = F)
#Arbol con hclust
hcrst<-hclust(RST)
plotTree(as.phylo(hcrst))
hcrst2<-hclust(RST2)
plotTree(as.phylo(hcrst2))

#RST_arl con Arlequin
RST_arl <- read.table("Archives/Arlequin_RST.txt", header = T, fill = T, dec = ".")
RST_arl <-as.matrix(RST_arl)
rownames(RST_arl) <- as.factor(selected_communities)
colnames(RST_arl) <- as.factor(selected_communities)
RST_arl <- as.dist(RST_arl)
RST_arl2 <- read.table("Archives/Arlequin_RST2.txt", header = T, fill = T, dec = ".")
RST_arl2 <-as.matrix(RST_arl2)
rownames(RST_arl2) <- as.factor(selected_communities2)
colnames(RST_arl2) <- as.factor(selected_communities2)
RST_arl2 <- as.dist(RST_arl2)
#Arbol upgma
phyRST_arl <- upgma(RST_arl)
phyRST_arl2 <- upgma(RST_arl2)
#phyRST_arl2<-root(phyRST_arl2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyRST_arl)
dend.rst.arl<-as.dendrogram(phyRST_arl)
plot.phylo(phyRST_arl)
is.rooted(phyRST_arl2)
dend.rst.arl2<-as.dendrogram(phyRST_arl2)
plot.phylo(phyRST_arl2)
#Arbol con Neighbor joining
njRST_arl<-nj(RST_arl)
njRST_arl2<-nj(RST_arl2)
#njRST_arl2<-root(njRST_arl2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njRST_arl)
njRST_arl$edge.length[njRST_arl$edge.length<0]<-0
plot.phylo(njRST_arl, use.edge.length = F)
njRST_arl2$edge.length[njRST_arl2$edge.length<0]<-0
plot.phylo(njRST_arl2, use.edge.length = F)
#Arbol con hclust
hcrst_arl<-hclust(RST_arl)
plot.phylo(as.phylo(hcrst_arl))
hcrst_arl2<-hclust(RST_arl2)
plot.phylo(as.phylo(hcrst_arl2))

#(Rst con otro nombre) con populations
ASD <-read.table("Archives/ASD.txt", header = FALSE, sep = "\t", dec = ".")
ASD <-as.data.frame(ASD)
rownames(ASD) <- as.factor(selected_communities)
ASD<-ASD[,-1]
colnames(ASD) <- as.factor(selected_communities)
ASD <- as.dist(ASD)
ASD2 <-read.table("Archives/ASD2.txt", header = FALSE, sep = "\t", dec = ".")
ASD2 <-as.data.frame(ASD2)
rownames(ASD2) <- as.factor(selected_communities2)
ASD2<-ASD2[,-1]
colnames(ASD2) <- as.factor(selected_communities2)
ASD2 <- as.dist(ASD2)
#Arbol con UPGMA
phyASD <- upgma(ASD)
phyASD2 <-upgma(ASD2)
#phyASD2<-root(phyASD2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyASD)
plot.phylo(phyASD)
dend.asd<-as.dendrogram(upgma(ASD))
is.rooted(phyASD2)
plot.phylo(phyASD2)
dend.asd<-as.dendrogram(upgma(ASD2))
#Arbol con Neighbor joining
njASD<-nj(ASD)
is.rooted(njASD)
njASD$edge.length[njASD$edge.length<0]<-0
plot.phylo(njASD, use.edge.length = F)
njASD2<-nj(ASD2)
#njASD2<-root(njASD2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njASD2)
njASD2$edge.length[njASD2$edge.length<0]<-0
plot.phylo(njASD2, use.edge.length = F)
#Arbol con hclust
hcASD<-hclust(ASD)
plot(hcASD)
hcASD<-hclust(ASD2)
plot(hcASD2)

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
plot.phylo(phyDSW)
plot.phylo(phyDSW2)
dend.DSW<-as.dendrogram(upgma(DSW))
dend.DSW2<-as.dendrogram(upgma(DSW2))
#Arbol con Neighbor Joining
njDSW<-nj(DSW)
is.rooted(njDSW)
njDSW$edge.length[njDSW$edge.length<0]<-0
plot.phylo(njDSW, use.edge.length = F)
njDSW2<-nj(DSW2)
#njDSW2<-root(njDSW2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njDSW2)
njDSW2$edge.length[njDSW2$edge.length<0]<-0
plot.phylo(njDSW2, use.edge.length = F)
#Arbol con hclust
hcDSW<-hclust(DSW)
plot(hcDSW)
hcDSW2<-hclust(DSW2)
plot(hcDSW2)

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

# Dmu2 con Arlequin 
Dmu2_arl <- read.table("Archives/Arlequin_dm2.txt", header = T, fill = T, dec = ".")
Dmu2_arl <-as.matrix(Dmu2_arl)
rownames(Dmu2_arl) <- as.factor(selected_communities)
colnames(Dmu2_arl) <- as.factor(selected_communities)
Dmu2_arl <- as.dist(Dmu2_arl)
Dmu2_arl2 <- read.table("Archives/Arlequin_dm22.txt", header = T, fill = T, dec = ".")
Dmu2_arl2 <-as.matrix(Dmu2_arl2)
rownames(Dmu2_arl2) <- as.factor(selected_communities2)
colnames(Dmu2_arl2) <- as.factor(selected_communities2)
Dmu2_arl2 <- as.dist(Dmu2_arl2)
#Arbol upgma
phyDmu2_arl <- upgma(Dmu2_arl)
phyDmu2_arl2 <- upgma(Dmu2_arl2)
#phyDmu2_arl2<-root(phyDmu2_arl2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyDmu2_arl)
dend.Dmu2.arl<-as.dendrogram(phyDmu2_arl)
plot.phylo(phyDmu2_arl)
is.rooted(phyDmu2_arl2)
dend.Dmu2.arl2<-as.dendrogram(phyDmu2_arl2)
plot.phylo(phyDmu2_arl2)
#Arbol con Neighbor joining
njDmu2_arl<-nj(Dmu2_arl)
njDmu2_arl2<-nj(Dmu2_arl2)
#njDmu2_arl2<-root(njDmu2_arl2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njDmu2_arl)
njDmu2_arl$edge.length[njDmu2_arl$edge.length<0]<-0
plot.phylo(njDmu2_arl, use.edge.length = F)
njDmu2_arl2$edge.length[njDmu2_arl2$edge.length<0]<-0
plot.phylo(njDmu2_arl2, use.edge.length = F)
#Arbol con hclust
hcDmu2_arl<-hclust(Dmu2_arl)
plot.phylo(as.phylo(hcDmu2_arl))
hcDmu2_arl2<-hclust(Dmu2_arl2)
plot.phylo(as.phylo(hcDmu2_arl2))

#FST
FST <- as.matrix(genet.dist(STR_hierfstat, diploid=T, method = "Fst"), labels=T)
colnames(FST) <- rownames(FST) <- levels(STR_hierfstat$pop)
FST <- as.dist(FST)
FST2 <- as.matrix(genet.dist(STR2_hierfstat, diploid=T, method = "Fst"), labels=T)
colnames(FST2) <- rownames(FST)2 <- levels(STR2_hierfstat$pop)
FST2 <- as.dist(FST2)
#Arbol con UPGMA
dend.Fst <- as.dendrogram(upgma(FST))
phyFST <- as.phylo.dendrogram(dend.Fst)
dend.Fst2 <- as.dendrogram(upgma(FST2))
phyFST2 <- as.phylo.dendrogram(dend.Fst2)
#phyFST2<-root(phyFST2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(phyFST)
plot(dend.Fst)
plot.phylo(phyFST)
is.rooted(phyFST2)
plot(dend.Fst2)
plot.phylo(phyFST2)
#Arbol con Neighbor joining
njFST<-nj(FST)
njFST2<-nj(FST2)
#njFST2<-root(njFST2, outgroup = "PUCLARO", resolve.root = TRUE)
is.rooted(njFST)
njFST$edge.length[njFST$edge.length<0]<-0
plot.phylo(njFST, use.edge.length = F)
is.rooted(njFST2)
njFST2$edge.length[njFST2$edge.length<0]<-0
plot.phylo(njFST2, use.edge.length = F)
#Arbol con hclust
hcFST<-as.phylo(hclust(FST))
dend.FST<- as.dendrogram(hcFST)
plot(hcFST)
hcFST2<-as.phylo(hclust(FST2))
dend.FST2<- as.dendrogram(hcFST2)
plot(hcFST2)

### SE UTILIZARÁ DSW POR EL MOMENTO ###
png("Figures/DSWtree.png")
plot.phylo(phyDSW)
dev.off()

write.dendrogram(dend.DSW, file = "Figures/treeDSW.phy", edges = FALSE)
write.nexus(phyDSW, file = "Figures/treeDSW.nex", translate = TRUE)

