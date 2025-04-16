##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################
### Reordenar este script. Lo necesito antes de los Dendroplot ###
#### OBJETIVO 5 ####
### Part 3 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

################################## K MEANS #####################################
<<<<<<< Updated upstream
<<<<<<< Updated upstream
### SE ELIMINA POR COMPLETO ESTE SCRIPT DEL PROYECTO FINAL
=======
## NOTA: Este script quedaría fuera de los finales
>>>>>>> Stashed changes
=======
## NOTA: Este script quedaría fuera de los finales
>>>>>>> Stashed changes
################ LIBERIAS #################
library(cluster)
library(dplyr)
library(factoextra)
library(fpc)
library(ggplot2)
library(NbClust)
library(phytools)



#### DATA MUESTREADA ####
#Árbol
df_tree <- as.data.frame(cophenetic.phylo(consensus_tree))
df_tree <- na.omit(df_tree)
df_tree <- as.data.frame(lapply(df_tree, as.numeric))
df_tree_scaled <- scale(df_tree)

#### N° of clusters ####
## WSS ##
fviz_nbclust(df_tree_scaled, kmeans, method = "wss") + 
  geom_vline(xintercept = 3, linetype = 2) # 2-4 con DPS
wss <- (nrow(df_tree_scaled)-1)*sum(apply(df_tree_scaled,2,var))
for (i in 1:14) wss[i] <- sum(kmeans(df_tree_scaled,
                                     centers=i)$withinss)
plot(1:14, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
wss #3?

## Silueta ##
fviz_nbclust(df_tree_scaled, kmeans, method = "silhouette") #2

## ASW ##
pamk(df_tree_scaled,krange=1:14,criterion="asw", usepam=TRUE,
     scaling=FALSE, alpha=0.001, diss=inherits(df_tree_scaled, "dist"),
     critout=FALSE, ns=10, seed=NULL) #2
## Gap stat ##
fviz_nbclust(df_tree_scaled, kmeans, method = "gap_stat") #4 


# Ejecutar el K-means con un número de clústeres seleccionado
km <- kmeans(df_tree_scaled, centers = 3, iter.max = 100, nstart = 100)
fviz_cluster(km, data = df_tree_scaled)
par(mfrow=c(1,1))
fviz_cluster(km, data = df_tree_scaled)
plot(consensus_tree)
dev.off()
## 3 clusters principales DPS ##


#### DATA TOTAL ####
#Árbol
df_tree_total <- as.data.frame(as.matrix(surname_matrix))
df_tree_total <- as.data.frame(lapply(df_tree_total, as.numeric))
df_tree_scaled_t <- scale(df_tree_total)

#### N° of clusters ####
## WSS ##
fviz_nbclust(df_tree_scaled_t, kmeans, method = "wss") + 
  geom_vline(xintercept = 3, linetype = 2) #3 o 4?

wss <- (nrow(df_tree_scaled_t)-1)*sum(apply(df_tree_scaled_t,2,var))
for (i in 1:14) wss[i] <- sum(kmeans(df_tree_scaled_t,
                                     centers=i)$withinss)
plot(1:14, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")  # 3 o 6
wss
## Silueta ##
fviz_nbclust(df_tree_scaled_t, kmeans, method = "silhouette") #2

## ASW ##
pamk_res <-pamk(df_tree_scaled_t,krange=1:14,criterion="asw", usepam=TRUE,
     scaling=FALSE, alpha=0.001, diss=inherits(df_tree_scaled_t, "dist"),
     critout=FALSE, ns=10, seed=NULL) 
pam_res <- pam(df_tree_scaled_t, 6)
fviz_cluster(pam_res, palette = c("#00AFBB", "#FC4E07"), ellipse.type = "t", repel = TRUE, ggtheme = theme_classic())

## Gap stat ##
fviz_nbclust(df_tree_scaled_t, kmeans, method = "gap_stat") #7

## NbClust ## No funcionó

#<-NbClust(diss = surname_matrix, min.nc=2, max.nc=7, method="kmeans", index="all")
mds_coords <- cmdscale(surname_matrix, k = 5)
mds_df <- as.data.frame(mds_coords)
resumclust <- NbClust(data = mds_df,diss = surname_matrix, distance = NULL, min.nc = 2,
                           max.nc = 7, method = "average", index = "alllong")
print(resumclust) # 2

# Ejecutar el K-means con un número de clústeres seleccionado
pam_res <- pam(df_tree_scaled_t, 3)
fviz_cluster(pam_res, palette = c("#00AFBB", "#FC4E07"), ellipse.type = "t", repel = F, ggtheme = theme_classic())

km <- kmeans(df_tree_scaled_t, centers = 2, iter.max = 100, nstart = 100)
fviz_cluster(km, data = df_tree_scaled_t)

## 2 clusters principales ##
