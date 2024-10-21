##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

################################## K MEANS #####################################

################ LIBERIAS #################

library(factoextra)
library(cluster)
library(phytools)


#### DATA ####
df_tree <- as.data.frame(cophenetic.phylo(consensus_tree))
df_tree <- na.omit(df_tree)
df_tree <- as.data.frame(lapply(df_tree, as.numeric))
df_tree_scaled <- scale(df_tree)

#### N° of clusters ####
## suma de cuadrados ##
num <- numeric(15)  # Inicializar el vector num
n_rows <- nrow(df_tree_scaled)

# Calcular la suma de cuadrados dentro de los grupos
for (i in 2:min(15, n_rows)) {
  num[i] <- sum(kmeans(df_tree_scaled, centers = i)$withinss)
}

# Completar el vector num si hay menos de 15 filas
if (n_rows < 15) {
  num[(n_rows + 1):15] <- NA
}

# Visualizar los resultados
plot(1:15, 
     num, 
     type = "b", 
     xlab = "Número de Clusters",  
     ylab = "Suma de cuadrados dentro de grupos",
     col = "red",
     lwd = 2)

## Otros métodos ##
fviz_nbclust(df_tree_scaled, kmeans, method = "wss") + 
  geom_vline(xintercept = 4, linetype = 2) #4 a 7 ?

fviz_nbclust(df_tree_scaled, kmeans, method = "silhouette") #2

fviz_nbclust(df_tree_scaled, kmeans, method = "gap_stat") #9

# Ejecutar el K-means con un número de clústeres seleccionado
km <- kmeans(df_tree_scaled, centers = 8, iter.max = 100, nstart = 100)
fviz_cluster(km, data = df_tree_scaled)

## 2 clusters principales ##

