##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 2 ####
### To estimate the traits of surnames' diversity, concentration of commoners' rights and inheritance's agnatic bias for each community based on the distributions of surnames within communities ###

## Cargar paquetes y librerias ##
library(cluster)
library(dplyr)
library(e1071)
library(factoextra)
library(fpc)
library(ggplot2)
library(Hmisc)
library(NbClust)
library(phytools)
library(REAT)

##Calcular traits
# Definir la funcion gini
gini <- function (x, weights = rep(1, length = length(x))) {
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox] / sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu / nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

# Definir la funcion principal traits
traits <- function(comuneros, group_by_cols = c("community","commune")) {
  # Asegurarse de que group_by_cols es un vector
  if (!is.vector(group_by_cols)) {
    group_by_cols <- as.vector(group_by_cols)
  }
  
  # Calcular los �ndices
  result <- comuneros %>%
    group_by(across(all_of(group_by_cols))) %>%
    summarise(
      N = n(),
      S = n_distinct(surname_father) / N,
      R = mean(rights, na.rm = TRUE),
      G = gini(shares),
      A = mean(sex == "M", na.rm = TRUE),
      M = sum(rights < 1, na.rm = TRUE) / N,
    )
  
  return(result)
}
 
result <-traits(comuneros) 

#Editar tabla
result[is.na(result)] <- 0
result <- as.data.frame(result)
head(result)

#Descripciones estadisticas de traits
Hmisc::describe(result) #Descripciones generales
summary(result) #Mediana, rango
var(result$N) #Varianza
skewness(result$N) #Asimetria
#Normalidad
#histogramas
png(filename = "Figures/Traits_hist.png")
par(mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
par(mfrow = c(2, 3)) 
exclude_cols <- c("community", "commune", "R")
include_cols <- setdiff(colnames(result), exclude_cols)
for (col in include_cols) {
  hist(result[[col]], breaks = 100, main = col, xlab = col)
}
dev.off()
par(mfrow=c(1,1))
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
shapiro.test(result$N)
shapiro.test(result$S)
shapiro.test(result$R)
shapiro.test(result$G)
shapiro.test(result$A)
shapiro.test(result$M)
#Ningún índice se distribuye de manera normal

#Crear cluster
#### G y M DATA TOTAL ####
GM_df <- as.data.frame(select(result,community, G,M))
rownames(GM_df) <- GM_df$community
GM_df <- GM_df[, -1]
# Escalar los datos (opcional)
GM_scaled <- scale(GM_df)

#### N° of clusters ####
## WSS ##
set.seed(123) 
fviz_nbclust(GM_df, kmeans, method = "wss") + 
  geom_vline(xintercept = 2, linetype = 2) #2 o 3

fviz_nbclust(GM_scaled, kmeans, method = "wss") + 
  geom_vline(xintercept = 2, linetype = 2) #2 o 3

wss_df <- (nrow(GM_df)-1)*sum(apply(GM_df,2,var))
for (i in 1:14) wss_df[i] <- sum(kmeans(GM_df,
                                        centers=i)$withinss)
plot(1:14, wss_df, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
wss_df  #2

wss_sc <- (nrow(GM_scaled)-1)*sum(apply(GM_scaled,2,var))
for (i in 1:14) wss_sc[i] <- sum(kmeans(GM_scaled,
                                        centers=i)$withinss)
plot(1:14, wss_sc, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
wss_sc  #2

## Silueta ##
fviz_nbclust(GM_df, kmeans, method = "silhouette") #3
fviz_nbclust(GM_scaled, kmeans, method = "silhouette") #3

## ASW ##
pamk_res <-pamk(GM_df,krange=1:14,criterion="asw", usepam=TRUE,scaling=FALSE, alpha=0.001, diss=inherits(GM_df, "dist"),critout=FALSE, ns=10, seed=NULL) 
pamk_res #4

pamk_res <-pamk(GM_scaled,krange=1:14,criterion="asw", usepam=TRUE,scaling=F, alpha=0.001, diss=inherits(GM_scaled, "dist"),critout=FALSE, ns=10, seed=NULL) 
pamk_res #2

## Gap stat ##
fviz_nbclust(GM_df, kmeans, method = "gap_stat") #5
fviz_nbclust(GM_scaled, kmeans, method = "gap_stat") #5

## NbClust ##
resumclust_df <- NbClust(data = GM_df,diss = NULL, min.nc = 2,max.nc = 7, method = "average", index = "alllong")
#print(resumclust_df) # 3
resumclust_sc <- NbClust(data = GM_scaled,diss = NULL, min.nc = 2,max.nc = 7, method = "average", index = "alllong")
#print(resumclust_sc) # 4

# Ejecutar el K-means con un número de clústeres seleccionado
pam_res <- pam(GM_df, 3)
fviz_cluster(pam_res, palette = "Set2", repel = F, ggtheme = theme_classic(),show.clust.cent = FALSE,labelsize = 0)

km <- kmeans(GM_df, centers = 3, iter.max = 100, nstart = 100)
fviz_cluster(km, data = GM_df,show.clust.cent = FALSE,labelsize = 0)


## 3 clusters principales si no se escalan los datos, 4 si se escalan##

set.seed(123) # Para reproducibilidad
kmeans_result <- kmeans(GM_df, centers = 3, nstart = 1000)


# Agregar los clusters al tibble original
GM_df <- GM_df %>%
  mutate(cluster = unname(kmeans_result$cluster))

# Visualización de los clusters
ggplot(GM_df, aes(x = G, y = M, color = factor(cluster))) +
  geom_point(size = 3) +
  labs(title = "Clusters de comunidades basados en G y M",
       color = "Cluster") +
  theme_minimal()

GM_df$community<-row.names(GM_df)

result_traits <- merge(result,GM_df, by =c("community","G","M"))
colnames(result_traits)
hist(result_traits$cluster)


###TEST DE KRUSTALL-WALLIS ENTRE CLUSTER DE G Y M
## Arbol muestra
#G
G_kw <-result_traits %>%
  filter(community %in% selected_communities) %>%
  select(community,G,cluster) %>%
  {kruskal.test(G~cluster,data=.)} #diff entre medianas
print(G_kw)
#M
M_kw <-result_traits %>%
  filter(community %in% selected_communities) %>%
  select(community,M,cluster) %>%
  {kruskal.test(M~cluster,data=.)} #diff entre medianas
print(M_kw)
## Arbol total
#G
kruskal.test(result_traits$G~result_traits$cluster) #diff entre medianas
#M
kruskal.test(result_traits$M~result_traits$cluster) #diff entre medianas

### Regresiones
pca_result <- prcomp(result_traits[, c("G", "M")], scale. = TRUE)
summary(pca_result)
result_traits_pca<-result_traits
result_traits_pca$x <- x

#TODO:Error in eval(ei, envir) : object 'x' not found

result_traits_pca$PC2 <- pca_result$x[, 2]
result_traits_pca$PC1 <- pca_result$x[, 1]
# Modelo de regresión con interacción
pcr_model <- lm(M ~ PC1 * cluster + PC2 * cluster, data = result_traits_pca)
summary(pcr_model)
kruskal_PC1 <- kruskal.test(PC1 ~ cluster, data = result_traits_pca)
print(kruskal_PC1)
kruskal_PC2 <- kruskal.test(PC2 ~ cluster, data = result_traits_pca)
print(kruskal_PC2)
