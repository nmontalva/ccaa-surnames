##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 5 ####
### Part 4 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

######## Crear gráfico Kmeans vs traits + gráfico de árbol #####################

############################## COMUNIDADES MUESTREADAS #########################

library(PlotTools)
library(phytools)
library(wesanderson)

########### Se usa un data.frame basado en el script de kmeans
df_tree_scaled_numeric <- as.data.frame(df_tree_scaled)  # Asegúrate de que sea un data frame de valores numéricos
kmeans_result <-kmeans(df_tree_scaled_numeric, centers = 3, iter.max = 100, nstart = 100)
Cluster <- as.factor(kmeans_result$cluster)
colnames(kmeans_result$centers)
pca <- prcomp(df_tree_scaled, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:2])  # Mantener solo las dos primeras componentes

# Función para estimar estados ancestrales
estimacion_estados_ancestrales <- function(tree, trait_vector) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  contMap(tree, sorted_trait_vector, plot = F)
}
cluster_symbols <- c(21, 22,24) 

############# S #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        S_trait = S_trait$S, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D

# Estimación de estados ancestrales para el árbol filogenético + colores
obj.S <- estimacion_estados_ancestrales(consensus_tree, sv1)
y.S <- c(obj.S$tree$maps)  # Obtener mapeo de colores del árbol
lims.S <- obj.S$lims  # Límites de valores del rasgo
cols.S <- obj.S$cols  # Paleta de colores utilizada
trans.S <- seq(lims.S[1], lims.S[2], length.out = length(cols.S))  # Rango interpolado
tip_values.S <- sv1[consensus_tree$tip.label]  # Valores del rasgo ordenados
tip_colors.S <- sapply(tip_values.S, function(value) {
  idx <- which.min(abs(trans.S - value))  # Índice del color correspondiente
  cols.S[as.character(idx)]  # Color asociado
})
names(tip_colors.S) <- names(tip_values.S)

# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/S_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.S,
     pch = as.numeric(plot_data$cluster),
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - S")

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.01), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),title = "Cluster", cex = 1)

# Gráfico de contMap
plot(obj.S, type = "phylogram", legend = 0.3 * max(nodeHeights(consensus_tree)), ftype = "i",  fsize = c(1), leg.txt = "S")

# Cerrar el PNG
dev.off()


############# R #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        R_trait = R_trait$R, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# Estimación de estados ancestrales para el árbol filogenético + colores
obj.r <- estimacion_estados_ancestrales(consensus_tree, rv1)
y.r <- c(obj.r$tree$maps) 
lims.r <- obj.r$lims 
cols.r <- obj.r$cols 
trans.r <- seq(lims.r[1], lims.r[2], length.out = length(cols.r))  
tip_values.r <- rv1[consensus_tree$tip.label]  
tip_colors.r <- sapply(tip_values.r, function(value) {
  idx <- which.min(abs(trans.r - value)) 
  cols.r[as.character(idx)] 
})
names(tip_colors.r) <- names(tip_values.r)

# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/R_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.r,
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - R")

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.01), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
plot(obj.r, type = "phylogram", legend = 0.3 * max(nodeHeights(consensus_tree)), ftype = "i",  fsize = c(1), leg.txt = "R")

# Cerrar el PNG
dev.off()

############# A #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        A_trait = A_trait$A, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# Estimación de estados ancestrales para el árbol filogenético + colores
obj.a <- estimacion_estados_ancestrales(consensus_tree, av1)
y.a <- c(obj.a$tree$maps) 
lims.a <- obj.a$lims 
cols.a <- obj.a$cols 
trans.a <- seq(lims.a[1], lims.a[2], length.out = length(cols.a))  
tip_values.a <- av1[consensus_tree$tip.label]  
tip_colors.a <- sapply(tip_values.a, function(value) {
  idx <- which.min(abs(trans.a - value)) 
  cols.a[as.character(idx)] 
})
names(tip_colors.a) <- names(tip_values.a)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/A_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.a,
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - A")

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.01), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
plot(obj.a, type = "phylogram", legend = 0.3 * max(nodeHeights(consensus_tree)), ftype = "i",  fsize = c(1), leg.txt = "A")

# Cerrar el PNG
dev.off()

############# G #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        G_trait = G_trait$G, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# Estimación de estados ancestrales para el árbol filogenético + colores
obj.g <- estimacion_estados_ancestrales(consensus_tree, gv1)
y.g <- c(obj.g$tree$maps) 
lims.g <- obj.g$lims 
cols.g <- obj.g$cols 
trans.g <- seq(lims.g[1], lims.g[2], length.out = length(cols.g))  
tip_values.g <- gv1[consensus_tree$tip.label]  
tip_colors.g <- sapply(tip_values.g, function(value) {
  idx <- which.min(abs(trans.g - value)) 
  cols.g[as.character(idx)] 
})
names(tip_colors.g) <- names(tip_values.g)

# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/G_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.g,
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - G")

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.01), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
plot(obj.g, type = "phylogram", legend = 0.3 * max(nodeHeights(consensus_tree)), ftype = "i",  fsize = c(1), leg.txt = "G")

# Cerrar el PNG
dev.off()

############# M #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        M_trait = M_trait$M, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# Estimación de estados ancestrales para el árbol filogenético + colores
obj.m <- estimacion_estados_ancestrales(consensus_tree, mv1)
y.m <- c(obj.m$tree$maps) 
lims.m <- obj.m$lims 
cols.m <- obj.m$cols 
trans.m <- seq(lims.m[1], lims.m[2], length.out = length(cols.m))  
tip_values.m <- mv1[consensus_tree$tip.label]  
tip_colors.m <- sapply(tip_values.m, function(value) {
  idx <- which.min(abs(trans.m - value)) 
  cols.m[as.character(idx)] 
})
names(tip_colors.m) <- names(tip_values.m)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/M_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.m,
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - M")

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.01), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
plot(obj.m, type = "phylogram", legend = 0.3 * max(nodeHeights(consensus_tree)), ftype = "i",  fsize = c(1), leg.txt = "M")

# Cerrar el PNG
dev.off()


############################## COMUNIDADES TOTALES #############################
########### Se usa un data.frame basado en el script de kmeans
df_tree_scaled_numeric <- as.data.frame(df_tree_scaled_t)
kmeans_result <-kmeans(df_tree_scaled_numeric, centers = 2, iter.max = 100, nstart = 100)
Cluster <- as.factor(kmeans_result$cluster)
colnames(kmeans_result$centers)
pca <- prcomp(df_tree_scaled_t, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:2])  # Mantener solo las dos primeras componentes

# Función para estimar estados ancestrales
estimacion_estados_ancestrales <- function(tree, trait_vector) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  contMap(tree, sorted_trait_vector, plot = F)
}
cluster_symbols <- c(21, 22) 
############# S #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        S_trait = S_trait2$S, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# Estimación de estados ancestrales para el árbol filogenético + colores
obj.S <- estimacion_estados_ancestrales(y_total, sv2)
y.S <- c(obj.S$tree$maps)  # Obtener mapeo de colores del árbol
lims.S <- obj.S$lims  # Límites de valores del rasgo
cols.S <- obj.S$cols  # Paleta de colores utilizada
trans.S <- seq(lims.S[1], lims.S[2], length.out = length(cols.S))  
tip_values.S <- sv2[y_total$tip.label]  # Valores del rasgo ordenados
tip_colors.S <- sapply(tip_values.S, function(value) {
  idx <- which.min(abs(trans.S - value))  # Índice del color correspondiente
  cols.S[as.character(idx)]  # Color asociado
})
names(tip_colors.S) <- names(tip_values.S)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/S_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
 # Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.S,  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - S")

# Gráfico de contMap
plot(obj.S, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "S")

# Cerrar el PNG
dev.off()

# Agregar leyenda para formas (Cluster)
legend(inset = c(-0.25,0.7), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1,y.intersp = 0,xpd = T)
par(xpd = FALSE)
############# R #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        R_trait = R_trait2$R, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
obj.r <- estimacion_estados_ancestrales(y_total, rv2)
y.r <- c(obj.r$tree$maps) 
lims.r <- obj.r$lims 
cols.r <- obj.r$cols 
trans.r <- seq(lims.r[1], lims.r[2], length.out = length(cols.r))  
tip_values.r <- rv2[y_total$tip.label]  
tip_colors.r <- sapply(tip_values.r, function(value) {
  idx <- which.min(abs(trans.r - value)) 
  cols.r[as.character(idx)] 
})
names(tip_colors.r) <- names(tip_values.r)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/R_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.r,  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - R")

# Gráfico de contMap
plot(obj.r, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "R")
dev.off()
############# A #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        A_trait = A_trait2$A, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
obj.a <- estimacion_estados_ancestrales(y_total, av2)
y.a <- c(obj.a$tree$maps) 
lims.a <- obj.a$lims 
cols.a <- obj.a$cols 
trans.a <- seq(lims.a[1], lims.a[2], length.out = length(cols.a))  
tip_values.a <- av2[y_total$tip.label]  
tip_colors.a <- sapply(tip_values.a, function(value) {
  idx <- which.min(abs(trans.a - value)) 
  cols.a[as.character(idx)] 
})
names(tip_colors.a) <- names(tip_values.a)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/A_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.a,  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - A")

# Gráfico de contMap
plot(obj.a, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "A")

# Cerrar el PNG
dev.off()

############# G #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        G_trait = G_trait2$G, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
obj.g <- estimacion_estados_ancestrales(y_total, gv2)
y.g <- c(obj.g$tree$maps) 
lims.g <- obj.g$lims 
cols.g <- obj.g$cols 
trans.g <- seq(lims.g[1], lims.g[2], length.out = length(cols.g))  
tip_values.g <- gv2[y_total$tip.label]  
tip_colors.g <- sapply(tip_values.g, function(value) {
  idx <- which.min(abs(trans.g - value)) 
  cols.g[as.character(idx)] 
})
names(tip_colors.g) <- names(tip_values.g)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/G_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.g,  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - G")

# Gráfico de contMap
plot(obj.g, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "G")

# Cerrar el PNG
dev.off()
############# M #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        M_trait = M_trait2$M, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
obj.m <- estimacion_estados_ancestrales(y_total, mv2)
y.m <- c(obj.m$tree$maps) 
lims.m <- obj.m$lims 
cols.m <- obj.m$cols 
trans.m <- seq(lims.m[1], lims.m[2], length.out = length(cols.m))  
tip_values.m <- mv2[y_total$tip.label]  
tip_colors.m <- sapply(tip_values.m, function(value) {
  idx <- which.min(abs(trans.m - value)) 
  cols.m[as.character(idx)] 
})
names(tip_colors.m) <- names(tip_values.m)
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/M_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = tip_colors.m,  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - M")

# Gráfico de contMap
plot(obj.m, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "M")

# Cerrar el PNG
dev.off()
