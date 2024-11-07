##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

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

# Definir colores usando la paleta rainbow
#colors <- wes_palette("Zissou1", 15, type = "continuous")
paleta <- colorRampPalette(c("#ff0500", "#ff7f00", "#ffff00", "#00ff00", "#00ffff", "#0000ff", "#3200ff"))
colors <- paleta(15)
# Función para estimar estados ancestrales
estimacion_estados_ancestrales <- function(tree, trait_vector, leg_txt) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = TRUE)
  plot(obj, type = "phylogram", legend = 0.3 * max(nodeHeights(tree)), ftype = "i",  fsize = c(0.3), leg.txt = leg_txt)
  return(list(anc=anc, obj = obj))
}
cluster_symbols <- c(21, 22,24) 

############# S #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        S_trait = S_trait$S, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/S_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(6, 2, 6, 2) + 0.1)  # Ajusta márgenes según sea necesario

plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$S_trait, breaks = length(colors)))]),  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - S")
#Labels
#text(jitter(plot_data$PC1,amount=0.1), jitter(plot_data$PC2,amount=0.1), labels = colnames(plot_data[1:15]), pos = 3, cex = 0.5, col = "grey")  # Cambia el tamaño de las etiquetas
#Leyenda para colores (S)
SpectrumLegend(
  x = "bottomright",
  inset = c(0.01, 0.1),
  palette = rev(colors),
  legend = round(c(min(plot_data$S_trait), max(plot_data$S_trait)), 3),
  lwd = 10,  # Ancho de la barra de color
  horiz = FALSE,  # Leyenda vertical
  title = "S",
  cex = 1,
  bty = "o",  # Cuadro alrededor de la leyenda
  seg.len = 1,  # Largo de los segmentos
  y.intersp = 2,
  adj = c(0, 0.5)
)

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.5), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(consensus_tree, sv1, "S")
plot(obj$obj, type = "phylogram", legend = 0.6 * max(nodeHeights(consensus_tree)), ftype = "i", fsize = c(0.6), leg.txt = "S")

# Cerrar el PNG
dev.off()


############# R #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        R_trait = R_trait$R, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/R_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(6, 2, 6, 2) + 0.1)  # Ajusta márgenes según sea necesario

plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$R_trait, breaks = length(colors)))]),  # Asignar colores según R_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - R")
#Labels
#text(jitter(plot_data$PC1,amount=0.1), jitter(plot_data$PC2,amount=0.1), labels = colnames(plot_data[1:15]), pos = 3, cex = 0.5, col = "grey")  # Cambia el tamaño de las etiquetas
#Leyenda para colores (R)
SpectrumLegend(
  x = "bottomright",
  inset = c(0.01, 0.1),
  palette = rev(colors),
  legend = round(c(min(plot_data$R_trait), max(plot_data$R_trait)), 2),
  lwd = 10,  # Ancho de la barra de color
  horiz = FALSE,  # Leyenda vertical
  title = "R",
  cex = 1,
  bty = "o",  # Cuadro alrededor de la leyenda
  seg.len = 1,  # Largo de los segmentos
  y.intersp = 2,
  adj = c(0, 0.5)
)

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.5), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(consensus_tree, rv1, "R")
plot(obj$obj, type = "phylogram", legend = 0.6 * max(nodeHeights(consensus_tree)), ftype = "i", fsize = c(0.6), leg.txt = "R")

# Cerrar el PNG
dev.off()

############# A #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        A_trait = A_trait$A, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/A_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(6, 2, 6, 2) + 0.1)  # Ajusta márgenes según sea necesario

plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$A_trait, breaks = length(colors)))]),  # Asignar colores según R_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - A")
#Labels
#text(jitter(plot_data$PC1,amount=0.1), jitter(plot_data$PC2,amount=0.1), labels = colnames(plot_data[1:15]), pos = 3, cex = 0.5, col = "grey")  # Cambia el tamaño de las etiquetas
#Leyenda para colores (A)
SpectrumLegend(
  x = "bottomright",
  inset = c(0.01, 0.1),
  palette = rev(colors),
  legend = round(c(min(plot_data$A_trait), max(plot_data$A_trait)), 2),
  lwd = 10,  # Ancho de la barra de color
  horiz = FALSE,  # Leyenda vertical
  title = "A",
  cex = 1,
  bty = "o",  # Cuadro alrededor de la leyenda
  seg.len = 1,  # Largo de los segmentos
  y.intersp = 2,
  adj = c(0, 0.5)
)

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.5), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(consensus_tree, av1, "A")
plot(obj$obj, type = "phylogram", legend = 0.6 * max(nodeHeights(consensus_tree)), ftype = "i", fsize = c(0.6), leg.txt = "A")

# Cerrar el PNG
dev.off()

############# G #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        G_trait = G_trait$G, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/G_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(6, 2, 6, 2) + 0.1)  # Ajusta márgenes según sea necesario

plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$G_trait, breaks = length(colors)))]),  # Asignar colores según R_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - G")
#Labels
#text(jitter(plot_data$PC1,amount=0.1), jitter(plot_data$PC2,amount=0.1), labels = colnames(plot_data[1:15]), pos = 3, cex = 0.5, col = "grey")  # Cambia el tamaño de las etiquetas
#Leyenda para colores (G)
SpectrumLegend(
  x = "bottomright",
  inset = c(0.01, 0.1),
  palette = rev(colors),
  legend = round(c(min(plot_data$G_trait), max(plot_data$G_trait)), 2),
  lwd = 10,  # Ancho de la barra de color
  horiz = FALSE,  # Leyenda vertical
  title = "G",
  cex = 1,
  bty = "o",  # Cuadro alrededor de la leyenda
  seg.len = 1,  # Largo de los segmentos
  y.intersp = 2,
  adj = c(0, 0.5)
)

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.5), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(consensus_tree, gv1, "G")
plot(obj$obj, type = "phylogram", legend = 0.6 * max(nodeHeights(consensus_tree)), ftype = "i", fsize = c(0.6), leg.txt = "G")

# Cerrar el PNG
dev.off()

############# M #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        M_trait = M_trait$M, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/M_trait_cluster.png", width = 5000, height = 2000, res = 300)

# Configuración de múltiples gráficos
par(mfrow = c(1, 2), mar = c(6, 2, 6, 2) + 0.1)  # Ajusta márgenes según sea necesario

plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$M_trait, breaks = length(colors)))]),  # Asignar colores según R_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - M")
#Labels
#text(jitter(plot_data$PC1,amount=0.1), jitter(plot_data$PC2,amount=0.1), labels = colnames(plot_data[1:15]), pos = 3, cex = 0.5, col = "grey")  # Cambia el tamaño de las etiquetas
#Leyenda para colores (M)
SpectrumLegend(
  x = "bottomright",
  inset = c(0.01, 0.1),
  palette = rev(colors),
  legend = round(c(min(plot_data$M_trait), max(plot_data$M_trait)), 2),
  lwd = 10,  # Ancho de la barra de color
  horiz = FALSE,  # Leyenda vertical
  title = "M",
  cex = 1,
  bty = "o",  # Cuadro alrededor de la leyenda
  seg.len = 1,  # Largo de los segmentos
  y.intersp = 2,
  adj = c(0, 0.5)
)

# Agregar leyenda para formas (Cluster)
legend("bottomright",inset = c(0.01, 0.5), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1)

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(consensus_tree, mv1, "M")
plot(obj$obj, type = "phylogram", legend = 0.6 * max(nodeHeights(consensus_tree)), ftype = "i", fsize = c(0.6), leg.txt = "M")

# Cerrar el PNG
dev.off()


############################## COMUNIDADES TOTALES #############################
########### Se usa un data.frame basado en el script de kmeans
df_tree_scaled_numeric <- as.data.frame(df_tree_scaled_t)  # Asegúrate de que sea un data frame de valores numéricos
kmeans_result <-kmeans(df_tree_scaled_numeric, centers = 2, iter.max = 100, nstart = 100)
Cluster <- as.factor(kmeans_result$cluster)
colnames(kmeans_result$centers)
pca <- prcomp(df_tree_scaled_t, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:2])  # Mantener solo las dos primeras componentes

# Definir colores usando la paleta rainbow
#colors <- wes_palette("Zissou1", 15, type = "continuous")
paleta <- colorRampPalette(c("#ff0500", "#ff7f00", "#ffff00", "#00ff00", "#00ffff", "#0000ff", "#3200ff"))
colors <- paleta(170)
# Función para estimar estados ancestrales
estimacion_estados_ancestrales <- function(tree, trait_vector) {
  sorted_trait_vector <- trait_vector[sort(tree$tip.label)]
  anc <- fastAnc(tree, sorted_trait_vector, vars = TRUE, CI = TRUE)
  obj <- contMap(tree, sorted_trait_vector, plot = F)
  return(list(anc=anc, obj = obj))
}
cluster_symbols <- c(21, 22) 
############# S #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        S_trait = S_trait2$S, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/S_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
 # Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$S_trait, breaks = length(colors)))]),  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - S")

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(y_total, sv2)
plot(obj$obj, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "S")

#Leyenda para colores (S) No funcionan las leyendas
SpectrumLegend(
  inset = c(0, ),
  palette = rev(colors),
  legend = round(c(min(plot_data$S_trait), max(plot_data$S_trait)), 3),
  lwd = 10,  # Ancho de la barra de color
  horiz = FALSE,  # Leyenda vertical
  title = "S",
  cex = 1,
  bty = "o",  # Cuadro alrededor de la leyenda
  seg.len = 1,  # Largo de los segmentos
  y.intersp = -2,
  adj = c(0, 0),
  xpd =T
)

# Agregar leyenda para formas (Cluster)
legend(inset = c(-0.25,0.7), legend = levels(plot_data$cluster), pch = 1:length(levels(plot_data$cluster)),
       title = "Cluster", cex = 1,y.intersp = 0,xpd = T)
par(xpd = FALSE)

# Cerrar el PNG
dev.off()

############# R #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        R_trait = R_trait2$R, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/R_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$R_trait, breaks = length(colors)))]),  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - R")

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(y_total, rv2)
plot(obj$obj, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "R")
dev.off()
############# A #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        A_trait = A_trait2$A, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/A_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$A_trait, breaks = length(colors)))]),  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - A")

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(y_total, av2)
plot(obj$obj, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "A")

# Cerrar el PNG
dev.off()

############# G #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        G_trait = G_trait2$G, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/G_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$G_trait, breaks = length(colors)))]),  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - G")

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(y_total, gv2)
plot(obj$obj, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "G")

# Cerrar el PNG
dev.off()
############# M #############
plot_data <- data.frame(df_tree_scaled_numeric, 
                        M_trait = M_trait2$M, 
                        cluster = as.factor(kmeans_result$cluster))
plot_data <- cbind(plot_data, pca_data)#Reducir a 2D
# ANCESTRAL STATES FOR SAMPLED COMMUNITIES
png("Figures/M_total_cluster.png", width = 6000, height = 2000, res = 300)

# Configuración de múltiples gráficos
# Ajusta márgenes según sea necesario
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
plot(plot_data$PC1, plot_data$PC2,
     col = rev(colors[as.numeric(cut(plot_data$M_trait, breaks = length(colors)))]),  # Asignar colores según S_trait
     pch = as.numeric(plot_data$cluster),  # Diferentes símbolos para cada cluster
     cex = 1.5,
     bg = Cluster,
     xlab = "Dim 1",
     ylab = "Dim 2",
     main = "Gráfico de Clústeres - M")

# Gráfico de contMap
obj <- estimacion_estados_ancestrales(y_total, mv2)
plot(obj$obj, type = "phylogram", legend = 0.5 * max(nodeHeights(y_total)), ftype = "off", leg.txt = "M")

# Cerrar el PNG
dev.off()