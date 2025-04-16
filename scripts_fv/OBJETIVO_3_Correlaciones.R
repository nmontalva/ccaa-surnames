##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

###  Correr luego de OBJETIVO_3_STR ###

## Cargar paquetes y librerias
library(Hmisc)

## Correlación entre matrices 
<<<<<<< Updated upstream
matrices <- list(RST = c(RST), ASD = c(ASD), GST = c(GST), Nei = c(Nei), cs = c(cs), DSW = c(DSW), Dmu2 = c(Dmu2), FST=c(FST),DPS=c(DPS))
=======
matrices <- list(RST = c(RST), RST2 = c(RST2), ASD = c(ASD), GST = c(GST), Nei = c(Nei), cs = c(cs), DSW = c(DSW), Dmu2 = c(Dmu2), FST=c(FST))

#=======
#TODO: REVISAR. El objeto RST no existe. 
# También probé corriendo 3.1 y 3.2. antes de correr este, pero tampoco.
# Así que no puedo probar este script.
# 
#=======
<<<<<<< HEAD
# REVISION: Este script es antiguo y se corre con una versión más antigua de STR que contenía todas las matrices de distancia que creamos. Ese script aún existe pero lo quité del repositorio para dejarlo más limpio. Quizás deberíamos dejar un espacio para almacenar ese tipo de scripts

>>>>>>> Stashed changes
=======

>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22
combined_data <- do.call(cbind, matrices)
rcorr_results <- rcorr(combined_data)

# Extraer las matrices de correlación y p-valores
cor_matrix <- rcorr_results$r
print(cor_matrix)
pval_matrix <- rcorr_results$P
print(pval_matrix)
# Crear una nueva matriz para combinar correlación y p-valores
combined_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(combined_matrix) <- colnames(combined_matrix) <- colnames(combined_data)

# Rellenar la parte inferior con correlaciones y la parte superior con p-valores
for (i in 1:ncol(combined_matrix)) {
  for (j in 1:ncol(combined_matrix)) {
    if (i == j) {
      combined_matrix[i, j] <- 1  # La diagonal siempre será 1
    } else if (i > j) {
      combined_matrix[i, j] <- round(cor_matrix[i, j], 4)  # Correlaciones en la parte inferior
    } else {
      combined_matrix[i, j] <- round(pval_matrix[i, j], 4)  # P-valores en la parte superior
    }
  }
}

# Mostrar la matriz de correlacion y p-valores
View(combined_matrix) ## Contiene las correlaciones en el triangulo inferior y los p-valores en el triangulo superior
