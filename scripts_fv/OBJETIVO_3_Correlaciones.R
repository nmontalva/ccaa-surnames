##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

###  Correr luego de OBJETIVO_3_STR ###

## Cargar paquetes y librerias
library(Hmisc)

## Correlación entre matrices 
matrices <- list(RST = c(RST), RST2 = c(RST2), ASD = c(ASD), GST = c(GST), Nei = c(Nei), cs = c(cs), DSW = c(DSW), Dmu2 = c(Dmu2))
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
