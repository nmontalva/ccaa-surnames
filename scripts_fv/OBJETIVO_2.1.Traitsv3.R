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
  
  # Calcular los indices
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
var(result$G) #Varianza
skewness(result$G) #Asimetria
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

#### G y M DATA TOTAL ####
GM_df <- as.data.frame(dplyr::select(result,community, G,M))
rownames(GM_df) <- GM_df$community
GM_df <- GM_df[, -1]

### Agregar datos logit (usamos una verión resumida y ajustada de prep_traits)
# Función que agrega columna logit a result
add_logit_column <- function(df, variable) {
  n <- nrow(df)
  adj <- (df[[variable]] * (n - 1) + 0.5) / n
  logit <- log(adj / (1 - adj))
  df[[paste0(variable, "_logit")]] <- logit
  return(df)
}
# Aplicamos la función a result
result <- result %>%
  add_logit_column("A") %>%
  add_logit_column("S") %>%
  add_logit_column("G") %>%
  add_logit_column("M")
# Verificación
head(result[, c("community", "A", "A_logit", "S", "S_logit", "G", "G_logit", "M", "M_logit")])

