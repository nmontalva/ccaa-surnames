##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 2 ####
### To estimate the traits of surnames' diversity, concentration of commoners' rights and inheritance's agnatic bias for each community based on the distributions of surnames within communities ###

## Cargar paquetes y librerias ##
library(dplyr)
library(e1071)
library(factoextra)
library(fpc)
library(ggplot2)
library(Hmisc)
library(NbClust)
library(phytools)
library(REAT)
library(stargazer)
library(grDevices)


conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)

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
  result_traits <- comuneros %>%
    group_by(across(all_of(group_by_cols))) %>%
    summarise(
      N = n(),
      S = n_distinct(surname_father) / N,
      R = mean(rights, na.rm = TRUE),
      G = gini(shares),
      A = mean(sex == "M", na.rm = TRUE),
      M = sum(rights < 1, na.rm = TRUE) / N,
    )
  
  return(result_traits)
}
 
result_traits <-traits(comuneros) 

#Editar tabla
result_traits[is.na(result_traits)] <- 0
result_traits <- as.data.frame(result_traits)
head(result_traits)

#Descripciones estadisticas de traits
get_stats_df <- function(data, vars) {
  data %>%
    select(all_of(vars)) %>%
    sapply(function(x) {
      c(
        Mean = mean(x, na.rm = TRUE),
        Median = median(x, na.rm = TRUE),
        SD = sd(x, na.rm = TRUE),
        Min = min(x, na.rm = TRUE),
        Max = max(x, na.rm = TRUE),
        Variance = var(x, na.rm = TRUE),
        Skewness = skewness(x, na.rm = TRUE),
        Kurtosis = kurtosis(x, na.rm = TRUE),
        n = sum(!is.na(x))
        )} )%>%
        t() %>%
        as.data.frame() %>%
        round(4)
}
get_stats_ms <- function(data, vars,selected_communities) {
  filtered_data <- data %>% 
    filter(community %in% selected_communities)
  stats <- filtered_data %>%
    select(all_of(vars)) %>%
    sapply(function(x) {
      c(
        Mean = mean(x, na.rm = TRUE),
        Median = median(x, na.rm = TRUE),
        SD = sd(x, na.rm = TRUE),
        Min = min(x, na.rm = TRUE),
        Max = max(x, na.rm = TRUE),
        Variance = var(x, na.rm = TRUE),
        Skewness = skewness(x, na.rm = TRUE),
        Kurtosis = kurtosis(x, na.rm = TRUE),
        n = sum(!is.na(x))
      )} )%>%
    t() %>%
    as.data.frame() %>%
    round(4)
  return(stats)
}
original_stats <- get_stats_df(result_traits, c("N","S", "R", "G", "A", "M"))
stargazer(original_stats,
          type = "latex",
          title = "Descriptive Statistics - Original Traits",
          summary = FALSE,
          rownames = TRUE,
          out = "outputs/Figures/original_traits_stats.tex")
#muestreadas
muestred_stats <- get_stats_ms(result_traits, c("N","S", "R", "G", "A", "M"),selected_communities)
stargazer(muestred_stats,
          type = "latex",
          title = "Descriptive Statistics - Original Muestred Traits",
          summary = FALSE,
          rownames = TRUE,
          out = "outputs/Figures/original_muestred_stats.tex")
#Normalidad
#histogramas
svg(filename = "outputs/Figures/Traits_hist.svg")
par(mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
par(mfrow = c(2, 3)) 
exclude_cols <- c("community", "commune", "R")
include_cols <- setdiff(colnames(result_traits), exclude_cols)
for (col in include_cols) {
  hist(result_traits[[col]], breaks = 100, main = col, xlab = col)
}
dev.off()
par(mfrow=c(1,1))
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
shapiro.test(result_traits$N)
shapiro.test(result_traits$S)
shapiro.test(result_traits$R)
shapiro.test(result_traits$G)
shapiro.test(result_traits$A)
shapiro.test(result_traits$M)
#Ningún índice se distribuye de manera normal

#### G y M DATA TOTAL ####
GM_df <- as.data.frame(dplyr::select(result_traits,community, G,M))
rownames(GM_df) <- GM_df$community
GM_df <- GM_df[, -1]

### Agregar datos logit (usamos una verión resumida y ajustada de prep_traits)
# Función que agrega columna logit a result_traits
add_logit_column <- function(df, variable) {
  n <- nrow(df)
  adj <- (df[[variable]] * (n - 1) + 0.5) / n
  logit <- log(adj / (1 - adj))
  df[[paste0(variable, "_logit")]] <- logit
  return(df)
}
# Aplicamos la función a result_traits
result_traits <- result_traits %>%
  add_logit_column("A") %>%
  add_logit_column("S") %>%
  add_logit_column("G") %>%
  add_logit_column("M")
# Verificación
head(result_traits[, c("community", "A", "A_logit", "S", "S_logit", "G", "G_logit", "M", "M_logit")])
# Descripciones estadísticas logit
logit_stats <- get_stats_df(result_traits, c("S_logit", "G_logit", "A_logit", "M_logit"))
stargazer(logit_stats,
          type = "latex",
          title = "Descriptive Statistics - Logit Transformed Traits",
          summary = FALSE,
          rownames = TRUE,
          out = "outputs/Figures/logit_traits_stats.tex")
#muestreadas
muestred_stats <- get_stats_ms(result_traits, c("S_logit", "G_logit", "A_logit", "M_logit"),selected_communities)
stargazer(muestred_stats,
          type = "latex",
          title = "Descriptive Statistics - Original Muestred Traits",
          summary = FALSE,
          rownames = TRUE,
          out = "outputs/Figures/logit_muestred_stats.tex")
#### G y M DATA TOTAL ####
GM_logit <- as.data.frame(dplyr::select(result_traits,community, G_logit,M_logit))
rownames(GM_logit) <- GM_logit$community
GM_logit <- GM_logit[, -1]
