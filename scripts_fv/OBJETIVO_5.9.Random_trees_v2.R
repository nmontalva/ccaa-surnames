## =============================================================================
## SCRIPT COMPLETO: ANÁLISIS DE SEÑAL FILOGENÉTICA (K) Y GRAFICADO (G y M)
## =============================================================================

# 1. Cargar librerías necesarias
library(dplyr)
library(geiger)
library(phytools)
library(picante)
library(ape)

# Verificación de seguridad: Asegurar que existen los datos
if (!exists("GM_logit") || !exists("y_total")) {
  stop("ERROR: Necesitas tener cargados los objetos 'GM_logit' (datos) e 'y_total' (árbol) antes de correr este script.")
}

## -----------------------------------------------------------------------------
## 1. PREPARACIÓN DE DATOS
## -----------------------------------------------------------------------------
# Extraer valores G y M
G_values <- GM_logit %>% dplyr::select(G_logit)
M_values <- GM_logit %>% dplyr::select(M_logit)

# Nombrar vectores con los nombres de las especies
trait_data_G <- setNames(G_values$G_logit, row.names(G_values))
trait_data_M <- setNames(M_values$M_logit, row.names(M_values))

# Definir número de iteraciones para las simulaciones
iter <- 1000  # Puedes bajarlo a 100 para pruebas rápidas

## -----------------------------------------------------------------------------
## 2. DEFINICIÓN DE FUNCIONES
## -----------------------------------------------------------------------------

# A) Función para generar árboles con tips aleatorizados
generate_random_trees <- function(phy, n=iter) {
  random_trees <- replicate(n, {
    rt <- phy
    rt$tip.label <- sample(phy$tip.label)  # Permutación de tips
    return(rt)
  }, simplify = FALSE)
  class(random_trees) <- "multiPhylo"
  return(random_trees)
}

# B) Función para cálculo de señal filogenética (Lambda y K)
phylo_signal_analysis <- function(phy, trait, tree_set) {
  # Señal observada
  obs_lambda <- phylosig(phy, trait, method="lambda", test=TRUE, nsim = iter)
  obs_K <- phylosig(phy, trait, method="K", test=TRUE, nsim = iter)
  
  # Señal en árboles randomizados
  rand_lambda <- sapply(tree_set, function(x) phylosig(x, trait, method="lambda", test=TRUE, nsim = iter))
  rand_K <- sapply(tree_set, function(x) phylosig(x, trait, method="K", test=TRUE, nsim = iter))
  
  return(list(obs_lambda=obs_lambda, obs_K=obs_K, 
              rand_lambda=rand_lambda, rand_K=rand_K))
}

# C) Función de GRAFICADO (Solo K, sin p-valor, SVG)
plot_K_only_GM_total <- function(res_G, res_M, archivo) {
  
  # Crear directorio si no existe
  dir_path <- dirname(archivo)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Configurar SVG
  svg(archivo, width = 10, height = 5) 
  
  # Layout: 1 fila, 2 columnas
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
  
  # --- GRÁFICO PARA G ---
  obs_K_G <- res_G$obs_K$K
  sim_K_G <- as.numeric(res_G$rand_K["K", ])
  
  hist(sim_K_G, breaks = 30, col = "lightblue",
       main = "K simulado - Trait G\n(Árbol Completo)",
       xlab = "K simulado", border = "white")
  abline(v = obs_K_G, col = "red", lwd = 2)
  legend("topright",
         legend = paste("K observado =", round(obs_K_G, 3)),
         col = "red", lwd = 2, bty = "n")
  
  # --- GRÁFICO PARA M ---
  obs_K_M <- res_M$obs_K$K
  sim_K_M <- as.numeric(res_M$rand_K["K", ])
  
  hist(sim_K_M, breaks = 30, col = "lightblue",
       main = "K simulado - Trait M\n(Árbol Completo)",
       xlab = "K simulado", border = "white")
  abline(v = obs_K_M, col = "red", lwd = 2)
  legend("topright",
         legend = paste("K observado =", round(obs_K_M, 3)),
         col = "red", lwd = 2, bty = "n")
  
  dev.off()
}

## -----------------------------------------------------------------------------
## 3. EJECUCIÓN DEL ANÁLISIS
## -----------------------------------------------------------------------------
random_trees_complete <- generate_random_trees(y_total, iter)

results_G_complete <- phylo_signal_analysis(y_total, trait_data_G, random_trees_complete)

results_M_complete <- phylo_signal_analysis(y_total, trait_data_M, random_trees_complete)

## -----------------------------------------------------------------------------
## 4. GENERAR GRÁFICO FINAL
## -----------------------------------------------------------------------------

plot_K_only_GM_total(
  res_G = results_G_complete, 
  res_M = results_M_complete, 
  archivo = "outputs/Figures/phylo_signal_K_GM_total.svg"
)
