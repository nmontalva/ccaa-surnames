##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################
##### PART 2 #####
##### CORRER PRIMERO CREACIÓN DE ESTADOS ANCESTRALES ####

### Extraer valores de estados ancestrales y número de clados acumulados

## Comunidades muestreadas
# Función para crear gráficos de dispersión
create_scatter_plot <- function(clade_numbers, anc_states, measure_name) {
  file_name <- paste0("Figures/Estdos ancestrales plot. Sampled communities -", measure_name, ".png")
  png(file_name, width = 800, height = 600)
  plot(clade_numbers, anc_states,
       type = "l",
       xlab = "Number of clade", 
       ylab = "Ancestral State", 
       main = paste("Sampled communities -", measure_name), 
       pch = 16)
  dev.off()
}

# Listas de estados ancestrales y nombres de medidas
anc_states_list <- list(svc$anc$ace, gvc$anc$ace, mvc$anc$ace, avc$anc$ace)
measure_names <- c("S", "G", "M", "A")

# Crear gráficos de dispersión para cada medida y guardarlos como PNG
for (i in seq_along(anc_states_list)) {
  clade_numbers <- 1:length(anc_states_list[[i]])
  create_scatter_plot(clade_numbers, anc_states_list[[i]], measure_names[i])
}


## Comunidades total
## Función para crear gráficos de dispersión Originales
create_scatter_plot <- function(clade_numbers, anc_states, measure_name, window_size = 10) {
  file_name <- paste0("Figures/Estados ancestrales plot. Total communities -", measure_name, ".png")
  png(file_name, width = 800, height = 600)
  plot(clade_numbers,anc_states,
       type = "l",
       xlab = "Number of clade", 
       ylab = "Ancestral State", 
       main = paste("Total communities -", measure_name), 
       pch = 16)
  dev.off()
}

# Listas de estados ancestrales y nombres de medidas
anc_states_list <- list(svt$anc$ace, gvt$anc$ace, mvt$anc$ace, avt$anc$ace, gs1vt$anc$ace, gs2vt$anc$ace, r2vt$anc$ace)
measure_names <- c("S", "G", "M", "A", "GS1", "GS2", "R2")

# Crear gráficos de dispersión para cada medida y guardarlos como PNG
for (i in seq_along(anc_states_list)) {
  clade_numbers <- 1:length(anc_states_list[[i]])
  create_scatter_plot(clade_numbers, anc_states_list[[i]], measure_names[i])
}




## Función para crear gráficos con datos suavizados ( Media Movil)
if (!requireNamespace("stats", quietly = TRUE)) install.packages("stats")
library(stats)
# Función para suavizar datos
smooth_data <- function(data, window_size) {
  filter <- rep(1 / window_size, window_size)
  smoothed_data <- stats::filter(data, filter, sides = 2)
  return(smoothed_data)
}
create_scatter_plot <- function(clade_numbers, anc_states, measure_name, window_size = 10) {
  file_name <- paste0("Figures/Estados ancestrales suave. Total communities -", measure_name, ".png")
  png(file_name, width = 800, height = 600)
  smoothed_anc_states <- smooth_data(anc_states, window_size)
  plot(clade_numbers,smoothed_anc_states,
       type = "l",
       xlab = "Number of clade", 
       ylab = "Ancestral State", 
       main = paste("Total communities -", measure_name), 
       pch = 16)
  dev.off()
}
# Listas de estados ancestrales y nombres de medidas
anc_states_list <- list(svt$anc$ace, gvt$anc$ace, mvt$anc$ace, avt$anc$ace, gs1vt$anc$ace, gs2vt$anc$ace, r2vt$anc$ace)
measure_names <- c("S", "G", "M", "A", "GS1", "GS2", "R2")
# Crear gráficos de dispersión para cada medida y guardarlos como PNG
for (i in seq_along(anc_states_list)) {
  clade_numbers <- 1:length(anc_states_list[[i]])
  create_scatter_plot(clade_numbers, anc_states_list[[i]], measure_names[i])
}
## Función para crear gráficos con datos suavizados (Savitzky-Golay filtering)
if (!requireNamespace("gsignal", quietly = TRUE)) install.packages("gsignal")
library(gsignal)
# Función para suavizar datos usando Savitzky-Golay
smooth_data_sg <- function(data, window_size, poly_order = 3) {
  smoothed_data <- sgolayfilt(data, p = poly_order, n = window_size)
  return(smoothed_data)
}
# Función para crear gráficos de dispersión con suavizado
create_scatter_plot <- function(clade_numbers, anc_states, measure_name, window_size = 11, poly_order = 3) {
  file_name <- paste0("Figures/Estados ancestrales suaves S-G. Total communities -", measure_name, ".png")
  png(file_name, width = 800, height = 600)
  smoothed_anc_states <- smooth_data_sg(anc_states, window_size, poly_order)
  plot(clade_numbers, smoothed_anc_states,
       type = "l",
       xlab = "Number of clade", 
       ylab = "Ancestral State", 
       main = paste("Total communities -", measure_name), 
       pch = 16)
  dev.off()
}
# Listas de estados ancestrales y nombres de medidas
anc_states_list <- list(svt$anc$ace, gvt$anc$ace, mvt$anc$ace, avt$anc$ace, gs1vt$anc$ace, gs2vt$anc$ace, r2vt$anc$ace)
measure_names <- c("S", "G", "M", "A", "GS1", "GS2", "R2")
# Crear gráficos de dispersión para cada medida y guardarlos como PNG
for (i in seq_along(anc_states_list)) {
  clade_numbers <- 1:length(anc_states_list[[i]])
  create_scatter_plot(clade_numbers, anc_states_list[[i]], measure_names[i])
}
