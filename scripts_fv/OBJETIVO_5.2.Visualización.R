##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 5 ####
### Part 2 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

## CARGAR LIBRERIAS ##
library(cowplot)
library(GGally)
library(png)
library(grid)


## Matriz de graficos ##
## Comunidades muestreadas
# Grafico 1
annotate_r_squared <- function(data1, mapping, ...) {
  x <- eval_data_col(data1, mapping$x)
  y <- eval_data_col(data1, mapping$y)
  
  result <- calc_r(x, y)
  r_squared <- result$r_squared
  p_value <- result$p_value
  
  if (p_value < 0.001) {
    stars <- "***"
  } else if (p_value < 0.01) {
    stars <- "**"
  } else if (p_value < 0.05) {
    stars <- "*"
  } else {
    stars <- ""
  }
  label <- paste("Corr = ", round(r_squared, 4),stars)
  
  ggplot() +
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 4, color = "black")
}

image_paths1 <- c(
  "Figures/S_A_muestra.png",
  "Figures/S_G_muestra.png",
  "Figures/A_G_muestra.png",
  "Figures/S_M_muestra.png",
  "Figures/A_M_muestra.png",
  "Figures/G_M_muestra.png")

plot_with_image <- function(image_path) {
  # Cargar la imagen usando cowplot
  image <- cowplot::ggdraw() + cowplot::draw_image(image_path, scale = 1)
  return(image)
}
image_index <- 1

# Definir la funci?n personalizada para ggally que utilice la imagen correspondiente
custom_image_plot <- function(data1, mapping) {
  image_path <- image_paths1[image_index]
  p <- plot_with_image(image_path)
  image_index <<- image_index + 1
  ggplot() + annotation_custom(ggplotGrob(p), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void()
}


png("Figures/Scatterplot_muestra_1.png",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)
ggpairs(data1,
        columnLabels = c("S","A","G","M"), 
        lower= list(continuous = custom_image_plot, image_path = image_paths1),
        upper = list(continuous =wrap(annotate_r_squared))
)
dev.off()

#=======
#TODO: REVISAR. Error: object 'data1' not found
#=======

# Grafico 2
annotate_r_squared <- function(data1, mapping, ...) {
  x <- eval_data_col(data1, mapping$x)
  y <- eval_data_col(data1, mapping$y)
  
  result <- calc_r(x, y)
  r_squared <- result$r_squared
  p_value <- result$p_value
  
  if (p_value < 0.001) {
    stars <- "***"
  } else if (p_value < 0.01) {
    stars <- "**"
  } else if (p_value < 0.05) {
    stars <- "*"
  } else {
    stars <- ""
  }
  label <- paste("Corr = ", round(r_squared, 4),stars)
  
  ggplot() +
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 4, color = "black")
}


png("Figures/Scatterplot_muestra_2.png",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)
ggpairs(data1,
        columnLabels = c("S", "A","G","M"), 
        lower= list(continuous =function(data1, mapping, ...) {
          ggplot(data = data1, mapping = mapping) +
            geom_point() +
            geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE)
        }),
        upper = list(continuous =wrap(annotate_r_squared))
)
dev.off()

#=======
#TODO: REVISAR. Error: object 'data1' not found
#=======
# REVISION: Si corre el objetivo 5.1 debería crearse el objeto data1

## Comunidades total
#Grafico 1
annotate_r_squared <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  result <- calc_r(x, y)
  r_squared <- result$r_squared
  p_value <- result$p_value
  
  if (p_value < 0.001) {
    stars <- "***"
  } else if (p_value < 0.01) {
    stars <- "**"
  } else if (p_value < 0.05) {
    stars <- "*"
  } else {
    stars <- ""
  }
  label <- paste("Corr = ", round(r_squared, 4),stars)
  
  ggplot() +
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 4, color = "black")
}

image_paths <- c(
  "Figures/S_A_total.png",
  "Figures/S_G_total.png",
  "Figures/A_G_total.png",
  "Figures/S_M_total.png",
  "Figures/A_M_total.png",
  "Figures/G_M_total.png"
)

plot_with_image <- function(image_path) {
  # Cargar la imagen usando cowplot
  image <- cowplot::ggdraw() + cowplot::draw_image(image_path, scale = 1)
  return(image)
}
image_index <- 1

# Definir la funci?n personalizada para ggally que utilice la imagen correspondiente
custom_image_plot <- function(data, mapping) {
  image_path <- image_paths[image_index]
  p <- plot_with_image(image_path)
  image_index <<- image_index + 1
  ggplot() + annotation_custom(ggplotGrob(p), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void()
}


png("Figures/Scatterplot_total_1.png",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)

ggpairs(data,
        columnLabels = c("S", "A","G","M"), 
        lower= list(continuous = custom_image_plot, image_path = image_paths),
        upper = list(continuous =wrap(annotate_r_squared)))

dev.off()

#=======
#TODO: REVISAR. Error in 1:dim(data)[2] : argument of length 0
#=======

# Gr?fico 2 Correlaciones
annotate_r_squared <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  result <- calc_r(x, y)
  r_squared <- result$r_squared
  p_value <- result$p_value
  
  if (p_value < 0.001) {
    stars <- "***"
  } else if (p_value < 0.01) {
    stars <- "**"
  } else if (p_value < 0.05) {
    stars <- "*"
  } else {
    stars <- ""
  }
  label <- paste("Corr = ", round(r_squared, 4),stars)
  
  ggplot() +
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 4, color = "black")
}

png("Figures/Scatterplot_total_2.png",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)
ggpairs(data,
        columnLabels = c("S", "A","G","M"), 
        lower= list(continuous =function(data, mapping, ...) {
          ggplot(data = data, mapping = mapping) +
            geom_point() +
            geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE)
        }),
        upper = list(continuous =wrap(annotate_r_squared))
        )
dev.off()

##TODO: 1: Package `magick` is required to draw images. Image not drawn.

#=======
#TODO: REVISAR. Error in 1:dim(data)[2] : argument of length 0
#=======
## REVISION: Deberian arreglarse todos los problemas si se corre completo el objetivo 5.1