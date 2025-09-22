##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 5 ####
### Part 2 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

## CARGAR LIBRERIAS ##
library(cowplot)
library(GGally)
library(grDevices)
library(grid)
library(magick)
library(stargazer)

## Matriz de graficos ##

calc_r <- function(x, y) {
  r <- cor.test(x, y)  # Calcular el coeficiente de correlaci?n de Pearson
  r_squared <- r$estimate
  p_value<-r$p.value
  return(list(r_squared = r_squared, p_value = p_value))
}

## Comunidades muestreadas
# Grafico 1
annotate_r_squared <- function(data1, mapping, ...) {
  x <- eval_data_col(data1, mapping$x)
  y <- eval_data_col(data1, mapping$y)
  
  result_r <- calc_r(x, y)
  r_squared <- result_r$r_squared
  p_value <- result_r$p_value
  
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
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 7, color = "black")
}

image_paths1 <- c(
  "outputs/Figures/S_A_muestra.svg",
  "outputs/Figures/S_G_muestra.svg",
  "outputs/Figures/A_G_muestra.svg",
  "outputs/Figures/S_M_muestra.svg",
  "outputs/Figures/A_M_muestra.svg",
  "outputs/Figures/G_M_muestra.svg")

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


svg("outputs/Figures/Scatterplot_muestra_1.svg",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)
ggpairs(data1,
        columnLabels = c("S","A","G","M"), 
        lower= list(continuous = custom_image_plot, image_path = image_paths1),
        upper = list(continuous =wrap(annotate_r_squared))
)
dev.off()

# Grafico 2
annotate_r_squared <- function(data1, mapping, ...) {
  x <- eval_data_col(data1, mapping$x)
  y <- eval_data_col(data1, mapping$y)
  
  result_r <- calc_r(x, y)
  r_squared <- result_r$r_squared
  p_value <- result_r$p_value
  
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
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 7, color = "black")
}

svg("outputs/Figures/Scatterplot_muestra_2.svg",width = 3000, height = 3000, res = 300)
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
#Gráfico 3
image_paths1 <- c(
  "outputs/Figures/S_A_muestra.svg",
  "outputs/Figures/S_G_muestra.svg",
  "outputs/Figures/A_G_muestra.svg",
  "outputs/Figures/S_M_muestra.svg",
  "outputs/Figures/A_M_muestra.svg",
  "outputs/Figures/G_M_muestra.svg")


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
svg("outputs/Figures/Scatterplot_muestra_3.svg",width = 3000, height = 3000, res = 300)

# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)
ggpairs(data1,
        columnLabels = c("S","A","G","M"), 
        lower= list(continuous = custom_image_plot, image_path = image_paths1),
        upper = list(continuous =function(data1, mapping, ...) {
          ggplot(data = data1, mapping = mapping) +
            geom_point() +
            geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE)
        })
)
dev.off()

# Extraer valores sin grafico
extract_r_values <- function(data, mapping) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  calc_r(x, y)  # Devuelve lista con r_squared y p_value
}

# 3. Función para procesar todos los pares
get_all_correlations <- function(data, vars) {
  combinations <- combn(vars, 2, simplify = FALSE)
  results <- list()
  
  for (i in seq_along(combinations)) {
    pair <- combinations[[i]]
    mapping <- aes_string(x = pair[1], y = pair[2])
    results[[i]] <- list(
      pair = paste(pair[1], "~", pair[2]),
      values = extract_r_values(data, mapping)
    )
  }
  
  do.call(rbind, lapply(results, as.data.frame))
}
corr_table <- get_all_correlations(
  data = data1,
  vars = c("S_pic1", "G_pic1", "A_pic1", "M_pic1")
)

# Mostrar tabla ordenada por significancia
corr_table[order(corr_table$values.p_value), ]
# Convierte a markdown
md_table <- corr_table %>%
  knitr::kable(format = "markdown")
# Copia al portapapeles
clipr::write_clip(md_table) # Ahora pega (Ctrl+V) directamente en Notion

## Comunidades total
#Grafico 1
annotate_r_squared <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  result_r <- calc_r(x, y)
  r_squared <- result_r$r_squared
  p_value <- result_r$p_value
  
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
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 7, color = "black")
}

image_paths <- c(
  "outputs/Figures/S_A_total.svg",
  "outputs/Figures/S_G_total.svg",
  "outputs/Figures/A_G_total.svg",
  "outputs/Figures/S_M_total.svg",
  "outputs/Figures/A_M_total.svg",
  "outputs/Figures/G_M_total.svg"
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


svg("outputs/Figures/Scatterplot_total_1.svg",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)

ggpairs(data,
        columnLabels = c("S", "A","G","M"), 
        lower= list(continuous = custom_image_plot, image_path = image_paths),
        upper = list(continuous =wrap(annotate_r_squared)))

dev.off()

# Gr?fico 2 Correlaciones
annotate_r_squared <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  result_r <- calc_r(x, y)
  r_squared <- result_r$r_squared
  p_value <- result_r$p_value
  
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
    annotate("text", x = mean(x), y = max(y), label = label, hjust = 0.5, vjust = 0.5, size = 7, color = "black")
}

svg("outputs/Figures/Scatterplot_total_2.svg",width = 3000, height = 3000, res = 300)
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
# Gráfico 3
image_paths <- c(
  "outputs/Figures/S_A_total.svg",
  "outputs/Figures/S_G_total.svg",
  "outputs/Figures/A_G_total.svg",
  "outputs/Figures/S_M_total.svg",
  "outputs/Figures/A_M_total.svg",
  "outputs/Figures/G_M_total.svg"
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
svg("outputs/Figures/Scatterplot_total_3.svg",width = 3000, height = 3000, res = 300)
# Ajustar margen y tama?o de texto para evitar colapso
par(mar = c(1, 1, 0.5, 0.5) + 0.1)

ggpairs(data,
        columnLabels = c("S", "A","G","M"), 
        lower= list(continuous = custom_image_plot, image_path = image_paths),
        upper = list(continuous =function(data, mapping, ...) {
          ggplot(data = data, mapping = mapping) +
            geom_point() +
            geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE)
        }))
dev.off()

#Tabla 
corr_total <- get_all_correlations(
  data = data,
  vars = c("S_pic", "G_pic", "A_pic", "M_pic")
)

# Mostrar tabla ordenada por significancia
corr_total[order(corr_total$values.p_value), ]
# Convierte a markdown
md_total <- corr_total %>%
  knitr::kable(format = "markdown")
# Copia al portapapeles
clipr::write_clip(md_total) # Ahora pega (Ctrl+V) directamente en Notion
