##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 5 ####
### Part 2 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

## CARGAR LIBRERIAS ##
library(ggplot2)
library(cowplot)
library(magick)
library(dplyr)

# 1. PREPARACIÓN DE DATOS (Solo G y M del dataset Total)
if(exists("data")) {
  # Ajustamos nombres. Si tus columnas se llaman distinto (ej: G_mean), ajustalo aquí.
  df_gm <- as.data.frame(data) %>%
    dplyr::select(G = G_pic, M = M_pic) %>% 
    na.omit()
} else {
  stop("Error: El objeto 'data' (Total) no está cargado.")
}

# 2. CARGAR EL ÁRBOL (Lado Derecho)
path_img <- "outputs/Figures/G_M_total.svg"

if(file.exists(path_img)) {
  img <- magick::image_read(path_img)
  img_raster <- as.raster(img)
  
  # Creamos un plot vacío que solo contiene la imagen
  plot_tree <- ggplot() +
    theme_void() +
    xlim(0, 1) + ylim(0, 1) +
    annotation_raster(img_raster, xmin=0, xmax=1, ymin=0, ymax=1) +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
} else {
  # Si no existe la imagen, creamos un cuadro de error
  plot_tree <- ggplot() + theme_void() + 
    annotate("text", x=0.5, y=0.5, label="Imagen no encontrada\n(G_M_total.svg)", color="red")
}

# 3. CREAR EL SCATTERPLOT
# Calculamos correlación primero para ponerla en el texto
test <- cor.test(df_gm$G, df_gm$M)
r_sq <- test$estimate^2
p_val <- test$p.value
stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "")))
label_text <- paste0("Corr = ", round(test$estimate, 4), stars, "\n(R² = ", round(r_sq, 4), ")")

plot_scatter <- ggplot(df_gm, aes(x = G, y = M)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
  annotate("text", x = mean(range(df_gm$G)), y = max(df_gm$M), 
           label = label_text, size = 6, fontface = "bold", vjust = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# 4. COMBINAR Y GUARDAR
svg("outputs/Figures/Scatterplot_G_M.svg", width = 14, height = 7)

# Unimos los dos gráficos lado a lado
final_plot <- cowplot::plot_grid(plot_tree,plot_scatter, ncol = 2, rel_widths = c(1, 1))
print(final_plot)

dev.off()
