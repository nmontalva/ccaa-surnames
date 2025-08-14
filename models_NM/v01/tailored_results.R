# Selected results
load("models_NM/sandbox/model_resullts.RData")

## Idea:
# 1. After running evolutionary analysis.R
# 2. Pick the resulting object (i.e. "results")
# 3. Use it to build tailored outputs.

## A table to compare models # LISTO
# Compare BM, single OU, multi OU
# For both G and M
# Include relevant numbers such as:
# - AIC
# - What else?
# Comparison G
Publish::publish(results$G$comparison)
# Comparison M
Publish::publish(results$M$comparison)

## Nicer trees
# Better colours
# I don't like the colours
# Titles lack the pointer to the trait name.
# Remove "NA" from legend
# Considering to add alpha besides theta
# Eliminate the numeric scale at the bottom 
# For some reason, the smaller display sinked to RStudio visual device looks better thant the sinked output. No idea why.
# View plots for a specific variable
plot_regimes <- function(result_obj, regime_colors, title_text, file_path) {
  library(dplyr)
  library(ggplot2)
  
  # Extraer datos únicos por régimen
  regime_info <- result_obj$plots$tree_plot$data %>%
    filter(!is.na(regime)) %>%
    distinct(regime, theta_original, alpha) %>%
    mutate(label_txt = sprintf("%s (θ=%.2f)", regime, theta_original))  # solo θ
  
  # Crear vector de etiquetas para la leyenda
  regime_labels <- setNames(regime_info$label_txt, regime_info$regime)
  
  # Asignar nombres a colores
  names(regime_colors) <- regime_info$regime
  
  # α único
  alpha_val <- unique(na.omit(regime_info$alpha))
  alpha_text <- sprintf("\u03B1 = %.2f", alpha_val)
  
  # Modificar el plot
  plot_out <- result_obj$plots$tree_plot +
    scale_color_manual(
      values = regime_colors,
      labels = regime_labels,
      na.translate = FALSE
    ) +
    ggtitle(title_text, subtitle = alpha_text) +
    theme(
      legend.position = "right",
      legend.justification = "top"
    )
  
  # Guardar
  png(file_path, width = 2200, height = 1900, res = 300)
  print(plot_out)
  dev.off()
  
  return(plot_out)
}
T_model_G <- plot_regimes(
  result_obj = results$G,
  regime_colors = c("#aaFF00", "#0066FF", "#CC00FF"),
  title_text = "Phylogenetic regimes for G",
  file_path = "outputs/Figures/T_model_g.png"
)

T_model_M <- plot_regimes(
  result_obj = results$M,
  regime_colors = c("#aaFF00", "#0066FF", "#CC00FF", "#FF0000"),
  title_text = "Phylogenetic regimes for M",
  file_path = "outputs/Figures/T_model_m.png"
)

print(T_model_G)
print(T_model_M)

# Gráfico combinado
library(cowplot)
library(ggplot2)

# Modificar los plots para que la leyenda esté en el borde
T_model_G2 <- T_model_G + theme(legend.position = "left",
                                plot.title = element_text(hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5) )+
  ggtitle("G")
T_model_M2 <- T_model_M +
  scale_x_reverse() + 
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  ggtitle("M")

# Combinar con cowplot
p_combined <- plot_grid(
  T_model_G2,
  T_model_M2,
  ncol = 2,
  align = "h",
  rel_widths = c(1, 1) # igual tamaño
)
final_plot <- plot_grid(
  ggdraw() + draw_label("Phylogenetic regimes comparison", fontface = 'bold', size = 18),
  p_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)  # altura del título y del contenido
)
print(final_plot)
# Guardar como PNG
png("outputs/Figures/T_model_G_M_combined.png", width = 4000, height = 2000, res = 300)
print(p_combined)
dev.off()

## Regimes
# Summaries already have thetha and logit(theta).
# I think we may also need alpha and other numbers. Maybe checkout what papers normally report.
# Curiously, it seems like ´$parameters´ is empty (return NAs)
# Read surface documentation to figure out what are does values stored at n_regimes: k, kprime, deltak, etc.
# Regimes dist (plot) is nice, but that should just go on a table (?)
print(results$G$summary$regimes)
print.table(results$G$surface$final_model$n_regimes)
# k (the number of regime shifts, counting the basal regime as 1), kprime, (the number of regimes, some of which may be reached by multiple shifts), deltak (k-kprime, a measure of convergence), c (the number of shifts to convergent regimes, another measure of convergence), kprime_conv (the number of convergent regimes shifted to multiple times), and kprime_nonconv (the number of nonconvergent regimes only shifted to once)
#      k         kprime         deltak              c    kprime_conv 
#9              3              6              8              2 
#kprime_nonconv 
#1
ggplot_build(results$G$plots$regime_dist)$data[[1]]
reg_summary_G <-as.data.frame(table(results$G$summary$data$regime))
colnames(reg_summary_G) <- c("regime", "count")
reg_summary_G$percent <- round(100 * reg_summary_G$count / sum(reg_summary_G$count), 1)
print(reg_summary_G)
plot(results$G$plots$regime_dist)
plot(results$G$plots$aic_plot)


print(results$M$summary$regimes)
print.table(results$M$surface$final_model$n_regimes)
gplot_build(results$M$plots$regime_dist)$data[[1]]
reg_summary_M <-as.data.frame(table(results$M$summary$data$regime))
colnames(reg_summary_M) <- c("regime", "count")
reg_summary_M$percent <- round(100 * reg_summary_M$count / sum(reg_summary_M$count), 1)
print(reg_summary_M)



# What's next?
# I need to list actual things to get and sort them.
# Resulting number I need
# Tables
# Plots (includin tree(s))
