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

library(cowplot)
library(ggplot2)

# Modificar los plots para que la leyenda esté en el borde
T_model_G2 <- T_model_G + theme(legend.position = "left",
                                plot.title = element_text(hjust = 0.5) )+
  ggtitle("G")
T_model_M2 <- T_model_M +
  scale_x_reverse() + 
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5) )+
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
library(ggplot2)

# Función para generar tablas compatibles con Markdown de Notion
library(ggplot2)

process_results_markdown <- function(results, type) {
  res <- results[[type]]
  
  # 1. Datos del plot
  plot_data <- ggplot_build(res$plots$regime_dist)$data[[1]]
  
  # 2. Resumen de regimes
  reg_summary <- as.data.frame(table(res$summary$data$regime))
  colnames(reg_summary) <- c("regime", "count")
  reg_summary$percent <- round(100 * reg_summary$count / sum(reg_summary$count), 1)
  
  # 3. Detalles de regimes
  if (!is.null(res$summary$regimes)) {
    reg_details <- res$summary$regimes
    reg_table <- merge(reg_summary, reg_details, by = "regime", all = TRUE)
  } else {
    reg_table <- reg_summary
  }
  
  # 4. Número de regimes en el modelo final
  final_model <- res$surface$final_model$n_regimes
  final_model_df <- data.frame(
    parameter = names(final_model),
    value = unlist(final_model, use.names = FALSE)
  )
  
  # Función para convertir a Markdown
  to_markdown <- function(df) {
    header <- paste("|", paste(colnames(df), collapse=" | "), "|")
    separator <- paste("|", paste(rep("---", ncol(df)), collapse=" | "), "|")
    rows <- apply(df, 1, function(x) paste("|", paste(x, collapse=" | "), "|"))
    md_table <- paste(c(header, separator, rows), collapse="\n")
    return(md_table)
  }
  
  # Convertir ambas tablas a Markdown
  reg_table_md <- to_markdown(reg_table)
  final_model_md <- to_markdown(final_model_df)
  
  # 5. Imprimir tablas Markdown
  cat("\n**Tabla de resumen y detalles de regimes para", type, "**\n")
  cat(reg_table_md, "\n\n")
  
  cat("**Tabla de número de regimes en el modelo final**\n")
  cat(final_model_md, "\n\n")
  
  # 6. Graficar
  print(res$plots$regime_dist)
  if (!is.null(res$plots$aic_plot)) {
    print(res$plots$aic_plot)
  }
  
  # 7. Retornar información
  return(list(plot_data = plot_data, reg_table_md = reg_table_md, final_model_md = final_model_md))
}

# Uso
res_G_md <- process_results_markdown(results, "G")
res_M_md <- process_results_markdown(results, "M")


# What's next?
# I need to list actual things to get and sort them.
# Resulting number I need
# Tables
# Plots (includin tree(s))
