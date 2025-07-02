library(ggplot2)

# Extraer datos del paso final
fit_step <- bwd[[length(bwd)]]
regimes <- fit_step$fit$G_logit@regimes$regs
names(regimes) <- tree$tip.label  # asegurarse que los nombres coincidan con los tips

# Crear un data.frame combinando datos y regímenes
df_plot <- as.data.frame(odata)
df_plot$regimen <- regimes[rownames(df_plot)]  # agregar columna con régimen

# Asegurar que los valores sean factores para ordenarlos en gráfico
df_plot$regimen <- factor(df_plot$regimen, levels = sort(unique(df_plot$regimen)))

# Función para graficar un rasgo por régimen
plot_trait_by_regime <- function(trait_name) {
  ggplot(df_plot, aes(x = regimen, y = .data[[trait_name]], color = regimen)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    theme_minimal() +
    labs(
      title = paste("Distribución de", trait_name, "por régimen"),
      x = "Régimen",
      y = trait_name
    )
}

# Graficar ambos rasgos
plot_trait_by_regime("G_logit")
plot_trait_by_regime("M_logit")
