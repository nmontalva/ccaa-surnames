# Prepara datos
rasgo <- "G_logit"  # O "M_logit"
fit_obj <- bwd[[length(bwd)]]$fit[[rasgo]]
theta_vals <- as.numeric(fit_obj@theta[[1]])
regimenes <- names(fit_obj@theta[[1]])
cols <- rainbow(length(regimenes))

# Divide layout en dos: árbol + leyenda
layout(matrix(1:2, 1, 2), widths = c(3, 1))

# Plot del árbol sin etiquetas largas
surfaceTreePlot(tree,
                bwd[[length(bwd)]],
                cols = cols,
                show.tip.label = FALSE,
                lwd = 2)

title(paste("Regímenes adaptativos para", rasgo), line = -1)

# Leyenda en panel derecho
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = paste0(regimenes, ": θ = ", round(theta_vals, 2)),
       fill = cols,
       cex = 0.8,
       bty = "n")

