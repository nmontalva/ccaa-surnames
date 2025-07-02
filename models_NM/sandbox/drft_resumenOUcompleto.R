## Resumen:

# ============================
# INSPECCIÓN AUTOMÁTICA DE RESULTADOS
# ============================

# Usar último paso del modelo
fit_step <- bwd[[length(bwd)]]

# Recorrer todos los rasgos
for (trait in names(fit_step$fit)) {
  cat("\n==== Rasgo:", trait, "====\n")
  
  hansen <- fit_step$fit[[trait]]
  regimes <- hansen@regimes$regs
  theta <- hansen@theta[[1]]
  alpha <- (hansen@sqrt.alpha)^2
  sigma2 <- (hansen@sigma)^2
  var_est <- sigma2 / (2 * alpha)
  
  # Tabla resumen
  tabla <- data.frame(
    Regimen = names(theta),
    Theta = as.numeric(theta),
    Alpha = alpha,
    Sigma2 = sigma2,
    VarEst = var_est
  )
  print(tabla)
  
  # Distribución de tips por régimen
  cat("\nRegímenes encontrados:", length(unique(regimes)), "\n")
  print(table(regimes))
}

# Mostrar shifts (de surfaceSummary)
cat("\n===== RESUMEN GENERAL =====\n")
print(surfaceSummary(bwd))
