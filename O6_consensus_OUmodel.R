# 1. Cargar librerías
library(ape)
library(ouch)
library(surface)

# 2. Asegurarse de que los nodos internos estén etiquetados
consensus_tree <- nameNodes(consensus_tree2)

# 3. Filtrar GM_df para que contenga solo las comunidades del árbol
tips <- consensus_tree$tip.label
GM_sub <- GM_df[rownames(GM_df) %in% tips, ]

# 4. Convertir el árbol y datos para surface
otree_data <- convertTreeData(consensus_tree, GM_sub[1:2])
otree <- otree_data[[1]]
odata <- otree_data[[2]]

# Inicializar lista de resultados
results <- list()

# 5. Correr surfaceForward (con árbol pequeño debería ser manejable)
fwd <- surfaceForward(
  otree, odata,
  aic_threshold = -2,  # acepta cualquier mejora
  exclude = 0.1,        # no excluye candidatos
  verbose = TRUE
)
results$forward <- fwd

# Paso 2: backward para colapsar regímenes convergentes
bwd <- surfaceBackward(otree, odata, fwd[[1]], aic_threshold = 0, max_steps = 300, verbose = TRUE)
results$backward <- bwd

# Mostrar resumen de los resultados
summary(bwd)

 surfaceTreePlot(consensus_tree, bwd[[1]], cols=c("black","red"))
 #Aquí no hay nada que graficar por que solo hay un régimen

 ### ejemplo de visualización con otros datos
 
 mod <- startingModel(otree, odata, shifts = c("3"="b", "5"="c"))
 surfaceTreePlot(consensus_tree, mod[[1]], cols=c("black","red","blue"))
surfaceTraitPlot(GM_sub[1:2], mod[[1]], cols=c("black","red","blue"))
 
 #Regimenes
  # Extraer etiquetas y regímenes
 labels <- mod[[1]]$fit$G@nodelabels
 regimenes <- mod[[1]]$fit$G@regimes$regs1
 
 # Filtrar solo tips (los que están en el árbol original)
 es_tip <- labels %in% consensus_tree$tip.label
 
 # Regímenes por tip (named vector)
 regimenes_por_tip <- setNames(regimenes[es_tip], labels[es_tip])
 print(regimenes_por_tip)
 
 # Número de regímenes distintos
 length(unique(regimenes_por_tip))
 
 # Agrupados por régimen
 split(names(regimenes_por_tip), regimenes_por_tip)
 
 # el óptimo adaptativo del rasgo bajo cada régimen.
 mod[[1]]$fit$G@theta
 
 # la "fuerza de atracción" hacia el óptimo
 (mod[[1]]$fit$G@sqrt.alpha)^2
 
 #variabilidad estocástica ("ruido")
 mod[[1]]$fit$G@sigma
 
 #Varianza estacionaria (cuánta dispersión esperamos alrededor del theta en equilibrio)
 alpha <-  (mod[[1]]$fit$G@sqrt.alpha)^2
 sigma_cuadrado <- (mod[[1]]$fit$G@sigma)^2
 varianza_estacionaria <- sigma_cuadrado / (2 * alpha)
 varianza_estacionaria
 
 # Extraer valores
 theta <- mod[[1]]$fit$G@theta$G
 alpha <- (mod[[1]]$fit$G@sqrt.alpha)^2
 sigma2 <- (mod[[1]]$fit$G@sigma)^2
 var_est <- sigma2 / (2 * alpha)
 
 # Tabla
 tabla_regimenes <- data.frame(
   Regimen = names(theta),
   Theta = as.numeric(theta),
   Alpha = alpha,
   Sigma2 = sigma2,
   VarianzaEstacionaria = var_est
 )
 print(tabla_regimenes)
 
 # Gráfico de theta ± sqrt(varianza estacionaria)
 library(ggplot2)
 
 ggplot(tabla_regimenes, aes(x = Regimen, y = Theta)) +
   geom_point(size = 3) +
   geom_errorbar(aes(ymin = Theta - sqrt(VarianzaEstacionaria),
                     ymax = Theta + sqrt(VarianzaEstacionaria)),
                 width = 0.2) +
   ylab("Theta ± sd estacionaria") +
   ggtitle("Óptimos adaptativos y su dispersión esperada") +
   theme_minimal()
 
 #Otro mejor gráfico
 # Extraer objeto hansen
 fit_obj <- mod[[1]]$fit$G
 
 # Extraer theta y nombres de regímenes
 theta_vals <- fit_obj@theta$G
 regimenes <- names(theta_vals)
 
 # Definir colores automáticamente (puedes usar otros si prefieres)
 cols <- rainbow(length(regimenes))
 
 # Dibujar árbol con colores asignados
 surfaceTreePlot(consensus_tree, mod[[1]], cols = cols)
 
 # Añadir leyenda adaptativa con colores y theta
 legend("bottomleft",
        legend = paste0(regimenes, ": θ = ", round(theta_vals, 3)),
        fill = cols,
        title = "Regímenes")
 