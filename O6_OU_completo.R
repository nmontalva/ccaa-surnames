# Preparación inicial ----
options(mc.cores = parallel::detectCores(logical = FALSE))  # Solo Linux, usa núcleos físicos

#Set timer
old <- Sys.time() # get start time

# Convertir árbol y datos con TOT

# Paso previo: escalado:
# Crear nuevas variables transformadas
GM_df$G_logit <- log((GM_df$G + 1e-6) / (1 - GM_df$G + 1e-6))
GM_df$M_logit <- log((GM_df$M + 1e-6) / (1 - GM_df$M + 1e-6))

# Opcional: escalado z-score
GM_df$G_z <- scale(GM_df$G)
GM_df$M_z <- scale(GM_df$M)

# Seleccionar los rasgos para análisis
GM_sub <- GM_df[tips, c("G_logit", "M_logit")]  # se puede cambiar por G_z, M_z

tree <- nameNodes(y_total)
tips <- tree$tip.label
otree_data <- convertTreeData(tree, GM_sub)
otree <- otree_data[[1]]
odata <- otree_data[[2]]

# Con Consensus
 # tree <- nameNodes(consensus_tree2)
 # tips <- tree$tip.label
 # GM_sub <- GM_df[tips, 1:2]  # Ordena y selecciona columnas numéricas
 # otree_data <- convertTreeData(tree, GM_sub)
 # otree <- otree_data[[1]]
 # odata <- otree_data[[2]]

# Otros datos para testear problemas más rápido. Comentar para no correr.
 # data(surfaceDemo)
 # tree <- surfaceDemo$tree
 # dat <- surfaceDemo$sim$dat
 # colnames(dat) <- c("G","M","GregM")
 # olist <- convertTreeData(tree,dat)
 # otree <- olist[[1]]
 # odata <- olist[[2]]
##

# Inicializar lista de resultados
results <- list()

# 5. Correr surfaceForward (con árbol pequeño debería ser manejable)
fwd <- surfaceForward(
  otree, odata,
  aic_threshold = 0,
  exclude = 0,
  verbose = TRUE
)

results$forward <- fwd

# Paso 2: backward para colapsar regímenes convergentes
bwd <- surfaceBackward(otree, odata, fwd[[length(fwd)]], aic_threshold = 0, verbose = TRUE)
results$backward <- bwd

# Mostrar resumen de los resultados
summary(bwd)

### ejemplo de visualización con otros datos

mod <- bwd

#Regimenes
# Extraer etiquetas y regímenes
labels <- mod[[length(bwd)]]$fit$G@nodelabels
regimenes <- mod[[length(bwd)]]$fit$G@regimes$regs

# Filtrar solo tips (los que están en el árbol original)
es_tip <- labels %in% tree$tip.label

# Regímenes por tip (named vector)
regimenes_por_tip <- setNames(regimenes[es_tip], labels[es_tip])
print(regimenes_por_tip)

# Número de regímenes distintos
length(unique(regimenes_por_tip))

# Agrupados por régimen
split(names(regimenes_por_tip), regimenes_por_tip)

# el óptimo adaptativo del rasgo bajo cada régimen.
mod[[length(bwd)]]$fit$G@theta

# la "fuerza de atracción" hacia el óptimo
(mod[[length(bwd)]]$fit$G@sqrt.alpha)^2

#variabilidad estocástica ("ruido")
mod[[length(bwd)]]$fit$G@sigma

#Varianza estacionaria (cuánta dispersión esperamos alrededor del theta en equilibrio)
alpha <-  (mod[[length(bwd)]]$fit$G@sqrt.alpha)^2
sigma_cuadrado <- (mod[[length(bwd)]]$fit$G@sigma)^2
varianza_estacionaria <- sigma_cuadrado / (2 * alpha)
varianza_estacionaria

# Extraer valores
theta <- mod[[length(bwd)]]$fit$G@theta
alpha <- (mod[[length(bwd)]]$fit$G@sqrt.alpha)^2
sigma2 <- (mod[[length(bwd)]]$fit$G@sigma)^2
var_est <- sigma2 / (2 * alpha)

# 
theta_vector <- mod[[length(bwd)]]$fit$G@theta[[1]]  # [[1]] accede al vector interno
# Obtiene SOLO los nombres del vector (a, b, d), ignorando el nombre de la lista padre
nombres_finales <- names(theta_vector)

# Tabla
tabla_regimenes <- data.frame(
  Regimen = nombres_finales,
  Theta = as.numeric(unlist(mod[[length(bwd)]]$fit$G@theta)),
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
fit_obj <- mod[[length(bwd)]]$fit$G

# Extraer theta y nombres de regímenes
theta_vals <- as.numeric(unlist(mod[[length(bwd)]]$fit$G@theta))
regimenes <- nombres_finales

# Definir colores automáticamente (puedes usar otros si prefieres)
cols <- rainbow(length(regimenes))

# Dibujar árbol con colores asignados
surfaceTreePlot(tree, mod[[length(mod)]], cols = cols)

# Añadir leyenda adaptativa con colores y theta
legend("bottomleft",
 legend = paste0(regimenes, ": θ = ", plogis(round(theta_vals, length(regimenes)))),
 fill = cols,
 title = "Regímenes")

# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

