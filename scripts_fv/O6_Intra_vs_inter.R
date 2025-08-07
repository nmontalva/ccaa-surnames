# 6.1. Varianza entre e intra clados

# ----------------------------------------------
# Análisis de varianza intra/inter clado en el árbol
# Objetivo: Evaluar cómo se distribuyen los valores de un rasgo (Gini, Multigeniture)
# a lo largo del árbol, evaluando en distintos niveles (alturas) del árbol si:
# - Las diferencias son entre clados (evolución hacia diferentes óptimos) o
# - Hay alta variabilidad dentro de los clados (intercalado)
# ----------------------------------------------

library(ape)       # Para trabajar con árboles filogenéticos
library(dplyr)     # Para manipulación de datos
library(ggplot2)   # Para visualización

# === INPUTS ===
tree <- y_total                 # Árbol filogenético de comunidades
trait <- GM_df$G               # Reemplazar por GM_df$M si analizas Multigeniture
names(trait) <- rownames(GM_df)

# Filtrar para que coincidan los tips
common <- intersect(tree$tip.label, names(trait))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common))
trait <- trait[common]

# === Convertir el objeto phylo a hclust ===
# Esto es necesario porque cutree() funciona solo con hclust
tree_hclust <- as.hclust.phylo(tree)

# === Definir niveles (alturas) en el árbol ===
heights <- seq(0.01, max(tree_hclust$height), length.out = 100)

# === Preparar estructura para resultados ===
results <- data.frame()

# === Bucle principal ===
for (h in heights) {
  # Cortar el árbol en clusters a la altura h
  clades <- cutree(tree_hclust, h = h)
  
  df <- data.frame(group = factor(clades), value = trait[names(clades)])
  
  # Solo proceder si hay más de un clado
  if (length(unique(df$group)) > 1) {
    kw <- kruskal.test(value ~ group, data = df)
    
    # Calcular varianza intra-clado
    var_within <- df %>%
      group_by(group) %>%
      dplyr::summarize(var = var(value), .groups = "drop") %>%
      pull(var) %>%
      mean(na.rm = TRUE)
    
    var_total <- var(df$value, na.rm = TRUE)
    var_between <- max(var_total - var_within, 0)
    
    results <- rbind(results, data.frame(
      height = h,
      p_value = kw$p.value,
      n_clades = length(unique(df$group)),
      var_within = var_within,
      var_between = var_between,
      var_total = var_total
    ))
  }
}

# === Visualización: líneas para varianza relativa y p-valor ===
ggplot(results, aes(x = height)) +
  geom_line(aes(y = var_within / var_total), color = "blue") +
  geom_line(aes(y = var_between / var_total), color = "red") +
  geom_line(aes(y = p_value), color = "black", linetype = "dashed") +
  labs(
    title = "Varianza intra/inter clado vs altura del árbol",
    x = "Altura del árbol",
    y = "Proporción / p-valor (Kruskal-Wallis)",
    caption = "Rojo = varianza inter clado | Azul = intra | Negro punteado = p-valor"
  ) +
  theme_minimal()

## Otras visualizaciones

library(ggplot2)
library(dplyr)
library(tidyr)

# Normalizar
results_long <- results %>%
  select(height, var_within, var_between) %>%
  pivot_longer(cols = c(var_within, var_between), names_to = "type", values_to = "variance") %>%
  mutate(type = dplyr::recode(type, var_within = "Intra-clado", var_between = "Inter-clado")) %>%
  drop_na()

ggplot(results_long, aes(x = height, y = variance, fill = type)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("Intra-clado" = "blue", "Inter-clado" = "red")) +
  labs(title = "Distribución de la varianza intra e inter clado",
       x = "Altura del árbol", y = "Varianza",
       fill = "Tipo de varianza") +
  theme_minimal()

# Ratio
results <- results %>%
  mutate(var_ratio = var_within / var_between)

ggplot(results, aes(x = height, y = var_ratio)) +
  geom_line(color = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  labs(title = "Razón entre varianza intra e inter clado",
       x = "Altura del árbol", y = "Varianza intra / inter") +
  theme_minimal()