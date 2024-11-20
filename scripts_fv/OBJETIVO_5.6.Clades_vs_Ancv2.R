##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################
#### OBJETIVO 5 ####
### Part 6 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###


### Extraer valores de estados ancestrales y numero de clados acumulados
#Librerías
library(ape)
library(dplyr)
library(ggplot2)
library(phytools)

################################################################################
######################## Comunidades muestreadas ###############################
################################################################################

# Cluster
df_tree <- as.data.frame(cophenetic.phylo(consensus_tree))
df_tree <- na.omit(df_tree)
df_tree <- as.data.frame(lapply(df_tree, as.numeric))
df_tree_scaled <- scale(df_tree)
df_tree_scaled_numeric <- as.data.frame(df_tree_scaled)  # Asegúrate de que sea un data frame de valores numéricos
kmeans_result <-kmeans(df_tree_scaled_numeric, centers = 3, iter.max = 100, nstart = 100)
Cluster <- as.factor(c(kmeans_result$cluster))
names(Cluster) <- colnames(kmeans_result$centers)

# Crear una tabla de clusters
clusters_table <- data.frame(tip.label = names(Cluster), cluster = Cluster)
clusters_table

################################## S ###########################################
# Tree
S_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.127973,LA_POLVADA:0.127973)[&CI={0.221892,0.387561},ancstate={0.304726}]:0.123261,EL_DIVISADERO:0.251235)[&CI={0.211899,0.418322},ancstate={0.315111}]:0.142584,RINCONADA_DE_PUNITAQUI:0.393819)[&CI={0.173239,0.40556},ancstate={0.289399}]:0.077413,EL_ALTAR:0.471233)[&CI={0.168678,0.414167},ancstate={0.291423}]:0.163097,CANELILLA_OVALLE:0.634329)[&CI={0.142831,0.417768},ancstate={0.280299}]:0.084733,MANQUEHUA:0.719062)[&CI={0.121324,0.412715},ancstate={0.26702}]:0.280938,((((((CANELA_BAJA:0.169277,CASTILLO_MAL_PASO:0.169277)[&CI={0.105892,0.25871},ancstate={0.182301}]:0.011386,PUNITAQUI:0.180663)[&CI={0.106218,0.258429},ancstate={0.182324}]:0.02304,MONTE_PATRIA:0.203702)[&CI={0.104629,0.267465},ancstate={0.186047}]:0.125087,BARRAZA:0.328789)[&CI={0.109971,0.322882},ancstate={0.216427}]:0.069566,LA_CALERA:0.398355)[&CI={0.101371,0.336418},ancstate={0.218895}]:0.215934,HUENTELAUQUEN:0.614289)[&CI={0.086222,0.386008},ancstate={0.236115}]:0.385711,EL_ESPINAL:1)[&CI={0.099883,0.43725},ancstate={0.268567}];")
par(mar=c(1,1,1,1))
plot(S_tree,font=2); nodelabels(bg="white",cex =0.5)
terminal_nodes <-S_tree$tip.label
terminal_nodes
#ancestral_nodes <- Ancestors(S_tree,1:15,"all")
paths <- nodepath(S_tree,1:15)
paths

# Extraer los estados ancestrales
S_anc_states <- svc$anc$ace
S_anc_states
#Extraer estados actuales
S_trait_vector <- setNames(S_trait$S, rownames(S_trait))
print(S_trait_vector)


#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster <- clusters_table[clusters_table$tip.label == tip, "cluster"]
  
  # Obtener el valor de S para el nodo terminal desde S_trait
  S_value_terminal <- S_trait_vector[tip]
  
  # Agregar el nodo terminal al dataframe con su valor de S y su cluster
  result_list <- append(result_list, list(data.frame(
    Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
    Cluster = cluster,
    S_value = S_value_terminal
  )))
  
  # Agregar cada nodo ancestral en el path (excluyendo el terminal)
  for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
    # Asignar el valor de S para el nodo ancestral
    S_value <- ifelse(as.character(node) %in% names(S_anc_states), S_anc_states[as.character(node)], NA)
    
    # Agregar el nodo ancestral, cluster y valor de S al resultado
    result_list <- append(result_list, list(data.frame(
      Node = node,
      Cluster = cluster,  # Asignamos el mismo cluster del terminal
      S_value = S_value
    )))
  }
}

# Combinar todos los dataframes en uno solo
result_df_S <- do.call(rbind, result_list)
result_df_S <- result_df_S[!duplicated(result_df_S$Node), ]  # Eliminar filas duplicadas de cada nodo

row.names(result_df_S)<-NULL
print(result_df_S)

# Expansión de result_df_S basado en paths (distancia desde la raíz)
expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  # Obtener el cluster y valores de S para cada nodo en el path
  cluster <- result_df_S$Cluster[match(path[length(path)], result_df_S$Node)]
  S_values <- result_df_S$S_value[match(path, result_df_S$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  # Crear un dataframe temporal para el path
  data.frame(Step = steps_from_root, Cluster = cluster, S_value = S_values, Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de S
  last_step <- max(df$Step)
  last_value <-tail(df$S_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    S_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))

# Graficar los valores de S a lo largo del camino de nodos desde la raíz para cada cluster
ggplot(expanded_df_extended, aes(x = Step, y = S_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de S a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de S",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

png("Figures/S_TRAIT.png",width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = S_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de S a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de S",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

dev.off()
  
################################## R ###########################################
# Cluster
clusters_table
# Tree
R_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.127973,LA_POLVADA:0.127973)[&CI={0.852732,1.124395},ancstate={0.988563}]:0.123261,EL_DIVISADERO:0.251235)[&CI={0.797287,1.135777},ancstate={0.966532}]:0.142584,RINCONADA_DE_PUNITAQUI:0.393819)[&CI={0.798782,1.179741},ancstate={0.989261}]:0.077413,EL_ALTAR:0.471233)[&CI={0.75163,1.15418},ancstate={0.952905}]:0.163097,CANELILLA_OVALLE:0.634329)[&CI={0.618735,1.069575},ancstate={0.844155}]:0.084733,MANQUEHUA:0.719062)[&CI={0.586604,1.064424},ancstate={0.825514}]:0.280938,((((((CANELA_BAJA:0.169277,CASTILLO_MAL_PASO:0.169277)[&CI={0.725941,0.976531},ancstate={0.851236}]:0.011386,PUNITAQUI:0.180663)[&CI={0.713754,0.963348},ancstate={0.838551}]:0.02304,MONTE_PATRIA:0.203702)[&CI={0.676076,0.943093},ancstate={0.809585}]:0.125087,BARRAZA:0.328789)[&CI={0.552073,0.901204},ancstate={0.726638}]:0.069566,LA_CALERA:0.398355)[&CI={0.496686,0.882114},ancstate={0.6894}]:0.215934,HUENTELAUQUEN:0.614289)[&CI={0.491837,0.983424},ancstate={0.73763}]:0.385711,EL_ESPINAL:1)[&CI={0.439061,0.992273},ancstate={0.715667}];
")
paths
# Extraer los estados ancestrales
R_anc_states <- rvc$anc$ace
R_anc_states
#Extraer estados actuales
R_trait_vector <- setNames(R_trait$R, rownames(R_trait))
print(R_trait_vector)

#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster <- clusters_table[clusters_table$tip.label == tip, "cluster"]
  
  # Obtener el valor de R para el nodo terminal desde R_trait
  R_value_terminal <- R_trait_vector[tip]
  
  # Agregar el nodo terminal al dataframe con su valor de R y su cluster
  result_list <- append(result_list, list(data.frame(
    Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
    Cluster = cluster,
    R_value = R_value_terminal
  )))
  
  # Agregar cada nodo ancestral en el path (excluyendo el terminal)
  for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
    # Asignar el valor de R para el nodo ancestral
    R_value <- ifelse(as.character(node) %in% names(R_anc_states), R_anc_states[as.character(node)], NA)
    
    # Agregar el nodo ancestral, cluster y valor de R al resultado
    result_list <- append(result_list, list(data.frame(
      Node = node,
      Cluster = cluster,  # Asignamos el mismo cluster del terminal
      R_value = R_value
    )))
  }
}

# Combinar todos los dataframes en uno solo
result_df_R <- do.call(rbind, result_list)
result_df_R <- result_df_R[!duplicated(result_df_R$Node), ]  # Eliminar filas duplicadas de cada nodo
row.names(result_df_R)<-NULL
print(result_df_R)

# Expansión de result_df_R basado en paths (distancia desde la raíz)
expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  # Obtener el cluster y valores de R para cada nodo en el path
  cluster <- result_df_R$Cluster[match(path[length(path)], result_df_R$Node)]
  R_values <- result_df_R$R_value[match(path, result_df_R$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  # Crear un dataframe temporal para el path
  data.frame(Step = steps_from_root, Cluster = cluster, R_value = R_values, Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de R
  last_step <- max(df$Step)
  last_value <-tail(df$R_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    R_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))

# Graficar los valores de M a lo largo del camino de nodos desde la raíz para cada cluster
ggplot(expanded_df_extended, aes(x = Step, y = R_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de R a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de R",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

png("Figures/R_TRAIT.png",width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = R_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de R a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de R",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

dev.off()
################################## A ###########################################
# Cluster
clusters_table
# Tree
A_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.127973,LA_POLVADA:0.127973)[&CI={0.516337,0.577262},ancstate={0.546799}]:0.123261,EL_DIVISADERO:0.251235)[&CI={0.515608,0.59152},ancstate={0.553564}]:0.142584,RINCONADA_DE_PUNITAQUI:0.393819)[&CI={0.526669,0.612105},ancstate={0.569387}]:0.077413,EL_ALTAR:0.471233)[&CI={0.535683,0.625962},ancstate={0.580822}]:0.163097,CANELILLA_OVALLE:0.634329)[&CI={0.522888,0.623997},ancstate={0.573443}]:0.084733,MANQUEHUA:0.719062)[&CI={0.515228,0.622387},ancstate={0.568808}]:0.280938,((((((CANELA_BAJA:0.169277,CASTILLO_MAL_PASO:0.169277)[&CI={0.531959,0.588158},ancstate={0.560059}]:0.011386,PUNITAQUI:0.180663)[&CI={0.528932,0.584907},ancstate={0.556919}]:0.02304,MONTE_PATRIA:0.203702)[&CI={0.5208,0.580683},ancstate={0.550741}]:0.125087,BARRAZA:0.328789)[&CI={0.513111,0.591409},ancstate={0.55226}]:0.069566,LA_CALERA:0.398355)[&CI={0.506295,0.592733},ancstate={0.549514}]:0.215934,HUENTELAUQUEN:0.614289)[&CI={0.492251,0.602498},ancstate={0.547375}]:0.385711,EL_ESPINAL:1)[&CI={0.498775,0.622842},ancstate={0.560808}];")
paths
# Extraer los estados ancestrales
A_anc_states <- avc$anc$ace
A_anc_states
#Extraer estados actuales
A_trait_vector <- setNames(A_trait$A, rownames(A_trait))
print(A_trait_vector)

#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster <- clusters_table[clusters_table$tip.label == tip, "cluster"]
  
  # Obtener el valor de S para el nodo terminal desde S_trait
  A_value_terminal <- A_trait_vector[tip]
  
  # Agregar el nodo terminal al dataframe con su valor de S y su cluster
  result_list <- append(result_list, list(data.frame(
    Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
    Cluster = cluster,
    A_value = A_value_terminal
  )))
  
  # Agregar cada nodo ancestral en el path (excluyendo el terminal)
  for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
    # Asignar el valor de S para el nodo ancestral
    A_value <- ifelse(as.character(node) %in% names(A_anc_states), A_anc_states[as.character(node)], NA)
    
    # Agregar el nodo ancestral, cluster y valor de S al resultado
    result_list <- append(result_list, list(data.frame(
      Node = node,
      Cluster = cluster,  # Asignamos el mismo cluster del terminal
      A_value = A_value
    )))
  }
}

# Combinar todos los dataframes en uno solo
result_df_A <- do.call(rbind, result_list)

# Para nodos ancestrales, asegurarse de asignar clusters
result_df_A <- result_df_A[!duplicated(result_df_A$Node), ]  # Eliminar filas duplicadas de cada nodo

# Ver el dataframe resultante
row.names(result_df_A)<-NULL
print(result_df_A)

# Expansión de result_df_A basado en paths (distancia desde la raíz)
expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  # Obtener el cluster y valores de A para cada nodo en el path
  cluster <- result_df_A$Cluster[match(path[length(path)], result_df_A$Node)]
  A_values <- result_df_A$A_value[match(path, result_df_A$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  # Crear un dataframe temporal para el path
  data.frame(Step = steps_from_root, Cluster = cluster, A_value = A_values, Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de A
  last_step <- max(df$Step)
  last_value <-tail(df$A_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    A_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))

# Graficar los valores de A a lo largo del camino de nodos desde la raíz para cada cluster
ggplot(expanded_df_extended, aes(x = Step, y = A_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de A a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de A",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

png("Figures/A_TRAIT.png",width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = A_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de A a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de A",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

dev.off()
################################## G ###########################################
# Cluster
clusters_table
# Tree
G_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.127973,LA_POLVADA:0.127973)[&CI={-0.0521,0.128952},ancstate={0.038426}]:0.123261,EL_DIVISADERO:0.251235)[&CI={-0.000347,0.225242},ancstate={0.112448}]:0.142584,RINCONADA_DE_PUNITAQUI:0.393819)[&CI={0.054363,0.308255},ancstate={0.181309}]:0.077413,EL_ALTAR:0.471233)[&CI={0.057315,0.325597},ancstate={0.191456}]:0.163097,CANELILLA_OVALLE:0.634329)[&CI={0.1019,0.402365},ancstate={0.252132}]:0.084733,MANQUEHUA:0.719062)[&CI={0.092656,0.411103},ancstate={0.251879}]:0.280938,((((((CANELA_BAJA:0.169277,CASTILLO_MAL_PASO:0.169277)[&CI={0.138494,0.305502},ancstate={0.221998}]:0.011386,PUNITAQUI:0.180663)[&CI={0.144077,0.31042},ancstate={0.227248}]:0.02304,MONTE_PATRIA:0.203702)[&CI={0.154424,0.33238},ancstate={0.243402}]:0.125087,BARRAZA:0.328789)[&CI={0.181639,0.41432},ancstate={0.29798}]:0.069566,LA_CALERA:0.398355)[&CI={0.190979,0.44785},ancstate={0.319415}]:0.215934,HUENTELAUQUEN:0.614289)[&CI={0.111258,0.43888},ancstate={0.275069}]:0.385711,EL_ESPINAL:1)[&CI={0.080061,0.448753},ancstate={0.264407}];
")
paths
# Extraer los estados ancestrales
G_anc_states <- gvc$anc$ace
G_anc_states
#Extraer estados actuales
G_trait_vector <- setNames(G_trait$G, rownames(G_trait))
print(G_trait_vector)

#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster <- clusters_table[clusters_table$tip.label == tip, "cluster"]
  
  # Obtener el valor de G para el nodo terminal desde G_trait
  G_value_terminal <- G_trait_vector[tip]
  
  # Agregar el nodo terminal al dataframe con su valor de G y su cluster
  result_list <- append(result_list, list(data.frame(
    Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
    Cluster = cluster,
    G_value = G_value_terminal
  )))
  
  # Agregar cada nodo ancestral en el path (excluyendo el terminal)
  for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
    # Asignar el valor de G para el nodo ancestral
    G_value <- ifelse(as.character(node) %in% names(G_anc_states), G_anc_states[as.character(node)], NA)
    
    # Agregar el nodo ancestral, cluster y valor de G al resultado
    result_list <- append(result_list, list(data.frame(
      Node = node,
      Cluster = cluster,  # Asignamos el mismo cluster del terminal
      G_value = G_value
    )))
  }
}

# Combinar todos los dataframes en uno solo
result_df_G <- do.call(rbind, result_list)

# Para nodos ancestrales, asegurarse de asignar clusters
result_df_G <- result_df_G[!duplicated(result_df_G$Node), ]  # Eliminar filas duplicadas de cada nodo

# Ver el dataframe resultante
row.names(result_df_G)<-NULL
print(result_df_G)

# Expansión de result_df_G basado en paths (distancia desde la raíz)
expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  # Obtener el cluster y valores de G para cada nodo en el path
  cluster <- result_df_G$Cluster[match(path[length(path)], result_df_G$Node)]
  G_values <- result_df_G$G_value[match(path, result_df_G$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  # Crear un dataframe temporal para el path
  data.frame(Step = steps_from_root, Cluster = cluster, G_value = G_values, Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de G
  last_step <- max(df$Step)
  last_value <-tail(df$G_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    G_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))

# Graficar los valores de G a lo largo del camino de nodos desde la raíz para cada cluster
ggplot(expanded_df_extended, aes(x = Step, y = G_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de G a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de G",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

png("Figures/G_TRAIT.png",width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = G_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de G a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de G",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

dev.off()

################################## M ###########################################
# Cluster
clusters_table
# Tree
M_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.127973,LA_POLVADA:0.127973)[&CI={-0.106009,0.185709},ancstate={0.03985}]:0.123261,EL_DIVISADERO:0.251235)[&CI={-0.065124,0.298355},ancstate={0.116616}]:0.142584,RINCONADA_DE_PUNITAQUI:0.393819)[&CI={-0.022553,0.386529},ancstate={0.181988}]:0.077413,EL_ALTAR:0.471233)[&CI={-0.012875,0.419393},ancstate={0.203259}]:0.163097,CANELILLA_OVALLE:0.634329)[&CI={0.065793,0.549916},ancstate={0.307854}]:0.084733,MANQUEHUA:0.719062)[&CI={0.066873,0.579967},ancstate={0.32342}]:0.280938,((((((CANELA_BAJA:0.169277,CASTILLO_MAL_PASO:0.169277)[&CI={0.101203,0.370293},ancstate={0.235748}]:0.011386,PUNITAQUI:0.180663)[&CI={0.114032,0.382052},ancstate={0.248042}]:0.02304,MONTE_PATRIA:0.203702)[&CI={0.135969,0.422699},ancstate={0.279334}]:0.125087,BARRAZA:0.328789)[&CI={0.180911,0.555816},ancstate={0.368363}]:0.069566,LA_CALERA:0.398355)[&CI={0.204241,0.618122},ancstate={0.411182}]:0.215934,HUENTELAUQUEN:0.614289)[&CI={0.105868,0.633746},ancstate={0.369807}]:0.385711,EL_ESPINAL:1)[&CI={0.137808,0.731859},ancstate={0.434834}];")
paths
# Extraer los estados ancestrales
M_anc_states <- mvc$anc$ace
M_anc_states
#Extraer estados actuales
M_trait_vector <- setNames(M_trait$M, rownames(M_trait))
print(M_trait_vector)

#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster <- clusters_table[clusters_table$tip.label == tip, "cluster"]
  
  # Obtener el valor de S para el nodo terminal desde M_trait
  M_value_terminal <- M_trait_vector[tip]
  
  # Agregar el nodo terminal al dataframe con su valor de M y su cluster
  result_list <- append(result_list, list(data.frame(
    Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
    Cluster = cluster,
    M_value = M_value_terminal
  )))
  
  # Agregar cada nodo ancestral en el path (excluyendo el terminal)
  for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
    # Asignar el valor de S para el nodo ancestral
    M_value <- ifelse(as.character(node) %in% names(M_anc_states), M_anc_states[as.character(node)], NA)
    
    # Agregar el nodo ancestral, cluster y valor de S al resultado
    result_list <- append(result_list, list(data.frame(
      Node = node,
      Cluster = cluster,  # Asignamos el mismo cluster del terminal
      M_value = M_value
    )))
  }
}

# Combinar todos los dataframes en uno solo
result_df_M <- do.call(rbind, result_list)

# Para nodos ancestrales, asegurarse de asignar clusters
result_df_M <- result_df_M[!duplicated(result_df_M$Node), ]  # Eliminar filas duplicadas de cada nodo

# Ver el dataframe resultante
row.names(result_df_M)<-NULL
print(result_df_M)

# Expansión de result_df_M basado en paths (distancia desde la raíz)
expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  # Obtener el cluster y valores de M para cada nodo en el path
  cluster <- result_df_M$Cluster[match(path[length(path)], result_df_M$Node)]
  M_values <- result_df_M$M_value[match(path, result_df_M$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  # Crear un dataframe temporal para el path
  data.frame(Step = steps_from_root, Cluster = cluster, M_value = M_values, Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de M
  last_step <- max(df$Step)
  last_value <-tail(df$M_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    M_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))

# Graficar los valores de M a lo largo del camino de nodos desde la raíz para cada cluster
ggplot(expanded_df_extended, aes(x = Step, y = M_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de M a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de M",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

png("Figures/M_TRAIT.png",width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = M_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 1) +
  geom_point(size = 2) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de M a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de M",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))  # Colores para cada cluster

dev.off()


################################################################################
########################## Todas las comunidades ###############################
################################################################################
# Cluster
df_tree_total <- as.data.frame(as.matrix(surname_matrix))
df_tree_total <- as.data.frame(lapply(df_tree_total, as.numeric))
df_tree_scaled_t <- scale(df_tree_total)
km <- kmeans(df_tree_scaled_t, centers = 2, iter.max = 100, nstart = 100)
Cluster <- as.factor(c(km$cluster))
names(Cluster) <- colnames(km$centers)

# Crear una tabla de clusters
clusters_table <- data.frame(tip.label = names(Cluster), cluster = Cluster)
clusters_table

################################## S ###########################################
# Obtener los nodos terminales
terminal_nodes <- svt$obj$tree$tip.label

# Obtener las rutas de nodos desde la raíz
paths <- nodepath(svt$obj$tree, 1:169)

# Extraer los estados ancestrales
S_anc_states <- svt$anc$ace

# Extraer estados actuales
S_trait_vector <- setNames(S_trait2$S, rownames(S_trait2))
print(S_trait_vector)

# Crear un data.frame para almacenar resultados
result_list <- list() 

# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <- clusters_table[clusters_table$tip.label == tip, ]
  cluster <- ifelse(nrow(cluster_row) > 0, cluster_row$cluster, NA)
  
  # Obtener el valor de S para el nodo terminal desde S_trait
  S_value_terminal <- ifelse(tip %in% names(S_trait_vector), S_trait_vector[tip], NA)
  
  # Verificar si hay un cluster y un valor de S válido
  if (!is.na(cluster) && !is.na(S_value_terminal)) {
    # Agregar el nodo terminal al dataframe con su valor de S y su cluster
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      S_value = S_value_terminal
    )
    
    # Agregar cada nodo ancestral en el path (excluyendo el terminal)
    for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
      # Asignar el valor de S para el nodo ancestral
      S_value <- ifelse(as.character(node) %in% names(S_anc_states), 
                        S_anc_states[as.character(node)], NA)
      
      # Asegurarse de que el nodo ancestral tenga un cluster asignado
      # Asignar el cluster del terminal (debe estar asignado ya)
      if (!is.na(cluster) && !is.na(S_value)) {
        result_list[[length(result_list) + 1]] <- data.frame(
          Node = node,
          Cluster = cluster,  # Asignamos el mismo cluster del terminal
          S_value = S_value
        )
      }
    }
  }
}

# Combinar todos los dataframes en uno solo
result_df_S <- do.call(rbind, result_list)

# Eliminar filas duplicadas de cada nodo (manteniendo el primero)
result_df_S <- result_df_S[!duplicated(result_df_S$Node), ]

# Mesetear los nombres de las filas
row.names(result_df_S) <- NULL

# Ver el dataframe resultante
print(result_df_S)

expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  cluster <- result_df_S$Cluster[match(path[length(path)], result_df_S$Node)]
  S_values <- result_df_S$S_value[match(path, result_df_S$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  data.frame(Step = steps_from_root, 
             Cluster = cluster, 
             S_value = S_values, 
             Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
max_steps
# Eliminar filas donde S_value es NA
expanded_df <- expanded_df[!is.na(expanded_df$S_value), ]

# Asegurarse de que el Cluster no tenga NA
expanded_df <- expanded_df[!is.na(expanded_df$Cluster), ]

expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de S
  last_step <- max(df$Step)
  last_value <-tail(df$S_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    S_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))
ggplot(expanded_df_extended, aes(x = Step, y = S_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de S a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de S",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster

# Graficar los valores de S a lo largo del camino de nodos desde la raíz para cada cluster
png("Figures/S_TRAIT_TOTAL.png", width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = S_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de S a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de S",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster
dev.off()

################################## R ###########################################
# Obtener los nodos terminales
terminal_nodes <- rvt$obj$tree$tip.label

# Obtener las rutas de nodos desde la raíz
paths <- nodepath(rvt$obj$tree, 1:169)

# Extraer los estados ancestrales
R_anc_states <- rvt$anc$ace

# Extraer estados actuales
R_trait_vector <- setNames(R_trait2$R, rownames(R_trait2))
print(R_trait_vector)

# Crear un data.frame para almacenar resultados
result_list <- list() 

# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <- clusters_table[clusters_table$tip.label == tip, ]
  cluster <- ifelse(nrow(cluster_row) > 0, cluster_row$cluster, NA)
  
  # Obtener el valor de R para el nodo terminal desde R_trait
  R_value_terminal <- ifelse(tip %in% names(R_trait_vector), R_trait_vector[tip], NA)
  
  # Verificar si hay un cluster y un valor de R válido
  if (!is.na(cluster) && !is.na(R_value_terminal)) {
    # Agregar el nodo terminal al dataframe con su valor de R y su cluster
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      R_value = R_value_terminal
    )
    
    # Agregar cada nodo ancestral en el path (excluyendo el terminal)
    for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
      # Asignar el valor de R para el nodo ancestral
      R_value <- ifelse(as.character(node) %in% names(R_anc_states), 
                        R_anc_states[as.character(node)], NA)
      
      # Asegurarse de que el nodo ancestral tenga un cluster asignado
      # Asignar el cluster del terminal (debe estar asignado ya)
      if (!is.na(cluster) && !is.na(R_value)) {
        result_list[[length(result_list) + 1]] <- data.frame(
          Node = node,
          Cluster = cluster,  # Asignamos el mismo cluster del terminal
          R_value = R_value
        )
      }
    }
  }
}

# Combinar todos los dataframes en uno solo
result_df_R <- do.call(rbind, result_list)

# Eliminar filas duplicadas de cada nodo (manteniendo el primero)
result_df_R <- result_df_R[!duplicated(result_df_R$Node), ]

# Resetear los nombres de las filas
row.names(result_df_R) <- NULL

# Ver el dataframe resultante
print(result_df_R)

expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  cluster <- result_df_R$Cluster[match(path[length(path)], result_df_R$Node)]
  R_values <- result_df_R$R_value[match(path, result_df_R$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  data.frame(Step = steps_from_root, 
             Cluster = cluster, 
             R_value = R_values, 
             Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
max_steps
# Eliminar filas donde R_value es NA
expanded_df <- expanded_df[!is.na(expanded_df$R_value), ]

# Asegurarse de que el Cluster no tenga NA
expanded_df <- expanded_df[!is.na(expanded_df$Cluster), ]

expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de R
  last_step <- max(df$Step)
  last_value <-tail(df$R_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    R_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))
ggplot(expanded_df_extended, aes(x = Step, y = R_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de R a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de R",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster

# Graficar los valores de R a lo largo del camino de nodos desde la raíz para cada cluster
png("Figures/R_TRAIT_TOTAL.png", width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = R_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de R a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de R",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster
dev.off()

################################## A ###########################################
# Obtener los nodos terminales
terminal_nodes <- avt$obj$tree$tip.label

# Obtener las rutas de nodos desde la raíz
paths <- nodepath(avt$obj$tree, 1:170)

# Extraer los estados ancestrales
A_anc_states <- avt$anc$ace

# Extraer estados actuales
A_trait_vector <- setNames(A_trait2$A, rownames(A_trait2))
print(A_trait_vector)

# Crear un data.frame para almacenar resultados
result_list <- list() 

# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <- clusters_table[clusters_table$tip.label == tip, ]
  cluster <- ifelse(nrow(cluster_row) > 0, cluster_row$cluster, NA)
  
  # Obtener el valor de A para el nodo terminal desde A_trait
  A_value_terminal <- ifelse(tip %in% names(A_trait_vector), A_trait_vector[tip], NA)
  
  # Verificar si hay un cluster y un valor de A válido
  if (!is.na(cluster) && !is.na(A_value_terminal)) {
    # Agregar el nodo terminal al dataframe con su valor de A y su cluster
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      A_value = A_value_terminal
    )
    
    # Agregar cada nodo ancestral en el path (excluyendo el terminal)
    for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
      # Asignar el valor de A para el nodo ancestral
      A_value <- ifelse(as.character(node) %in% names(A_anc_states), 
                        A_anc_states[as.character(node)], NA)
      
      # Asegurarse de que el nodo ancestral tenga un cluster asignado
      # Asignar el cluster del terminal (debe estar asignado ya)
      if (!is.na(cluster) && !is.na(A_value)) {
        result_list[[length(result_list) + 1]] <- data.frame(
          Node = node,
          Cluster = cluster,  # Asignamos el mismo cluster del terminal
          A_value = A_value
        )
      }
    }
  }
}

# Combinar todos los dataframes en uno solo
result_df_A <- do.call(rbind, result_list)

# Eliminar filas duplicadas de cada nodo (manteniendo el primero)
result_df_A <- result_df_A[!duplicated(result_df_A$Node), ]

# Resetear los nombres de las filas
row.names(result_df_A) <- NULL

# Ver el dataframe resultante
print(result_df_A)

expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  cluster <- result_df_A$Cluster[match(path[length(path)], result_df_A$Node)]
  A_values <- result_df_A$A_value[match(path, result_df_A$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  data.frame(Step = steps_from_root, 
             Cluster = cluster, 
             A_value = A_values, 
             Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
max_steps
# Eliminar filas donde A_value es NA
expanded_df <- expanded_df[!is.na(expanded_df$A_value), ]

# Asegurarse de que el Cluster no tenga NA
expanded_df <- expanded_df[!is.na(expanded_df$Cluster), ]

expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de A
  last_step <- max(df$Step)
  last_value <-tail(df$A_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    A_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))
ggplot(expanded_df_extended, aes(x = Step, y = A_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de A a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de A",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster

# Graficar los valores de A a lo largo del camino de nodos desde la raíz para cada cluster
png("Figures/A_TRAIT_TOTAL.png", width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = A_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de A a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de A",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster
dev.off()

################################## G ###########################################
# Obtener los nodos terminales
terminal_nodes <- gvt$obj$tree$tip.label

# Obtener las rutas de nodos desde la raíz
paths <- nodepath(gvt$obj$tree, 1:170)

# Extraer los estados ancestrales
G_anc_states <- gvt$anc$ace

# Extraer estados actuales
G_trait_vector <- setNames(G_trait2$G, rownames(G_trait2))
print(G_trait_vector)

# Crear un data.frame para almacenar resultados
result_list <- list() 

# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <- clusters_table[clusters_table$tip.label == tip, ]
  cluster <- ifelse(nrow(cluster_row) > 0, cluster_row$cluster, NA)
  
  # Obtener el valor de G para el nodo terminal desde G_trait
  G_value_terminal <- ifelse(tip %in% names(G_trait_vector), G_trait_vector[tip], NA)
  
  # Verificar si hay un cluster y un valor de G válido
  if (!is.na(cluster) && !is.na(G_value_terminal)) {
    # Agregar el nodo terminal al dataframe con su valor de G y su cluster
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      G_value = G_value_terminal
    )
    
    # Agregar cada nodo ancestral en el path (excluyendo el terminal)
    for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
      # Asignar el valor de G para el nodo ancestral
      G_value <- ifelse(as.character(node) %in% names(G_anc_states), 
                        G_anc_states[as.character(node)], NA)
      
      # Asegurarse de que el nodo ancestral tenga un cluster asignado
      # Asignar el cluster del terminal (debe estar asignado ya)
      if (!is.na(cluster) && !is.na(G_value)) {
        result_list[[length(result_list) + 1]] <- data.frame(
          Node = node,
          Cluster = cluster,  # Asignamos el mismo cluster del terminal
          G_value = G_value
        )
      }
    }
  }
}

# Combinar todos los dataframes en uno solo
result_df_G <- do.call(rbind, result_list)

# Eliminar filas duplicadas de cada nodo (manteniendo el primero)
result_df_G <- result_df_G[!duplicated(result_df_G$Node), ]

# Resetear los nombres de las filas
row.names(result_df_G) <- NULL

# Ver el dataframe resultante
print(result_df_G)

expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  cluster <- result_df_G$Cluster[match(path[length(path)], result_df_G$Node)]
  G_values <- result_df_G$G_value[match(path, result_df_G$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  data.frame(Step = steps_from_root, 
             Cluster = cluster, 
             G_value = G_values, 
             Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
max_steps
# Eliminar filas donde G_value es NA
expanded_df <- expanded_df[!is.na(expanded_df$G_value), ]

# Asegurarse de que el Cluster no tenga NA
expanded_df <- expanded_df[!is.na(expanded_df$Cluster), ]

expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de G
  last_step <- max(df$Step)
  last_value <-tail(df$G_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    G_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))
ggplot(expanded_df_extended, aes(x = Step, y = G_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de G a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de G",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster

# Graficar los valores de G a lo largo del camino de nodos desde la raíz para cada cluster
png("Figures/G_TRAIT_TOTAL.png", width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = G_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de G a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de G",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster
dev.off()

################################## M ###########################################
# Obtener los nodos terminales
terminal_nodes <- mvt$obj$tree$tip.label

# Obtener las rutas de nodos desde la raíz
paths <- nodepath(mvt$obj$tree, 1:170)

# Extraer los estados ancestrales
M_anc_states <- mvt$anc$ace

# Extraer estados actuales
M_trait_vector <- setNames(M_trait2$M, rownames(M_trait2))
print(M_trait_vector)

# Crear un data.frame para almacenar resultados
result_list <- list() 

# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <- clusters_table[clusters_table$tip.label == tip, ]
  cluster <- ifelse(nrow(cluster_row) > 0, cluster_row$cluster, NA)
  
  # Obtener el valor de M para el nodo terminal desde M_trait
  M_value_terminal <- ifelse(tip %in% names(M_trait_vector), M_trait_vector[tip], NA)
  
  # Verificar si hay un cluster y un valor de M válido
  if (!is.na(cluster) && !is.na(M_value_terminal)) {
    # Agregar el nodo terminal al dataframe con su valor de M y su cluster
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      M_value = M_value_terminal
    )
    
    # Agregar cada nodo ancestral en el path (excluyendo el terminal)
    for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
      # Asignar el valor de M para el nodo ancestral
      M_value <- ifelse(as.character(node) %in% names(M_anc_states), 
                        M_anc_states[as.character(node)], NA)
      
      # Asegurarse de que el nodo ancestral tenga un cluster asignado
      # Asignar el cluster del terminal (debe estar asignado ya)
      if (!is.na(cluster) && !is.na(M_value)) {
        result_list[[length(result_list) + 1]] <- data.frame(
          Node = node,
          Cluster = cluster,  # Asignamos el mismo cluster del terminal
          M_value = M_value
        )
      }
    }
  }
}

# Combinar todos los dataframes en uno solo
result_df_M <- do.call(rbind, result_list)

# Eliminar filas duplicadas de cada nodo (manteniendo el primero)
result_df_M <- result_df_M[!duplicated(result_df_M$Node), ]

# Mesetear los nombres de las filas
row.names(result_df_M) <- NULL

# Ver el dataframe resultante
print(result_df_M)

expanded_df <- do.call(rbind, lapply(1:length(paths), function(i) {
  path <- paths[[i]]
  cluster <- result_df_M$Cluster[match(path[length(path)], result_df_M$Node)]
  M_values <- result_df_M$M_value[match(path, result_df_M$Node)]
  # Crear una columna de pasos desde la raíz (distancia en pasos desde el nodo 16)
  steps_from_root <- seq(0, length(path)-1)
  data.frame(Step = steps_from_root, 
             Cluster = cluster, 
             M_value = M_values, 
             Terminal_Node = path[length(path)])
}))
max_steps <- max((expanded_df$Step))
max_steps
# Eliminar filas donde M_value es NA
expanded_df <- expanded_df[!is.na(expanded_df$M_value), ]

# Asegurarse de que el Cluster no tenga NA
expanded_df <- expanded_df[!is.na(expanded_df$Cluster), ]

expanded_df_extended <- do.call(rbind, lapply(split(expanded_df, expanded_df$Terminal_Node), function(df) {
  # Calcular valores adicionales de M
  last_step <- max(df$Step)
  last_value <-tail(df$M_value, 1)
  cluster <- unique(df$Cluster)
  
  # Generar filas adicionales hasta el máximo de pasos deseado
  additional_steps <- data.frame(
    Step = seq(last_step, max_steps),
    Cluster = cluster,
    M_value = last_value,  # Usamos el último valor
    Terminal_Node = df$Terminal_Node[1]
  )
  
  # Combinar el dataframe original con los pasos adicionales
  rbind(df, additional_steps)
}))
ggplot(expanded_df_extended, aes(x = Step, y = M_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de M a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de M",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster

# Graficar los valores de M a lo largo del camino de nodos desde la raíz para cada cluster
png("Figures/M_TRAIT_TOTAL.png", width = 2000, height = 750, res = 300)
ggplot(expanded_df_extended, aes(x = Step, y = M_value, color = as.factor(Cluster), group = interaction(Cluster, Terminal_Node))) +
  geom_line(size = 0.1, alpha = 0.6) +
  geom_point(size = 0.1) +  # Añadir puntos en cada nodo
  geom_smooth(aes(group = Cluster), method = "lm", se = FALSE, size = 1.2, linetype = "dashed") +  # Línea de tendencia
  labs(
    title = "Valores de M a lo largo de los nodos",
    x = "Paso desde la raíz",
    y = "Valor de M",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green"))  # Colores para cada cluster
dev.off()

