##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################
##### PAMT 2 #####

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
S_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.128179,LA_POLVADA:0.128179)[&CI={0.211964,0.375293},ancstate={0.293629}]:0.13116,MINCONADA_DE_PUNITAQUI:0.259339)[&CI={0.179189,0.384707},ancstate={0.281948}]:0.148585,EL_DIVISADEMO:0.407924)[&CI={0.197171,0.424892},ancstate={0.311031}]:0.057821,EL_ALTAM:0.465745)[&CI={0.194064,0.430633},ancstate={0.312349}]:0.17539,CANELILLA_OVALLE:0.641135)[&CI={0.17478,0.440191},ancstate={0.307485}]:0.040564,EL_ESPINAL:0.681698)[&CI={0.167812,0.441244},ancstate={0.304528}]:0.318302,(((((PUNITAQUI:0.042605,CANELA_BAJA:0.042605)[&CI={0.110537,0.208574},ancstate={0.159556}]:0.066623,CASTILLO_MAL_PASO:0.109227)[&CI={0.107083,0.256815},ancstate={0.181949}]:0.242948,MONTE_PATMIA:0.352175)[&CI={0.076137,0.313116},ancstate={0.194627}]:0.131883,(BAMMAZA:0.264236,LA_CALEMA:0.264236)[&CI={0.117594,0.344242},ancstate={0.230918}]:0.219821)[&CI={0.084848,0.336995},ancstate={0.210921}]:0.219526,MANQUEHUA:0.703583)[&CI={0.067191,0.368958},ancstate={0.218074}]:0.296417,HUENTELAUQUEN:1)[&CI={0.093097,0.419418},ancstate={0.256257}];
")
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

# Para nodos ancestrales, asegurarse de asignar clusters
result_df_S <- result_df_S[!duplicated(result_df_S$Node), ]  # Eliminar filas duplicadas de cada nodo

# Ver el dataframe resultante
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

png("Figures/S_TMAIT.png",width = 2000, height = 750, res = 300)
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
  
################################## M ###########################################
# Cluster
clusters_table
# Tree
M_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.128179,LA_POLVADA:0.128179)[&CI={0.871058,1.147971},ancstate={1.009515}]:0.13116,MINCONADA_DE_PUNITAQUI:0.259339)[&CI={0.854767,1.203206},ancstate={1.028986}]:0.148585,EL_DIVISADEMO:0.407924)[&CI={0.738828,1.124912},ancstate={0.93187}]:0.057821,EL_ALTAM:0.465745)[&CI={0.700664,1.101748},ancstate={0.901206}]:0.17539,CANELILLA_OVALLE:0.641135)[&CI={0.528749,0.978733},ancstate={0.753741}]:0.040564,EL_ESPINAL:0.681698)[&CI={0.500054,0.963638},ancstate={0.731846}]:0.318302,(((((PUNITAQUI:0.042605,CANELA_BAJA:0.042605)[&CI={0.810824,0.977038},ancstate={0.893931}]:0.066623,CASTILLO_MAL_PASO:0.109227)[&CI={0.774211,1.028071},ancstate={0.901141}]:0.242948,MONTE_PATMIA:0.352175)[&CI={0.568442,0.970221},ancstate={0.769332}]:0.131883,(BAMMAZA:0.264236,LA_CALEMA:0.264236)[&CI={0.415893,0.800157},ancstate={0.608025}]:0.219821)[&CI={0.514279,0.941775},ancstate={0.728027}]:0.219526,MANQUEHUA:0.703583)[&CI={0.523303,1.034927},ancstate={0.779115}]:0.296417,HUENTELAUQUEN:1)[&CI={0.50012,1.053371},ancstate={0.776746}];")
#ancestral_nodes <- Ancestors(M_tree,1:15,"all")
paths
# Extraer los estados ancestrales
M_anc_states <- rvc$anc$ace
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
  
  # Obtener el valor de M para el nodo terminal desde M_trait
  M_value_terminal <- M_trait_vector[tip]
  
  # Agregar el nodo terminal al dataframe con su valor de M y su cluster
  result_list <- append(result_list, list(data.frame(
    Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
    Cluster = cluster,
    M_value = M_value_terminal
  )))
  
  # Agregar cada nodo ancestral en el path (excluyendo el terminal)
  for (node in head(path, -1)) {  # head(path, -1) excluye el último nodo terminal
    # Asignar el valor de M para el nodo ancestral
    M_value <- ifelse(as.character(node) %in% names(M_anc_states), M_anc_states[as.character(node)], NA)
    
    # Agregar el nodo ancestral, cluster y valor de M al resultado
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

# Expansión de result_df_S basado en paths (distancia desde la raíz)
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

png("Figures/M_TMAIT.png",width = 2000, height = 750, res = 300)
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
################################## A ###########################################
# Cluster
clusters_table
# Tree
A_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.128179,LA_POLVADA:0.128179)[&CI={0.514536,0.581409},ancstate={0.547973}]:0.13116,MINCONADA_DE_PUNITAQUI:0.259339)[&CI={0.515488,0.599635},ancstate={0.557561}]:0.148585,EL_DIVISADEMO:0.407924)[&CI={0.523321,0.616559},ancstate={0.56994}]:0.057821,EL_ALTAM:0.465745)[&CI={0.530646,0.627506},ancstate={0.579076}]:0.17539,CANELILLA_OVALLE:0.641135)[&CI={0.517552,0.626221},ancstate={0.571886}]:0.040564,EL_ESPINAL:0.681698)[&CI={0.513769,0.625723},ancstate={0.569746}]:0.318302,(((((PUNITAQUI:0.042605,CANELA_BAJA:0.042605)[&CI={0.557544,0.597684},ancstate={0.577614}]:0.066623,CASTILLO_MAL_PASO:0.109227)[&CI={0.536124,0.59743},ancstate={0.566777}]:0.242948,MONTE_PATMIA:0.352175)[&CI={0.493407,0.590435},ancstate={0.541921}]:0.131883,(BAMMAZA:0.264236,LA_CALEMA:0.264236)[&CI={0.504465,0.597264},ancstate={0.550864}]:0.219821)[&CI={0.494888,0.598126},ancstate={0.546507}]:0.219526,MANQUEHUA:0.703583)[&CI={0.488011,0.611567},ancstate={0.549789}]:0.296417,HUENTELAUQUEN:1)[&CI={0.48735,0.620958},ancstate={0.554154}];")
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

png("Figures/A_TMAIT.png",width = 2000, height = 750, res = 300)
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
G_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.128179,LA_POLVADA:0.128179)[&CI={-0.046512,0.145575},ancstate={0.049532}]:0.13116,MINCONADA_DE_PUNITAQUI:0.259339)[&CI={0.030048,0.271751},ancstate={0.150899}]:0.148585,EL_DIVISADEMO:0.407924)[&CI={0.03501,0.302825},ancstate={0.168917}]:0.057821,EL_ALTAM:0.465745)[&CI={0.040635,0.318857},ancstate={0.179746}]:0.17539,CANELILLA_OVALLE:0.641135)[&CI={0.09487,0.407012},ancstate={0.250941}]:0.040564,EL_ESPINAL:0.681698)[&CI={0.091494,0.413069},ancstate={0.252282}]:0.318302,(((((PUNITAQUI:0.042605,CANELA_BAJA:0.042605)[&CI={0.105225,0.220523},ancstate={0.162874}]:0.066623,CASTILLO_MAL_PASO:0.109227)[&CI={0.120508,0.296603},ancstate={0.208556}]:0.242948,MONTE_PATMIA:0.352175)[&CI={0.136283,0.414986},ancstate={0.275634}]:0.131883,(BAMMAZA:0.264236,LA_CALEMA:0.264236)[&CI={0.25067,0.517224},ancstate={0.383947}]:0.219821)[&CI={0.155646,0.452188},ancstate={0.303917}]:0.219526,MANQUEHUA:0.703583)[&CI={0.093622,0.448522},ancstate={0.271072}]:0.296417,HUENTELAUQUEN:1)[&CI={0.057334,0.44111},ancstate={0.249222}];")
#ancestral_nodes <- Ancestors(G_tree,1:15,"all")
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

png("Figures/G_TMAIT.png",width = 2000, height = 750, res = 300)
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
M_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.128179,LA_POLVADA:0.128179)[&CI={-0.103305,0.202191},ancstate={0.049443}]:0.13116,MINCONADA_DE_PUNITAQUI:0.259339)[&CI={-0.041574,0.342832},ancstate={0.150629}]:0.148585,EL_DIVISADEMO:0.407924)[&CI={-0.007128,0.418807},ancstate={0.20584}]:0.057821,EL_ALTAM:0.465745)[&CI={0.012878,0.455363},ancstate={0.234121}]:0.17539,CANELILLA_OVALLE:0.641135)[&CI={0.148357,0.644789},ancstate={0.396573}]:0.040564,EL_ESPINAL:0.681698)[&CI={0.165674,0.677109},ancstate={0.421392}]:0.318302,(((((PUNITAQUI:0.042605,CANELA_BAJA:0.042605)[&CI={0.07256,0.255932},ancstate={0.164246}]:0.066623,CASTILLO_MAL_PASO:0.109227)[&CI={0.052299,0.332363},ancstate={0.192331}]:0.242948,MONTE_PATMIA:0.352175)[&CI={0.109608,0.552859},ancstate={0.331234}]:0.131883,(BAMMAZA:0.264236,LA_CALEMA:0.264236)[&CI={0.283215,0.707143},ancstate={0.495179}]:0.219821)[&CI={0.140948,0.61257},ancstate={0.376759}]:0.219526,MANQUEHUA:0.703583)[&CI={0.052061,0.616495},ancstate={0.334278}]:0.296417,HUENTELAUQUEN:1)[&CI={0.040801,0.651159},ancstate={0.34598}];")
#ancestral_nodes <- Ancestors(S_tree,1:15,"all")
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

png("Figures/M_TMAIT.png",width = 2000, height = 750, res = 300)
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
png("Figures/S_TMAIT_TOTAL.png", width = 2000, height = 750, res = 300)
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

################################## M ###########################################
# Obtener los nodos terminales
terminal_nodes <- rvt$obj$tree$tip.label

# Obtener las rutas de nodos desde la raíz
paths <- nodepath(rvt$obj$tree, 1:169)

# Extraer los estados ancestrales
M_anc_states <- rvt$anc$ace

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
png("Figures/M_TMAIT_TOTAL.png", width = 2000, height = 750, res = 300)
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

# Mesetear los nombres de las filas
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
png("Figures/A_TMAIT_TOTAL.png", width = 2000, height = 750, res = 300)
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

# Mesetear los nombres de las filas
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
png("Figures/G_TMAIT_TOTAL.png", width = 2000, height = 750, res = 300)
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
