##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 5 ####
### Part 5 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###
##IMPORTANTE: CORRER LOS SCRIPTS DE LOS OBJETIVOS 2 y 4 ##

#### Cargar paquetes y librerias ####
library(ape)
library(dplyr)
library(dendextend)
library(phytools)
library(ggplot2)

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
## Necesito que en el eje x estén los clados promediados ##
## En el eje y1 deben estar los valores de G y en ele eje y2 los valores de M ##
################################## G ###########################################
# Cluster
clusters_table
# Tree
G_tree <- read.tree(text = "(((((((GUALLIGUAICA:0.127973,LA_POLVADA:0.127973)[&CI={-0.0521,0.128952},ancstate={0.038426}]:0.123261,EL_DIVISADERO:0.251235)[&CI={-0.000347,0.225242},ancstate={0.112448}]:0.142584,RINCONADA_DE_PUNITAQUI:0.393819)[&CI={0.054363,0.308255},ancstate={0.181309}]:0.077413,EL_ALTAR:0.471233)[&CI={0.057315,0.325597},ancstate={0.191456}]:0.163097,CANELILLA_OVALLE:0.634329)[&CI={0.1019,0.402365},ancstate={0.252132}]:0.084733,MANQUEHUA:0.719062)[&CI={0.092656,0.411103},ancstate={0.251879}]:0.280938,((((((CANELA_BAJA:0.169277,CASTILLO_MAL_PASO:0.169277)[&CI={0.138494,0.305502},ancstate={0.221998}]:0.011386,PUNITAQUI:0.180663)[&CI={0.144077,0.31042},ancstate={0.227248}]:0.02304,MONTE_PATRIA:0.203702)[&CI={0.154424,0.33238},ancstate={0.243402}]:0.125087,BARRAZA:0.328789)[&CI={0.181639,0.41432},ancstate={0.29798}]:0.069566,LA_CALERA:0.398355)[&CI={0.190979,0.44785},ancstate={0.319415}]:0.215934,HUENTELAUQUEN:0.614289)[&CI={0.111258,0.43888},ancstate={0.275069}]:0.385711,EL_ESPINAL:1)[&CI={0.080061,0.448753},ancstate={0.264407}];")
terminal_nodes <-G_tree$tip.label
terminal_nodes

#Path
paths <- nodepath(G_tree,1:15)
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

################################## M ###########################################
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
      Cluster = cluster,
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

################################## G-M #########################################
#Dataframe de ambos
result_df_GM<-merge(result_df_G, 
                    result_df_M,
                    all.x = T,
                    by = c("Node", "Cluster"))
result_df_GM
max_path_length<-max(sapply(paths,length))
max_path_length
extended_paths <- lapply(paths, function(path){
  if (length(path)<max_path_length){
    c(path,rep(tail(path,1),max_path_length - length(path)))
  }else{
    path
      }
})
extended_paths

calc_means<-function(nodos,df) {
  #Extraer valores de G y M
  G_values <- df$G_value[df$Node %in% nodos]
  M_values <- df$M_value[df$Node %in% nodos]
  #Calcular promedio- NA == 0
  G_mean<-median(G_values)
  M_mean<-median(M_values)
  G_mean<-ifelse(is.na(G_values),0,G_values)
  M_mean<-ifelse(is.na(M_values),0,M_values)
  return(data.frame(G_mean=G_mean,M_mean=M_mean))
}

levels_means<-list()
for(level in 1:max_path_length){
  nodos_at_level<-sapply(extended_paths, function(path)path[level])
  mean_2<-calc_means(nodos_at_level,result_df_GM)
  levels_means[[level]]<-data.frame(Level=level,G_mean=mean_2$G_mean,M_mean=mean_2$M_mean)
}
level_means_df<- do.call(rbind,levels_means)
print(level_means_df)
averages <- level_means_df %>%
  group_by(Level) %>%
  summarise(across(everything(),median))
averages

hist_data<-ggplot_build(ggplot(averages, aes(x = Level)) + 
  geom_histogram(binwidth = 0.5))$data[[1]]
y_max <-max(hist_data$count)

# Crear el gráfico
ggplot(averages, aes(x = Level)) + 
  geom_point(aes(y=G_mean), color="blue") + # Puntos para G
  geom_line(aes(y=G_mean), color="blue",linetype="dashed")+ #Línea para G
  geom_point(aes(y=M_mean), color="red") + # Puntos para M
  geom_line(aes(y=M_mean), color="red",linetype="dashed")+ #Línea para M
  scale_y_continuous(
    name = "G",
    sec.axis = sec_axis(~.*2,name = "M")
    ) +
  annotate("text",x=mean(averages$Level),y=mean(averages$G_mean),label="G trait",hjust=0.2,vjust=-4.2,size=5,color="blue")+
  annotate("text",x=mean(averages$Level),y=mean(averages$M_mean),label="M trait",hjust=0.1,vjust=-7,size=5,color="red")+
  labs(
    x="Number of clades",
    title = "The G and M traits of the tree of sampled communities",
  ) +
  theme_classic() 

########################## Líneas por cluster ##################################
# Filtrar el nodo 16 (raiz)
node_16_G <- result_df_G[result_df_G$Node == 16, ]
node_16_M <- result_df_M[result_df_M$Node == 16, ]
# Crear tres copias cambiando el valor del Cluster
node_16_cluster_1G <- node_16_G
node_16_cluster_1G$Cluster <- 1

node_16_cluster_2G <- node_16_G
node_16_cluster_2G$Cluster <- 2

node_16_cluster_3G <- node_16_G
node_16_cluster_3G$Cluster <- 3

node_16_cluster_1M <- node_16_M
node_16_cluster_1M$Cluster <- 1

node_16_cluster_2M <- node_16_M
node_16_cluster_2M$Cluster <- 2

node_16_cluster_3M <- node_16_M
node_16_cluster_3M$Cluster <- 3

# Unir las filas modificadas al dataframe original
result_df_G <- rbind(result_df_G, node_16_cluster_1G, node_16_cluster_2G, node_16_cluster_3G)
result_df_G
result_df_M <- rbind(result_df_M, node_16_cluster_1M,node_16_cluster_2M, node_16_cluster_3M)
result_df_M
#Eliminar duplicados
result_df_G <- result_df_G %>% distinct()
result_df_M <- result_df_M %>% distinct()

#Dataframe de ambos
result_df_GM <- merge(result_df_G, 
                      result_df_M, 
                      all.x = TRUE, 
                      by = c("Node", "Cluster"))

result_df_GM
max_path_length<-max(sapply(paths,length))
max_path_length
extended_paths <- lapply(paths, function(path){
  if (length(path)<max_path_length){
    c(path,rep(tail(path,1),max_path_length - length(path)))
  }else{
    path
  }
})
extended_paths

calc_means <- function(nodos, df) {
  # Filtrar los datos para los nodos en el nivel actual
  df_filtered <- df %>% filter(Node %in% nodos)
  
  # Calcular promedios por cluster
  G_means <- df_filtered %>%
    group_by(Cluster) %>%
    summarise(G_mean = median(G_value, na.rm = TRUE), .groups = 'drop')
  
  M_means <- df_filtered %>%
    group_by(Cluster) %>%
    summarise(M_mean = median(M_value, na.rm = TRUE), .groups = 'drop')
  
  # Combinar promedios de G y M para cada cluster
  cluster_means <- merge(G_means, M_means, by = "Cluster", all = TRUE)
  
  # Llenar NAs con 0
  cluster_means[is.na(cluster_means)] <- 0
  
  return(cluster_means)
}

levels_means <- list()
for (level in 1:max_path_length) {
  # Obtener nodos en el nivel actual
  nodos_at_level <- sapply(extended_paths, function(path) path[level])
  
  # Calcular los promedios para cada cluster en este nivel
  mean_data <- calc_means(nodos_at_level, result_df_GM)
  
  # Agregar la columna de nivel
  mean_data$Level <- level
  
  # Almacenar resultados en la lista
  levels_means[[level]] <- mean_data
}

# Convertir lista en dataframe
level_means_df <- do.call(rbind, levels_means)
print(level_means_df)

# Calcular promedio general por nivel
averages_cluster <- level_means_df %>%
  group_by(Level, Cluster) %>%
  summarise(across(starts_with("G_mean"), median, na.rm = TRUE),
            across(starts_with("M_mean"), median, na.rm = TRUE),
            .groups = 'drop')
averages_cluster

# Crear el gráfico con ggplot2
ggplot(averages_cluster) +
  # Línea continua para G_mean
  geom_line(aes(x = Level, y = G_mean, color = factor(Cluster), linetype = "G_mean", group = Cluster), 
            size = 1) +
  # Línea punteada para M_mean
  geom_line(aes(x = Level, y = M_mean, color = factor(Cluster), linetype = "M_mean", group = Cluster), 
            size = 1) +
  # Configuración de colores
  scale_color_manual(values = c("blue", "red", "green")) +
  # Configuración de tipos de línea
  scale_linetype_manual(values = c("G_mean" = "solid", "M_mean" = "dotted")) +
  labs(
    title = "The G and M traits of the tree of sampled communities",
    x = "Number of clades",
    y = "G mean",
    color = "Cluster",
    linetype = "Trait"
  ) +
  # Agregar eje secundario para M
  scale_y_continuous(
    sec.axis = sec_axis(~ ., name = "M mean")
  ) +
  theme_classic()
 
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

# Tree
terminal_nodes <-gvt$obj$tree$tip.label
terminal_nodes
#Path
paths <- nodepath(gvt$obj$tree,1:170)
paths
# Extraer los estados ancestrales
G_anc_states <- gvt$anc$ace
G_anc_states
#Extraer estados actuales
G_trait_vector <- setNames(G_trait2$G, rownames(G_trait2))
print(G_trait_vector)

#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <-clusters_table[clusters_table$tip.label==tip,]
  cluster <- ifelse(nrow(cluster_row)>0,cluster_row$cluster,NA)
  # Obtener el valor de M para el nodo terminal desde M_trait
  M_value_terminal <- ifelse(tip %in% names(M_trait_vector),M_trait_vector[tip],NA)
  if(!is.na(cluster)&&!is.na(G_value_terminal)){
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      G_value = G_value_terminal)
  for(node in head(path,-1)){
    G_value<-ifelse(as.character(node)%in% names(G_anc_states),
                    G_anc_states[as.character(node)],NA)
    if(!is.na(cluster)&&!is.na(G_value)) {
      result_list[[length(result_list)+1]]<-data.frame(
        Node =node,
        Cluster = cluster,
        G_value = G_value
        )
      }
    }  
  }
}


  
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


# Combinar todos los dataframes en uno solo
result_df_G <- do.call(rbind, result_list)

# Para nodos ancestrales, asegurarse de asignar clusters
result_df_G <- result_df_G[!duplicated(result_df_G$Node), ]  # Eliminar filas duplicadas de cada nodo

# Ver el dataframe resultante
row.names(result_df_G)<-NULL
print(result_df_G)
################################## M ###########################################
paths
# Extraer los estados ancestrales
M_anc_states <- mvt$anc$ace
M_anc_states
#Extraer estados actuales
M_trait_vector <- setNames(M_trait2$M, rownames(M_trait2))
print(M_trait_vector)

#Crear un data.frame
result_list <- list() 
# Iterar sobre cada nodo terminal y su ruta desde la raíz
for (i in seq_along(paths)) {
  # Nodo terminal y ruta de nodos desde la raíz hasta el nodo terminal
  tip <- terminal_nodes[i]
  path <- paths[[i]]
  
  # Obtener el cluster del nodo terminal desde clusters_table
  cluster_row <-clusters_table[clusters_table$tip.label==tip,]
  cluster <- ifelse(nrow(cluster_row)>0,cluster_row$cluster,NA)
  # Obtener el valor de M para el nodo terminal desde M_trait
  M_value_terminal <- ifelse(tip %in% names(M_trait_vector),M_trait_vector[tip],NA)
  if(!is.na(cluster)&&!is.na(M_value_terminal)){
    result_list[[i]] <- data.frame(
      Node = path[length(path)],   # Último nodo en la ruta, que es el terminal
      Cluster = cluster,
      M_value = M_value_terminal)
    for(node in head(path,-1)){
      M_value<-ifelse(as.character(node)%in% names(M_anc_states),
                      M_anc_states[as.character(node)],NA)
      if(!is.na(cluster)&&!is.na(M_value)) {
        result_list[[length(result_list)+1]]<-data.frame(
          Node =node,
          Cluster = cluster,
          M_value = M_value
        )
      }
    }  
  }
}



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


# Combinar todos los dataframes en uno solo
result_df_M <- do.call(rbind, result_list)

# Para nodos ancestrales, asegurarse de asignar clusters
result_df_M <- result_df_M[!duplicated(result_df_M$Node), ]  # Eliminar filas duplicadas de cada nodo

# Ver el dataframe resultante
row.names(result_df_M)<-NULL
print(result_df_M)
################################## G-M #########################################
#Dataframe de ambos
result_df_GM<-merge(result_df_G, 
                    result_df_M,
                    all.x = T,
                    by = c("Node","Cluster"))
result_df_GM
max_path_length<-max(sapply(paths,length))
max_path_length
extended_paths <- lapply(paths, function(path){
  if (length(path)<max_path_length){
    c(path,rep(tail(path,1),max_path_length - length(path)))
  }else{
    path
  }
})
extended_paths

calc_means<-function(nodos,df) {
  #Extraer valores de G y M
  G_values <- df$G_value[df$Node %in% nodos]
  M_values <- df$M_value[df$Node %in% nodos]
  #Calcular promedio- NA == 0
  G_mean<-median(G_values)
  M_mean<-median(M_values)
  G_mean<-ifelse(is.na(G_values),0,G_values)
  M_mean<-ifelse(is.na(M_values),0,M_values)
  return(data.frame(G_mean=G_mean,M_mean=M_mean))
}


levels_means<-list()
for(level in 1:max_path_length){
  nodos_at_level<-sapply(extended_paths, function(path)path[level])
  mean_2<-calc_means(nodos_at_level,result_df_GM)
  levels_means[[level]]<-data.frame(Level=level,G_mean=mean_2$G_mean,M_mean=mean_2$M_mean)
}
level_means_df<- do.call(rbind,levels_means)
print(level_means_df)

averages_total <- level_means_df %>%
  group_by(Level) %>%
  summarise(across(everything(),median))
averages_total

# Crear el gráfico
ggplot(averages_total, aes(x = Level)) + 
  geom_point(aes(y=G_mean), color="blue") + # Puntos para G
  geom_line(aes(y=G_mean), color="blue",linetype="dashed")+ #Línea para G
  geom_point(aes(y=M_mean), color="red") + # Puntos para M
  geom_line(aes(y=M_mean), color="red",linetype="dashed")+ #Línea para M
  scale_y_continuous(
    name = "G",
    sec.axis = sec_axis(~.*2,name = "M")
  ) +
  annotate("text",x=mean(averages_total$Level),y=mean(averages_total$G_mean),label="G trait",hjust=0.7,vjust=-0.6,size=4,color="blue")+
  annotate("text",x=mean(averages_total$Level),y=mean(averages_total$M_mean),label="M trait",hjust=0.25,vjust=-0.6,size=4,color="red")+
  labs(
    x="Number of clades",
    title = "The G and M traits for the tree of total communities",
  ) +
  theme_classic()

########################## Líneas por cluster ##################################
# Filtrar el nodo 16 (raiz)
node_171_G <- result_df_G[result_df_G$Node == 171, ]
node_171_M <- result_df_M[result_df_M$Node == 171, ]
# Crear tres copias cambiando el valor del Cluster
node_171_cluster_1G <- node_171_G
node_171_cluster_1G$Cluster <- 1

node_171_cluster_2G <- node_171_G
node_171_cluster_2G$Cluster <- 2


node_171_cluster_1M <- node_171_M
node_171_cluster_1M$Cluster <- 1

node_171_cluster_2M <- node_171_M
node_171_cluster_2M$Cluster <- 2


# Unir las filas modificadas al dataframe original
result_df_G <- rbind(result_df_G, node_171_cluster_1G, node_171_cluster_2G)
tail(result_df_G)
result_df_M <- rbind(result_df_M, node_171_cluster_1M,node_171_cluster_2M)
tail(result_df_M)
#Eliminar duplicados
result_df_G <- result_df_G %>% distinct()
result_df_M <- result_df_M %>% distinct()

#Dataframe de ambos
result_df_GM <- merge(result_df_G, 
                      result_df_M, 
                      all.x = TRUE, 
                      by = c("Node", "Cluster"))

result_df_GM
max_path_length<-max(sapply(paths,length))
max_path_length
extended_paths <- lapply(paths, function(path){
  if (length(path)<max_path_length){
    c(path,rep(tail(path,1),max_path_length - length(path)))
  }else{
    path
  }
})
extended_paths

calc_means <- function(nodos, df) {
  # Filtrar los datos para los nodos en el nivel actual
  df_filtered <- df %>% filter(Node %in% nodos)
  
  # Calcular promedios por cluster
  G_means <- df_filtered %>%
    group_by(Cluster) %>%
    summarise(G_mean = median(G_value, na.rm = TRUE), .groups = 'drop')
  
  M_means <- df_filtered %>%
    group_by(Cluster) %>%
    summarise(M_mean = median(M_value, na.rm = TRUE), .groups = 'drop')
  
  # Combinar promedios de G y M para cada cluster
  cluster_means <- merge(G_means, M_means, by = "Cluster", all = TRUE)
  
  # Llenar NAs con 0
  cluster_means[is.na(cluster_means)] <- 0
  
  return(cluster_means)
}

levels_means <- list()
for (level in 1:max_path_length) {
  # Obtener nodos en el nivel actual
  nodos_at_level <- sapply(extended_paths, function(path) path[level])
  
  # Calcular los promedios para cada cluster en este nivel
  mean_data <- calc_means(nodos_at_level, result_df_GM)
  
  # Agregar la columna de nivel
  mean_data$Level <- level
  
  # Almacenar resultados en la lista
  levels_means[[level]] <- mean_data
}

# Convertir lista en dataframe
level_means_df <- do.call(rbind, levels_means)
print(level_means_df)

# Calcular promedio general por nivel
averages_t_cluster <- level_means_df %>%
  group_by(Level, Cluster) %>%
  summarise(across(starts_with("G_mean"), median, na.rm = TRUE),
            across(starts_with("M_mean"), median, na.rm = TRUE),
            .groups = 'drop')
averages_t_cluster

# Crear el gráfico con ggplot2
ggplot(averages_t_cluster) +
  # Línea continua para G_mean
  geom_line(aes(x = Level, y = G_mean, color = factor(Cluster), linetype = "G_mean", group = Cluster), 
            size = 1) +
  # Línea punteada para M_mean
  geom_line(aes(x = Level, y = M_mean, color = factor(Cluster), linetype = "M_mean", group = Cluster), 
            size = 1) +
  # Configuración de colores
  scale_color_manual(values = c("blue", "red")) +
  # Configuración de tipos de línea
  scale_linetype_manual(values = c("G_mean" = "solid", "M_mean" = "dotted")) +
  labs(
    title = "The G and M traits of the tree of sampled communities",
    x = "Number of clades",
    y = "G mean",
    color = "Cluster",
    linetype = "Trait"
  ) +
  # Agregar eje secundario para M
  scale_y_continuous(
    sec.axis = sec_axis(~ ., name = "M mean")
  ) +
  theme_classic()
