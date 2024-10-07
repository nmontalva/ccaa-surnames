##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 4 ####
### To compare trees builds with different data and their deviations from the consensus tree ###

## Cargar paquetes y librerias ##
library(dendextend)
library(dplyr)
library(geosphere)
library(vegan)

### Cargar DATOS geograficos ###
shape.data <- sf::st_read("scripts_fv/Datos/comunidades_reprojected.shx")
shape.data$Name <-  gsub("\\s*\\(\\d+\\)", "", shape.data$Name)
shape.data$Name <- toupper(shape.data$Name)
shape.data$Name <- gsub(" ", "_", shape.data$Name)
shape.data$Name[grepl("LA_RINCONADA_DE_PUNITAQUI" , shape.data$Name)] <- "RINCONADA_DE_PUNITAQUI"
rownames(shape.data) <- shape.data$Name
my_points_v<- dplyr::select(shape.data,Name,geometry)
rownames(my_points_v)<-shape.data$Name

shape.data2<-shape.data
shape.data2$Name[shape.data2$Name %in% c("GUALLIGUAICA")] <- "PUCLARO"
rownames(shape.data2) <- shape.data2$Name
my_points_w<- dplyr::select(shape.data2,Name,geometry)
rownames(my_points_w)<-shape.data2$Name


##Test de Mantel##
#3. Comunidades muestreadas: ?rbol de consenso #
#Crear una matriz de distancia con datos de coordenadas de comunidades seleccionadas
my_points_v <- my_points_v %>% filter(row.names(shape.data) %in% selected_communities)
coords <- sf::st_coordinates(my_points_v)# Estraer coordenadas
geo_muestra <- sf::st_distance(my_points_v)# Calcular matriz de distancias
rownames(geo_muestra) <- as.factor(rownames(my_points_v))
colnames(geo_muestra)<- as.factor(rownames(my_points_v))
geo_muestra

my_points_w <- my_points_w %>% filter(row.names(shape.data2) %in% selected_communities2)
coords <- sf::st_coordinates(my_points_w)# Estraer coordenadas
geo_muestra2 <- sf::st_distance(my_points_w)# Calcular matriz de distancias
rownames(geo_muestra2) <- as.factor(rownames(my_points_w))
colnames(geo_muestra2)<- as.factor(rownames(my_points_w))
geo_muestra2

#Crar una matriz de distancia con arbol de consenso
con1 <-as.matrix(cophenetic.phylo(consensus_tree1)) #hy primero, s/PUCLARO
con2 <-as.matrix(cophenetic.phylo(consensus_tree2)) #PhyDSW primero, s/PUCLARO
con3 <-as.matrix(cophenetic.phylo(consensus_tree3)) #hy primero, c/PUCLARO
con4 <-as.matrix(cophenetic.phylo(consensus_tree4)) #PhyDSW primero, c/PUCLARO

mantel_function <- function(t1,t2) {
  #Conversir en matrix
  t1 <- as.matrix(t1)
  t2 <- as.matrix(t2)
  # Encontrar los row.names en com?n
  common_rows <- intersect(row.names(t1), row.names(t2))
  # Filtrar las matrices para que sol# Filtrar las matrices para que sol# Filtrar las matrices para que solo contengan las filas y columnas con row.names en com?n
  mat1_filtered <- t1[common_rows, common_rows]
  mat2_filtered <- t2[common_rows, common_rows]
  # Identificar las filas/columnas que quedaron fuera en mat1
  rows_outside_mat1 <- setdiff(row.names(t1), common_rows)
  mat1_outside <- t1[rows_outside_mat1, rows_outside_mat1]
  # Identificar las filas/columnas que quedaron fuera en mat2
  rows_outside_mat2 <- setdiff(row.names(t2), common_rows)
  mat2_outside <- t2[rows_outside_mat2, rows_outside_mat2]
  # Mostrar las filas/columnas que quedaron fuera
  if(length(rows_outside_mat1) > 0){
    print("Filas/columnas que quedaron fuera de la matriz gen?tica:")
    print(t1[rows_outside_mat1, rows_outside_mat1])
  }
  
  if(length(rows_outside_mat2) > 0){
    print("Filas/columnas que quedaron fuera de la matriz geogr?fica:")
    print(t2[rows_outside_mat2, rows_outside_mat2])
  }
  #Distancias
  mat2_filtered<-as.dist(mat2_filtered)
  mat1_filtered <- as.dist(mat1_filtered)
  #Generar mantel
  mantel.test <- mantel(
    mat2_filtered, # La primera de las dos matrices de disimilitud
    mat1_filtered, # La segunda matriz
    method = "pearson", # M?todo, el m?s com?n es pearson puede ser ( "pearson", "spearman", o "kendall")
    permutations = 999, # N?mero de permutaciones, usar how() para determinarlas
    strata = NULL, # Valor integrer o factor que identifique estratificaci?n para las permutaciones
    na.rm = FALSE, # Remover NS 
    parallel = getOption("mc.cores") # Opciones para procesos paralelos
  )
  
  #Visualizar mantel
  return(print(mantel.test))
}

mantel_con1 <- mantel_function(con1, geo_muestra)
mantel_con2 <- mantel_function(con2, geo_muestra)
mantel_con3 <- mantel_function(con3, geo_muestra2)
mantel_con4 <- mantel_function(con4,geo_muestra2)

mantel_con1_ap <- mantel_function(con1, surname_matrix_muestra)
mantel_con2_ap <- mantel_function(con2, surname_matrix_muestra)
mantel_con3_ap <- mantel_function(con3, surname_matrix_muestra2)
mantel_con4_ap <- mantel_function(con4, surname_matrix_muestra2)

mantel_con1_str <- mantel_function(con1, as.matrix(DSW))
mantel_con2_str <- mantel_function(con2, as.matrix(DSW))
mantel_con3_str <- mantel_function(con3, as.matrix(DSW2))
mantel_con4_str <- mantel_function(con4, as.matrix(DSW2))
