##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

### Cargar paquetes y librerias ###
library(dendextend)
library(dplyr)
library(geosphere)
library(vegan)

### Cargar DATOS geogr?ficos ###
shape.data <- sf::st_read("scripts_fv/Datos/comunidades_reprojected.shx")
shape.data$Name <-  gsub("\\s*\\(\\d+\\)", "", shape.data$Name)
shape.data$Name <- toupper(shape.data$Name)
shape.data$Name <- gsub(" ", "_", shape.data$Name)
shape.data$Name[grepl("LA_RINCONADA_DE_PUNITAQUI" , shape.data$Name)] <- "RINCONADA_DE_PUNITAQUI"
#shape.data$Name[shape.data$Name %in% c("GUALLIGUAICA")] <- "PUCLARO"
rownames(shape.data) <- shape.data$Name
my_points_v<- dplyr::select(shape.data,Name,geometry)
rownames(my_points_v)<-shape.data$Name

###Test de Mantel###
###3. Comunidades muestreadas: STRs ###
##Crear una matriz de distancia con datos de coordenadas de comunidades seleccionadas
my_points_v <- my_points_v %>% filter(row.names(shape.data) %in% selected_communities)
coords <- sf::st_coordinates(my_points_v)# Estraer coordenadas
geo_muestra <- sf::st_distance(my_points_v)# Calcular matriz de distancias
rownames(geo_muestra) <- as.factor(rownames(my_points_v))
colnames(geo_muestra)<- as.factor(rownames(my_points_v))
geo_muestra

##Matriz de distancia con datos de apellidos
surname_matrix_muestra
surname_matrix_muestra2

###Intersecci?n entre distintas matrices
##Mantel function
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


##Resultados de mantel tests con geo_muestra
mantel_GST <- mantel_function(GST, geo_muestra)
mantel_Nei <- mantel_function(Nei, geo_muestra)
mantel_CS <- mantel_function(cs, geo_muestra)
mantel_RST <- mantel_function(RST,geo_muestra)
mantel_RST2 <- mantel_function(RST2,geo_muestra)
mantel_ASD <- mantel_function(ASD, geo_muestra)
mantel_DSW <- mantel_function(DSW, geo_muestra)
mantel_Dmu2 <- mantel_function(Dmu2, geo_muestra)
mantel_FST <- mantel_function(FST, geo_muestra)

##Resultados de mantel tests con apellidos
mantel_GST_Ap <- mantel_function(GST,surname_matrix_muestra)
mantel_Nei_Ap <- mantel_function(Nei, surname_matrix_muestra)
mantel_CS_Ap <- mantel_function(cs, surname_matrix_muestra)
mantel_RST_Ap <- mantel_function(RST,surname_matrix_muestra)
mantel_RST2_Ap <- mantel_function(RST2,surname_matrix_muestra)
mantel_ASD_Ap <- mantel_function(ASD, surname_matrix_muestra)
mantel_DSW_Ap <- mantel_function(DSW, surname_matrix_muestra)
mantel_Dmu2_Ap <- mantel_function(Dmu2, surname_matrix_muestra)
mantel_FST_Ap <- mantel_function(FST, surname_matrix_muestra)

#Revisando matrices
Nei_cs<-mantel_function(Nei,cs)
RST_ASD<-mantel_function(RST,ASD)
RST_RST2<-mantel_function(RST,RST2) #Son muy parecidas pero no idénticas
GST_GST<-mantel_function(GST,GST2)
