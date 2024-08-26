##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### ESPACIO DE TRABAJO ####
setwd("C:/Users/Kibif/Desktop/Proyecto desigualdad agropastores/Directorio_proyecto")

#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

### Cargar paquetes y librerias ###
library(dendextend)
library(dplyr)
library(geosphere)
library(vegan)

### Cargar DATOS geográficos ###
coordenadas <- read.csv("Datos/coordenadas.csv", header = T)
coordenadas$ï..community <- gsub(" ", "_", coordenadas$ï..community)
coordenadas$ï..community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$ï..community)] <- "RINCONADA_DE_PUNITAQUI"

###Administración de datos geográficos ###
rownames(coordenadas) <- coordenadas$ï..community
colnames(coordenadas)<- c("community","lon","lat","org_name")
my_points_t <- dplyr::select(coordenadas, lon, lat)
rownames(my_points_t) <- coordenadas$community
common_communities <- unique(comuneros$community)

###Test de Mantel###
###3. Comunidades muestreadas: STRs ###
##Crear una matriz de distancia con datos de coordenadas de comunidades seleccionadas
selected_communities <- unique(STR$pop)
my_points_t <- my_points_t %>% filter(row.names(my_points_t) %in% selected_communities)
geo_muestra <- distm (my_points_t, fun = distGeo )
rownames(geo_muestra) <- as.factor(rownames(my_points_t))
colnames(geo_muestra) <- as.factor(rownames(my_points_t))

##Matriz de distancia con datos de apellidos
surname_distance_muestra

###Intersección entre distintas matrices
##Mantel function
mantel_function <- function(t1,t2) {
  #Conversir en matrix
  t1 <- as.matrix(t1)
  t2 <- as.matrix(t2)
  # Encontrar los row.names en común
  common_rows <- intersect(row.names(t1), row.names(t2))
  # Filtrar las matrices para que sol# Filtrar las matrices para que sol# Filtrar las matrices para que solo contengan las filas y columnas con row.names en común
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
    print("Filas/columnas que quedaron fuera de la matriz genética:")
    print(genetic_matrix[rows_outside_mat1, rows_outside_mat1])
  }
  
  if(length(rows_outside_mat2) > 0){
    print("Filas/columnas que quedaron fuera de la matriz geográfica:")
    print(geo_matrix[rows_outside_mat2, rows_outside_mat2])
  }
  #Distancias
  mat2_filtered<-as.dist(mat2_filtered)
  mat1_filtered <- as.dist(mat1_filtered)
  #Generar mantel
  mantel.test <- mantel(
    mat2_filtered, # La primera de las dos matrices de disimilitud
    mat1_filtered, # La segunda matriz
    method = "pearson", # Método, el más común es pearson puede ser ( "pearson", "spearman", o "kendall")
    permutations = 999, # Número de permutaciones, usar how() para determinarlas
    strata = NULL, # Valor integrer o factor que identifique estratificación para las permutaciones
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
mantel_ASD <- mantel_function(ASD, geo_muestra)

##Resultados de mantel tests con apellidos
mantel_GST_Ap <- mantel_function(GST,surname_matrix_muestra)
mantel_Nei_Ap <- mantel_function(Nei, surname_matrix_muestra)
mantel_CS_Ap <- mantel_function(cs, surname_matrix_muestra)
mantel_RST_Ap <- mantel_function(RST,surname_matrix_muestra)
mantel_ASD_Ap <- mantel_function(ASD, surname_matrix_muestra)

#Revisando matrices
Nei_Nei<-mantel_function(Nei,Nei2)
RST_ASD<-mantel_function(RST,ASD)
GST_GST<-mantel_function(GST,GST3)
FST_FST<-mantel_function(FST,FST)
