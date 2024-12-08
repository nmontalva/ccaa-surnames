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
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T)
coordenadas$ï..community <- gsub(" ", "_", coordenadas$ï..community)
coordenadas$ï..community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$ï..community)] <- "RINCONADA_DE_PUNITAQUI"

##Administraci?n de datos ##
rownames(coordenadas) <- coordenadas$ï..community
colnames(coordenadas)<- c("community","lon","lat","org_name")
my_points_t <- dplyr::select(coordenadas, lon, lat)
rownames(my_points_t) <- coordenadas$community
common_communities <- unique(comuneros$community)

##Test de Mantel##
#3. Comunidades muestreadas: ?rbol de consenso #
#Crear una matriz de distancia con datos de coordenadas de comunidades seleccionadas
selected_communities <- unique(STR$pop)
my_points_t <- my_points_t %>% filter(row.names(my_points_t) %in% selected_communities)
geo_muestra <- distm (my_points_t, fun = distGeo )
rownames(geo_muestra) <- as.factor(rownames(my_points_t))
colnames(geo_muestra) <- as.factor(rownames(my_points_t))
geo_muestra <- as.matrix(geo_muestra)
geo_muestra

# como arbol
Geo_tree <- upgma(as.dist(geo_muestra),method="average")
plot.phylo(Geo_tree)

##Matriz deGeo_tree##Matriz de distancia con datos de apellidos
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
 
  #Distancias
  mat2_filtered<-as.dist(mat2_filtered)
  mat1_filtered <- as.dist(mat1_filtered)
  set.seed(152)
  #Generar mantel
  mantel.test <- mantel.rtest(
    mat2_filtered, # La primera de las dos matrices de disimilitud
    mat1_filtered, # La segunda matriz
    nrepet = 1000000
  )
  
  #Visualizar mantel
  return(print(mantel.test))
}

##Resultados de mantel tests con geo_muestra
mantel_DPS<- mantel_function(DPS,geo_muestra)

##Resultados de mantel tests con apellidos
mantel_DPS_Ap<- mantel_function(DPS,surname_matrix_muestra)

