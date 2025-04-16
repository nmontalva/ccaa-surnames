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
library(phangorn)

### Cargar DATOS geogr?ficos ###
<<<<<<< Updated upstream
<<<<<<< Updated upstream
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T)
<<<<<<< HEAD
coordenadas$ï..community <- gsub(" ", "_", coordenadas$ï..community)
coordenadas$ï..community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$ï..community)] <- "RINCONADA_DE_PUNITAQUI"
=======
=======
>>>>>>> Stashed changes
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T, fileEncoding = "UTF-8-BOM")
=======
>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22
coordenadas$community <- gsub(" ", "_", coordenadas$community)
coordenadas$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$community)] <- "RINCONADA_DE_PUNITAQUI"

#=======
#TODO: REVISAR. Al igual que en el objetivo_1.2. las lineas 15,16 y 17 arrojan error por un problema de codificación de caracteres.
# Arreglé manualmente, pero seguro que se va a revertir cuando se abra desde el equipo con el problema.
# Vamos a tener que resolverlo, o seguirá pasando.
<<<<<<< HEAD
#=======
## REVISION: Agregué la misma solución que en 1.2
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
# 
#=======
>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22

##Administraci?n de datos ##
rownames(coordenadas) <- coordenadas$community
colnames(coordenadas)<- c("community","lon","lat","org_name")
my_points_t <- dplyr::select(coordenadas, lon, lat)
rownames(my_points_t) <- coordenadas$community
common_communities <- unique(comuneros$community)

##Test de Mantel##
#3. Comunidades muestreadas: ?rbol de consenso #
#Crear una matriz de distancia con datos de coordenadas de comunidades seleccionadas
selected_communities <- unique(STR$pop)
<<<<<<< HEAD
<<<<<<< Updated upstream
my_points_t <- my_points_t %>% filter(row.names(my_points_t) %in% selected_communities)
=======
my_points_t <- my_points_t %>% dplyr::filter(row.names(my_points_t) %in% selected_communities)

<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
my_points_t <- my_points_t %>% dplyr::filter(row.names(my_points_t) %in% selected_communities)

#=======
#TODO: REVISAR. El paquete conflicted nos obliga a definir cuál proveedor de la función "filter estamos usando.
# Asumí que es la de dplyr:
# my_points_t <- my_points_t %>% dplyr::filter(row.names(my_points_t) %in% selected_communities)
# Hay que revisar si es correcto.
# 
#=======

>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22
geo_muestra <- distm (my_points_t, fun = distGeo )

#=======
#TODO: REVISAR. Aquí aparecío el objeto "geo_muestra".
# ¿Es ese el qie faltaba para correr O1.2?
# De ser así, habría que incluir este código allá en vez de acá.
# O bien, correr esto primero y O1.2. después 
# (e inidicar en alguna parte que deben correrse en ese orden)
# 
#=======

rownames(geo_muestra) <- as.factor(rownames(my_points_t))
colnames(geo_muestra) <- as.factor(rownames(my_points_t))
geo_muestra <- as.matrix(geo_muestra)
geo_muestra

# como arbol
<<<<<<< Updated upstream
<<<<<<< Updated upstream
Geo_tree <- upgma(as.dist(geo_muestra),method="average")
<<<<<<< HEAD
=======
Geo_tree <- phangorn::upgma(as.dist(geo_muestra),method="average")
>>>>>>> Stashed changes
=======
Geo_tree <- phangorn::upgma(as.dist(geo_muestra),method="average")
>>>>>>> Stashed changes
=======

#=======
#TODO: REVISAR. Ninguna de las librerías cargadas hasta ahora provee la función UPGMA.
# Podría asumir que es upgma de Phangorn,  upgma de numbat, u otro.
# Asumí phangorn:
# library(phangorn)
# Revisar si es correcto, y si lo es, añadir al preambulo.
# 
#=======

>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22
plot.phylo(Geo_tree)

##Matriz de Geo_tree##Matriz de distancia con datos de apellidos
surname_matrix_muestra
<<<<<<< HEAD
<<<<<<< Updated upstream
<<<<<<< Updated upstream
surname_matrix_muestra2
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
#surname_matrix_muestra2

#=======
#TODO: REVISAR. El objeto surname_matrix_muestra2 no existe.
# ¿Tal vez se creaba en "Objetivo_1.2.Manteltest.R" que no corrió completamente?
#  
#=======
>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22

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

