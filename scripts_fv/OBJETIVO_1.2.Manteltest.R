##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 1 ####
###  To build a phylogenetic tree showing relationships between communities based on the distributions of surnames within and between communities. ###

## Cargar paquetes y librerias ##
library(dendextend)
library(dplyr)
library(geosphere)
library(vegan)

### Cargar DATOS geogr?ficos ###
<<<<<<< Updated upstream
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T)
coordenadas$ï..community <- gsub(" ", "_", coordenadas$ï..community)
coordenadas$ï..community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$ï..community)] <- "RINCONADA_DE_PUNITAQUI"

=======
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T, fileEncoding = "UTF-8-BOM")
coordenadas$community <- gsub(" ", "_", coordenadas$community)
coordenadas$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$community)] <- "RINCONADA_DE_PUNITAQUI"

#=======
#TODO: REVISAR. Las lineas 15,16 y 17 arrojan error por un problema de codificación de caracteres.
# Arreglé manualmente, pero seguro que se va a revertir cuando se abra desde el equipo con el problema.
# Vamos a tener que resolverlo, o seguirá pasando.
# 
#=======
## REVISION: Encontré el error! Era un problema de marca de orden de bytes (BOM) en sistemas Windows. Lo agregué a la lectura del archivo csv, no debería seguir ocurriendo.
>>>>>>> Stashed changes
##Test de Mantel##
#1. Todas las comunidades # #Revisar la matriz de distancia de apellidos
rownames(coordenadas) <- coordenadas$ï..community
colnames(coordenadas)<- c("community","lon","lat","org_name")
my_points_t <- dplyr::select(coordenadas, lon, lat)
rownames(my_points_t) <- coordenadas$community
common_communities <- unique(comuneros$community)
geo_total <- distm (my_points_t, fun = distGeo )
rownames(geo_total) <- as.factor(rownames(my_points_t))
colnames(geo_total) <- as.factor(rownames(my_points_t))
geo_total <- as.matrix(geo_total)
geo_total
#Matriz de distancia con datos de apellidos (si se corrio el script "OBJETIVO_1_Surnames.R" es la matriz surname_matrix)
surname_matrix <- as.matrix(surname_matrix)

#Intersecci?n entre ambas matrices
# Encontrar los row.names en com?n
common_rows <- intersect(row.names(surname_matrix), row.names(geo_total))
<<<<<<< Updated upstream
=======
#=======
#TODO: REVISAR. No existe el objeto 'geo_total'.
# Error: object 'geo_total' not found
# Como esto me impide seguir revisando, voy a correr la siguiente linea, que dejaré comentada:
# common_rows <- row.names(surname_matrix)
# Obviamente esto es "trampa". Hay que volver a revisar todo después.
#=======
## REVISIÓN: Agregué la creación de Geo_total. Este error se dió porque creé primero la versión 3.2 de Manteltest y no me dí cuenta que no estaba ese archivo creado en este script.
>>>>>>> Stashed changes
common_rows

# Filtrar las matrices para que solo contengan las filas y columnas con row.names en com?n
mat1_filtered <- surname_matrix[common_rows, common_rows]
mat2_filtered <- geo_total[common_rows, common_rows]

# Identificar las filas/columnas que quedaron fuera en mat1
rows_outside_mat1 <- setdiff(row.names(surname_matrix), common_rows)
mat1_outside <- surname_matrix[rows_outside_mat1, rows_outside_mat1]

# Identificar las filas/columnas que quedaron fuera en mat2
rows_outside_mat2 <- setdiff(row.names(geo_total), common_rows)
mat2_outside <- geo_total[rows_outside_mat2, rows_outside_mat2]

# Mostrar las filas/columnas que quedaron fuera
print("Filas/columnas que quedaron fuera de la matriz de apellidos:")
print(mat1_outside) ## Quedan fuera Huascoaltinos y totoral
print("Filas/columnas que quedaron fuera de la matriz geográfica:")
print(mat2_outside)

#Matrices de Distancias
mat2_filtered<-as.dist(mat2_filtered)
mat1_filtered <- as.dist(mat1_filtered)

#Generar mantel
set.seed(152)
mantel.rtest(
    mat2_filtered, # La primera de las dos matrices de disimilitud
    mat1_filtered, # La segunda matriz
    nrepet = 1000000
  )

#2. Comunidades muestreadas #
#Administración de datos
selected_communities <- unique(STR$pop)
my_points_t <- my_points_t %>% dplyr::filter(row.names(my_points_t) %in% selected_communities)
geo_muestra <- distm (my_points_t, fun = distGeo )
rownames(geo_muestra) <- as.factor(rownames(my_points_t))
colnames(geo_muestra) <- as.factor(rownames(my_points_t))
geo_muestra <- as.matrix(geo_muestra)
geo_muestra
#Crear una matriz de distancia con datos de apellidos de las comunidades muestreadas 
surname_matrix_muestra
surname_matrix_muestra <- as.matrix(surname_matrix_muestra)

#Intersecci?n entre ambas matrices
# Encontrar los row.names en com?n
common_rows <- intersect(row.names(surname_matrix_muestra), row.names(geo_muestra))

# Filtrar las matrices para que solo contengan las filas y columnas con row.names en com?n
mat1_filtered <- surname_matrix_muestra[common_rows, common_rows]
mat2_filtered <- geo_muestra[common_rows, common_rows]
# Identificar las filas/columnas que quedaron fuera en mat1
rows_outside_mat1 <- setdiff(row.names(surname_matrix_muestra), common_rows)
mat1_outside <- surname_matrix_muestra[rows_outside_mat1, rows_outside_mat1]

# Identificar las filas/columnas que quedaron fuera en mat2
rows_outside_mat2 <- setdiff(row.names(geo_muestra), common_rows)
mat2_outside <- geo_muestra[rows_outside_mat2, rows_outside_mat2]

# Mostrar las filas/columnas que quedaron fuera #Esto deber?a generar una matriz de extensi?n 0x0
print("Filas/columnas que quedaron fuera de mat1:")
print(mat1_outside)
print("Filas/columnas que quedaron fuera de mat2:")
print(mat2_outside)

#configurar como matrices de Distancias
mat2_filtered<-as.dist(mat2_filtered)
mat1_filtered <- as.dist(mat1_filtered)
set.seed(152)
#Generar mantel
mantel.rtest(
  mat2_filtered, # La primera de las dos matrices de disimilitud
  mat1_filtered, # La segunda matriz
  nrepet = 1000000
)