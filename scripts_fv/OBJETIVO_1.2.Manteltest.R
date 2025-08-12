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

### Cargar DATOS geograficos ###

coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T, fileEncoding = "UTF-8-BOM")
coordenadas$community <- gsub(" ", "_", coordenadas$community)
coordenadas$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$community)] <- "RINCONADA_DE_PUNITAQUI"

##Test de Mantel##
#1. Todas las comunidades # #Revisar la matriz de distancia de apellidos
rownames(coordenadas) <- coordenadas$community
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

# Filtrar las matrices para que solo contengan las filas y columnas con row.names en com?n
surname_filtered <- surname_matrix[common_rows, common_rows]
geo_filtered <- geo_total[common_rows, common_rows]


# Identificar las filas/columnas que quedaron fuera en mat1
rows_outside_mat1 <- setdiff(row.names(surname_matrix), common_rows)
surname_outside <- surname_matrix[rows_outside_mat1, rows_outside_mat1]

# Identificar las filas/columnas que quedaron fuera en mat2
rows_outside_mat2 <- setdiff(row.names(geo_total), common_rows)
geo_outside <- geo_total[rows_outside_mat2, rows_outside_mat2]

# Mostrar las filas/columnas que quedaron fuera
print("Filas/columnas que quedaron fuera de la matriz de apellidos:")
print(surname_outside) ## Quedan fuera Huascoaltinos y totoral
print("Filas/columnas que quedaron fuera de la matriz geográfica:")
print(geo_outside)

#Matrices de Distancias
surname_filtered<-as.dist(surname_filtered)
geo_filtered<-as.dist(geo_filtered)

#Generar mantel
mantel.rtest(
  geo_filtered, # La primera de las dos matrices de disimilitud
  surname_filtered, # La segunda matriz
    nrepet = iter
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
surname_muestra <- surname_matrix_muestra[common_rows, common_rows]
geo_muestra <- geo_muestra[common_rows, common_rows]
# Identificar las filas/columnas que quedaron fuera en mat1
rows_outside_mat1 <- setdiff(row.names(surname_matrix_muestra), common_rows)
surname_muestra_outside <- surname_matrix_muestra[rows_outside_mat1, rows_outside_mat1]

# Identificar las filas/columnas que quedaron fuera en mat2
rows_outside_mat2 <- setdiff(row.names(geo_muestra), common_rows)
geo_muestra_outside <- geo_muestra[rows_outside_mat2, rows_outside_mat2]

# Mostrar las filas/columnas que quedaron fuera #Esto deber?a generar una matriz de extensi?n 0x0
print("Filas/columnas que quedaron fuera de mat1:")
print(surname_muestra_outside)
print("Filas/columnas que quedaron fuera de mat2:")
print(geo_muestra_outside)

#configurar como matrices de Distancias
mat2_filtered<-as.dist(mat2_filtered)
mat1_filtered <- as.dist(mat1_filtered)
surname_muestra<-as.dist(surname_muestra)
geo_muestra <- as.dist(geo_muestra)
#Generar mantel
mantel.rtest(
  geo_muestra, # La primera de las dos matrices de disimilitud
  surname_muestra, # La segunda matriz
  nrepet = iter
)

