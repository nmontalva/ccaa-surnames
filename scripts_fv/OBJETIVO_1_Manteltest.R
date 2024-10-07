##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### ESPACIO DE TRABAJO ####
#getwd()
#setwd("C:/Users/Kibif/Desktop/Proyecto desigualdad agropastores/Directorio_proyecto")

#### OBJETIVO 1 ####
###  To build a phylogenetic tree showing relationships between communities based on the distributions of surnames within and between communities. ###


## Cargar paquetes y librerias ##
library(dendextend)
library(dplyr)
library(geosphere)
library(vegan)

## Cargar DATOS ##
shape.data <- sf::st_read("scripts_fv/Datos/comunidades_reprojected.shx")
shape.data$Name <-  gsub("\\s*\\(\\d+\\)", "", shape.data$Name)
shape.data$Name <- toupper(shape.data$Name)
shape.data$Name <- gsub(" ", "_", shape.data$Name)
shape.data$Name[grepl("LA_RINCONADA_DE_PUNITAQUI" , shape.data$Name)] <- "RINCONADA_DE_PUNITAQUI"


##Administraci?n de datos ##
rownames(shape.data) <- shape.data$Name
my_points_v<- dplyr::select(shape.data,Name,geometry)
rownames(my_points_v)<-shape.data$Name

##Test de Mantel##
#1. Todas las comunidades # #Revisar la matriz de distancia de apellidos
# Crear una matriz de distancia con datos de coordenadas 
coords <- sf::st_coordinates(my_points_v)# Estraer coordenadas
geo_total <- sf::st_distance(my_points_v)# Calcular matriz de distancias
rownames(geo_total) <- as.factor(rownames(my_points_v))
colnames(geo_total)<- as.factor(rownames(my_points_v))

#Matriz de distancia con datos de apellidos (si se corri? el script "OBJETIVO_1_Surnames.R" es la matriz surname_matrix)
surname_matrix <- as.matrix(surname_matrix)

#Intersecci?n entre ambas matrices
# Encontrar los row.names en com?n
common_rows <- intersect(row.names(surname_matrix), row.names(geo_total))
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
mantel(
    mat2_filtered, # La primera de las dos matrices de disimilitud
    mat1_filtered, # La segunda matriz
    method = "pearson", # M?todo, el m?s com?n es pearson puede ser ( "pearson", "spearman", o "kendall")
    permutations = 999, # N?mero de permutaciones, usar how() para determinarlas
    strata = NULL, # Valor integrer o factor que identifique estratificaci?n para las permutaciones
    na.rm = FALSE, # Remover NS 
    parallel = getOption("mc.cores") # Opciones para procesos paralelos
  )

#2. Comunidades muestreadas #
#Administración de datos
shape.data$Name[shape.data$Name %in% c("GUALLIGUAICA")] <- "PUCLARO"
rownames(shape.data) <- shape.data$Name
my_points_v<- dplyr::select(shape.data,Name,geometry)
rownames(my_points_v)<-shape.data$Name

#Crear una matriz de distancia con datos de coordenadas de comunidades muestreadas
selected_communities <- unique(STR$pop)
selected_communities
my_points_v <- my_points_v %>% filter(row.names(shape.data) %in% selected_communities)
coords <- sf::st_coordinates(my_points_v)# Estraer coordenadas
geo_muestra <- sf::st_distance(my_points_v)# Calcular matriz de distancias
rownames(geo_muestra) <- as.factor(rownames(my_points_v))
colnames(geo_muestra)<- as.factor(rownames(my_points_v))
View(geo_muestra)
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

#Generar mantel
mantel(
  mat2_filtered, # La primera de las dos matrices de disimilitud
  mat1_filtered, # La segunda matriz
  method = "pearson", # M?todo, el m?s com?n es pearson puede ser ( "pearson", "spearman", o "kendall")
  permutations = 999, # N?mero de permutaciones, usar how() para determinarlas
  strata = NULL, # Valor integrer o factor que identifique estratificaci?n para las permutaciones
  na.rm = FALSE, # Remover NS 
  parallel = getOption("mc.cores") # Opciones para procesos paralelos
)
