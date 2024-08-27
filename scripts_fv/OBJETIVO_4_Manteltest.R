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

## Cargar DATOS ##
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T)
coordenadas$community <- gsub(" ", "_", coordenadas$community)
coordenadas$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$community)] <- "RINCONADA_DE_PUNITAQUI"

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
my_points_t <- my_points_t %>% dplyr::filter(row.names(my_points_t) %in% selected_communities)
geo_muestra <- distm (my_points_t, fun = distGeo )
rownames(geo_muestra) <- as.factor(rownames(my_points_t))
colnames(geo_muestra) <- as.factor(rownames(my_points_t))
geo_muestra <- as.matrix(geo_muestra)
geo_muestra

#Crar una matriz de distancia con arbol de consenso
consenso_matrix <- as.matrix(cophenetic(consensus_tree))

#Intersecci?n entre ambas matrices
# Encontrar los row.names en com?n
common_rows <- intersect(row.names(consenso_matrix), row.names(geo_muestra))
# Filtrar las matrices para que solo contengan las filas y columnas con row.names en com?n
mat1_filtered <- consenso_matrix[common_rows, common_rows]
mat2_filtered <- geo_muestra[common_rows, common_rows]
# Identificar las filas/columnas que quedaron fuera en mat1
rows_outside_mat1 <- setdiff(row.names(consenso_matrix), common_rows)
mat1_outside <- consenso_matrix[rows_outside_mat1, rows_outside_mat1]

# Identificar las filas/columnas que quedaron fuera en mat2
rows_outside_mat2 <- setdiff(row.names(geo_muestra), common_rows)
mat2_outside <- geo_muestra[rows_outside_mat2, rows_outside_mat2]

# Mostrar las filas/columnas que quedaron fuera
print("Filas/columnas que quedaron fuera de mat1:")
print(mat1_outside)
print("Filas/columnas que quedaron fuera de mat2:")
print(mat2_outside)

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
# Explorar el objeto
mantel.test
