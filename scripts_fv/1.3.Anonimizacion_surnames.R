# anonimizacion de surnames
##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#Anonimizacion de apellidos
library(Biodem)
library(digest)

## Cargar DATOS ##
comuneros <- read.csv("scripts_fv/Datos/commoners.csv", fileEncoding = "UTF-8-BOM")
comuneros$community <- gsub(" ", "_", comuneros$community)
comuneros$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , comuneros$community)] <- "RINCONADA_DE_PUNITAQUI"


# Hacer copia del dataframe original
anon_comuneros <- comuneros

# Anonimizar columnas de apellidos (puedes agregar o quitar según necesites)
cols_apellidos <- c("firstname","firstname1","firstname2","surname", "surname_father", "surname_mother")

for (col in cols_apellidos) {
  anon_comuneros[[col]] <- sapply(
    tolower(trimws(anon_comuneros[[col]])),
    function(x) ifelse(is.na(x) || x == "", NA, substr(digest(x, algo = "sha256"), 1, 10))
  )
}

# Distancia con apellidos reales
surname_distance_matrix <- function(comuneros,
                                    group_by_col="community") {
  # cross tabulate
  surnames_freq <- table(comuneros$surname_father,
                         comuneros[[group_by_col]])
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # hedrick returns values of similarity
  # transform them into values of dissimilarity (distance)
  as.dist(1-hedkin)
  #as.dist(1-hedkin) #IMPORTANT: This is likely wrong. We figured out a better way.
  #as.dist(Biodem::Fst(hedkin, sum(colSums(surnames_freq) )))
}
orig <- surname_distance_matrix(comuneros)
# Hacer distancia anonimizada
anon_dist <- surname_distance_matrix(anon_comuneros)
# Comparar
cor(as.vector(orig), as.vector(anon_dist))  # debería dar 1

# Guardar en CSV (sin modificar el archivo original)
write.csv(anon_comuneros, "scripts_fv/Datos/comuneros_anon.csv", row.names = FALSE, fileEncoding = "UTF-8")

#Eliminar objetos
rm(anon,anon_comuneros,anon_dist,col,cols_apellidos,comuneros,orig,surname_distance_matrix)
