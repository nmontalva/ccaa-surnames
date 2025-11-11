# anonimizacion de surnames
##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#Anonimizacion de apellidos

library(digest)

## Cargar DATOS ##
comuneros <- read.csv("scripts_fv/Datos/commoners.csv", fileEncoding = "UTF-8-BOM")
comuneros$community <- gsub(" ", "_", comuneros$community)
comuneros$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , comuneros$community)] <- "RINCONADA_DE_PUNITAQUI"


# Hacer copia del dataframe original
anon_comuneros <- comuneros

# Anonimizar columnas de apellidos (puedes agregar o quitar según necesites)
cols_apellidos <- c("surname", "surname_father", "surname_mother")

for (col in cols_apellidos) {
  anon_comuneros[[col]] <- sapply(
    tolower(trimws(anon_comuneros[[col]])),
    function(x) ifelse(is.na(x) || x == "", NA, substr(digest(x, algo = "sha256"), 1, 10))
  )
}

# Guardar en CSV (sin modificar el archivo original)
write.csv(anon_comuneros, "scripts_fv/Datos/comuneros_anon.csv", row.names = FALSE, fileEncoding = "UTF-8")
