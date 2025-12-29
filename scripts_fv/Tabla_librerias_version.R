info_sesion <- sessionInfo()
char_packages <- names(info_sesion$otherPkgs)

# Función para extraer la cita en texto plano
get_cite <- function(pkg){
  cit <- citation(pkg)
  # Intentamos sacar la primera cita (a veces hay varias)
  if(length(cit) > 0){
    format(cit[[1]], style = "text") 
  } else {
    "R Core Team"
  }
}

# Tabla con versión Y cita
tabla_completa <- data.frame(
  Library = char_packages,
  Version = sapply(char_packages, function(x) as.character(packageVersion(x))),
  Citation = sapply(char_packages, get_cite), # Agrega la cita real
  row.names = NULL
)

#alphabetic format
tabla_completa <- tabla_completa[order(tabla_completa$Library), ]