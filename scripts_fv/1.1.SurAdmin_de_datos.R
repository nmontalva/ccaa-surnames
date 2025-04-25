######################Limpieza de la base de datos########################
##########################################################################

#Instalación de paquetes
library(dplyr)
library(ade4)
library(tidyr)


#Abrir y explorar bases de datos
primera<- read.csv("scripts_fv/Datos/all_tabulated.csv", sep=",", header=TRUE, fill = T)
segunda<- read.csv("scripts_fv/Datos/Muestras Fondecyt 11160402 Alelos - RACK 1.csv", sep =",", header = TRUE, fill = TRUE)

# Renombrar la primera columna y crear variable de alelos
primera <- primera %>%
  rename(Sample.Name = X) %>% 
  mutate(Alleles = rep(1:2, length.out = nrow(primera)))

# Eliminar columna X.1 y reemplazar -9 y datos vacíos por NA
primera <- primera %>%
  dplyr::select(-X.1) %>%
  mutate(across(c(D3S1358, TH01, D21S11, D18S51, Penta.E, D5S818, D13S317, D7S820, D16S539, CSF1PO, Penta.D, vWA, D8S1179, TPOX, FGA), ~ na_if(., -9)))

segunda <- segunda %>%
  mutate(across(c(Allele.1, Allele.2), ~ na_if(., "")))

################Unificar las bases de datos######################
#Transformar columnas en filas

X <- segunda %>%
  pivot_wider(names_from = Marker, values_from = c(Allele.1, Allele.2))

x1 <- X[1:17] %>% mutate(Alleles = 1)
x2 <- X[18:33]
x3 <- X[1]
x4 <- cbind(x2, x3) %>% mutate(Alleles = 2)

# Renombrar columnas de alelos de manera eficiente
new_names <- c("D3S1358", "TH01", "D21S11", "D18S51", "Penta.E", "D5S818", "D13S317", "D7S820", 
               "D16S539", "CSF1PO", "Penta.D", "AMEL", "vWA", "D8S1179", "TPOX", "FGA")

colnames(x1)[2:17] <- new_names
colnames(x4)[1:16] <- new_names

# Unir las bases de datos x1 y x4
union_df <- bind_rows(x1, x4)
colnames(union_df)
# Valores raros (OL)
# Cambiar "OL" a NA en todas las columnas de marcadores
union_df <- union_df %>%
  mutate(across(c("D3S1358", "TH01", "D21S11", "D18S51", "Penta.E", 
                  "D5S818", "D13S317", "D7S820", "D16S539", "CSF1PO", 
                  "Penta.D", "vWA", "D8S1179", "TPOX", "FGA"), ~ ifelse(. == "OL", NA, .)))
# Cambiar tipos de datos a factor y numérico donde corresponda
union_df <- union_df %>%
  mutate(across(c("D3S1358", "TH01", "D21S11", "D18S51", "Penta.E", 
                  "D5S818", "D13S317", "D7S820", "D16S539", "CSF1PO", 
                  "Penta.D", "vWA", "D8S1179", "TPOX", "FGA"), ~ as.numeric(as.character(.))))

#TODO: REVISAR. Me salen 2 warnings al correr la línea 50:
## RESOLUCION: Para evitar esto, mejor generar los NA que corresponan antes. CORREGIDO

# Fusionar ambas bases de datos
mi.final <- merge(union_df, primera, all = TRUE)

# Eliminar columna AMEL
mi.final <- mi.final %>%
  dplyr::select(-AMEL)

#Borrar otras bases 
primera <- NULL
segunda <- NULL 
X <- NULL 
x1 <- NULL 
x2 <- NULL
x3 <- NULL
x4 <- NULL
union_df <- NULL

rm(primera, segunda, union_df, X, x1, x2, x3, x4)
# Convertir columnas a clase numérica, excepto "Sample.Name" y "Alleles"
mi.final[ , 2:16] <- lapply(mi.final[ , 2:16], function(x) as.numeric(as.character(x)))

# Agregar variable 'community' basado en los patrones de 'Sample.Name'
mi.final <- mi.final %>%
  mutate(community = case_when(
    grepl("^DV", Sample.Name) ~ "EL_Divisadero",
    grepl("^BZ", Sample.Name) ~ "BARRAZA",
    grepl("^EA", Sample.Name) ~ "EL_ALTAR",
    grepl("^MQ", Sample.Name) ~ "MANQUEHUA",
    grepl("^PQ", Sample.Name) ~ "PUNITAQUI",
    grepl("^RP", Sample.Name) ~ "RINCONADA_DE_PUNITAQUI",
    grepl("^CA", Sample.Name) ~ "CANELILLA_OVALLE",
    grepl("^CB", Sample.Name) ~ "CANELA_BAJA",
    grepl("^(CL|CR)", Sample.Name) ~ "LA_CALERA",
    grepl("^ES", Sample.Name) ~ "EL_ESPINAL",
    grepl("^GG", Sample.Name) ~ "GUALLIGUAICA",
    grepl("^HT", Sample.Name) ~ "HUENTELAUQUEN",
    grepl("^LG", Sample.Name) ~ "CASTILLO_MAL_PASO",
    grepl("^MP", Sample.Name) ~ "MONTE_PATRIA",
    grepl("^VP", Sample.Name) ~ "LA_POLVADA",
    TRUE ~ "OTHER"
  ))

# Filtrar los casos donde 'community' no es 'OTHER' ni NA
mi.final <- mi.final %>%
  dplyr::filter(!community %in% c('OTHER', NA))

# Renombrar 'Sample.Name' a 'Ind' y 'community' a 'Pop'
mi.final <- mi.final %>%
  rename(Ind = Sample.Name, Pop = community)
mi.final<-data.frame(mi.final)
# Verificar el resultado final
View(mi.final)

# Transformación de datos
mi.final_wide <- mi.final %>%
  mutate(row_num = Alleles) %>%
  pivot_wider(
    names_from = row_num,
    values_from = -c(Ind, Pop),
    names_sep = "_"
  )
colnames(mi.final_wide)

#Quitar las últimas 4 columnas
STRv2 <-data.frame(mi.final_wide%>% dplyr::select(-tail(names(.), 4)))
colnames(STRv2)
STRv2
#Guardar en STRv.2
names(STRv2) <- gsub("_1|_2", "", names(STRv2))
write.csv(STRv2, "scripts_fv/Datos/STRV2.csv", row.names = FALSE)
