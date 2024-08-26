######################Limpieza de la base de datos########################
##########################################################################

#Instalación de paquetes

library(dplyr)
library(ade4)
library(tidyr)


#Abrir y explorar bases de datos

primera<- read.csv("scripts_fv/Datos/all_tabulated.csv", sep=",", header=TRUE, fill = T)
segunda<- read.csv("scripts_fv/Datos/Muestras Fondecyt 11160402 Alelos - RACK 1.csv", sep =",", header = TRUE, fill = TRUE)
data.frame(names(primera))
data.frame(names(segunda))

#Cambiar nombres de variable X

names(primera)[1]<-"Sample.Name"

#Crear variable que muestre los dos alelos.

primera$Alleles <- rep(1:2)

#eliminar variable X.1

primera$X.1 <- NULL


#Cambiar -9 y datos blancos a NA

for (i in -9) {
  primera$D3S1358[primera$D3S1358 == i] <- NA
  primera$TH01 [primera$TH01 == i] <- NA
  primera$D21S11 [primera$D21S11 == i] <- NA
  primera$D18S51 [primera$D18S51 == i] <- NA
  primera$Penta.E [primera$Penta.E == i] <- NA
  primera$D5S818 [primera$D5S818 == i] <- NA
  primera$D13S317 [primera$D13S317 == i] <- NA
  primera$D7S820 [primera$D7S820 == i] <- NA
  primera$D16S539 [primera$D16S539 == i] <- NA
  primera$CSF1PO [primera$CSF1PO == i] <- NA
  primera$Penta.D [primera$Penta.D == i] <- NA
  primera$AMEL [primera$AMEL == i] <- NA
  primera$vWA [primera$vWA == i] <- NA
  primera$D8S1179 [primera$D8S1179 == i] <- NA
  primera$TPOX [primera$TPOX == i] <- NA
  primera$FGA [primera$FGA == i] <- NA
  
}

for (i in "") {
  segunda$Allele.1 [segunda$Allele.1 == i] <- NA
  segunda$Allele.2 [segunda$Allele.2 == i] <- NA
  
}



################Unificar las bases de datos######################
#Transformar columnas en filas

segunda %>%
  pivot_wider(names_from= Marker, values_from= c(Allele.1,Allele.2)) %>%
  as.data.frame() -> X

#Separar los alelos
x1<-X[1:17]
x1$Alleles <- (1)

x2 <- X[18:33]
x3 <-X[1]
x4 <- cbind(x2,x3)
x4$Alleles <-(2)

# cambiar nombres de variables

for (i in 1) {
  
  names(x1)[i+1] <- "D3S1358"
  names(x1)[i+2] <- "TH01"
  names(x1)[i+3] <- "D21S11"
  names(x1)[i+4] <- "D18S51"
  names(x1)[i+5] <- "Penta.E"
  names(x1)[i+6] <- "D5S818"
  names(x1)[i+7] <- "D13S317"
  names(x1)[i+8] <- "D7S820"
  names(x1)[i+9] <- "D16S539"
  names(x1)[i+10] <- "CSF1PO"
  names(x1)[i+11] <- "Penta.D"
  names(x1)[i+12] <- "AMEL"
  names(x1)[i+13] <- "vWA"
  names(x1)[i+14] <- "D8S1179"
  names(x1)[i+15] <- "TPOX"
  names(x1)[i+16] <- "FGA"
  
  names(x4)[i] <- "D3S1358"
  names(x4)[i+1] <- "TH01"
  names(x4)[i+2] <- "D21S11"
  names(x4)[i+3] <- "D18S51"
  names(x4)[i+4] <- "Penta.E"
  names(x4)[i+5] <- "D5S818"
  names(x4)[i+6] <- "D13S317"
  names(x4)[i+7] <- "D7S820"
  names(x4)[i+8] <- "D16S539"
  names(x4)[i+9] <- "CSF1PO"
  names(x4)[i+10] <- "Penta.D"
  names(x4)[i+11] <- "AMEL"
  names(x4)[i+12] <- "vWA"
  names(x4)[i+13] <- "D8S1179"
  names(x4)[i+14] <- "TPOX"
  names(x4)[i+15] <- "FGA"
  
}


#Unir bases de datos en una sola

union<- union_all(x1,x4)

#Cambiar variables a integrales

primera %>%
  str()

union %>%
  str()

as.factor(union$Sample.Name)
as.integer(union$D3S1358)
as.numeric(union$TH01)
as.numeric(union$D21S11)
as.integer(union$D18S51)
as.integer(union$`Penta.E`)
as.integer(union$D5S818)
as.integer(union$D13S317)
as.integer(union$D7S820)
as.integer(union$D16S539)
as.integer(union$CSF1PO)
as.numeric(union$`Penta.D`)
as.factor(union$AMEL)
as.integer(union$vWA)
as.integer(union$D8S1179)
as.integer(union$TPOX)
as.numeric(union$FGA)

#Base final 

mi.final <-merge(union, primera, all=T)

#Borrar otras bases 


primera <- NULL
segunda <- NULL 
X <- NULL 
x1 <- NULL 
x2 <- NULL
x3 <- NULL
x4 <- NULL
union <- NULL

# eliminar columna "AMEL"

mi.final <- subset(mi.final, select = -AMEL)

# convertir columnas (exceptuando "Sample.name" y "Alleles") a clase numérica

i <- c(2:16)
mi.final[ ,i] <- apply(mi.final[ , i], 2,
                       function(x) as.numeric(as.character(x)))
sapply(mi.final, class) # comprobar clase de cada columna

# truncar columnas que incluyen decimales (decimal corresponde a número de bases de una repetición incompleta) 

mi.final$TH01 = trunc(mi.final$TH01)
mi.final$D21S11 = trunc(mi.final$D21S11)
mi.final$Penta.E = trunc(mi.final$Penta.E)
mi.final$Penta.D = trunc(mi.final$Penta.D)
mi.final$FGA = trunc(mi.final$FGA)


rm(primera, segunda, union, X, x1, x2, x3, x4)

## Agregar community
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
## Quitar los que no tienen community O tienen NA
mi.final <- mi.final[!(apply(mi.final, 1, function(row) any(row %in% c('OTHER', NA)))), ]
mi.final

## Sample.Name to Ind & Community to Pop
STR_ll <- mi.final %>%
  rename(Ind = Sample.Name,
         Pop = community)

