##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################


#### OBJETIVO 3 ####
### To build a phylogenetic tree showing relationships between communities based on the distributions of genetic microsatellite markers within and between communities.###

## Cargar paquetes y librerias
library(adegenet)
library(readr)
library(tidyverse)

## Cargar bases de datos ##
STR <- read.csv("scripts_fv/Datos/STR.csv", sep = ",")
STR$pop <- gsub(" ", "_", STR$pop)
colnames(STR) <- gsub("\\.1$", " ", colnames(STR))
write.table(STR, file = "scripts_fv/Datos/STR1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
STR [, 3:ncol(STR)] <- lapply(STR[, 3:ncol(STR)], function(x) sprintf("%02d", x))
locis <- colnames(STR[,3:32])
locis_unicos <- unique(trimws(locis))

pop<-STR$pop

STR_alelos_slash  <- read_delim("scripts_fv/Datos/STR_ceros.csv",";", escape_double = FALSE, trim_ws = TRUE)

STR_alelos_slash <- STR_alelos_slash %>% 
  unite( D3S1358, D3S1358.1, D3S1358.2, sep = "/", remove = TRUE) %>% 
  unite(TH01, TH01.1, TH01.2, sep = "/", remove = TRUE) %>% 
  unite(D21S11, D21S11.1, D21S11.2, sep = "/", remove = TRUE) %>% 
  unite(D18S51, D18S51.1, D18S51.2, sep = "/", remove = TRUE) %>% 
  unite(Penta_E, Penta_E.1, Penta_E.2, sep = "/", remove = TRUE) %>%
  unite(D5S818, D5S818.1, D5S818.2, sep = "/", remove = TRUE) %>% 
  unite(D13S317, D13S317.1, D13S317.2, sep = "/", remove = TRUE) %>% 
  unite(D7S820, D7S820.1, D7S820.2, sep = "/", remove = TRUE) %>% 
  unite(D16S539, D16S539.1, D16S539.2, sep = "/", remove = TRUE) %>% 
  unite(CSF1PO, CSF1PO.1, CSF1PO.2, sep = "/", remove = TRUE) %>% 
  unite(Penta_D, Penta_D.1, Penta_D.2, sep = "/", remove = TRUE) %>% 
  unite(vWA, vWA.1, vWA.2, sep = "/", remove = TRUE) %>% 
  unite(D8S1179, D8S1179.1, D8S1179.2, sep = "/", remove = TRUE) %>% 
  unite(TPOX, TPOX.1, TPOX.2, sep = "/", remove = TRUE) %>% 
  unite(FGA, FGA.1, FGA.2, sep = "/", remove = TRUE)

STR_alelos_slash <- STR_alelos_slash[ !(STR_alelos_slash$ID %in% c("AS", "DS", "JG", "MR", "NM", "VG", "WW182", "pos_ctrl")), ]
STR_alelos_slash<-as.data.frame(STR_alelos_slash)
row.names(STR_alelos_slash) <- STR_alelos_slash$ID 
STR_alelos_slash[1] <- NULL


STR_genind <-df2genind(STR_alelos_slash,sep = "/",ncode = 3,ind.names = row.names(STR_alelos_slash),loc.names = colnames(STR_alelos_slash),pop = pop,ploidy = 2,NA.char = "000")

<<<<<<< HEAD
<<<<<<< Updated upstream
=======
=======
>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22
#=======
#TODO: REVISAR. La última línea me arroja un Warning:
#Warning message:
#  In df2genind(STR_alelos_slash, sep = "/", ncode = 3, ind.names = row.names(STR_alelos_slash),  :
#                 Individuals with no scored loci have been removed
<<<<<<< HEAD
#=======
<<<<<<< Updated upstream
## REVISIÓN: Aquí se remueve el individuo BZ023 el cual no tiene 000/000 en todos los loci
>>>>>>> Stashed changes
=======
## REVISIÓN: Aquí se remueve el individuo BZ023 el cual no tiene 000/000 en todos los loci
>>>>>>> Stashed changes
=======
#=======
>>>>>>> 841c4755a13e22ff3c2cbd31b954c62774cf7b22
