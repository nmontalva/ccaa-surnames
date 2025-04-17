##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################
#### OBJETIVO 5 ####
### Part 8 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

#### LIBRERÍAS
library(spdep)
library(dplyr)
library(geosphere)
library(SDPDmod)

#Datos geográficos
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T)
coordenadas$community <- gsub(" ", "_", coordenadas$community)
coordenadas$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$community)] <- "RINCONADA_DE_PUNITAQUI"
rownames(coordenadas) <- coordenadas$ï..community
colnames(coordenadas)<- c("community","lon","lat","org_name")
coordenadas<-select(coordenadas,community,lon,lat)

#Data community
View(GM_df)
# Ajustar el modelo de regresión de G ~ M
modelo <- lm(G ~ M, data = GM_df)
G.reg.M <- predict(modelo)

#Crear df
# Combinar tablas por una clave común (por ejemplo, "community")
GM_df$community <- row.names(GM_df)
community_data <- GM_df %>%
  left_join(coordenadas, by = "community" ) %>%
  mutate(
    G.reg.M = predict(lm(M ~ G, data = .))  # Calcular valores ajustados (G~M)
  ) %>%
  select(community, G, M, lon, lat, G.reg.M)
community_data <- na.omit(community_data)
# Verificar el resultado
colnames(community_data)
head(community_data)

# Create spatial points from coordinates
community_data$G.reg.M
coords <- community_data[, c("lon", "lat")]

# Generate distance-based neighbors (e.g., within a distance of 1.5 units)
# Units are somewhat arbitrary, ew may need to find a number by trial an error
neighbors <- dnearneigh(coords, 0, 1.5)

# Convert to spatial weights
weights <- nb2listw(neighbors, style = "W") 

# K-neighbours

# Generate k-nearest neighbors (determining k is also a bit trial an error)
n <- nrow(community_data)
k_est <- sqrt(n)/2
k_est_2 <- n^(1/2)
neighbors_knn <- knn2nb(knearneigh(coords, k = k_est_2)) #Or k_est_2 ## NO FUNCIONA

#TODO Warning message:
#In matrix(z$nn, np, k, byrow = TRUE) :
#  data length [2177] is not a sub-multiple or multiple of the number of rows [168]

# Convert to spatial weights
weights_knn <- nb2listw(neighbors_knn, style = "W")

# Moran's I for G,M
moran_G <- moran.test(community_data$G, weights)
print(moran_G)
# I>0, Positive autocorrelation (similar values cluster).
# I=0, Random distribution.
# I<0, Negative autocorrelation (dissimilar values cluster). *** ESTE ES EL CASO

moran_M <- moran.test(community_data$M, weights)
print(moran_M)
# I>0, Positive autocorrelation (similar values cluster).
# I=0, Random distribution.
# I<0, Negative autocorrelation (dissimilar values cluster). *** ESTE ES EL CASO

moran_GM <- moran.test(community_data$G.reg.M, weights)
print(moran_GM)
# I>0, Positive autocorrelation (similar values cluster).
# I=0, Random distribution.
# I<0, Negative autocorrelation (dissimilar values cluster). *** ESTE ES EL CASO

