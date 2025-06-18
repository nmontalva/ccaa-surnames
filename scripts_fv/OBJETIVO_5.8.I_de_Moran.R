##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################
#### OBJETIVO 5 ####
### Part 8 ###
### To plot over tips of the resulting trees the estimated values of each trait, and compute the most likely values for each at the internal nodes ###

#### LIBRERÍAS
library(dplyr)
library(geosphere)
library(SDPDmod)
library(sf)
library(spdep)

#Datos geográficos
coordenadas <- read.csv("scripts_fv/Datos/coordenadas.csv", header = T, fileEncoding = "UTF-8-BOM")
coordenadas$community <- gsub(" ", "_", coordenadas$community)
coordenadas$community[grepl("LA_RINCONADA_DE_PUNITAQUI" , coordenadas$community)] <- "RINCONADA_DE_PUNITAQUI"
rownames(coordenadas) <- coordenadas$community
colnames(coordenadas)<- c("community","lon","lat","org_name")
coordenadas<-select(coordenadas,community,lon,lat)

#Data community
#View(GM_logit)
<<<<<<< HEAD
# Ajustar el modelo de regresión de M ~ G
modelo <- lm(M_logit ~ G_logit, data = GM_logit) # Modelo M en función de G
GM_logit$G.reg.M <- predict(modelo)  # Guardar valores predichos
=======
# Ajustar el modelo de regresión de G ~ M
modelo <- lm(G_logit ~ M_logit, data = GM_logit)
G.reg.M <- predict(modelo)
>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01

#Crear df
# Combinar tablas por una clave común (por ejemplo, "community")
GM_logit$community <- row.names(GM_logit)
<<<<<<< HEAD
# Crear dataframe combinado
community_data <- GM_logit %>%
  mutate(community = row.names(.)) %>%  # Asegurar que tenemos nombres de comunidad
  left_join(coordenadas, by = "community") %>%
  select(community, G_logit, M_logit, lon, lat, G.reg.M) %>%
  na.omit()  # Eliminar filas con NA

=======
community_data <- GM_logit %>%
  left_join(coordenadas, by = "community" ) %>%
  mutate(
    G.reg.M = predict(lm(M_logit ~ G_logit, data = .))  # Calcular valores ajustados (G~M)
  ) %>%
  select(community, G_logit, M_logit, lon, lat, G.reg.M)
community_data <- na.omit(community_data)
>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01
# Verificar el resultado
colnames(community_data)
head(community_data)

# Create spatial points from coordinates
#coords <- community_data[, c("lon", "lat")]
points_sf <- st_as_sf(community_data, coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(32619)  # Project to UTM for accurate distance calculations

# Generate distance-based neighbors (e.g., within a distance of 1.5 units)
# Define neighbors - using k-nearest neighbors (k=3)
coords <- st_coordinates(points_sf)
k <- 3 # Number of neighbors
neighbors <- knn2nb(knearneigh(coords, k = k))

# Convert to spatial weights
<<<<<<< HEAD
weights <- nb2listw(neighbors, style = "W") # Estilo de estandarización por filas
weights

# Moran's I for G,M
#moran_G <- moran.test(community_data$G_logit, weights)
#print(moran_G)
=======
weights <- nb2listw(neighbors, style = "W") 

# K-neighbours

# Generate k-nearest neighbors (determining k is also a bit trial an error)
n <- nrow(community_data)
k_est <- sqrt(n)/2
k_est_2 <- n^(1/2)
#neighbors_knn <- knn2nb(knearneigh(coords, k = k_est_2)) #Or k_est_2 ## NO FUNCIONA
# Convert to spatial weights
#weights_knn <- nb2listw(neighbors_knn, style = "W")

# Moran's I for G,M
moran_G <- moran.test(community_data$G_logit, weights)
print(moran_G)
>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01
# I>0, Positive autocorrelation (similar values cluster).
# I=0, Random distribution.
# I<0, Negative autocorrelation (dissimilar values cluster). *** ESTE ES EL CASO

<<<<<<< HEAD
#moran_M <- moran.test(community_data$M_logit, weights)
#print(moran_M)
=======
moran_M <- moran.test(community_data$M_logit, weights)
print(moran_M)
>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01
# I>0, Positive autocorrelation (similar values cluster).
# I=0, Random distribution.
# I<0, Negative autocorrelation (dissimilar values cluster). *** ESTE ES EL CASO

moran_GM <- moran.test(community_data$G.reg.M, weights)
print(moran_GM)
# I>0, Positive autocorrelation (similar values cluster). *** ESTE ES EL CASO
# I=0, Random distribution.
<<<<<<< HEAD
# I<0, Negative autocorrelation (dissimilar values cluster). 

#Monte Carlo simulation for p-value
MC <- moran.mc(community_data$G.reg.M, weights, nsim = 999, alternative = "greater")
print(MC)
# I>0, Positive autocorrelation (similar values cluster). *** ESTE ES EL CASO
# I=0, Random distribution.
# I<0, Negative autocorrelation (dissimilar values cluster). 

=======
# I<0, Negative autocorrelation (dissimilar values cluster). *** ESTE ES EL CASO
>>>>>>> c02a8b857ca6e13ebf759e5063a5bab8e3cdea01
