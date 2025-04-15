##########################################################################################################################################################
#########Project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.############
##########################################################################################################################################################

#### OBJETIVO 4 ####
### To compare trees builds with different data and their deviations from the consensus tree ###

#### Otras formas de evaluar similitud entre ?rboles ####
##librer?as
library(dendextend)
library(TreeDist)

##cor.dendlist
cor.dendlist(dendlist(d1 = sort(dend.DPS), d2 = sort(hd)), method = "FM_index",k=9) #ERROR, revisar

##cor_cophenetic
cor_cophenetic(dend.DPS,hd,method_coef = "kendall")

##Bk_plot
Bk(phyNei,hy,k = nleaves(phyCS)-1)
Bk_plot(
  phyCS,
  hy
)
Bk_permutations(phyRST,hy)
cor_FM_index(dend.rst,hd, k=5)

#trees ultrametric and rooted
GST_u <- force.ultrametric(njGST, method = c("extend"))
GST_u$root.edge <- 0
Nei_u <-force.ultrametric(njNei,method = c("extend"))
Nei_u$root.edge <- 0
CS_u <-force.ultrametric(njCS,method = c("extend"))
CS_u$root.edge <- 0
RST_u <-force.ultrametric(njRST,method = c("extend"))
RST_u$root.edge <- 0
ASD_u <-force.ultrametric(njASD,method = c("extend"))
ASD_u$root.edge <- 0
RST2_u <-force.ultrametric(njRST2,method = c("extend"))
RST2_u$root.edge <- 0
DSW_u <-force.ultrametric(njDSW,method = c("extend"))
DSW_u$root.edge <- 0
Dmu2_u <-force.ultrametric(njDmu2,method = c("extend"))
Dmu2_u$root.edge <- 0
FST_u <- force.ultrametric(njFST, method = c("extend"))
FST_u$root.edge <- 0

Geo_tree <- upgma(as.dist(geo_muestra),method="average")
plotTree(Geo_tree)
Geo_tree_nj <- nj(as.dist(geo_muestra))
plotTree(Geo_tree_nj)
Geo_u <-force.ultrametric(Geo_tree_nj,method = c("extend"))
Geo_u$root.edge <- 0
hc_u <-force.ultrametric(hc_nj, method=c("extend"))
hc_u$root.edge <- 0
plotTree(hc_u)

# Encontrar p-valor
set.seed(10000)
the_cor <- cor_bakers_gamma(hy,hy)
the_cor2 <- cor_bakers_gamma(as.dendrogram(phyDPS), as.dendrogram(hy))
R <- 1000
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- hd
for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = F)
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(hd, dend_mixed)
}
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft", legend = c("cor", "cor2"), fill = c(2,1))
sum(the_cor2 < cor_bakers_gamma_results)/ R
the_cor2
set.seed(NULL)

##dist.dendlist
dist.dendlist(dendlist(d1 = as.dendrogram(hcrst), d2 = hd))

## Compare phylo
comparePhylo(ape::as.phylo(hcrst), hy, plot = TRUE, force.rooted = TRUE) # SMM es el ?nico ?rbol que comparte split con el de apellidos

## TreeDist
distance <- TreeDistance(phyRST,hy)
distance
dis_info<-ClusteringInfoDistance(phyRST, hy, reportMatching = TRUE)
dis_info

visual<-VisualizeMatching(ClusteringInfoDistance, phyRST,hy)  ## Sí corre

## Tanglegrama
e <- dendlist(as.dendrogram(phyRST), as.dendrogram(hy)) %>%
  dendextend::untangle(method = "step1side") %>% # Find the best alignment layout
  entanglement()                     # Alignment?quality
e

# Align and plot two dendrograms side by side
dendlist(as.dendrogram(hcrst), as.dendrogram(hc)) %>%
  dendextend::untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(highlight_distinct_edges = FALSE, # Turn-off dashed lines,
             main = paste("entanglement =", round(e,2)),
             main_left = "STRs",
             main_right = "Surnames",
             sort = TRUE,
             color_lines = TRUE,
             intersecting = FALSE,
             k_branches = 4,
             rank_branches = TRUE 
  )                       # Draw the?two?dendrograms


dend_list <- dendlist(dend.rst, hd)
library(corrplot)
corrplot(cor.dendlist(dend_list, method = "baker"),"pie","lower") #ERROR
cor.dendlist(dend_list, method = "baker")
#??cor.denlist

par(mfrow=c(1,1))

## NOTA: Este era un script desordenado en dónde probaba distintos tipos de comparaciones filogenéticas. De todas estas comparaciones se terminó por usar Baker-Gamma, deberiamos dejarlo en el objetivo 4.1 y quitar este sript del github