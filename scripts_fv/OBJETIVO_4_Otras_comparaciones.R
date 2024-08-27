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
cor.dendlist(dendlist(d1 = sort(dend.rst), d2 = sort(hd)), method = "FM_index",k=9) #ERROR, revisar

##cor_cophenetic
cor_cophenetic(dend.rst,hd,method_coef = "kendall")

##Bk_plot
Bk(phyNei,hy,k = nleaves(phyCS)-1)
Bk_plot(
  phyCS,
  hy
)
Bk_permutations(phyRST,hy)
cor_FM_index(dend.rst,hd, k=5)

#cor_bakers_gamma
cor_bakers_gamma(phyGST,hy)
cor_bakers_gamma(phyNei,hy)
cor_bakers_gamma(phyCS,hy)
cor_bakers_gamma(ape::as.phylo(hcrst),hy)
cor_bakers_gamma(phyASD,hy)

# Encontrar p-valor
set.seed(10000)
set.seed(NULL)
the_cor <- cor_bakers_gamma(hy,hy)
the_cor2 <- cor_bakers_gamma(ape::as.phylo(hcrst), hy)

R <- 1000
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- hd
for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = F)
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(as.dendrogram(hcrst), dend_mixed)
}
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft", legend = c("cor", "cor2"), fill = c(2,1))
sum(the_cor2 < cor_bakers_gamma_results)/ R

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

