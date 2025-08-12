# Alternative names for some objects
## === INTRO === ##

# Define some naming conventions for objects:
# snake_case https://en.wikipedia.org/wiki/Snake_case

# D* = Distance matrices

# Da, Ds  # surname distances (1 - Hedrick’s H)
# Dg_*            # genetic distances (pick your chosen, e.g., D_ps)
# Da_geo, D_geo # geography distances (km)

# T* = Trees

# Ta  = surnames tree (all tips, UPGMA, ultrametric)
# Ts  = surnames tree (16 tips subset, UPGMA)
# Tg  = genetic tree (16 tips, UPGMA on Dps)  # final choice per Methods
# Tc  = consensus tree (15/16 tips — verify)

## === PROCEDURE === ##
# We don't need to do it all over again.
# We only need to add the renaming script at the end of the wrapper

# To be filled and uncomented

## == MATRICES == #

## All communities:
Da <-  surname_matrix # Surnames distances for all communities
Da_geo <- geo_filtered # Geographic distances for all communities

## Sampled (16) communities (15?).
Ds <- surname_matrix_muestra# Surnames distances for sampled communities
D_geo <- geo_muestra # Geographic distances for sampled communities

## Plus many genetic distances:

 Dg_dps <- DPS
 Dg <- Dg_dps # Since this is the one we choose to use in the end
# D_ch <- cs #Object cs not found
 # D_sw <- DSW #Object cs not found
 # G_st <- GST #Not found
 # D_s <- Nei # Nei object not found
# D_m <- 
 #F_st <- FST #FST not found
# theta_w <- 
# D_r <- 
# C_p <- 
# D_a <- 
# X2 <- 
# R_st <- RST
# phi <- 
# D_st <- 
# delta_mu2 <- Dmu2

##== CONSENSUS ==##
# I don't know if this was made from a Tree or from matrices. @VasEstay

## == TREES == #

## Ta (All communities)
Ta <- y_total

## Ts (Surnames in sampled communities)
Ts <- hy

## Tg (Genetic (Dps) in sampled communties)
Tg <- phyDPS

## Tc (Consensus (Dps x Surnames) in sampled communities)
# Tc_ori <- 
# Tc_logit <- 

## T_geo (Geographic distance matrix converted to tree)
# T_geo <- 
Tc <- consensus_tree

## == TRAITS == #

result_traits # Traits in all communities

GM_df # Only G and M traits in all communities
