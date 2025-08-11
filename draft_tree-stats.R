# Load libraries
library(ape)
library(phytools)
library(treestats)

# Read your tree (replace with your actual file)
# Example: Ta_tree <- read.tree("Ta_tree.nwk")
Ta_tree <- y_total

# Basic plot to visualize tree shape
#plot(Ta_tree, cex = 0.4)
#title("Visual inspection of tree shape")

# Colless' index
colless_index <- colless(Ta_tree, normalization = "pda")
cat("Colless' index (normalized):", colless_index, "\n")

# Sackin's index
sackin_index <- sackin(Ta_tree, normalization = "pda")
cat("Sackin's index (normalized):", sackin_index, "\n")

# Additional tip: You can simulate a balanced tree for comparison
#set.seed(123)
#balanced_tree <- stree(n = length(Ta_tree$tip.label), type = "balanced")
#cat("Colless index of perfectly balanced tree:", colless(balanced_tree, normalization = "pda"), "\n")

# For comparison

# Colless' index
colless_index <- colless(consensus_tree, normalization = "pda")
cat("Colless' index (normalized):", colless_index, "\n")

# Sackin's index
sackin_index <- sackin(consensus_tree, normalization = "pda")
cat("Sackin's index (normalized):", sackin_index, "\n")

# Additional tip: You can simulate a balanced tree for comparison
#set.seed(123)
#balanced_tree <- stree(n = length(consensus_tree$tip.label), type = "balanced")
#cat("Colless index of perfectly balanced tree:", colless(balanced_tree, normalization = "pda"), "\n")

# A simple Counter:
sum(GM_df$G == 0)   # Number of communities with G = 0
sum(GM_df$G == 1)   # Number of communities with G = 0
sum(GM_df$M == 0)   # Number with M = 0
sum(GM_df$M == 1)   # Number with M = 1