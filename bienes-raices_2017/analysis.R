library(stringr)
library(tidyverse)
library(Biodem)

commoners_csv <- "commoners.csv"
col_names <- c("community", "right_id",
               "firstname", "surname",
               "shares", "commune",
               "firstname1", "firstname2",
               "surname_father", "surname_mother")
read_commoners_csv <- function(filename=commoners_csv) {
  read_csv(commoners_csv,
           skip=1,
           col_names=col_names,
           col_types="ciccdccccc")
}

# frequencies <- function(df, col="surname_father") {
#   df %>%
#     # group_by_(col, "community") %>%
#     group_by_("community", col) %>%
#     summarise(freq=n()) %>%
#     arrange(community)
# }

dendrogram <- function(df) {
  # cross tabulate
  surnames_freq <- table(df$surname_father, df$community)
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # a very simple method. see also as.dist(), dist()
  hedkin.dist <- as.dist(1-hedkin) 
  # hierarchical clustering of the distance matrix
  clust_hedkin <- hclust(hedkin.dist)
  # plot the dendrogram
  plot(clust_hedkin)
}

surname_diversity <- function(surnames_freq) {
  data_surnames <- as.data.frame(surnames)
  # disaggreagete frequencies
  desag_surnames <- data.frame(
    rep(data_surnames$Cognome, data_surnames$Freq),
    rep(data_surnames$Population, data_surnames$Freq))
  # put readable columns' names
  colnames(desag_surnames) <- c("Cognome","Population")
  
  # a self-made function to compute "s"
  fun_s <- function(x) length(unique(x$Cognome))
  # get "s" of each population
  ss <- by(desag_surnames, desag_surnames$Population, fun_s)
  # get "n" of each population
  ns <- by(desag_surnames, desag_surnames$Population, nrow)
  
  # this will make a table of n and s by population
  S_table <- cbind(sapply(ns, I), sapply(ss, I))
  colnames(S_table) <- c("n", "s")
  
  # get S
  S_table <- as.data.frame(S_table)
  S_table$S <- S_table$s / S_table$n
  S_table
}

