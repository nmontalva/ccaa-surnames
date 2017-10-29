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

dendroplot <- function(df, hc) {
  # generate dendrogram from hclust data
  hcd <- dendro_data(hc, type="rectangle")
  # get rid of those factors
  hcd$labels$label <- as.character(hcd$labels$label)
  # make palette
  communities <- df %>%
    group_by(region, province, commune, community) %>%
    summarise()
  communes <- unique(communities$commune)
  colours <- colorRampPalette(c("yellow", "brown"))(
    length(communes))
  names(colours) <- communes
  commune <- communities$commune
  names(commune) <- communities$community
  # output to pdf
  pdf("dendro.pdf", width=12, height=40)
  # plot
  p <- ggplot() +
    geom_segment(data=segment(hcd),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(hcd),
              aes(x=x, y=y, label=label, hjust=0),
              nudge_y=0.01,
              size=3) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=commune[label], hjust=0,
                  colour=colours[commune[label]]),
              nudge_y=0.7,
              size=3,
              show.legend=FALSE) +
    coord_flip() + scale_y_reverse(expand=c(0, 0.6)) +
    theme(axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank()) +
    labs(y="")
  # save plot
  # ggsave("dendro.pdf", width=12, height=60, units="cm")
  print(p)
  dev.off()
}

dendrogram <- function(df) {
  # cross tabulate
  surnames_freq <- table(df$surname_father, df$community)
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # a very simple method. see also as.dist(), dist()
  hedkin_dist <- as.dist(1-hedkin)
  # hierarchical clustering of the distance matrix
  clust_hedkin <- hclust(hedkin_dist)
  # plot the dendrogram
  dendroplot(df, clust_hedkin)
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

