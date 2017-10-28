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

dendroplot <- function(hc) {
  # hcd = as.dendrogram(clust_hedkin)
  # # plot(clust_hedkin, hang=-1, sub="", xlab="", ylab="")
  # # xlab="Communities", ylab="Time")
  # # cut dendrogram in 4 clusters
  # # hcd2 = cutree(hcd, 4)
  # # function to get color labels
  # colour_labels <- function(n) {
  #   if (is.leaf(n)) {
  #     a <- attributes(n)
  #     attr(n, "nodePar") <-
  #       c(a$nodePar, lab.col=colour_community(a$label))
  #   }
  #   n
  # }
  # # hcd2 = dendrapply(hcd, colour_labels)
  # # pdf("dendro.pdf", width=12, height=60)
  # # plot(hcd, sub="", xlab="", ylab="", horiz=TRUE)
  # pdf("dendro.pdf", width=60, height=12)
  # plot(hcd, sub="", xlab="", ylab="", cex=30)
  # # pdf("dendro.pdf", width=60, height=12)
  # # plot(clust_hedkin, hang=-1, sub="", xlab="", ylab="")
  # generate dendrogram from hclust data
  hcd <- dendro_data(hc, type="rectangle")
  # make palette
  communes <- read_commune_csv()
  colours <- colorRampPalette(c("yellow", "brown"))(
    nrow(communes))
  names(colours) <- communes$commune
  # output to pdf
  pdf("dendro.pdf", width=12, height=40)
  # plot
  p <- ggplot() +
    geom_segment(data=segment(hcd),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(hcd),
              aes(x=x, y=y, label=label, hjust=0,
                  colour=colours[label]),
              nudge_y=0.01,
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
  dendroplot(clust_hedkin)
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

