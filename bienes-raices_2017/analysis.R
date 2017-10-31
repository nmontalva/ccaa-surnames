library(stringr)
library(tidyverse)
library(Biodem)
library(ggdendro)
library(reldist)
# library(ade4)

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

stats_communities <- function(commoners) {
  commoners %>%
    group_by(region, province, commune, community) %>%
    summarise(N=n(),
              S=n_distinct(surname_father)/N,
              R=mean(shares),
              G=gini(shares),
              A=mean(sex == "M"))
}

dendroplot <- function(hc, sc, save_as=NULL) {
  # generate dendrogram from hclust data
  hcd <- dendro_data(hc, type="rectangle")
  # get rid of those factors
  hcd$labels$label <- as.character(hcd$labels$label)
  # vectors to obtain commune, S, R, G, A of communities
  vector_of <- function(col) {
    v <- sc[[col]]
    names(v) <- sc$community
    v
  }
  # province <- vector_of("province")
  commune <- vector_of("commune")
  N <- vector_of("N")
  S <- vector_of("S")
  R <- vector_of("R")
  G <- vector_of("G")
  A <- vector_of("A")
  # output to pdf
  if (!is.null(save_as))
    pdf(save_as, width=12, height=40)
  # plot
  p <- ggplot() +
    geom_segment(data=segment(hcd),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(hcd),
              aes(x=x, y=y, label=label, hjust=0),
              nudge_y=0.01,
              size=3) +
    annotate("text", y=-0.30, x=172.5,
             label="Community", fontface="bold", size=4) +
    annotate("text", y=-1.30, x=172.5,
             label="Commune", fontface="bold", size=4) +
    annotate("text", y=-1.64, x=172.5,
             label="#", fontface="bold", size=4) +
    annotate("text", y=-1.86, x=172.5,
             label="S", fontface="bold", size=4) +
    annotate("text", y=-2.06, x=172.5,
             label="R", fontface="bold", size=4) +
    annotate("text", y=-2.26, x=172.5,
             label="G", fontface="bold", size=4) +
    annotate("text", y=-2.46, x=172.5,
             label="A", fontface="bold", size=4) +
    geom_text(data=label(hcd),
              aes(x=x, y=y, hjust=0,
                  label=commune[label],
                  colour=commune[label]),
              nudge_y=1.1,
              size=3,
              show.legend=FALSE) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=N[label],
                  size=N[label]*1.3/max(N[label])+2.7),
              nudge_y=1.6,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(S[label], digits=3),
                  size=S[label]*1.3/max(S[label])+2.7),
              nudge_y=1.8,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(R[label], digits=3),
                  size=R[label]*1.3/max(R[label])+2.7),
              nudge_y=2.0,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(G[label], digits=3),
                  size=G[label]*1.3/max(G[label])+2.7),
              nudge_y=2.2,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(A[label], digits=3),
                  size=A[label]*1.3/max(A[label])+2.7),
              nudge_y=2.4,
              hjust=0) +
    scale_size_identity() +
    # scale_colour_manual(values=colours) +
    # scale_colour_brewer(palette="Paired") +
    # scale_colour_discrete(l=50) +
    coord_flip() + scale_y_reverse(expand=c(0, 0.6)) +
    theme(axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank()) +
    labs(y="")
  # save plot
  if (!is.null(save_as)) {
    # ggsave(save_as, width=12, height=60, units="cm")
    print(p)
    dev.off()
  } else p
}

surname_distance_matrix <- function(commoners) {
  # cross tabulate
  surnames_freq <- table(commoners$surname_father,
                         commoners$community)
  # generates Hedrick (1971) kinship matrix
  # there are other methods (i.e. lasker, uri)
  hedkin <- hedrick(surnames_freq)
  # hedrick returns values of similarity
  # transform them into values of dissimilarity (distance)
  as.dist(1-hedkin)
}

hclust_default_method <- "complete"
surname_clustering <- function(commoners,
                               method=hclust_default_method) {
  surname_distance_matrix(commoners) %>%
    # hierarchical clustering of the distance matrix
    hclust(method=method)
}

surname_dendrogram <- function(commoners, save_as=NULL,
                               method=hclust_default_method) {
  hc <- surname_clustering(commoners, method)
  sc <- stats_communities(commoners)
  dendroplot(hc, sc, save_as=save_as)
}

k_histogram <- function(commoners,
                        k,
                        save_as=NULL,
                        hclust_method=hclust_default_method) {
  sc <- stats_communities(commoners)
  sur_clust <- surname_clustering(commoners, hclust_method)
  kgroups <- cutree(sur_clust, k=k)
  sc <- tibble(community=names(kgroups),
               branch=as.character(kgroups)) %>%
    right_join(sc, by="community")
  if (!is.null(save_as))
    pdf(save_as, width=5, height=4)
  p <- ggplot(data=sc) +
    geom_histogram(aes(x=G, y=..density.., fill=branch),
                   bins=10, alpha=1/k,
                   position="identity", boundary=0)
  if (!is.null(save_as)) {
    # ggsave(save_as, width=12, height=60, units="cm")
    print(p)
    dev.off()
  } else p
}

# compare surnames matrix with G matrix
# compare_dist <- function(commoners) {
#   sc <- stats_communities(commoners)
#   G_dist <- dist(sc$G)
#   hedkin_dist <- surname_distance_matrix(commoners)
#   # R: I don't understand why one would use quasieuclid or not
#   mantel.randtest(hedkin_dist, # quasieuclid(hedkin_dist),
#                   G_dist, nrepet=9999)
# }

