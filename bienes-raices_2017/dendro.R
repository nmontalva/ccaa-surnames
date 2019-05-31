library(ggdendro)
# library("ggpubr")

dendroplot <- function(hc, save_as=NULL,
                       group_by_col="community") {
  # generate dendrogram from hclust data
  hcd <- dendro_data(hc, type="rectangle")
  # get rid of those factors
  hcd$labels$label <- as.character(hcd$labels$label)
  # traits
  container <- if (group_by_col == "community") "commune"
  else if (group_by_col == "commune") "province"
  else if (group_by_col == "province") "region"
  tc <- traits(commoners, c(group_by_col, container))
  # vectors to obtain commune, S, R, G, A of communities
  vector_of <- function(target_col) {
    v <- tc[[target_col]]
    names(v) <- tc[[group_by_col]]
    v
  }
  location <- if (!is.null(container))
    vector_of(container)
  N <- vector_of("N")
  S <- vector_of("S")
  R <- vector_of("R")
  G <- vector_of("G")
  A <- vector_of("A")
  # useful coordinates
  lastrow <- nrow(hcd$labels)
  x0 <- hcd$labels$x[[lastrow]]
  y0 <- hcd$labels$y[[lastrow]]
  x1 <- x0 + 1 + 0.5 * lastrow / 170
  ydiff <- if (is.null(container)) 0.4 else 0
  # output to pdf
  if (!is.null(save_as)) {
    pdf(save_as,
        width=8 + 4 * lastrow / 170,
        height=1 + 40 * lastrow / 170)
  }
  size <- function(xs) {
    xs * 1.3 / max(xs) + (1.7 + lastrow / 170)
  }
  # plot
  p <- ggplot() +
    geom_segment(data=segment(hcd),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    # because of the coord_flip at the end, x and y are flipped
    geom_text(data=label(hcd),
              aes(x=x, y=y, label=label, hjust=0),
              nudge_y=0.01,
              size=3) +
    annotate("text", x=x1, y=y0-0.215, # y=y0-0.30,
             label=str_to_title(group_by_col),
             fontface="bold", size=4) +
    (if (!is.null(container))
      annotate("text", x=x1, y=y0-1.265, # y=y0-1.30,
               label=str_to_title(container),
               fontface="bold", size=4)) +
    annotate("text", x=x1, y=y0-1.64+ydiff,
             label="#", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-1.86+ydiff,
             label="S", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-2.06+ydiff,
             label="R", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-2.26+ydiff,
             label="G", fontface="bold", size=4) +
    annotate("text", x=x1, y=y0-2.46+ydiff,
             label="A", fontface="bold", size=4) +
    (if (!is.null(container))
      geom_text(data=label(hcd),
                aes(x=x, y=y, hjust=0,
                    label=location[label],
                    colour=location[label]),
                nudge_y=1.1,
                size=3,
                show.legend=FALSE)) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=N[label],
                  size=size(N[label])),
              nudge_y=1.6-ydiff,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(S[label], digits=3),
                  size=size(S[label])),
              nudge_y=1.8-ydiff,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(R[label], digits=3),
                  size=size(R[label])),
              nudge_y=2.0-ydiff,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(G[label], digits=3),
                  size=size(G[label])),
              nudge_y=2.2-ydiff,
              hjust=0) +
    geom_text(data=label(hcd),
              aes(x=x, y=y,
                  label=prettyNum(A[label], digits=3),
                  size=size(A[label])),
              nudge_y=2.4-ydiff,
              hjust=0) +
    scale_size_identity() +
    coord_flip() +
    scale_y_reverse(expand=c(0, 0.6)) +
    theme_dendro()
  # save plot
  if (!is.null(save_as)) {
    # ggsave(save_as, width=12, height=60, units="cm")
    print(p)
    dev.off()
  } else p
}

surname_dendrogram <- function(commoners, save_as=NULL,
                               hclust_method=hclust_default_method,
                               group_by_col="community") {
  hc <- surname_clustering(commoners, hclust_method, group_by_col)
  dendroplot(hc, save_as, group_by_col)
}

