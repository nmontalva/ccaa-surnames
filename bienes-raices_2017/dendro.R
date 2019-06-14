library(ggplot2)
library(ggdendro)
library(glue)

base_dendrogram <- function(dd) {
  # hjust = 0.5 for centered text
  # hjust = 0 for left-aligned text
  ggplot(data = label(dd), aes(x = x, y = y, hjust = 0)) +
    geom_segment(data = segment(dd),
                 aes(xend = xend, yend = yend)) +
    scale_size_identity() +
    coord_flip() +
    # TODO: what is the expand parameter for and why 0.6?
    scale_y_reverse(expand = c(0, 0.6)) +
    theme_dendro() # legend.position="top")
}

add_leaves <- function(plot,
                       # supra_division,
                       # colour_by_var=NULL,
                       title="",
                       x_title = 0,
                       y_title = 0,
                       y_offset = 0.01) {
  # sd_var <- sym(supra_division)
  plot +
    annotate("text", x = x_title,
             y = y_title - y_offset, # + 0.05,
             label = str_to_title(title),
             fontface = "bold", size = 4) +
    # TODO: try geom_label
    # TODO: conditionally add colours depending on whether
    # colour_by_var is null or not
    geom_text(aes(label = label,
                  # colour = !!sd_var),
                  size = 3),
              # x and y are flipped because of coord_flip
              nudge_y = y_offset) # +
    # guides(fill=guide_legend(title.position="top"))
}

# adm_division: administrative division
get_supra_division <- function(adm_division) {
  if (adm_division == "community") "commune"
  else if (adm_division == "commune") "province"
  else if (adm_division == "province") "region"
  else stop(glue("administrative division ",
                 "'{adm_division}' not recognised"))
}

# TODO: colour scale could relate to some sort of geographic ordering
add_supra_division <- function(plot,
                               supra_division,
                               x_title = 0,
                               y_title = 0,
                               y_offset = 0) {
  sd_var <- sym(supra_division)
  plot +
    annotate("text", x = x_title,
             y = y_title - y_offset, # + 0.05,
             label = str_to_title(supra_division),
             fontface = "bold", size = 4) +
    geom_text(aes(label = !!sd_var,
                  colour = !!sd_var,
                  size = 3),
              nudge_y = y_offset,
              show.legend = FALSE)
}

add_trait <- function(plot,
                      trait_name,
                      x_title = 0,
                      y_title = 0,
                      y_offset = 0,
                      size_fn = 3) {
  if (is.numeric(size_fn))
    size_fn <- function(x) size_fn
  else if (!is_function(size_fn))
    stop("parameter 'size_fn' not numeric nor a function")
  trait_var <- sym(trait_name)
  plot +
    # TODO: is there some parameter like TikZ's anchor=...?
    annotate("text", x = x_title,
             y = y_title - y_offset, # + 0.05,
             label = str_to_title(trait_name),
             fontface = "bold", size = 4) +
    geom_text(aes(label = prettyNum(!!trait_var, digits=3),
                  size = size_fn(!!trait_var)),
              nudge_y = y_offset)
}

surname_dendrogram <- function(commoners,
                               save_as=NULL,
                               hclust_method=hclust_default_method,
                               group_by="community") {
  hc <- surname_clustering(commoners, hclust_method, group_by)
  # generate dendrogram from hclust data
  # type="rectangle" will draw rectangular lines
  dd <- dendro_data(hc, type="rectangle")
  # get rid of those factors
  # apart from getting rid of the warning in the left_join,
  # is this useful? does it have unintended consequences?
  dd$labels$label <- as.character(dd$labels$label)
  # traits
  supra_division <- get_supra_division(group_by)
  ts <- traits(commoners, c(group_by, supra_division))
  # add traits to dd
  dd$labels <- dd$labels %>%
    left_join(ts, by=c("label" = group_by))
  # traits
  # supra_division <- get_supra_division(group_by)
  # ts <- traits(commoners, c(group_by, supra_division))
  # hc$order[i] = j means row j moves to row i (counterintuitevely)
  # reorder `ts` rows using hc$order and then invert it
  # ts <- ts[hc$order,][nrow(ts):1,]
  #
  # vectors for commune, S, R, G, A of communities
  # location <- if (!is.null(supra_division))
  #   column_to_vector(ts, supra_division, group_by)
  # N <- column_to_vector(ts, "N", group_by)
  # S <- column_to_vector(ts, "S", group_by)
  # R <- column_to_vector(ts, "R", group_by)
  # G <- column_to_vector(ts, "G", group_by)
  # A <- column_to_vector(ts, "A", group_by)
  # some useful coordinates, lengths and sizes
  lastrow <- nrow(dd$labels)
  x0 <- dd$labels$x[[lastrow]]
  y0 <- dd$labels$y[[lastrow]]
  x1 <- x0 + 1 + 0.5 * lastrow / 170
  y1 <- y0 - 0.215 # 0.30
  # y_offset <- if (is.null(supra_division)) 0.4 else 0
  # why this function?
  size_of_label <- function(xs) {
    xs * 1.3 / max(xs) + (1.7 + lastrow / 170)
  }
  saving_plot_as(
    save_as,
    width = 8 + 4 * lastrow / 170,
    height = 1 + 40 * lastrow / 170,
    plot_fn = function() {
      base_dendrogram(dd) %>%
        add_leaves(group_by, x1, y1, 0.01) %>%
        add_supra_division(
          supra_division, x1, y0-0.165, 1.1) %>%
        add_trait("N", x1, y0-0.06, 1.6, # 1.475, # y_offset,
                  size_fn = size_of_label) %>%
        add_trait("S", x1, y0-0.06, 1.8,
                  size_fn = size_of_label) %>%
        add_trait("R", x1, y0-0.06, 2.0,
                  size_fn = size_of_label) %>%
        add_trait("G", x1, y0-0.06, 2.2,
                  size_fn = size_of_label) %>%
        add_trait("A", x1, y0-0.06, 2.4,
                  size_fn = size_of_label)
    })
    # plot_fn=function() {
    #   ggplot() +
    #     geom_segment(data=segment(dd),
    #                  aes(x=x, y=y, xend=xend, yend=yend)) +
    #     # because of the coord_flip at the end, x and y are flipped
    #     geom_text(data=label(dd),
    #               aes(x=x, y=y, label=label, hjust=0),
    #               nudge_y=0.01,
    #               size=3) +
    #     annotate("text", x=x1, y=y0-0.215, # y=y0-0.30,
    #              label=str_to_title(group_by),
    #              fontface="bold", size=4) +
    #     (if (!is.null(supra_division))
    #       annotate("text", x=x1, y=y0-1.265, # y=y0-1.30,
    #                label=str_to_title(supra_division),
    #                fontface="bold", size=4)) +
    #     annotate("text", x=x1, y=y0-1.64+ydiff,
    #              label="#", fontface="bold", size=4) +
    #     annotate("text", x=x1, y=y0-1.86+ydiff,
    #              label="S", fontface="bold", size=4) +
    #     annotate("text", x=x1, y=y0-2.06+ydiff,
    #              label="R", fontface="bold", size=4) +
    #     annotate("text", x=x1, y=y0-2.26+ydiff,
    #              label="G", fontface="bold", size=4) +
    #     annotate("text", x=x1, y=y0-2.46+ydiff,
    #              label="A", fontface="bold", size=4) +
    #     (if (!is.null(supra_division))
    #       geom_text(data=label(dd),
    #                 aes(x=x, y=y, hjust=0,
    #                     label=location[label],
    #                     colour=location[label]),
    #                 nudge_y=1.1,
    #                 size=3,
    #                 show.legend=FALSE)) +
    #     geom_text(data=label(dd),
    #               aes(x=x, y=y,
    #                   label=N[label],
    #                   size=size_of_label(N[label])),
    #               nudge_y=1.6-ydiff,
    #               hjust=0) +
    #     geom_text(data=label(dd),
    #               aes(x=x, y=y,
    #                   label=prettyNum(S[label], digits=3),
    #                   size=size_of_label(S[label])),
    #               nudge_y=1.8-ydiff,
    #               hjust=0) +
    #     geom_text(data=label(dd),
    #               aes(x=x, y=y,
    #                   label=prettyNum(R[label], digits=3),
    #                   size=size_of_label(R[label])),
    #               nudge_y=2.0-ydiff,
    #               hjust=0) +
    #     geom_text(data=label(dd),
    #               aes(x=x, y=y,
    #                   label=prettyNum(G[label], digits=3),
    #                   size=size_of_label(G[label])),
    #               nudge_y=2.2-ydiff,
    #               hjust=0) +
    #     geom_text(data=label(dd),
    #               aes(x=x, y=y,
    #                   label=prettyNum(A[label], digits=3),
    #                   size=size_of_label(A[label])),
    #               nudge_y=2.4-ydiff,
    #               hjust=0) +
    #     scale_size_identity() +
    #     coord_flip() +
    #     scale_y_reverse(expand=c(0, 0.6)) +
    #     theme_dendro()
    # })
}

