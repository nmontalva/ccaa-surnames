library(ggplot2)
library(ggdendro)
library(glue)
library(reshape2)

base_dendrogram <- function(dd) {
  # hjust = 0.5 for centered text
  # hjust = 0 for left-aligned text
  ggplot(data = label(dd), aes(x = x, y = y, hjust = 0)) +
    geom_segment(data = segment(dd),
                 aes(xend = xend, yend = yend)) +
    scale_size_identity() +
    coord_flip() +
    # the expand parameter makes room for the labels
    scale_y_reverse(expand = c(0, 0.6)) +
    theme_dendro() # legend.position="top")
}

add_leaves <- function(plot,
                       title="",
                       x_title = 0,
                       y_title = 0,
                       y_offset = 0.01,
                       colour_by_var = NULL,
                       size = 3) {
  if (!(is.numeric(size) || is_function(size)))
    stop("parameter 'size' not numeric nor a function")
  plot +
    annotate("text", x = x_title,
             y = y_title - y_offset,
             label = str_to_title(title),
             fontface = "bold", size = 4) +
    if (is_null(colour_by_var)) {
      geom_text(aes(label = label,
                    size = if (is.numeric(size)) size
                    else size(label)),
                # x and y are flipped because of coord_flip
                nudge_y = y_offset)
    } else {
      geom_text(aes(label = label,
                    colour = !!sym(colour_by_var),
                    size = if (is.numeric(size)) size
                    else size(label)),
                nudge_y = y_offset) +
      guides(fill=guide_legend(title.position="top"))
    }
}

# adm_division: administrative division
get_supra_division <- function(adm_division) {
  if (adm_division == "community") "commune"
  else if (adm_division == "commune") "province"
  else if (adm_division == "province") "region"
  else stop(glue("administrative division ",
                 "'{adm_division}' not recognised"))
}

add_trait <- function(plot,
                      trait_name,
                      x_title = 0,
                      y_title = 0,
                      y_offset = 0,
                      colour_by_var = NULL,
                      size = 3) {
  if (!(is.numeric(size) || is_function(size)))
    stop("parameter 'size' not numeric nor a function")
  trait_var <- sym(trait_name)
  plot +
    # TODO: is there some parameter like TikZ's anchor=...? hjust=0?
    annotate("text", x = x_title,
             y = y_title - y_offset,
             label = str_to_title(trait_name),
             fontface = "bold", size = 4) +
    if (is_null(colour_by_var)) {
      geom_text(aes(label = prettyNum(!!trait_var, digits=3),
                    size = if (is.numeric(size)) size
                    else size(!!trait_var)),
                nudge_y = y_offset)
    } else {
      geom_text(aes(label = prettyNum(!!trait_var, digits=3),
                    colour = !!sym(colour_by_var),
                    size = if (is.numeric(size)) size
                    else size(!!trait_var)),
                nudge_y = y_offset,
                show.legend = FALSE)
    }
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
  # some useful coordinates, lengths and sizes
  lastrow <- nrow(dd$labels)
  x0 <- dd$labels$x[[lastrow]]
  y0 <- dd$labels$y[[lastrow]]
  x1 <- x0 + 1 + 0.5 * lastrow / 170
  # y_offset <- if (is_null(supra_division)) 0.4 else 0
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
        add_leaves(group_by, x1, y0-0.215, 0.01) %>%
        # TODO: colour scale could relate to
        # some sort of geographic ordering
        add_trait(supra_division, x1, y0-0.165, 1.1,
                  colour_by_var = supra_division) %>%
        add_trait("N", x1, y0-0.06, 1.6,
                  size = size_of_label) %>%
        add_trait("S", x1, y0-0.06, 1.8,
                  size = size_of_label) %>%
        add_trait("R", x1, y0-0.06, 2.0,
                  size = size_of_label) %>%
        add_trait("G", x1, y0-0.06, 2.2,
                  size = size_of_label) %>%
        add_trait("A", x1, y0-0.06, 2.4,
                  size = size_of_label)
    })
}

# the inspiration for this plot comes from
# https://stackoverflow.com/questions/44646488/
dendro <- function(commoners,
                   pvalues,
                   save_as=NULL,
                   hclust_method=hclust_default_method,
                   group_by="community") {
  hc <- surname_clustering(commoners, hclust_method, group_by)
  # generate dendrogram from hclust data
  # type="rectangle" will draw rectangular lines
  dd <- dendro_data(hc, type="rectangle")
  lvls <- levels(dd$labels$label)
  # traits
  supra_division <- get_supra_division(group_by)
  ts <- traits(commoners, c(group_by, supra_division))
  ts[[group_by]] <- factor(ts[[group_by]], levels=lvls)
  # num_groups <- nrow(ts)
  # add traits to dd
  dd$labels <- dd$labels %>%
    left_join(ts, by=c("label" = group_by)) %>%
    melt(id.vars=c("label", supra_division, "x", "y"))
  # plotting
  saving_plot_as(
    save_as,
    width = 22,
    height = 41,
    plot_fn = function() {
      dendr <- ggplot(data = label(dd), aes(x = x, y = y)) +
        geom_segment(data = segment(dd),
                     aes(xend = xend,
                         yend = yend)) +
        # the expand parameter makes room for the labels
        scale_y_reverse(expand = c(0, 0.05)) +
        theme_dendro() +
        theme(legend.position="left") +
        geom_tile(aes(y = -0.3,
                      fill = !!sym(supra_division)),
                  height = 0.5) +
        # geom_label(aes(label = label,
        #                fill = !!sym(supra_division)),
        #            nudge_y = 0.01,
        #            size = 3) +
        guides(fill=guide_legend(title = NULL)) +
        coord_flip()
      bars <- ggplot(data = label(dd)) +
        geom_col(aes(x = label,
                     y = value,
                     fill = variable),
                 show.legend = FALSE) +
        facet_grid(. ~ variable, scales = "free_x") +
        scale_x_discrete("") +
        scale_y_continuous("") +
        coord_flip() +
        # theme_transparent()
        theme(rect = element_rect(fill = "transparent"),
              text = element_text(face="bold"))
      colourbar <- colourbar_pvalues(pvalues)
      # using cowplot as suggested in
      # https://stackoverflow.com/questions/44646488/
      cowplot::ggdraw() +
        cowplot::draw_plot(dendr, 0, 0, 0.6, 1) +
        cowplot::draw_plot(bars, 0.42, 0.0405, 0.56, 0.921) +
        cowplot::draw_plot(colourbar, 0.087, 0.01, 0.364, 0.02)
      # another approach from
      # https://stackoverflow.com/questions/6673162/
      # grid.newpage()
      # print(p1, vp=viewport(0.8, 0.8, x=0.4, y=0.4))
      # print(p2, vp=viewport(0.52, 0.2, x=0.45, y=0.9))
    })
}

