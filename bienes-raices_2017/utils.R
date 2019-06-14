library('ggplot2')

# if rownames_col is not null,
# we use that column for the rownames of the vector
column_to_vector <- function(dfr,
                             col,
                             rownames_col=NULL) {
  vec <- dfr[[col]]
  if (!is.null(rownames_col))
    names(vec) <- dfr[[rownames_col]]
  vec
}

saving_plot_as <- function(filename,
                           width=NA,
                           height=NA,
                           ratio=NA, # width/height
                           units=c("in", "cm", "mm"),
                           plot_fn=NULL) {
  if (is.null(plot_fn))
    stop("plot_fn must be specified")
  if (is.na(height) && !any(is.na(c(ratio, width)))) {
    height <- width/ratio
  } else if (is.na(width) && !any(is.na(c(ratio, height)))) {
    width <- height*ratio
  }
  # output to pdf
  # if (!is.null(filename))
  #   pdf(filename, width=width, height=height)
  p <- plot_fn()
  # save plot
  if (!is.null(filename)) {
    ggsave(filename, plot=p, width=width, height=height, units=units)
    # print(p)
    # dev.off()
  } else p
}
