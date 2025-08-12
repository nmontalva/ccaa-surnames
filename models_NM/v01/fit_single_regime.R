#' Fit single regime models (BM and OU)
#' @param prepped Output from prep_traits
#' @return List with BM and OU model fits
fit_single_regime <- function(prepped) {
  require(geiger)
  
  bm <- fitContinuous(prepped$tree, prepped$var_vec, model = "BM")
  ou <- fitContinuous(prepped$tree, prepped$var_vec, model = "OU")
  
  return(list(
    BM = bm,
    OU = ou,
    AIC = data.frame(
      Model = c("BM", "OU"),
      AIC = c(bm$opt$aic, ou$opt$aic),
      stringsAsFactors = FALSE
    )
  ))
}
