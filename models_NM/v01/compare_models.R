compare_models <- function(var_summary) {
  require(dplyr)
  single_aicc <- var_summary$single_models$AICc
  multi_aicc <- tryCatch({
    val <- var_summary$surface$final_model$aic
    if (length(val) > 1) val <- val[1]  # asegurar único valor
    val
  }, error = function(e) NA)
  
  alpha_val <- var_summary$parameters$alpha
  sigma2_val <- var_summary$parameters$sigma2
  
  multi_df <- if (!is.na(multi_aicc) && !is.na(multi_aicc)) {
    data.frame(
      Model = "OU Multi",
      AICc = multi_aicc,
      Alpha = alpha_val,
      Sigma2 = sigma2_val,
      stringsAsFactors = FALSE
    )
  } else NULL
  
  results_df <- dplyr::bind_rows(single_aicc, multi_df)
  
  if (any(is.na(results_df$AICc))) {
    warning("Algunos AICc son NA, no se puede calcular Delta_AICc")
    results_df$Delta_AICc <- NA
  } else {
    results_df <- results_df %>% 
      dplyr::arrange(AICc) %>% 
      dplyr::mutate(Delta_AICc = AICc - min(AICc))
  }
  
  return(results_df)
}
