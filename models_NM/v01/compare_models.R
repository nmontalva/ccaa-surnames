compare_models <- function(var_summary) {
  require(dplyr)
  single_aic <- var_summary$single_models$AIC
  multi_aic <- tryCatch({
    val <- var_summary$surface$final_model$aic
    if (length(val) > 1) val <- val[1]  # asegurar único valor
    val
  }, error = function(e) NA)
  
  alpha_val <- var_summary$parameters$alpha
  sigma2_val <- var_summary$parameters$sigma2
  
  multi_df <- if (!is.na(multi_aic)) {
    data.frame(
      Model = "OU Multi",
      AIC = multi_aic,
      Alpha = alpha_val,
      Sigma2 = sigma2_val,
      stringsAsFactors = FALSE
    )
  } else NULL
  
  results_df <- dplyr::bind_rows(single_aic, multi_df)
  
  if (any(is.na(results_df$AIC))) {
    warning("Algunos AIC son NA, no se puede calcular Delta_AIC")
    results_df$Delta_AIC <- NA
  } else {
    results_df <- results_df %>% 
      dplyr::arrange(AIC) %>% 
      dplyr::mutate(Delta_AIC = AIC - min(AIC))
  }
  
  return(results_df)
}
