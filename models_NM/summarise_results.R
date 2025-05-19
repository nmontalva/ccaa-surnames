#' Summarize SURFACE results
#' @param prepped Output from prep_traits
#' @param surface_res Output from run_surface
#' @param single_regime Output from fit_single_regime
#' @return Comprehensive results summary
summarise_results <- function(prepped, surface_res, single_regime) {
  # Initialize default single regime results
  results <- list(
    data = data.frame(
      community = names(prepped$var_vec),
      original_value = prepped$data[[prepped$original_var]],
      transformed_value = prepped$var_vec,
      regime = "single_regime",
      stringsAsFactors = FALSE
    ),
    single_models = single_regime,
    regimes = data.frame(
      regime = "single_regime",
      theta = mean(prepped$var_vec),
      theta_original = if (prepped$analysis_var == paste0(prepped$original_var, "_logit")) {
        plogis(mean(prepped$var_vec))
      } else {
        mean(prepped$var_vec)
      },
      stringsAsFactors = FALSE
    ),
    parameters = data.frame(
      alpha = NA,
      sigma2 = NA,
      stat_var = NA,
      n_shifts = 0,
      delta_aic_vs_OU = NA
    )
  )
  
  # Try to extract surface results if available
  if (!is.null(surface_res$backward) && length(surface_res$backward) > 0) {
    final_model <- surface_res$backward[[length(surface_res$backward)]]
    
    if (length(final_model$fit) > 0 && prepped$analysis_var %in% names(final_model$fit)) {
      hmod <- final_model$fit[[prepped$analysis_var]]
      
      # Get regime assignments
      regimes <- as.character(hmod@regimes$regs[1:length(prepped$var_vec)])
      results$data$regime <- regimes
      
      # Get theta values
      theta <- if (is.matrix(hmod@theta)) hmod@theta[,1] else hmod@theta
      results$regimes <- data.frame(
        regime = names(theta),
        theta = theta,
        theta_original = if (prepped$analysis_var == paste0(prepped$original_var, "_logit")) {
          plogis(theta)
        } else {
          theta
        },
        stringsAsFactors = FALSE
      )
      
      # Get parameters
      results$parameters <- data.frame(
        alpha = hmod@alpha,
        sigma2 = hmod@sigma.squared,
        stat_var = hmod@sigma.squared/(2*hmod@alpha),
        n_shifts = length(unique(regimes)) - 1,
        delta_aic_vs_OU = final_model$aic - single_regime$AIC$AIC[2]
      )
    }
  }
  
  return(results)
}