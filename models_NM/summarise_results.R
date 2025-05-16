#' Summarize SURFACE results
#' @param prepped Output from prep_traits
#' @param surface_res Output from run_surface
#' @param single_regime Output from fit_single_regime
#' @return Comprehensive results summary
summarise_results <- function(prepped, surface_res, single_regime) {
  final_model <- surface_res$final_model
  
  # Initialize default single regime case
  regimes <- rep("single_regime", length(prepped$tree$tip.label))
  theta <- setNames(mean(prepped$var_vec), "single_regime")
  alpha <- sigma2 <- NA
  
  # Check if regimes were found
  if (length(final_model$fit) > 0 && prepped$analysis_var %in% names(final_model$fit)) {
    hmod <- final_model$fit[[prepped$analysis_var]]
    if (inherits(hmod, "hansentree")) {
      regimes <- hmod@regimes$regs[1:length(prepped$var_vec)]
      theta <- tryCatch({
        if (is.matrix(hmod@theta)) hmod@theta[,1] else hmod@theta
      }, error = function(e) {
        setNames(rep(mean(prepped$var_vec), length(unique(regimes))), 
                paste0("regime_", seq_along(unique(regimes))))
      })
      alpha <- tryCatch(hmod@alpha, error = function(e) NA)
      sigma2 <- tryCatch(hmod@sigma.squared, error = function(e) NA)
    }
  }
  
  # Create comprehensive results
  list(
    data = data.frame(
      community = names(prepped$var_vec),
      original_value = prepped$data[[prepped$original_var]],
      transformed_value = prepped$var_vec,
      regime = regimes,
      stringsAsFactors = FALSE
    ),
    single_models = single_regime,
    parameters = data.frame(
      alpha = alpha,
      sigma2 = sigma2,
      stat_var = if (!is.na(alpha) && !is.na(sigma2)) sigma2/(2*alpha) else NA
    ),
    regimes = data.frame(
      regime = names(theta),
      theta = theta,
      theta_original = if (prepped$analysis_var == paste0(prepped$original_var, "_logit")) {
        plogis(theta)
      } else {
        theta
      },
      stringsAsFactors = FALSE
    ),
    surface = surface_res
  )
}