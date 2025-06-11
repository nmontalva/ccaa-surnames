summarise_results <- function(prepped, surface_res, single_regime) {
  final_model <- surface_res$final_model
  
  # Initialize default values
  regimes <- rep("single_regime", length(prepped$tree$tip.label))
  theta <- setNames(mean(prepped$var_vec, na.rm = TRUE), "single_regime")
  alpha <- sigma2 <- NA
  
  # Check if regimes were found and extract properly
  if (length(final_model$fit) > 0 && prepped$analysis_var %in% names(final_model$fit)) {
    hmod <- final_model$fit[[prepped$analysis_var]]
    
    if (inherits(hmod, "hansentree")) {
      # Get all regime assignments (nodes + tips)
      all_regimes <- hmod@regimes$regs
      
      # Extract just the tip regimes (first n tips)
      n_tips <- length(prepped$tree$tip.label)
      regimes <- as.character(all_regimes[1:n_tips])
      
      # Extract theta values - handling different SURFACE output formats
      if (is.list(hmod@theta)) {
        theta <- hmod@theta[[1]]  # For list format
      } else if (is.matrix(hmod@theta)) {
        theta <- hmod@theta[,1]   # For matrix format
      } else {
        theta <- hmod@theta       # For vector format
      }
      
      # Extract parameters safely
      alpha <- tryCatch(hmod@alpha, error = function(e) NA)
      sigma2 <- tryCatch(hmod@sigma.squared, error = function(e) NA)
    }
  }
  
  # Create results data frame - ensuring proper ordering
  result_data <- data.frame(
    community = names(prepped$var_vec),
    original_value = prepped$data[[prepped$original_var]],
    transformed_value = prepped$var_vec,
    regime = regimes[match(names(prepped$var_vec), prepped$tree$tip.label)],
    stringsAsFactors = FALSE
  )
  
  # Create regime summary
  regime_summary <- data.frame(
    regime = names(theta),
    theta = theta,
    theta_original = if (prepped$analysis_var == paste0(prepped$original_var, "_logit")) {
      plogis(theta)
    } else {
      theta
    },
    stringsAsFactors = FALSE
  )
  
  # Return comprehensive results
  list(
    data = result_data,
    single_models = single_regime,
    parameters = data.frame(
      alpha = alpha,
      sigma2 = sigma2,
      stat_var = if (!is.na(alpha) && !is.na(sigma2)) sigma2/(2*alpha) else NA
    ),
    regimes = regime_summary,
    surface = surface_res
  )
}