run_evolutionary_cl <- function(data, variable, tree, logit_transform = TRUE, 
                                      aic_threshold = 0, verbose = FALSE) {
  # Required packages
  required_packages <- c("geiger", "surface", "dplyr", "ggplot2", "phytools", "ggtree")
  
  # Check and load required packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is needed but not installed. Please install it with install.packages('", pkg, "')"))
    }
    library(pkg, character.only = TRUE)
  }
  
  # Data validation
  if (!variable %in% colnames(data)) {
    stop(paste("Variable", variable, "not found in the dataset"))
  }
  
  if (!"community" %in% colnames(data)) {
    stop("Dataset must contain a 'community' column for identifying tree tips")
  }
  
  # Extract the variable data
  var_data <- data[[variable]]
  communities <- data$community
  
  # Step 1: Transform the variable if required
  if (logit_transform) {
    n <- nrow(data)
    eps <- 0.5 / n
    var_adj <- (var_data * (n - 1) + 0.5) / n
    var_logit <- log(var_adj / (1 - var_adj))
    
    # Add transformed variables to the data
    data[[paste0(variable, "_adj")]] <- var_adj
    data[[paste0(variable, "_logit")]] <- var_logit
    
    # Use logit transformed variable for analysis
    analysis_var <- var_logit
  } else {
    # Use original variable
    analysis_var <- var_data
  }
  
  # Create named vector
  var_vec <- setNames(analysis_var, communities)
  
  # Step 2: Fit BM and OU models
  bm_fit <- geiger::fitContinuous(tree, var_vec, model = "BM")
  ou_fit <- geiger::fitContinuous(tree, var_vec, model = "OU")
  
  # Compare AICc values
  bm_aicc <- bm_fit$opt$aicc
  ou_aicc <- ou_fit$opt$aicc
  
  # Print AICc comparison
  cat("Model comparison:\n")
  cat("BM AICc:", bm_aicc, "\n")
  cat("OU AICc:", ou_aicc, "\n")
  cat("AIC difference (BM - OU):", bm_aicc - ou_aicc, "\n\n")
  
  # Step 3: SURFACE analysis
  # Ensure unique node names
  named_tree <- surface::nameNodes(tree)
  
  # Prepare data frame for SURFACE
  dat_var <- data.frame(value = analysis_var, row.names = communities)
  colnames(dat_var)[1] <- variable # Rename column to match variable name
  
  # Convert to OUCH objects
  conv_result <- surface::convertTreeData(named_tree, dat_var)
  otree <- conv_result[[1]]
  odata <- conv_result[[2]]
  
  # Run SURFACE forward phase
  cat("Running SURFACE forward phase...\n")
  fwd_result <- surface::surfaceForward(
    otree, odata,
    aic_threshold = aic_threshold,
    verbose = verbose
  )
  
  # Run SURFACE backward phase
  cat("Running SURFACE backward phase...\n")
  bwd_result <- surface::surfaceBackward(
    otree, 
    odata, 
    fwd_result[[length(fwd_result)]], 
    aic_threshold = aic_threshold, 
    verbose = verbose
  )
  
  # Get SURFACE summary
  surf_summary <- surface::surfaceSummary(bwd_result)
  
  # Parameter table
  param_tab <- data.frame(
    alpha = surf_summary$alpha[variable],
    sigma2 = surf_summary$sigma_squared[variable]
  )
  param_tab$stat_var <- param_tab$sigma2 / (2 * param_tab$alpha)
  rownames(param_tab) <- variable
  
  cat("\nModel parameters", if(logit_transform) "(logit scale)" else "", ":\n")
  print(param_tab)
  
  # Regime optima table
  theta_values <- surf_summary$theta[, 1]
  theta_tab <- data.frame(
    regime = names(theta_values),
    theta = theta_values
  )
  
  # If logit transformed, add back-transformed values
  if (logit_transform) {
    theta_tab[[paste0("theta_", variable)]] <- round(plogis(theta_values), 4)
  }
  
  cat("\nRegime optima:\n")
  print(theta_tab)
  
  # Tip assignments
  hmod <- bwd_result[[length(bwd_result)]]$fit[[variable]]
  regs_fac <- hmod@regimes$regs
  tip_regs <- as.character(regs_fac[1:length(communities)])
  names(tip_regs) <- communities
  
  # Build tip data frame
  tip_df <- data.frame(
    community = names(tip_regs),
    regime = tip_regs,
    stringsAsFactors = FALSE
  )
  
  # Summary of regimes
  reg_summary <- tip_df %>%
    dplyr::count(regime) %>%
    dplyr::mutate(pct = round(100 * n / sum(n), 1))
  
  cat("\nRegime distribution:\n")
  print(reg_summary)
  
  # Create color palette
  pal <- setNames(rainbow(nrow(theta_tab)), theta_tab$regime)
  
  # Plots
  plots <- list()
  
  # 1. Bar chart
  bar_plot <- ggplot(reg_summary, aes(x = regime, y = n, fill = regime)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = pal) +
    labs(
      title = paste("Number of Communities per Regime (", variable, ")", sep = ""),
      x = "Regime",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  plots$bar_plot <- bar_plot
  
  # 2. Tree plot with ggtree
  # Prepare theta labels
  if (logit_transform) {
    theta_labels <- setNames(
      paste0(theta_tab$regime, ": θ=", format(round(theta_tab[[paste0("theta_", variable)]], 3), nsmall=3)),
      theta_tab$regime
    )
  } else {
    theta_labels <- setNames(
      paste0(theta_tab$regime, ": θ=", format(round(theta_tab$theta, 3), nsmall=3)),
      theta_tab$regime
    )
  }
  
  # Clean palette (remove NA entries)
  keep <- !is.na(names(pal))
  pal_clean <- pal[keep]
  theta_labels_clean <- theta_labels[keep]
  
  # Tree plot with ggtree
  tree_plot <- ggtree(named_tree, size = 0.8) %<+% tip_df +
    geom_tree(aes(color = regime), size = 0.8) +
    scale_color_manual(
      values = pal_clean,
      labels = theta_labels_clean,
      name = "Regime (θ)",
      na.translate = FALSE
    ) +
    theme_tree2() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.margin = margin(5, 100, 5, 5)  # Add margin for legend
    ) +
    ggtitle(paste("SURFACE Regimes on Tree (", variable, ")", sep = ""))
  
  plots$tree_plot <- tree_plot
  
  # Return results
  results <- list(
    data = data,
    variable = variable,
    logit_transform = logit_transform,
    bm_fit = bm_fit,
    ou_fit = ou_fit,
    surface_summary = surf_summary,
    parameter_table = param_tab,
    regime_table = theta_tab,
    tip_regimes = tip_df,
    regime_summary = reg_summary,
    plots = plots,
    bwd_result = bwd_result
  )
  
  return(results)
}

# Example usage:
 result <- run_evolutionary_cl(GM_df, "M", consensus_tree)
# 
# # Access results
 result$parameter_table
 result$regime_table
# 
# # Show plots
 print(result$plots$bar_plot)
 print(result$plots$tree_plot)

# Function to run analysis for multiple variables
run_multi_var_analysis <- function(data, variables, tree, ...) {
  results <- list()
  
  for (var in variables) {
    cat("\n=====================================\n")
    cat("Analyzing variable:", var, "\n")
    cat("=====================================\n\n")
    
    results[[var]] <- run_evolutionary_analysis(data, var, tree, ...)
  }
  
  return(results)
}

# Example usage for multiple variables:
# results <- run_multi_var_analysis(GM_df, c("M", "V", "R"), y_total)
#
# # Access results for a specific variable
# results$M$regime_table
# print(results$M$plots$tree_plot)