run_evolutionary_analysis <- function(data, variable, tree, verbose = TRUE) {
  # Load required packages
  require(geiger)
  require(surface)
  require(ggplot2)
  require(ggtree)
  
  # Start timing
  start_time <- Sys.time()
  if (verbose) message("Starting analysis of ", variable, " on tree with ", length(tree$tip.label), " tips")
  
  # 1. Match data to tree tips
  if (!"community" %in% colnames(data)) stop("Data must contain 'community' column")
  data <- data[data$community %in% tree$tip.label, ]
  if (nrow(data) == 0) stop("No matching tips between data and tree")
  
  # 2. Apply logit transform if variable is proportion (0-1)
  if (all(data[[variable]] >= 0 & data[[variable]] <= 1, na.rm = TRUE)) {
    if (verbose) message("Applying logit transform to ", variable)
    n <- nrow(data)
    data[[paste0(variable, "_adj")]] <- (data[[variable]] * (n - 1) + 0.5) / n
    data[[paste0(variable, "_logit")]] <- log(data[[paste0(variable, "_adj")]] / (1 - data[[paste0(variable, "_adj")]]))
    analysis_var <- paste0(variable, "_logit")
  } else {
    analysis_var <- variable
  }
  
  # 3. Create named vector
  var_vec <- setNames(data[[analysis_var]], data$community)
  
  # 4. Fit BM and OU models
  if (verbose) message("Fitting BM and OU models")
  bm <- try(fitContinuous(tree, var_vec, model = "BM"))
  ou <- try(fitContinuous(tree, var_vec, model = "OU"))
  
  # 5. Prepare SURFACE analysis
  named_tree <- nameNodes(tree)
  dat_var <- data.frame(trait = var_vec, row.names = names(var_vec))
  colnames(dat_var) <- analysis_var
  
  conv <- convertTreeData(named_tree, dat_var)
  otree <- conv[[1]]
  odata <- conv[[2]]
  
  # 6. Run SURFACE
  if (verbose) message("Running SURFACE forward phase")
  fwd <- surfaceForward(otree, odata, verbose = verbose)
  
  if (verbose) message("Running SURFACE backward phase")
  bwd <- surfaceBackward(otree, odata, fwd[[length(fwd)]], verbose = verbose)
  
  # 7. Extract results (with robust error handling)
  final_model <- bwd[[length(bwd)]]
  
  # Initialize default values for single regime case
  regimes <- rep("single_regime", length(otree@nodelabels))
  theta <- setNames(mean(var_vec, na.rm = TRUE), "single_regime")
  alpha <- sigma2 <- NA
  
  # Check if regimes were found
  if (length(final_model$fit) > 0 && analysis_var %in% names(final_model$fit)) {
    hmod <- final_model$fit[[analysis_var]]
    if (inherits(hmod, "hansentree")) {
      regimes <- hmod@regimes$regs
      theta <- tryCatch({
        if (is.matrix(hmod@theta)) hmod@theta[,1] else hmod@theta
      }, error = function(e) {
        setNames(rep(mean(var_vec, na.rm = TRUE), length(unique(regimes))), 
                 paste0("regime_", seq_along(unique(regimes))))
      })
      alpha <- tryCatch(hmod@alpha, error = function(e) NA)
      sigma2 <- tryCatch(hmod@sigma.squared, error = function(e) NA)
    }
  }
  
  # Get tip regimes
  tip_regs <- as.character(regimes[1:length(var_vec)])
  names(tip_regs) <- names(var_vec)
  
  # 8. Compile results
  results <- list(
    data = {
      df <- data.frame(
        community = names(var_vec),
        original_value = data[[variable]],
        transformed_value = var_vec,
        regime = tip_regs,
        stringsAsFactors = FALSE
      )
      if (analysis_var != variable) {
        df$adjusted_value <- data[[paste0(variable, "_adj")]]
      }
      df
    },
    models = list(
      BM = bm,
      OU = ou,
      AICc = data.frame(
        Model = c("BM", "OU"),
        AICc = c(bm$opt$aicc, ou$opt$aicc),
        stringsAsFactors = FALSE
      )
    ),
    parameters = data.frame(
      alpha = alpha,
      sigma2 = sigma2,
      stat_var = if (!is.na(alpha) && !is.na(sigma2)) sigma2/(2*alpha) else NA,
      row.names = NULL
    ),
    regimes = data.frame(
      regime = names(theta),
      theta = theta,
      theta_original = if (analysis_var == paste0(variable, "_logit")) plogis(theta) else theta,
      stringsAsFactors = FALSE
    ),
    surface = list(
      forward = fwd,
      backward = bwd,
      final_model = final_model
    ),
    timing = Sys.time() - start_time
  )
  
  # 9. Create plots (only for trees with <= 50 tips)
  if (length(tree$tip.label) <= 50) {
    if (verbose) message("Creating visualizations")
    
    # Color palette
    n_regimes <- length(unique(results$regimes$regime))
    pal <- if (n_regimes == 1) {
      setNames("#1f77b4", results$regimes$regime[1])
    } else {
      setNames(rainbow(n_regimes), results$regimes$regime)
    }
    
    # AICc plot
    results$plots$aicc_plot <- ggplot(results$models$AICc, aes(x = Model, y = AICc, fill = Model)) +
      geom_col() +
      geom_text(aes(label = round(AICc, 1)), vjust = -0.5) +
      labs(title = paste("AICc Comparison for", variable)) +
      theme_minimal()
    
    # Tree plot
    tip_data <- data.frame(
      label = results$data$community,
      regime = results$data$regime,
      stringsAsFactors = FALSE
    )
    
    theta_labels <- paste0(
      results$regimes$regime, 
      ": θ=", 
      if (analysis_var == paste0(variable, "_logit")) {
        round(results$regimes$theta_original, 3)
      } else {
        round(results$regimes$theta, 3)
      }
    )
    
    results$plots$tree_plot <- ggtree(named_tree) %<+% tip_data +
      geom_tree(aes(color = regime), size = 0.8) +
      scale_color_manual(
        values = pal,
        labels = theta_labels,
        name = "Regime (θ)"
      ) +
      theme_tree2() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 12, face = "bold")
      ) +
      ggtitle(paste("SURFACE Regimes for", variable))
  }
  
  if (verbose) {
    message("Analysis completed successfully")
    message("Total time: ", round(results$timing, 1), " ", attr(results$timing, "units"))
  }
  
  return(results)
}