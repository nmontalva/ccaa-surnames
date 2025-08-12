evolutionary_analysis <- function(data, variables, tree, steps = 1:6, verbose = TRUE) {
  # Initialize results list
  results <- list()
  
  # Ensure variables is a character vector
  if (!is.character(variables)) {
    variables <- as.character(variables)
  }
  
  # Process each variable sequentially
  for (var in variables) {
    if (verbose) message("\nProcessing variable: ", var)
    
    # Initialize variable-specific results
    var_results <- list()
    var_clean <- gsub("[^[:alnum:]]", "_", var)  # Clean variable name
    
    tryCatch({
      # Step 1: Prepare data
      if (1 %in% steps) {
        if (verbose) message("Step 1: Preparing trait data")
        var_results$prepped <- prep_traits(data, var, tree)
      }
      
      # Step 2: Fit single regime models
      if (2 %in% steps && !is.null(var_results$prepped)) {
        if (verbose) message("Step 2: Fitting single regime models")
        var_results$single_regime <- fit_single_regime(var_results$prepped)
      }
      
      # Step 3: Run SURFACE
      if (3 %in% steps && !is.null(var_results$prepped)) {
        if (verbose) message("Step 3: Running SURFACE analysis")
        var_results$surface <- run_surface(var_results$prepped, verbose)
      }
      
      # Step 4: Summarize results
      if (4 %in% steps && !is.null(var_results$prepped) && !is.null(var_results$surface)) {
        if (verbose) message("Step 4: Summarizing results")
        var_results$summary <- summarise_results(
          var_results$prepped, 
          var_results$surface, 
          var_results$single_regime
        )
      }
      # Step5: Compare models
      if (5 %in% steps && !is.null(var_results$summary) && !is.null(var_results$single_regime)) {
        if (verbose) message ("Step 5: Comparing models")
        var_results$comparison <- tryCatch({
          compare_models(var_results$summary)
        }, error = function(e){
          warning("Model comparison failed for ", var, ": ", e$message)
          NULL
        })
      }
      # Step 6: Generate plots
      if (6 %in% steps && !is.null(var_results$summary)) {
        if (verbose) message("Step 6: Generating plots")
        var_results$plots <- plot_result(var_results$summary)
      }
      
      # Store successful results
      results[[var_clean]] <- var_results
      
    }, error = function(e) {
      warning("Failed to process variable '", var, "': ", e$message)
      results[[var_clean]] <- list(error = e$message)
    })
  }
  
  # Add metadata about the analysis
  attr(results, "tree") <- deparse(substitute(tree))
  attr(results, "variables") <- variables
  attr(results, "steps_completed") <- steps
  
  return(results)
}
