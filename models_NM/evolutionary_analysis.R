#' Wrapper function for evolutionary analysis
#' @param data Input dataframe
#' @param variable Trait variable to analyze
#' @param tree Phylogenetic tree
#' @param steps Which steps to run (1-5)
#' @param verbose Whether to show progress messages
#' @return Analysis results
evolutionary_analysis <- function(data, variable, tree, steps = 1:5, verbose = TRUE) {
  results <- list()
  
  # Step 1: Prepare data
  if (1 %in% steps) {
    if (verbose) message("Step 1: Preparing trait data")
    results$prepped <- prep_traits(data, variable, tree)
  }
  
  # Step 2: Fit single regime models
  if (2 %in% steps && !is.null(results$prepped)) {
    if (verbose) message("Step 2: Fitting single regime models")
    results$single_regime <- fit_single_regime(results$prepped)
  }
  
  # Step 3: Run SURFACE
  if (3 %in% steps && !is.null(results$prepped)) {
    if (verbose) message("Step 3: Running SURFACE analysis")
    results$surface <- run_surface(results$prepped, verbose)
  }
  
  # Step 4: Summarize results
  if (4 %in% steps && !is.null(results$prepped) && !is.null(results$surface)) {
    if (verbose) message("Step 4: Summarizing results")
    results$summary <- summarise_results(
      results$prepped, 
      results$surface, 
      results$single_regime
    )
  }
  
  # Step 5: Generate plots
  if (5 %in% steps && !is.null(results$summary)) {
    if (verbose) message("Step 5: Generating plots")
    results$plots <- plot_result(results$summary)
  }
  
  return(results)
}