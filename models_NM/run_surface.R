#' Run SURFACE analysis
#' @param prepped Output from prep_traits
#' @param verbose Whether to show progress messages
#' @return SURFACE results object
run_surface <- function(prepped, verbose = TRUE) {
  require(surface)
  
  named_tree <- nameNodes(prepped$tree)
  dat_var <- data.frame(trait = prepped$var_vec, row.names = names(prepped$var_vec))
  colnames(dat_var) <- prepped$analysis_var
  
  conv <- convertTreeData(named_tree, dat_var)
  
  if (verbose) message("Running SURFACE forward phase")
  fwd <- surfaceForward(conv[[1]], conv[[2]], verbose = verbose)
  
  if (verbose) message("Running SURFACE backward phase")
  bwd <- surfaceBackward(conv[[1]], conv[[2]], fwd[[length(fwd)]], verbose = verbose)
  
  return(list(
    forward = fwd,
    backward = bwd,
    final_model = bwd[[length(bwd)]],
    named_tree = named_tree
  ))
}