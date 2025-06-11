# surface_analysis.R
# — Iteration 1: only prep_traits()

#’ Transforms 0–1 traits into adjusted (no exact 0/1) + logit
#’ @param data   data.frame with one column “community” and one or more trait columns
#’ @param traits character vector of column‐names in data to transform
#’ @return       data.frame with added <trait>_adj and <trait>_logit columns
prep_traits <- function(data, traits) {
  for (tr in traits) {
    if (! tr %in% names(data)) {
      stop("Trait '", tr, "' not in data frame")
    }
    n <- nrow(data)
    adj <- (data[[tr]] * (n - 1) + 0.5) / n
    data[[paste0(tr, "_adj")]]   <- adj
    data[[paste0(tr, "_logit")]] <- log(adj / (1 - adj))
  }
  data
}

# 2. Single-regime BM vs OU comparison
#’ @param trait_vec named numeric vector of logit-transformed trait values
#’ @param tree      phylogenetic tree (of class “phylo”)
#’ @return named vector of AICs: c(BM = <…>, OU = <…>)
fit_single_regime <- function(trait_vec, tree) {
  bm <- fitContinuous(tree, trait_vec, model = "BM")
  ou <- fitContinuous(tree, trait_vec, model = "OU")
  c(BM = bm$opt$aic, OU = ou$opt$aic)
}

# 3. SURFACE pipeline: forward then backward phase
#’ @param trait_vec named numeric vector of logit-transformed trait values
#’ @param tree      phylo object (must have unique node names)
#’ @return list with elements:
#’         $fwd a list of forward-phase SURFACE models,
#’         $bwd a list of backward-phase SURFACE models
run_surface <- function(trait_vec, tree) {
  # ensure unique node labels
  named_tree <- nameNodes(tree)
  # convert to OUCH data
  conv <- convertTreeData(named_tree,
                          data.frame(val = trait_vec, row.names = names(trait_vec)))
  otree <- conv[[1]]
  odat  <- conv[[2]]
  # forward
  fwd <- surfaceForward(otree, odat, aic_threshold = 0, verbose = TRUE)
  # backward
  bwd <- surfaceBackward(otree, odat, fwd[[length(fwd)]],
                         aic_threshold = 0, verbose = TRUE)
  list(fwd = fwd, bwd = bwd)
}

