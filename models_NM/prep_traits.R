#' Prepare trait data for analysis
#' @param data Dataframe containing trait data
#' @param variable Name of the trait variable to analyze
#' @param tree Phylogenetic tree
#' @return List containing prepared data and metadata
prep_traits <- function(data, variable, tree) {
  # Match data to tree tips
  if (!"community" %in% colnames(data)) {
    stop("Data must contain 'community' column matching tree tips")
  }
  
  data <- data[data$community %in% tree$tip.label, ]
  if (nrow(data) == 0) stop("No matching tips between data and tree")
  
  # Apply logit transform if needed (for proportions 0-1)
  if (all(data[[variable]] >= 0 & data[[variable]] <= 1, na.rm = TRUE)) {
    n <- nrow(data)
    data[[paste0(variable, "_adj")]] <- (data[[variable]] * (n - 1) + 0.5) / n
    data[[paste0(variable, "_logit")]] <- log(data[[paste0(variable, "_adj")]] / 
                                         (1 - data[[paste0(variable, "_adj")]]))
    analysis_var <- paste0(variable, "_logit")
  } else {
    analysis_var <- variable
  }
  
  # Create named vector
  var_vec <- setNames(data[[analysis_var]], data$community)
  
  return(list(
    data = data,
    var_vec = var_vec,
    analysis_var = analysis_var,
    original_var = variable,
    tree = tree
  ))
}