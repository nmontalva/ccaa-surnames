# Version 3: Install and log required R packages with explicit success check

# Define packages
packages <- c(
  "dplyr", "ade4", "tidyr", "adegenet", "readr", "tidyverse", "Biodem", "reldist", "ape",
  "ggdendro", "ggpubr", "phylogram", "phytools", "dendextend", "geosphere", "vegan", "e1071",
  "Hmisc", "REAT", "ggplot2", "gridExtra", "stringr", "conflicted", "graph4lg", "TreeDist",
  "corrplot", "geiger", "car", "caper", "nlme", "Publish", "treeio", "geomorph", "devtools",
  "cowplot", "GGally", "png", "grid", "factoextra", "cluster", "NbClust", "fpc", "PlotTools",
  "wesanderson"
)

# Log files
log_skipped <- "skipped_packages.log"
log_installed <- "installed_packages.log"
log_failed <- "failed_packages.log"
log_errors <- "error_messages.log"

# Clear old logs
file.remove(log_skipped, log_installed, log_failed, log_errors)

# Helper function
install_if_needed <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    write(pkg, file = log_skipped, append = TRUE)
  } else {
    message("Installing: ", pkg)
    tryCatch({
      install.packages(pkg)
      if (requireNamespace(pkg, quietly = TRUE)) {
        write(pkg, file = log_installed, append = TRUE)
      } else {
        write(pkg, file = log_failed, append = TRUE)
      }
    }, error = function(e) {
      write(pkg, file = log_failed, append = TRUE)
      write(paste(pkg, e$message), file = log_errors, append = TRUE)
    })
  }
}

# Install all
invisible(lapply(packages, install_if_needed))

# Test if packages are loading
for (pkg in packages) {
  message("Trying to load: ", pkg)
  success <- tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) {
    message("❌ Failed to load ", pkg, ": ", e$message)
    FALSE
  })
}
