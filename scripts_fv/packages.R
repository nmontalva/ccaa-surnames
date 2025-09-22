# =====================================
# Install missing R packages safely
# =====================================

# ----- 1. Normal CRAN packages -----

packages <- c(
  "Biodem",
  "GGally",
  "Hmisc",
  "NbClust",
  "PlotTools",
  "Publish",
  "RColorBrewer",
  "REAT",
  "SDPDmod",
  "TreeDist",
  "ade4",
  "adegenet",
  "ape",
  "bayesplot",
  "brms",
  "caper",
  "car",
  "cluster",
  "conflicted",
  "corrplot",
  "cowplot",
  "dendextend",
  "dplyr",
  "e1071",
  "factoextra",
  "fpc",
  "fs",
  "geiger",
  "geomorph",
  "geosphere",
  "ggdendro",
  "ggplot2",
  "ggpubr",
  "googledrive",
  "graph4lg",
  "grDevices",
  "grid",
  "gridExtra",
  "igraph",
  "magick",
  "nlme",
  "phylogram",
  "phytools",
  "picante",
  "readr",
  "reldist",
  "rr2",
  "sf",
  "stargazer",
  "stringr",
  "surface",
  "tidyr",
  "tidyverse",
  "treestats",
  "vegan",
  "wesanderson"
)

# Excluded manually:
# - "BiocManager" → Installer utility.

# Log files
packages_log_dir <- "outputs/logs/packages_logs"
dir.create(packages_log_dir, showWarnings = FALSE)

log_skipped <- file.path(packages_log_dir, "skipped_packages.log")
log_installed <- file.path(packages_log_dir, "installed_packages.log")
log_failed <- file.path(packages_log_dir, "failed_packages.log")
log_errors <- file.path(packages_log_dir, "error_messages.log")

# Safely clear old logs
safe_remove <- function(file) if (file.exists(file)) file.remove(file)
invisible(lapply(c(log_skipped, log_installed, log_failed, log_errors), safe_remove))

# Installed packages snapshot
installed <- rownames(installed.packages())

# Helper function to install missing packages only
install_if_needed <- function(pkg) {
  if (pkg %in% installed) {
    write(pkg, file = log_skipped, append = TRUE)
  } else {
    message("Installing: ", pkg)
    tryCatch({
      install.packages(pkg) # Force binary where available
      installed <<- rownames(installed.packages()) # Refresh full list
      if (pkg %in% installed) {
        write(pkg, file = log_installed, append = TRUE)
      } else {
        write(pkg, file = log_failed, append = TRUE)
      }
    }, error = function(e) {
      write(pkg, file = log_failed, append = TRUE)
      write(paste(pkg, ":", e$message), file = log_errors, append = TRUE)
    })
  }
}

# Install missing CRAN packages
invisible(lapply(packages, install_if_needed))

# Test if packages load
invisible(lapply(packages, function(pkg) {
  message("Trying to load: ", pkg)
  tryCatch({
    library(pkg, character.only = TRUE)
  }, error = function(e) {
    message("Failed to load ", pkg, ": ", e$message)
  })
}))

# ---- Install treeio separately (requires BiocManager) ----

# Refresh installed packages snapshot
installed <- rownames(installed.packages())

if (!"treeio" %in% installed) {
  message("Installing treeio via BiocManager...")
  if (!"BiocManager" %in% installed) {
    install.packages("BiocManager")
  }
  BiocManager::install("treeio")
}

message("Trying to load: treeio")
tryCatch({
  library(treeio)
}, error = function(e) {
  message("Failed to load treeio: ", e$message)
})

# ---- Install ggtree separately (requires BiocManager) ----

# Refresh installed packages snapshot
installed <- rownames(installed.packages())

if (!"treeio" %in% installed) {
  message("Installing treeio via BiocManager...")
  if (!"BiocManager" %in% installed) {
    install.packages("BiocManager")
  }
  BiocManager::install("ggtree")
}

message("Trying to load: ggtree")
tryCatch({
  library(ggtree)
}, error = function(e) {
  message("Failed to load ggtree: ", e$message)
})

#-------------- CiteR---------------------
citeR <- function(packages) {
  Rbibs <- ""
  
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message("Skipping missing package: ", package)
      next
    }
    
    Rbib <- capture.output(print(citation(package), bibtex = TRUE))    
    Rbib <- mapply(function(x, y) Rbib[x:y], 
                   grep("  @.+[{]", Rbib), 
                   which(Rbib == "  }"))
    
    if (is.matrix(Rbib)) {
      Rbib[1, 1] <- gsub(",", paste0(package, ","), Rbib[1, 1])
      Rbib <- paste0(Rbib, collapse = "\n")
    } else {
      Rbib <- unlist(lapply(Rbib, function(x) {
        x[1] <- gsub(",", paste0(package, ","), x[1])
        paste0(unlist(x), collapse = "\n")
      }))
    }
    
    if (length(Rbib) > 1) {
      if (any(grepl("@Manual", Rbib))) {
        Rbib <- Rbib[grep("@Manual", Rbib)][1]
      } else {
        Rbib <- Rbib[1]
      }
    }
    
    Rbibs <- paste(Rbibs, Rbib, sep = "\n\n")
  }
  
  writeBin(charToRaw(utf8::as_utf8(Rbibs)), "packages.bib")
}

packages <- unlist(lapply(as.list(match.call()), deparse))[-1]
citeR(packages)
