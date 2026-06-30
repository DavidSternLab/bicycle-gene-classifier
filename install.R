#!/usr/bin/env Rscript
# install.R — install bicycle_classifier's R dependencies.
#
# Usage:
#   Rscript install.R
#
# Installs (idempotent: skips packages already present):
#   CRAN: dplyr, tidyr, ggplot2, optparse
#   Bioconductor: rtracklayer
#
# Tested on R 4.2.3. Should work on R >= 4.0.

cran_pkgs <- c("dplyr", "tidyr", "ggplot2", "optparse")
bioc_pkgs <- c("rtracklayer")

cat("R version:", R.version.string, "\n\n")

# CRAN ---------------------------------------------------------------------
to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install,
                   repos = "https://cloud.r-project.org",
                   dependencies = TRUE)
} else {
  cat("All CRAN packages already installed.\n")
}

# Bioconductor -------------------------------------------------------------
bioc_to_install <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(bioc_to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  cat("Installing Bioconductor packages:", paste(bioc_to_install, collapse = ", "), "\n")
  BiocManager::install(bioc_to_install, ask = FALSE, update = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

# Verify ------------------------------------------------------------------
cat("\nVerifying imports:\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
for (p in all_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-12s %s\n", p, if (ok) "OK" else "FAIL"))
}

cat("\nDone. You can now run:\n")
cat("  bin/bicycle_classifier -g data/example.gff3 -d /tmp/bicycle_test\n")
