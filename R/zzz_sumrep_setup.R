
#' Install and setup sumrep (official workflow)
#'
#' This function follows the official sumrep installation instructions:
#' 1) Install CRAN dependencies
#' 2) Install Bioconductor dependencies
#' 3) Clone sumrep from GitHub (if needed)
#' 4) Load or install sumrep
#'
#' @param sumrep_dir Local directory to clone sumrep into.
#' @param github_repo sumrep GitHub repo.
#' @param quiet Logical.
#'
#' @export
#'
shield_setup_sumrep <- function(
    sumrep_dir   = file.path(tools::R_user_dir("ShieldAIRR", "data"), "sumrep"),
    github_repo  = "https://github.com/matsengrp/sumrep.git",
    quiet        = FALSE
) {
  message("[ShieldAIRR] Checking sumrep availability...")
  # If already available, done
  if (requireNamespace("sumrep", quietly = TRUE)) {
    message("[ShieldAIRR] sumrep already installed.")
    return(invisible(TRUE))
  }
  # ------------------------------------------------------------------
  # 1. Install CRAN dependencies (official list)
  # ------------------------------------------------------------------
  cran_pkgs <- c(
    "alakazam", "ape", "CollessLike", "data.table", "dplyr",
    "entropy", "jsonlite", "magrittr", "Peptides",
    "RecordLinkage", "shazam", "seqinr", "stringdist",
    "stringr", "testthat", "textmineR", "yaml"
  )
  message("[ShieldAIRR] Installing CRAN dependencies for sumrep...")
  missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran, repos = "https://cran.rstudio.com")
  }
  # ------------------------------------------------------------------
  # 2. Install Bioconductor dependencies
  # ------------------------------------------------------------------
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.rstudio.com")
  }
  bioc_pkgs <- c("Biostrings", "GenomicAlignments")
  message("[ShieldAIRR] Installing Bioconductor dependencies for sumrep...")
  missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
  # ------------------------------------------------------------------
  # 3. Clone sumrep if needed
  # ------------------------------------------------------------------
  if (!dir.exists(sumrep_dir)) {
    message("[ShieldAIRR] Cloning sumrep from GitHub...")
    if (!requireNamespace("git2r", quietly = TRUE)) {
      install.packages("git2r", repos = "https://cran.rstudio.com")
    }
    git2r::clone(github_repo, sumrep_dir)
  } else {
    message("[ShieldAIRR] Using existing sumrep at: ", sumrep_dir)
  }
  # ------------------------------------------------------------------
  # 4. Load sumrep
  # ------------------------------------------------------------------
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.rstudio.com")
  }
  message("[ShieldAIRR] Loading sumrep via devtools::load_all()")
  devtools::load_all(sumrep_dir, quiet = quiet)
  message("[ShieldAIRR] sumrep setup complete.")
  invisible(TRUE)
}










