# Entry point: Rscript run_example.R [--prepare-only]
# Also supports: result <- source("run_example.R")$value

# Resolve this file before entering the function, supporting both Rscript and source().
.source_file <- if (sys.nframe() > 0L) sys.frame(1)$ofile else NULL
.is_sourced <- !is.null(.source_file)
if (is.null(.source_file)) {
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.file_arg) != 1L) stop("Run with Rscript or source() this file.")
  .source_file <- sub("^--file=", "", .file_arg)
}
.project_dir <- dirname(normalizePath(.source_file, mustWork = TRUE))

run_ckmr_example <- function(project_dir, prepare_only = FALSE) {
  # Package attachment order determines unqualified function lookup.
  required <- c("ggplot2", "sf", "dplyr", "igraph", "tidyr", "optimx", "stringr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install missing R packages first: ", paste(missing, collapse = ", "))
  }
  for (package in required) {
    suppressPackageStartupMessages(library(package, character.only = TRUE))
  }

  # Likelihood functions and prepared sampling objects share this environment.
  analysis <- new.env(parent = globalenv())
  analysis$project_dir <- normalizePath(project_dir, mustWork = TRUE)
  analysis$time1 <- Sys.time()
  stages <- c(
    "config.R", "R/helpers.R", "R/hessian_ci.R", "R/likelihood.R",
    "R/prepare_data.R"
  )
  for (stage in stages) {
    sys.source(file.path(analysis$project_dir, stage), envir = analysis)
    if (stage == "config.R") {
      inputs <- c(analysis$sample_file, analysis$po_file, analysis$fs_file)
      if (any(!file.exists(inputs))) {
        stop("Missing input file(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "))
      }
    }
  }
  message(sprintf(
    "Prepared %d adults, %d sampled patches, %d total nodes; %d PO and %d FS pairs.",
    analysis$numSampledAdults, analysis$numSampledPatches,
    nrow(analysis$distMat), nrow(analysis$kinship_PO_barrier),
    nrow(analysis$kinship_FS_barrier)
  ))
  if (prepare_only) {
    message("Preparation complete. Optimization and Hessian calculations were skipped.")
    return(invisible(analysis))
  }

  dir.create(analysis$output_dir, recursive = TRUE, showWarnings = FALSE)
  sys.source(file.path(analysis$project_dir, "R/fit_model.R"), envir = analysis)
  writeLines(capture.output(sessionInfo()), file.path(analysis$output_dir, "sessionInfo.txt"))
  invisible(analysis)
}

.args <- commandArgs(trailingOnly = TRUE)
# When sourced from another R script, ignore that script's command-line arguments.
if (!.is_sourced) {
  if (any(!.args %in% "--prepare-only")) stop("Usage: Rscript run_example.R [--prepare-only]")
} else {
  .args <- character()
}
run_ckmr_example(.project_dir, prepare_only = "--prepare-only" %in% .args)
