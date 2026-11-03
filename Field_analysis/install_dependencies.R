# One-time dependency installation; the analysis itself never installs packages.
required <- c("ggplot2", "sf", "dplyr", "igraph", "tidyr", "optimx", "stringr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
