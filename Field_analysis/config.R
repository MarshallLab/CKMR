# Configuration for the season-2 dispersal and barrier example.
# Paths are anchored by the runner, so no setwd() or machine-specific paths are needed.
sample_file <- file.path(project_dir, "data", "season2_filtered_v2.csv")
po_file <- file.path(project_dir, "data", "Kinship_PO_v3.csv")
fs_file <- file.path(project_dir, "data", "Kinship_FS_v3.csv")
output_dir <- file.path(project_dir, "output")
fit_summary_file <- file.path(output_dir, "barrier_hessian_pop5_ver2.csv")

# Select the study blocks; block 836 is excluded.
selected_blocks <- c(943:953, 835, 837:850)

# Use observed floors 2, 7, 12 as sampled layers 1, 3, 5.
# Layers 2 and 4 remain in the movement matrix as unsampled intermediate layers.
raw_sampled_floors <- c(2, 7, 12)
model_sampled_layers <- c(1, 3, 5)
floors <- 1:5
floor_height <- 7 # metres between model layers
node_population_size <- 5 # fixed adults per node (not an estimated parameter)
grid_size <- 45 # metres
rotation_angle <- -50 # degrees, applied by rotate_points()

# Parameter order: horizontal rate, vertical rate, barrier strength.
lambda <- 1 / 15.3
lambda_lower <- 0.001
lambda_upper <- 0.999
lambda_ver <- 1 / 5
lambda_ver_lower <- 0.001
lambda_ver_upper <- 0.999
barrier_strength <- 0.75
barrier_strength_lower <- 0.001
barrier_strength_upper <- 0.999

# Life-history constants stay in logLike_all() in R/likelihood.R.
# The barrier index remains 98 in BOTH likelihood functions. See docs/MODEL_NOTES.md
# before adapting this example to a different spatial layout or dataset.
