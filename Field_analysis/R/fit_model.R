# Fit the joint likelihood with a single-start, bounded optimx/nlminb model.
# Calculate parameter estimates, Hessian intervals, and the output summary.
# Choose initial values for the dispersal optimization algorithm:
# Starting values and bounds are defined in config.R.
DispersalPars <- c(lambda, lambda_ver, barrier_strength)
DispersalPars_lower <- c(lambda_lower, lambda_ver_lower, barrier_strength_lower)
DispersalPars_upper <- c(lambda_upper, lambda_ver_upper, barrier_strength_upper)

DispersalPars_fit <- optimx(
  par = DispersalPars, fn = logLike_all, method = "nlminb",
  control = list(reltol = 1e-5, starttests = FALSE, kkt = FALSE),
  upper = DispersalPars_upper, lower = DispersalPars_lower
)

# The estimated values of the adult parameters given the kinship data are:
N_A_est <- node_population_size * nrow(distMat)
lambda_est <- DispersalPars_fit$p1
lambda_ver_est <- DispersalPars_fit$p2
barrier_strength_est <- DispersalPars_fit$p3
theta_est <- c(lambda_est, lambda_ver_est, barrier_strength_est)
hessian_result <- estimate_hessian_ci(
  theta_hat = theta_est,
  lower = DispersalPars_lower,
  upper = DispersalPars_upper,
  fn = logLike_all,
  par_names = c("lambda", "lambda_ver", "barrier_strength")
)
fit_convcode <- if ("convcode" %in% names(DispersalPars_fit)) DispersalPars_fit$convcode[[1]] else NA_integer_
fit_method <- if ("method" %in% names(DispersalPars_fit)) as.character(DispersalPars_fit$method[[1]]) else "nlminb"

fit_summary <- data.frame(
  grid_size = grid_size,
  rotation_angle = rotation_angle,
  node_population_size = node_population_size,
  patches_per_floor = patches_per_floor,
  num_layers = length(floors),
  num_patches_total = nrow(distMat),
  num_sampled_adults = numSampledAdults,
  num_sampled_adult_females = numSampledAdultFemales,
  t_sampling_start = as.character(tSamplingStart),
  t_sampling_end = as.character(tSamplingEnd),
  N_A_fixed = N_A_est,
  fit_method = fit_method,
  fit_convcode = fit_convcode,
  lambda_est = lambda_est,
  lambda_se = get_ci_value(hessian_result, "lambda", "se"),
  lambda_ci_lower = get_ci_value(hessian_result, "lambda", "ci_lower"),
  lambda_ci_upper = get_ci_value(hessian_result, "lambda", "ci_upper"),
  lambda_ci_lower_clipped = get_ci_value(hessian_result, "lambda", "ci_lower_clipped"),
  lambda_ci_upper_clipped = get_ci_value(hessian_result, "lambda", "ci_upper_clipped"),
  lambda_ver_est = lambda_ver_est,
  lambda_ver_se = get_ci_value(hessian_result, "lambda_ver", "se"),
  lambda_ver_ci_lower = get_ci_value(hessian_result, "lambda_ver", "ci_lower"),
  lambda_ver_ci_upper = get_ci_value(hessian_result, "lambda_ver", "ci_upper"),
  lambda_ver_ci_lower_clipped = get_ci_value(hessian_result, "lambda_ver", "ci_lower_clipped"),
  lambda_ver_ci_upper_clipped = get_ci_value(hessian_result, "lambda_ver", "ci_upper_clipped"),
  barrier_strength_est = barrier_strength_est,
  barrier_strength_se = get_ci_value(hessian_result, "barrier_strength", "se"),
  barrier_strength_ci_lower = get_ci_value(hessian_result, "barrier_strength", "ci_lower"),
  barrier_strength_ci_upper = get_ci_value(hessian_result, "barrier_strength", "ci_upper"),
  barrier_strength_ci_lower_clipped = get_ci_value(hessian_result, "barrier_strength", "ci_lower_clipped"),
  barrier_strength_ci_upper_clipped = get_ci_value(hessian_result, "barrier_strength", "ci_upper_clipped"),
  barrier_strength_boundary_flag = get_ci_value(hessian_result, "barrier_strength", "boundary_flag"),
  hessian_ok = get_ci_value(hessian_result, "barrier_strength", "hessian_ok"),
  hessian_condition_number = get_ci_value(hessian_result, "barrier_strength", "hessian_condition_number"),
  hessian_positive_definite = get_ci_value(hessian_result, "barrier_strength", "hessian_positive_definite"),
  ci_method = get_ci_value(hessian_result, "barrier_strength", "ci_method"),
  neg_loglike = as.numeric(DispersalPars_fit$value),
  elapsed_seconds = as.numeric(difftime(Sys.time(), time1, units = "secs"))
)

write.csv(
  fit_summary,
  file = fit_summary_file,
  row.names = FALSE
)

time2 <- Sys.time()
print("updated FL_barrier, fixed node population, lambda/lambda_ver/barrier estimated, block level")
print(time2 - time1)
print('N_A_fixed, lambda, lambda_ver, barrier_strength:')
print(N_A_est)
print(lambda_est)
print(lambda_ver_est)
print(barrier_strength_est)
print("saved fit summary:")
print(fit_summary_file)
