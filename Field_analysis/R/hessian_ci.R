# Confidence intervals for bounded parameters.
# Evaluate the Hessian on the bounded-logit scale, invert it, then apply the
# delta method on the parameter scale. Both raw and clipped Wald limits are kept.
# Diagnostic flags are reported separately and do not suppress the intervals.
bounded_logit <- function(theta, lower, upper, eps = 1e-6) {
  scaled <- (theta - lower) / (upper - lower)
  scaled <- pmin(pmax(scaled, eps), 1 - eps)
  qlogis(scaled)
}

bounded_inv_logit <- function(z, lower, upper) {
  lower + (upper - lower) * plogis(z)
}

estimate_hessian_ci <- function(theta_hat, lower, upper, fn, par_names, alpha = 0.05) {
  theta_hat <- as.numeric(theta_hat)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  z_hat <- bounded_logit(theta_hat, lower, upper)

  negloglike_z <- function(z) {
    theta <- bounded_inv_logit(z, lower, upper)
    value <- fn(theta)
    if (!is.finite(value)) return(.Machine$double.xmax / 100)
    as.numeric(value)
  }

  hessian_z <- tryCatch(
    stats::optimHess(par = z_hat, fn = negloglike_z),
    error = function(e) e
  )

  result <- data.frame(
    parameter = par_names,
    estimate = theta_hat,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    ci_lower_clipped = NA_real_,
    ci_upper_clipped = NA_real_,
    boundary_flag = theta_hat <= lower + 1e-5 | theta_hat >= upper - 1e-5,
    hessian_ok = FALSE,
    hessian_condition_number = NA_real_,
    hessian_positive_definite = FALSE,
    ci_method = "wald_delta_from_logit_hessian",
    stringsAsFactors = FALSE
  )

  if (inherits(hessian_z, "error")) {
    return(list(summary = result))
  }

  hessian_z <- as.matrix(hessian_z)
  eig <- tryCatch(eigen(hessian_z, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  positive_definite <- all(is.finite(eig)) && min(eig) > 0
  condition_number <- tryCatch(kappa(hessian_z, exact = TRUE), error = function(e) NA_real_)
  cov_z <- tryCatch(solve(hessian_z), error = function(e) NULL)

  if (is.null(cov_z)) {
    result$hessian_condition_number <- condition_number
    result$hessian_positive_definite <- positive_definite
    return(list(summary = result))
  }

  scaled_hat <- plogis(z_hat)
  jacobian_diag <- (upper - lower) * scaled_hat * (1 - scaled_hat)
  jacobian <- diag(jacobian_diag, nrow = length(theta_hat))
  cov_theta <- jacobian %*% cov_z %*% jacobian

  se <- sqrt(pmax(diag(cov_theta), 0))
  z_crit <- stats::qnorm(1 - alpha / 2)
  ci_lower <- theta_hat - z_crit * se
  ci_upper <- theta_hat + z_crit * se

  result$se <- se
  result$ci_lower <- ci_lower
  result$ci_upper <- ci_upper
  result$ci_lower_clipped <- pmax(ci_lower, lower)
  result$ci_upper_clipped <- pmin(ci_upper, upper)
  result$hessian_ok <- all(is.finite(se)) && all(is.finite(cov_theta))
  result$hessian_condition_number <- condition_number
  result$hessian_positive_definite <- positive_definite

  list(summary = result)
}

get_ci_value <- function(hessian_result, parameter, column) {
  row <- hessian_result$summary[hessian_result$summary$parameter == parameter, , drop = FALSE]
  if (nrow(row) != 1) return(NA_real_)
  row[[column]]
}
