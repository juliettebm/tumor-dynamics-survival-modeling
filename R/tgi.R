# Reusable tumor-growth-inhibition model helpers.

fit_tgi_patient <- function(df) {
  required <- c("time", "SLD_obs")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  tryCatch({
    fit <- minpack.lm::nlsLM(
      SLD_obs ~ SLD0_est * (exp(-g_est * time) + exp(d_est * time) - 1),
      data = df,
      start = list(SLD0_est = df$SLD_obs[df$time == 0][1], g_est = 0.02, d_est = 0.01),
      lower = c(SLD0_est = 1, g_est = 0.0001, d_est = 0.0001),
      upper = c(SLD0_est = 500, g_est = 1, d_est = 1),
      control = minpack.lm::nls.lm.control(maxiter = 200)
    )
    coefs <- stats::coef(fit)
    tibble::tibble(g_hat = coefs["g_est"], d_hat = coefs["d_est"], convergence = TRUE)
  }, error = function(e) {
    tibble::tibble(g_hat = NA_real_, d_hat = NA_real_, convergence = FALSE)
  })
}
