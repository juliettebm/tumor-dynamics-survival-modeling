# Uncertainty helpers for continuous association metrics.

correlation_ci <- function(correlation, n, confidence = 0.95) {
  if (n <= 3L) stop("n must be greater than 3", call. = FALSE)
  if (!is.finite(correlation) || abs(correlation) >= 1) {
    stop("correlation must be finite and strictly between -1 and 1", call. = FALSE)
  }
  if (confidence <= 0 || confidence >= 1) {
    stop("confidence must be between 0 and 1", call. = FALSE)
  }

  z <- atanh(correlation)
  critical <- stats::qnorm(1 - (1 - confidence) / 2)
  interval <- tanh(z + c(-1, 1) * critical / sqrt(n - 3))
  names(interval) <- c("lower", "upper")
  interval
}
