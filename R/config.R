# Reproducibility settings shared by all notebooks and automated tests.
# Change this value deliberately: it changes every stochastic simulation result.
PROJECT_SEED <- 42L

set_project_seed <- function() {
  set.seed(PROJECT_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion")
  invisible(PROJECT_SEED)
}

