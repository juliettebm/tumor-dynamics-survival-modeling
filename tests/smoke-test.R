required_packages <- c(
  "rmarkdown", "knitr", "minpack.lm", "survival", "survminer",
  "JM", "nlme", "dplyr", "ggplot2"
)

source(file.path("R", "config.R"))

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

source_notebook <- normalizePath(
  file.path("notebooks", "01_simulation_and_joint_model.Rmd"),
  mustWork = TRUE
)
project_root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(TUMOR_PROJECT_ROOT = project_root)
on.exit(Sys.unsetenv("TUMOR_PROJECT_ROOT"), add = TRUE)
work_dir <- tempfile("tumor-dynamics-smoke-")
dir.create(work_dir)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
notebook <- file.path(work_dir, "01_simulation_and_joint_model.Rmd")
if (!file.copy(source_notebook, notebook, overwrite = TRUE)) {
  stop("Could not copy notebook into the isolated working directory", call. = FALSE)
}

message("Rendering notebook 01 in ", work_dir)
render_environment <- new.env(parent = globalenv())
rendered_file <- rmarkdown::render(
  input = notebook,
  output_file = "01_simulation_and_joint_model.html",
  output_dir = work_dir,
  intermediates_dir = work_dir,
  knit_root_dir = work_dir,
  quiet = TRUE,
  envir = render_environment
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert(file.exists(rendered_file), "Notebook render did not produce an HTML file")
assert(file.info(rendered_file)$size > 10000, "Rendered HTML is unexpectedly small")

patients_path <- file.path(work_dir, "patients_final.rds")
longitudinal_path <- file.path(work_dir, "longitudinal_data_final.rds")
assert(file.exists(patients_path), "patients_final.rds was not generated")
assert(file.exists(longitudinal_path), "longitudinal_data_final.rds was not generated")

patients <- readRDS(patients_path)
longitudinal <- readRDS(longitudinal_path)

assert(identical(PROJECT_SEED, 42L), "The documented project seed changed unexpectedly")

patient_columns <- c("patient_id", "g", "d", "SLD0", "os_time", "os_event")
longitudinal_columns <- c("patient_id", "time", "SLD_true", "SLD_obs")

assert(nrow(patients) == 200L, "Expected exactly 200 simulated patients")
assert(all(patient_columns %in% names(patients)), "Patient output schema changed")
assert(all(longitudinal_columns %in% names(longitudinal)), "Longitudinal output schema changed")
assert(length(unique(longitudinal$patient_id)) == 200L, "Longitudinal data do not cover all patients")
assert(nrow(longitudinal) >= 4000L, "Too few longitudinal observations were generated")
assert(all(is.finite(patients$g) & patients$g > 0), "Invalid tumor regression rates")
assert(all(is.finite(patients$d) & patients$d > 0), "Invalid tumor regrowth rates")
assert(all(patients$os_event %in% c(0L, 1L)), "Survival event indicator is not binary")
assert(all(is.finite(longitudinal$SLD_obs) & longitudinal$SLD_obs > 0), "Invalid observed tumor sizes")

# Deterministic reference checks: broad enough to tolerate harmless numerical
# differences between platforms, narrow enough to detect model or seed drift.
assert(exists("analysis_data", envir = render_environment), "Two-stage analysis output is missing")
assert(exists("cox_model", envir = render_environment), "Cox model output is missing")
analysis_data <- get("analysis_data", envir = render_environment)
cox_model <- get("cox_model", envir = render_environment)
cor_g <- cor(analysis_data$g, analysis_data$g_hat)
cor_d <- cor(analysis_data$d, analysis_data$d_hat)
cox_p_g <- summary(cox_model)$coefficients["g_hat", "Pr(>|z|)"]

message(sprintf(
  "Computed references: cor(g, g_hat)=%.6f, cor(d, d_hat)=%.6f, Cox p(g_hat)=%.6g",
  cor_g, cor_d, cox_p_g
))

assert(dplyr::between(cor_g, 0.89, 0.92), "Reference correlation cor(g, g_hat) drifted")
assert(dplyr::between(cor_d, 0.98, 1.00), "Reference correlation cor(d, d_hat) drifted")
assert(dplyr::between(cox_p_g, 0.015, 0.035), "Reference Cox p-value for g_hat drifted")

message("Smoke test passed: notebook rendered and generated outputs are valid")
