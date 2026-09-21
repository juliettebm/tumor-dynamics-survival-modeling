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

# Notebook 04 regression test: this catches unit mistakes and future-information
# leakage in the corrected repeated-dose model.
source_notebook_04 <- normalizePath(
  file.path("notebooks", "04_ADC_multi_lesion_PFS_OS.Rmd"),
  mustWork = TRUE
)
notebook_04 <- file.path(work_dir, "04_ADC_multi_lesion_PFS_OS.Rmd")
if (!file.copy(source_notebook_04, notebook_04, overwrite = TRUE)) {
  stop("Could not copy notebook 04 into the isolated working directory", call. = FALSE)
}

message("Rendering notebook 04 in ", work_dir)
render_environment_04 <- new.env(parent = globalenv())
rendered_file_04 <- rmarkdown::render(
  input = notebook_04,
  output_file = "04_ADC_multi_lesion_PFS_OS.html",
  output_dir = work_dir,
  intermediates_dir = work_dir,
  knit_root_dir = work_dir,
  quiet = TRUE,
  envir = render_environment_04
)

assert(file.exists(rendered_file_04), "Notebook 04 did not produce an HTML file")
assert(identical(get("Emax", render_environment_04), 0.003), "Notebook 04 repeated-dose Emax changed unexpectedly")
assert(exists("d_daily", get("lesions", render_environment_04)), "Notebook 04 does not use daily intrinsic regrowth")

time_grid_04 <- get("time_grid", render_environment_04)
trajectories_04 <- get("lesion_trajectories", render_environment_04)
os_data_04 <- get("os_data", render_environment_04)
prog_adc_04 <- get("prog_adc", render_environment_04)
total_sld_for_os_04 <- get("total_sld_for_os", render_environment_04)

assert(all(diff(time_grid_04) == 7), "Notebook 04 is expected to use a seven-day grid")
assert(all(is.finite(trajectories_04$SLD_adc) & trajectories_04$SLD_adc >= 0), "Invalid ADC trajectories")
assert(all(is.finite(trajectories_04$SLD_no_drug) & trajectories_04$SLD_no_drug >= 0), "Invalid no-drug trajectories")
assert(nrow(os_data_04) == 2L * nrow(patients), "OS output must contain both arms for every patient")
assert(nrow(total_sld_for_os_04) == length(time_grid_04) * nrow(patients), "Post-progression OS trajectories are incomplete")
assert(any(prog_adc_04$progression_event == 1L), "Stop-at-progression logic was not exercised")
assert(all(os_data_04$os_time %in% time_grid_04), "OS events must fall on the prospective time grid")
assert(!exists("final_sld", envir = render_environment_04), "OS must not depend on future final tumor burden")
assert(is.finite(get("pfs_p_value", render_environment_04)), "PFS log-rank result is invalid")
assert(is.finite(get("os_p_value", render_environment_04)), "OS log-rank result is invalid")

message("Smoke test passed: notebook 04 uses coherent units and prospective OS simulation")

# Notebook 04: the new-lesion process must be exercised and shared by both arms.
assert("progression_cause" %in% names(prog_adc_04), "Notebook 04 lost the progression cause column")
assert(any(prog_adc_04$progression_cause == "new_lesion"), "Notebook 04 never detected a new lesion")

# Notebook 03 regression test: new-lesion progression, PFS composite and loss to follow-up.
source_notebook_03 <- normalizePath(
  file.path("notebooks", "03_multi_lesion_RECIST_PFS.Rmd"),
  mustWork = TRUE
)
notebook_03 <- file.path(work_dir, "03_multi_lesion_RECIST_PFS.Rmd")
if (!file.copy(source_notebook_03, notebook_03, overwrite = TRUE)) {
  stop("Could not copy notebook 03 into the isolated working directory", call. = FALSE)
}

message("Rendering notebook 03 in ", work_dir)
render_environment_03 <- new.env(parent = globalenv())
rendered_file_03 <- rmarkdown::render(
  input = notebook_03,
  output_file = "03_multi_lesion_RECIST_PFS.html",
  output_dir = work_dir,
  intermediates_dir = work_dir,
  knit_root_dir = work_dir,
  quiet = TRUE,
  envir = render_environment_03
)

assert(file.exists(rendered_file_03), "Notebook 03 did not produce an HTML file")

progression_03 <- get("progression_times", render_environment_03)
pfs_03 <- get("pfs_data", render_environment_03)
comparison_03 <- get("comparison_data", render_environment_03)

assert(nrow(progression_03) == nrow(patients), "Notebook 03 must have one progression record per patient")
assert(all(c("target", "new_lesion") %in% progression_03$progression_cause), "Notebook 03 must exercise target and new-lesion progression")
assert(all(progression_03$progression_event == as.integer(progression_03$progression_cause != "none")), "Progression event and cause disagree")
assert(all(pfs_03$pfs_time <= pfs_03$pfs_time_true), "Loss to follow-up can only shorten PFS time")
assert(all(pfs_03$pfs_event <= pfs_03$pfs_event_true), "Loss to follow-up can only remove PFS events")
assert(all(pfs_03$pfs_time <= pfs_03$os_time), "PFS cannot exceed OS time")
assert(all(is.finite(comparison_03$time) & comparison_03$time >= 0), "Invalid PFS/OS times in notebook 03")

message("Smoke test passed: notebook 03 combines target and new-lesion progression coherently")
