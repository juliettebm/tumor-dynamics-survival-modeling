source(file.path("R", "recist.R"))
source(file.path("R", "pk_pd.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "tgi.R"))

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

trajectory <- data.frame(
  patient_id = c(1L, 1L, 1L, 2L, 2L),
  time = c(0, 6, 12, 0, 6),
  sld = c(40, 30, 40, 25, 24)
)
recist <- detect_progression(trajectory, "sld")
progression <- extract_progression_times(recist)
assert(identical(recist$is_progression, c(FALSE, FALSE, TRUE, FALSE, FALSE)), "RECIST rule drifted")
assert(progression$progression_event[progression$patient_id == 1] == 1L, "Progression event missing")
assert(progression$progression_time[progression$patient_id == 1] == 12, "Wrong progression time")
assert(progression$progression_event[progression$patient_id == 2] == 0L, "False progression event")

new_lesions <- data.frame(patient_id = c(1L, 2L), new_lesion_time = c(6, 6))
combined <- add_new_lesion_progression(progression, new_lesions)
assert(combined$progression_time[combined$patient_id == 1] == 6, "Earlier new lesion must win")
assert(combined$progression_cause[combined$patient_id == 1] == "new_lesion", "Wrong cause for patient 1")
assert(combined$progression_event[combined$patient_id == 2] == 1L, "New lesion alone must be a progression")
assert(combined$target_progression_time[combined$patient_id == 1] == 12, "Target time must be preserved")
none <- add_new_lesion_progression(progression, data.frame(patient_id = c(1L, 2L), new_lesion_time = NA_real_))
assert(none$progression_event[none$patient_id == 2] == 0L && none$progression_cause[none$patient_id == 2] == "none", "No-event case broken")
set.seed(1)
nl <- simulate_new_lesion_times(1:500, 0.01, seq(0, 104, by = 6))
assert(all(nl$new_lesion_time[!is.na(nl$new_lesion_time)] %in% seq(0, 104, by = 6)), "New lesions must be detected at scan times")
assert(all(is.na(simulate_new_lesion_times(1:20, 0, seq(0, 104, by = 6))$new_lesion_time)), "Zero rate must give no new lesions")

single <- single_dose_concentration(c(0, 1), 0.1, 5, 2)
assert(abs(single[1] - 2.5) < 1e-12, "Single-dose concentration at t=0 is wrong")
repeated <- repeated_dose_concentration(c(0, 1, 2), 0.1, c(0, 2), 5, 2)
assert(repeated[3] > single[2], "Repeated dosing did not accumulate concentration")
assert(all(emax_effect(c(0, 1, 10), 0.5, 1) >= 0), "Emax effect is invalid")

interval <- correlation_ci(0.83, 200)
assert(interval[["lower"]] < 0.83 && interval[["upper"]] > 0.83, "Correlation CI excludes its estimate")
invalid_fit <- try(fit_tgi_patient(data.frame(time = 0, wrong = 1)), silent = TRUE)
assert(inherits(invalid_fit, "try-error"), "Invalid TGI input should fail clearly")
message("Reusable-function tests passed")
