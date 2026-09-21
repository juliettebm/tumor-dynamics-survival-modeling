source(file.path("R", "recist.R"))
source(file.path("R", "pk_pd.R"))
source(file.path("R", "metrics.R"))

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

single <- single_dose_concentration(c(0, 1), 0.1, 5, 2)
assert(abs(single[1] - 2.5) < 1e-12, "Single-dose concentration at t=0 is wrong")
repeated <- repeated_dose_concentration(c(0, 1, 2), 0.1, c(0, 2), 5, 2)
assert(repeated[3] > single[2], "Repeated dosing did not accumulate concentration")
assert(all(emax_effect(c(0, 1, 10), 0.5, 1) >= 0), "Emax effect is invalid")

interval <- correlation_ci(0.83, 200)
assert(interval[["lower"]] < 0.83 && interval[["upper"]] > 0.83, "Correlation CI excludes its estimate")
message("Reusable-function tests passed")
