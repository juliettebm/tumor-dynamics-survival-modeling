# Reusable one-compartment concentration helpers used by the PK/PD reports.

single_dose_concentration <- function(time, elimination_rate, dose, volume) {
  if (elimination_rate <= 0 || dose <= 0 || volume <= 0) {
    stop("elimination_rate, dose and volume must be positive", call. = FALSE)
  }
  (dose / volume) * exp(-elimination_rate * time)
}

repeated_dose_concentration <- function(time, elimination_rate, dosing_times, dose, volume) {
  vapply(time, function(current_time) {
    active_doses <- dosing_times[dosing_times <= current_time]
    if (length(active_doses) == 0L) return(0)
    sum(single_dose_concentration(
      current_time - active_doses,
      elimination_rate = elimination_rate,
      dose = dose,
      volume = volume
    ))
  }, numeric(1))
}

emax_effect <- function(concentration, emax, ec50) {
  if (emax < 0 || ec50 <= 0) stop("emax must be non-negative and ec50 positive", call. = FALSE)
  emax * concentration / (ec50 + concentration)
}
