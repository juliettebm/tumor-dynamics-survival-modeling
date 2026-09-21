# Reusable helpers for the simplified target-lesion component of RECIST 1.1.

detect_progression <- function(data, sld_column) {
  required <- c("patient_id", "time", sld_column)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data |>
    dplyr::group_by(.data$patient_id) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::mutate(
      nadir = cummin(.data[[sld_column]]),
      relative_increase = (.data[[sld_column]] - .data$nadir) / .data$nadir,
      absolute_increase = .data[[sld_column]] - .data$nadir,
      is_progression = .data$relative_increase >= 0.20 & .data$absolute_increase >= 5
    ) |>
    dplyr::ungroup()
}

extract_progression_times <- function(recist_data) {
  required <- c("patient_id", "time", "is_progression")
  missing <- setdiff(required, names(recist_data))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  recist_data |>
    dplyr::group_by(.data$patient_id) |>
    dplyr::summarise(
      progression_time = if (any(.data$is_progression)) min(.data$time[.data$is_progression]) else NA_real_,
      progression_event = as.integer(any(.data$is_progression)),
      last_observed_time = max(.data$time),
      .groups = "drop"
    )
}

# New-lesion progression: an unequivocal new lesion is a RECIST 1.1 progression on its
# own, independent of target-lesion growth. Appearance is a homogeneous Poisson process
# (constant rate per time unit) and is only observed at the next scan.
simulate_new_lesion_times <- function(patient_ids, rate, scan_times) {
  if (rate < 0) stop("rate must be non-negative", call. = FALSE)
  appearance <- if (rate == 0) rep(Inf, length(patient_ids)) else stats::rexp(length(patient_ids), rate = rate)
  detected <- vapply(appearance, function(a) {
    scans <- scan_times[scan_times >= a]
    if (length(scans) == 0L) NA_real_ else min(scans)
  }, numeric(1))
  data.frame(patient_id = patient_ids, new_lesion_time = detected)
}

# Combine target-lesion progression with new-lesion detection: the earliest wins.
add_new_lesion_progression <- function(progression_times, new_lesions) {
  required <- c("patient_id", "progression_time", "progression_event")
  missing <- setdiff(required, names(progression_times))
  if (length(missing) > 0L || !all(c("patient_id", "new_lesion_time") %in% names(new_lesions))) {
    stop("Missing columns for new-lesion progression", call. = FALSE)
  }

  progression_times |>
    dplyr::left_join(new_lesions, by = "patient_id") |>
    dplyr::mutate(
      target_progression_time = .data$progression_time,
      progression_time = suppressWarnings(
        pmin(.data$target_progression_time, .data$new_lesion_time, na.rm = TRUE)
      ),
      progression_time = dplyr::if_else(is.finite(.data$progression_time),
                                        .data$progression_time, NA_real_),
      progression_event = as.integer(!is.na(.data$progression_time)),
      progression_cause = dplyr::case_when(
        is.na(.data$progression_time) ~ "none",
        is.na(.data$new_lesion_time) ~ "target",
        is.na(.data$target_progression_time) ~ "new_lesion",
        .data$new_lesion_time < .data$target_progression_time ~ "new_lesion",
        TRUE ~ "target"
      )
    )
}
