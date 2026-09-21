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
