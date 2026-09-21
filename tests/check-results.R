manifest_path <- file.path("results", "key_results.psv")
results <- utils::read.delim(manifest_path, sep = "|", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
readme <- paste(readLines("README.md", warn = FALSE, encoding = "UTF-8"), collapse = "\n")

for (index in seq_len(nrow(results))) {
  source_path <- results$source_file[[index]]
  if (!file.exists(source_path)) stop("Missing result source: ", source_path, call. = FALSE)
  source <- paste(readLines(source_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!grepl(results$source_claim[[index]], source, fixed = TRUE)) {
    stop("Source claim drifted in ", source_path, ": ", results$source_claim[[index]], call. = FALSE)
  }
  if (!grepl(results$readme_claim[[index]], readme, fixed = TRUE)) {
    stop("README claim drifted: ", results$readme_claim[[index]], call. = FALSE)
  }
}

message("README/result contract passed for ", nrow(results), " key claims")
