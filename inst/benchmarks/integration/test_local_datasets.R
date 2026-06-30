# User-supplied local datasets (optional local_paths.R)

cfg <- .integration_state$config
local_paths <- cfg$local_paths

if (!length(local_paths)) {
  integration_record(
    "local:user:none",
    "No local_paths.R datasets configured",
    "skip",
    list(reason = "copy local_paths.R.example to local_paths.R to enable")
  )
} else {
  for (nm in names(local_paths)) {
    path <- local_paths[[nm]]
    integration_run_test(
      paste0("local:", nm),
      paste("User local dataset:", nm),
      skip_if = if (!nzchar(path) || !file.exists(path)) {
        paste("file not found:", path)
      } else NULL,
      expr = {
        obj <- readRDS(path)
        scores <- PhenoMap(
          expression = obj,
          reference = "precog",
          cancer_type = "Pancreatic",
          verbose = FALSE
        )
        sc <- scores[[grep("Pancreatic", colnames(scores), value = TRUE)[1]]]
        integration_assert(all(is.finite(sc)), "non-finite scores")
        list(n = length(sc), path = path)
      }
    )
  }
}
