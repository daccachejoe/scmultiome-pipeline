# Loads config/pipeline.config into the R environment via Sys.setenv(), so
# scripts work whether invoked through run/runmultiome (which already
# exports these as shell env vars) or run standalone in an interactive
# session. Values already set in the environment are never overwritten.
load_pipeline_config <- function(path = "config/pipeline.config") {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
  for (line in lines) {
    kv <- regmatches(line, regexec("^([A-Za-z_][A-Za-z0-9_]*)=(.*)$", line))[[1]]
    if (length(kv) != 3) next
    key <- kv[2]
    value <- trimws(kv[3])
    value <- gsub("^['\"]|['\"]$", "", value)
    if (identical(Sys.getenv(key), "")) {
      do.call(Sys.setenv, setNames(list(value), key))
    }
  }
  invisible(NULL)
}
