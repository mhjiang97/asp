asp_default_options <- list(
  asp_tools = c(
    "rmats" = "rMATS",
    "cash" = "CASH",
    "leafcutter" = "LeafCutter",
    "majiq" = "MAJIQ",
    "spladder" = "SplAdder",
    "bandits" = "BANDITS",
    "suppa" = "SUPPA",
    "psichomics" = "psichomics"
  )
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = asp_default_options) %in% names(x = op))
  if (any(toset)) options(asp_default_options[toset])
  invisible()
}




