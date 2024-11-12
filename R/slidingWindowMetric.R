slidingWindowMetric <- function(x, k, technique, pseudo_count = 0.1, ...) {

  # TODO Add Roxygen comments to this function

  checkmate::assert(
    checkmate::checkClass(x, c("SimpleRleList")),
    checkmate::checkCount(k, positive = T),
    checkmate::checkNumber(k, lower = 1),
    checkmate::checkSubset(technique, c("sum", "mean", "mean_exclusive", "geometric_mean")),
    combine = "and"
  )

  if (technique == "sum") {
    S4Vectors::runsum(x,k,...)
  } else if (technique == "mean") {
    S4Vectors::runmean(x,k,...)
  } else if (technique == "mean_exclusive") {
    (S4Vectors::runmean(x, k = k, ...) - (x/k)) * k/(k-1)
  } else if (technique == "geometric_mean") {
    (x + pseudo_count) |>
      log() |>
      S4Vectors::runmean(x = _, k = k, ...) |>
      exp()
  }
}
