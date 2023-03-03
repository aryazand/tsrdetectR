
find_tss <- function(five_prime_grange, half_window_size, technique, pseudo_count = 1, threshold = NULL) {

  # Check input parameters -----------
  checkmate::assert(
    checkmate::checkClass(five_prime_grange, "GRanges"),
    checkmate::checkNumeric(GenomicRanges::width(five_prime_grange), lower = 1, upper = 1),
    checkmate::checkNumeric(BiocGenerics::score(five_prime_grange), lower = 0.01, any.missing = F),
    checkmate::checkCount(half_window_size, positive = T),
    checkmate::checkSubset(technique, c("mean_inclusive", "mean_exclusive", "geometric_mean")),
    checkmate::checkCount(pseudo_count, positive = T),
    combine = "and"
  )

  if (!is.null(threshold)) {
    checkmate::checkNumeric(threshold)
  }

  if (any(BiocGenerics::strand(five_prime_grange) == "*")) {
    stop("all intervals in five_prime_grange must be strand-specific")
  }

  ####

  # Add score-to-local-background ratio to GRanges ------------------

  # Calculate local background
  cvg <- GenomicRanges::coverage(five_prime_grange)
  w = 1 + half_window_size*2
  local_bg <- get_local_background(cvg, k = w, technique = technique, pseudo_count = 0.01, endrule = "constant")
  local_bg_grange <- GenomicRanges::GRanges(local_bg)

  # Add to local background to GRanges
  overlaps <- GenomicRanges::findOverlaps(five_prime_grange, local_bg_grange)
  five_prime_grange$bg_score <- BiocGenerics::score(local_bg_grange)[overlaps@to]

  # Calculate score-to-local-background ratio
  five_prime_grange$score_to_bg_ratio <- (five_prime_grange$score + pseudo_count)/(five_prime_grange$bg_score + pseudo_count)

  ####

  # Apply threshold --------------
  if (is.null(threshold)) {
    threshold <- get_threshold_value(five_prime_grange$score_to_bg_ratio)
  }

  message(paste("Thresholding on a value of", threshold))
  five_prime_grange$TSS <- five_prime_grange$score_to_bg_ratio > threshold

  return(five_prime_grange)
}

get_local_background <- function(x, k, technique, pseudo_count = 0.1, ...) {

  checkmate::assert(
    checkmate::checkClass(x, c("SimpleRleList")),
    checkmate::checkCount(k, positive = T),
    checkmate::checkNumber(k, lower = 2),
    checkmate::checkSubset(technique, c("mean_inclusive", "mean_exclusive", "geometric_mean")),
    combine = "and"
  )

  if (technique == "mean_inclusive") {
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

make_fiveprimeends_grange <- function(input_bed) {

  checkmate::assert(
    checkmate::checkClass(input_bed, "GRanges")
  )

  fiveprime <- input_bed |> GenomicRanges::resize(width = 1, fix = "start")
  fiveprime_unique <- fiveprime |> BiocGenerics::unique()
  fiveprime_unique_counts <- GenomicRanges::countOverlaps(fiveprime_unique, fiveprime)
  BiocGenerics::score(fiveprime_unique) <- fiveprime_unique_counts

  return(fiveprime_unique)
}

get_threshold_value <- function(y) {
  y <- sort(y)
  y_smooth <- zoo::rollmean(y, k = length(y)*0.1, fill = NA) |> log2()
  x = seq(0, 1, length.out = length(y_smooth))
  spline_function <- splinefun(x,y_smooth)
  threshold = which(spline_function(x, deriv = 1) > 1 & spline_function(x) > median(y_smooth, na.rm = T)) |> min()
  threshold = 2^y_smooth[threshold]
  return(threshold)
}

