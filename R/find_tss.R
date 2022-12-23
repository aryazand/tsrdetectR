
find_tss <- function(input_bed, half_window_size, technique, pseudo_count = 1, threshold = NULL) {

  # Check input parameters -----------
  checkmate::assert(
    checkmate::checkClass(input_bed, "GRanges"),
    checkmate::checkNumeric(GenomicRanges::score(input_bed), lower = 0, any.missing = F),
    checkmate::checkCount(half_window_size, positive = T),
    checkmate::checkSubset(technique, c("mean_inclusive", "mean_exclusive", "geometric_mean")),
    checkmate::checkCount(pseudo_count, positive = T),
    combine = "and"
  )

  if (!is.null(threshold)) {
    checkmate::checkNumeric(threshold)
  }

  if (any(strand(input_bed) == "*")) {
    stop("all intervals in input_bed must be strand-specific")
  }

  ####

  # Prepare input_bed -----------

  # Get 5' ends
  fiveprime_ends <- input_bed |> GenomicRanges::resize(width = 1, fix = "start")
  fiveprime_ends.split <- split(fiveprime_ends, strand(fiveprime_ends)) |>
    as.list() |>
    within(expr = rm(`*`)) #remove automatically produced non-strand specific element

  ####

  # Calculate local background ------------------
  cvg <- coverage(fiveprime_ends)
  w = 1 + half_window_size*2
  local_bg <- get_local_background(cvg, k = w, technique = technique, pseudo_count = 0.01, endrule = "constant")

  ####

  # Calculate Score-to-background-ratio ------------------
  score_to_bg_ratio <-
    purrr::map(fiveprime_ends.split, coverage) |>
    purrr::map(~ (.x + pseudo_count)/(local_bg + pseudo_count))

  ####

  # Convert into GRanges object ----------------

  score_to_bg_ratio.grange <-
    purrr::map2(.x = score_to_bg_ratio, .y = names(x = score_to_bg_ratio), ~GRanges(.x, strand = .y)) |>
    as("GRangesList") |>
    unlist()

  # Rename column
  names(mcols(score_to_bg_ratio.grange))[names(mcols(score_to_bg_ratio.grange)) == "score"] <- "score_to_bg"
  score_to_bg_ratio.grange <- unname(score_to_bg_ratio.grange)

  ####

  # Apply threshold --------------
  if (!is.null(threshold)) {
    score_to_bg_ratio.grange[score(score_to_bg_ratio.grange) > threshold]
  } else {
    message("Please provide a threshold value. Automatic thresholding is a work-in-progress. Returning unthresholded GRange object")
    score_to_bg_ratio.grange
  }
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

