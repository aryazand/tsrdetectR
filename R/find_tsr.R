
#' Identify Transcription Start Regions
#'
#' @param five_prime_grange A GRanges object with all elements of width 1
#'   represent the 5' location of PRO-seq reads
#' @param half_scoring_window_size The half window to score signal
#' @param half_background_window_size The window around each potential TSS to
#'   assess the local background signal
#' @param score.technique A string identifying the technique to assess local
#'   background signal
#' @param pseudo_count A pseudo-count when evaluating the ratio between signal
#'   at a potential TSS and local-background signal
#' @param threshold.technique method for determining the threshold
#'   score-to-background ratio that distinguishes a TSR from background
#' @param threshold.value Threshold value of signal-to-local background ratio to
#'   call a TSS
#'
#' @return The five_prime_grange with two column added or replaced:
#'   score_to_bg_ratio (a numeric) and TSS (a logical)
#' @export
#'
#' @examples
#' x = fragStartPosGRange(sample_NT2_data)
#' find_tsr(x, half_scoring_window_size = 50, half_background_window_size = 100,
#' score.technique = "mean", threshold.technique = "z-score", threshold.value = 1)
find_tsr <- function(five_prime_grange,
                     half_scoring_window_size,
                     half_background_window_size,
                     score.technique,
                     pseudo_count = 1,
                     threshold.technique = NULL,
                     threshold.value = NULL) {

  # Check input parameters -----------
  checkmate::assert(
    checkmate::checkClass(five_prime_grange, "GRanges"),
    checkmate::checkNumeric(GenomicRanges::width(five_prime_grange), lower = 1, upper = 1),
    checkmate::checkNumeric(BiocGenerics::score(five_prime_grange), lower = 0, any.missing = F),
    checkmate::checkNumber(five_prime_grange@seqinfo@seqlengths),
    checkmate::checkNumber(half_scoring_window_size, lower = 0),
    checkmate::checkNumber(half_background_window_size, lower = 0),
    checkmate::checkChoice(score.technique, c("mean", "mean_exclusive", "geometric_mean")),
    checkmate::checkCount(pseudo_count, positive = T),
    checkmate::checkChoice(threshold.technique, c("defined", "distribution", "z-score", "poisson", "negative binomial")),
    checkmate::checkNumeric(threshold.value),
    combine = "and"
  )

  if(threshold.technique %in% c("z-score", "poisson", "negative binomial") & !(score.technique %in% c("mean", "mean_exclusive"))) {
    stop("if threshold.technique is z-score, poisson, or negative binomial then score.technique must be mean or mean_exclusive")
  }

  if (any(BiocGenerics::strand(five_prime_grange) == "*")) {
    stop("all intervals in five_prime_grange must be strand-specific")
  }

  ####

  # create grange covering each position in genome ----------------------------------

  tss.grange.pos <- five_prime_grange[GenomicRanges::strand(five_prime_grange) == "+"] |> GenomicRanges::coverage() |> GenomicRanges::GRanges()
  tss.grange.pos <- GenomicRanges::GRanges(
    S4Vectors::Rle(tss.grange.pos@seqinfo@seqnames, tss.grange.pos@seqinfo@seqlengths),
    IRanges::IRanges(start = seq(tss.grange.pos@seqinfo@seqlengths), width = 1),
    strand = "+",
    score = rep(tss.grange.pos$score, GenomicRanges::width(tss.grange.pos@ranges))
  )

  tss.grange.neg <- five_prime_grange[GenomicRanges::strand(five_prime_grange) == "-"] |> GenomicRanges::coverage() |> GenomicRanges::GRanges()
  tss.grange.neg <- GenomicRanges::GRanges(
    S4Vectors::Rle(tss.grange.neg@seqinfo@seqnames, tss.grange.neg@seqinfo@seqlengths),
    IRanges::IRanges(start = seq(tss.grange.neg@seqinfo@seqlengths), width = 1),
    strand = "-",
    score = rep(tss.grange.neg$score, GenomicRanges::width(tss.grange.neg@ranges))
  )

  tss.grange <- c(tss.grange.pos, tss.grange.neg)
  GenomeInfoDb::seqlengths(tss.grange) <- GenomeInfoDb::seqlengths(five_prime_grange)

  # Calculate sliding window of scores ------------------

  # Calculate scores
  scoring_window <- 1 + half_scoring_window_size*2
  window_scores.grange <- slidingWindowMetricGrange(five_prime_grange, scoring_window, score.technique)

  # Add scores to fiveprime-end grange object
  overlaps_scores <- GenomicRanges::findOverlaps(tss.grange, window_scores.grange)
  tss.grange$window_scores <- BiocGenerics::score(window_scores.grange)[overlaps_scores@to]

  # Calculate sliding window of background scores
  background_window <- 1 + half_background_window_size*2
  background_scores.grange <- slidingWindowMetricGrange(five_prime_grange, background_window, score.technique = "mean")

  # Add scores to fiveprime-end grange object
  overlaps_scores <- GenomicRanges::findOverlaps(tss.grange, background_scores.grange)
  tss.grange$background_scores <- BiocGenerics::score(background_scores.grange)[overlaps_scores@to]

  # Calculate score-to-local-background ratio
  tss.grange$score_to_bg_ratio <- (tss.grange$window_scores + pseudo_count)/(tss.grange$background_scores + pseudo_count)

  ####

  # Apply threshold --------------
  if (threshold.technique == "defined") {
    message(paste("Thresholding on a value of", threshold))
    tss.grange$sig.score <- tss.grange$score_to_bg_ratio > threshold.value
  } else if (threshold.technique == "distribution") {
    threshold <- threshold_value_by_slope(tss.grange$score_to_bg_ratio, threshold.value)
    message(paste("Thresholding on a value of", threshold))
    tss.grange$sig.score <- tss.grange$score_to_bg_ratio > threshold
  } else if (threshold.technique == "z-score") {
    sd.pos <- zoo::rollapply(BiocGenerics::score(tss.grange[GenomicRanges::strand(tss.grange) == "+"]), width=background_window, FUN = stats::sd, fill = NA)
    sd.neg <- zoo::rollapply(BiocGenerics::score(tss.grange[GenomicRanges::strand(tss.grange) == "-"]), width=background_window, FUN = stats::sd, fill = NA)
    tss.grange$background_sd <- c(sd.pos, sd.neg)
    tss.grange$zscore <- (tss.grange$score - tss.grange$background_scores)/tss.grange$background_sd
    threshold <- threshold_value_by_slope(tss.grange$zscore, threshold.value)
    tss.grange$sig.score <- tss.grange$zscore > threshold
  }

  return(tss.grange)
}
