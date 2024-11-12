slidingWindowMetricGrange <- function(my_grange, scoring_window, score.technique, pseudo_count = 0.01, endrule = "constant")  {

  # TODO handle multichromosome granges
  # TODO add Roxygen comments to this function
  # TODO add assertions statements to ensure appropriate inputs
  # TODO test if this function can handle grange that doesn't include every position

  # Scores on positive strand
  cvg.pos <- GenomicRanges::coverage(my_grange[GenomicRanges::strand(my_grange) == "+"])
  window_scores.pos <- slidingWindowMetric(x = cvg.pos, k = scoring_window, technique = score.technique, pseudo_count = pseudo_count, endrule = endrule)
  window_scores.pos.grange <- GenomicRanges::GRanges(window_scores.pos)
  GenomicRanges::strand(window_scores.pos.grange) = "+"

  # Scores on negative strand
  cvg.neg <- GenomicRanges::coverage(my_grange[GenomicRanges::strand(my_grange) == "-"])
  window_scores.neg <- slidingWindowMetric(x = cvg.neg, k = scoring_window, technique = score.technique, pseudo_count = pseudo_count, endrule = endrule)
  window_scores.neg.grange <- GenomicRanges::GRanges(window_scores.neg)
  GenomicRanges::strand(window_scores.neg.grange) = "-"

  # combine window scores
  window_scores.grange <- c(window_scores.pos.grange, window_scores.neg.grange)

  return(window_scores.grange)
}
