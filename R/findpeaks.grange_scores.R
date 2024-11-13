findpeaks.grange_scores <- function(x, minpeakdistance) {

  # TODO Roxygen documentation

  checkmate::assert(
    checkmate::checkClass(x, "GRanges"),
    checkmate::checkNumeric(BiocGenerics::score(x), lower = 0, any.missing = F),
    checkmate::checkNumber(minpeakdistance),
    combine = "and"
  )

  peaks <- pracma::findpeaks(x$score, minpeakdistance = minpeakdistance)
  x$peak <- seq(length(x)) %in% peaks[,2]

  peaks.df <- dplyr::left_join(data.frame(V2 = seq(length(x))), as.data.frame(peaks)[,2:4])
  peaks.df <- peaks.df[,-1]
  names(peaks.df) <- c("peak.start", "peak.end")
  S4Vectors::mcols(x) <- cbind(S4Vectors::mcols(x), peaks.df)

  return(x)

}

