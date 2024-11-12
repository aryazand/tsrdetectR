#' Identify 5' end of all ranges in a GRanges object
#'
#' @param input_bed GRanges object
#'
#' @return GRanges object with all ranges of width 1
#' @export
#'
#' @examples
#' fragStartPosGRange(sample_NT2_data)
fragStartPosGRange <- function(input_bed) {

  checkmate::assert(
    checkmate::checkClass(input_bed, "GRanges")
  )

  fiveprime <- input_bed |> GenomicRanges::resize(width = 1, fix = "start")

  return(fiveprime)
}
