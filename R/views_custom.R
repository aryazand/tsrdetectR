views_custom <- function(x, center, half_window, upstream = half_window, downstream = half_window, seq_length, circular = F) {

  # ToDO: add roxygen documentation
  # TODO: add assertion for function inputs

  ##########################
  # calculate start and stop
  ##########################

  v_start <- center - upstream
  v_end <- center + downstream

  if(!circular) {
    # ensure no negative
    v_start[v_start < 0] <- 0
    v_end[v_end > seq_length] <- seq_length
  }

  ##########################
  # create views
  ##########################

  v <- IRanges::Views(x, start = v_start, end = v_end)

  return(v)
}
