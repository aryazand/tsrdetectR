normalize <- function(x) {
  (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
}

threshold_value_by_slope <- function(y0, deriv.thershold = 1) {

  # TODO add roxygen comments

  y <- sort(y0)
  y_smooth <- zoo::rollmean(y, k = length(y)*0.1, fill = NA) |> log2()
  y_smooth_normalized <- normalize(y_smooth)
  x = seq(0, 1, length.out = length(y_smooth_normalized))
  spline_function <- stats::splinefun(x,y_smooth_normalized)
  threshold = which(spline_function(x, deriv = 1) > deriv.thershold & spline_function(x) > stats::median(y_smooth_normalized, na.rm = T)) |> min()
  threshold = 2^y_smooth[threshold]
  return(threshold)
}
