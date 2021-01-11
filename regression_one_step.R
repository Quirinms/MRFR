#' One Step Forecast with Regression
#'
#' This function computes a one step forecast using the Regression approach.
#'
#' @param data One dimensional array of signal values
#' @param ccps 1D array containing the number of coefficients chosen per each
#' scale. It accounts for all wavelet levels and the last smooth coefficient
#' level. The last number of the array is the number of coefficients for the
#' smooth level. Therefore the length of this array is larger by one than the
#' parameter "scales".
#' @return List of parameter with float forecast.
#'
#' @examples
#' data(AirPassengers)
#' len_data = length(array(AirPassengers))
#' result = regression_one_step(array(AirPassengers)[1:(len_data-1)], c(1,1,1))
#' forecast = result$forecast
#' @export
regression_one_step <- function(data, ccps, agg_per_lvl){
  if(is.null(agg_per_lvl)){
    stop("Parameter agg_per_lvl is not defined")
  }
  if((length(ccps)-1) != length(agg_per_lvl)){
    stop("Length of ccps must be longer than that of agg_per_level by one.
    Parameter ccps needs one coefficient per wavelet level and one
         coefficient for the final smooth level.
         Parameter agg_per_level defines the number of levels of the decomposition.")
  }
  scales = length(ccps)-1
  # Decomposition
  dec_res <- decomposition(data, agg_per_lvl)
  # Training
  trs_res <- training_scheme(dec_res$data, dec_res$wmatrix, dec_res$rhwtCoeff,
                             dec_res$scales, ccps, agg_per_lvl)
  arr_future_points = trs_res$points_in_future
  matrix = trs_res$lsmatrix
  # Optimization method
  weights = regression_lsm_optimization(arr_future_points, matrix)
  # Forecast
  forecast = regression_prediction_scheme(weights, dec_res$wmatrix,
                                          dec_res$rhwtCoeff,
                                          ccps, agg_per_lvl)
  forecast = forecast$future_point
  result = list("forecast" = forecast)
  return(result)
}



