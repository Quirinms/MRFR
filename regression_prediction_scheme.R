#' One Step Forecast with Regression
#'
#' This function delivers the required coefficients from the decomposition,
#' which take part in forecasting the next point of the time series for the
#' auto regressive approach.
#'
#' @param weights 1D array containing regression weights
#' @param wmatrix 2D matrix with wavelet coefficients
#' @param rhwtCoeff 2D matrix with smooth coefficients
#' @param scales Number of wavelet scales (smooth scale excluded!)
#' @param ccps 1D array containing the number of coefficients chosen per each
#' scale. It accounts for all wavelet levels and the last smooth coefficient
#' level. The last number of the array is the number of coefficients for the
#' smooth level. Therefore the length of this array is larger by one than the
#' parameter "scales".
#' @return List of parameter float forecast.
#' @examples
#' data(AirPassengers)
#' len_data = length(array(AirPassengers))
#' ccps = c(1,1,1)
#' scales = length(ccps) - 1
#' # Decomposition
#' dec_res <- decomposition(array(AirPassengers), scales = scales)
#' # Training
#' trs_res <- training_scheme(dec_res$data, dec_res$wmatrix,
#' dec_res$rhwtCoeff, dec_res$scales, ccps)
#' arr_future_points = trs_res$points_in_future
#' matrix = trs_res$lsmatrix
#' # Optimization method
#' weights = regression_lsm_optimization(arr_future_points, matrix)
#' # Forecast
#' forecast = regression_prediction_scheme(weights, dec_res$wmatrix,
#' dec_res$rhwtCoeff, dec_res$scales, ccps)
#' forecast = forecast$future_point
#'@export
regression_prediction_scheme <- function(weights, wmatrix, rhwtCoeff, ccps,
                                         agg_per_lvl = NULL){
  time = dim(wmatrix)[2]
  future_point = 0
  counter = 1
  scales = length(agg_per_lvl)
  for(s in 0:scales){
    for(k in 0:(ccps[s+1]-1)){
      if(s != scales){
        index = time - k*agg_per_lvl[s+1]
        future_point = future_point + ((weights[[1]][counter])) * wmatrix[s+1, index]
      }
      else{
        index = time - k*agg_per_lvl[s]
        future_point = future_point + ((weights[[1]][counter])) * rhwtCoeff[s, index]
      }
      counter = counter + 1
    }
  }
  result = list("future_point" = future_point)
  return(result)
}



regression_prediction_scheme_deprecated <- function(weights, wmatrix, rhwtCoeff, scales,
                                         ccps, dyadic = TRUE,
                                         agg_per_lvl = NULL){
  time = dim(wmatrix)[2]
  future_point = 0
  counter = 1
  for(s in 0:scales){
    for(k in 0:(ccps[s+1]-1)){
      if(s != scales){
        index = get_index(s, scales, time, k, dyadic, agg_per_lvl)
        future_point = future_point + ((weights[[1]][counter])) * wmatrix[s+1, index]
      }
      else{
        index = get_index(s, scales, time, k, dyadic, agg_per_lvl)
        future_point = future_point + ((weights[[1]][counter])) * rhwtCoeff[s, index]
      }
      counter = counter + 1
    }
  }
  result = list("future_point" = future_point)
  return(result)
}

get_index_deprecated <- function(s, scales, time, k, dyadic, agg_per_lvl = NULL){
  index = 0
  level = s
  if(s != scales){
    level = s + 1
  }
  if(dyadic == TRUE){
    index = time - (k)*(2**level)
  }else{
    translation = 1
    for(i in 1:level){
      translation = agg_per_lvl[i] * translation
    }
    index = (time - (k * translation))
  }
  return(index)
}

#
#
#
#
#
