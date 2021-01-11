#' Rolling Window for Multiresolution Forecasts
#'
#' This function creates a rolling window forecast using a multi layer perceptron
#' with one hidden Layer or a regression.
#'
#' @param data One dimensional array of signal values
#' @param i 1D array containing the number of coefficients chosen per scales
#' @param ccps 1D array containing the number of coefficients chosen per scales
#' @param int_total_length 1D array containing the number of coefficients chosen per scales
#' @param horizon Integer indicating the number of steps to forecast ahead
#' @param window_size Integer indicating the number of days, for which a forecast
#' should be evaluated. The procedure uses automatically the last part of the time
#' series
#' @param method Character indicating Regression ("r") or the Neural Network ("n")
#' @return List of parameter with a 2D matrix of the forecast error.
#'
rolling_window_single <- function(i, data, ccps, agg_per_lvl, horizon = 14,
                                  window_size = 365, method = "r"){
  int_total_length  = length(data)                        # Length time series
  int_CFCP = int_total_length - window_size - horizon + i # Current Forecast Position
  dfTrain  = data[1:int_CFCP]
  dfTest   = data[int_CFCP+1:horizon]
  forecast = multi_step(data = dfTrain, steps = horizon, ccps = ccps,
                        method = method, agg_per_lvl = agg_per_lvl)
  arr_Error = as.numeric(forecast$forecast) - dfTest
  return(list("Error" = arr_Error))
}
