#' Multi Step Forecast
#'
#' This function creates a multi step forecast using a multi layer perceptron
#' with one hidden Layer or an auto regressive approach.
#' Several steps ahead are computed recursively.
#'
#' @param data One dimensional array of signal values
#' @param steps Integer indicating the number of steps to forecast ahead (forecast horizon)
#' @param scales Integer Number of wavelet scales to use
#' @param ccps 1D array containing the number of coefficients chosen per scales
#' (includes the number of coefficients for the smooth coefficient)
#' @param method Character indicating Regression ("r") or the Neural Network ("n")
#' @return List of parameter with array carrying "steps" many forecasts.
#' @examples
#' data(AirPassengers)
#' len_data = length(array(AirPassengers))
#' result = multi_step(array(AirPassengers)[1:(len_data-1)], 2, 2, c(1,1,1), method="r")
#' forecast = result$forecast
#' @export
multi_step <- function(data, steps, ccps, agg_per_lvl, method = "r"){
  multistep = c()
  for(i in 1:steps){
    forecast = 0
    if(method == "r"){
      forecast = regression_one_step(data, ccps, agg_per_lvl)
      forecast = forecast$forecast
    }
    if(method == "nn"){
      forecast = neuralnet_one_step(data, ccps, agg_per_lvl)
      forecast = forecast$forecast
    }
    if(method == "xgboost"){
      forecast = xgboost_one_step(data, ccps, agg_per_lvl)
      forecast = forecast$forecast
    }
    if(method == "svm"){
      forecast = svm_one_step(data, ccps, agg_per_lvl)
      forecast = forecast$forecast
    }
    
    multistep = c(multistep, forecast)
    data = array(unlist(c(data, forecast)))
  }
  multistep = array(unlist(multistep))
  # Cap forecasts exceeding certain limits
  max_val = max(data)
  upper_limit = 0
  if(max_val > 0){
    upper_limit = 1.3*max(data)
  }else{
    upper_limit = 0.7*max(data)
  }
  min_val = min(data)
  lower_limit = 0
  if(min_val > 0){
    lower_limit = 0.7*min(data)
  }else{
    lower_limit = 1.3*min(data)
  }
  for(i in 1:steps){
    forecast = multistep[i]
    if(forecast > upper_limit){
      forecast = upper_limit
    }
    if(forecast < lower_limit){
      forecast = lower_limit
    }
    multistep[i] = forecast
  }
  result = list("forecast" = multistep)    # Return
  return(result)
}


#
#
#
#
#
