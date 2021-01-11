#' Rolling Window for Multiresolution Forecasts
#'
#' This function creates a rolling window forecast using a multi layer perceptron
#' with one hidden Layer or a regression.
#'
#' @param data One dimensional array of signal values
#' @param ccps 1D array containing the number of coefficients chosen per scales
#' @param horizon Integer indicating the number of steps to forecast ahead
#' @param window_size Integer indicating the number of days, for which a forecast
#' should be evaluated. The procedure uses automatically the last part of the time
#' series
#' @param method Character indicating Regression ("r") or the Neural Network ("n")
#' @param str_save String, if not empty, then execute a save to the given
#' destination
#' @return List of parameter with a 2D matrix of the forecast error.
#' @export
rolling_window_parallel <- function(data, ccps, agg_per_lvl, horizon = 14,
                                    window_size = 365, method = "r"){
  cores = parallel::detectCores()
  cl <- parallel::makeCluster(cores[1]-1) #not to overload your computer
  parallel::clusterEvalQ(cl, source("decomposition.R"))
  parallel::clusterEvalQ(cl, source("training_scheme.R"))
  parallel::clusterEvalQ(cl, source("multi_step.R"))
  parallel::clusterEvalQ(cl, source("svm_one_step.R"))
  parallel::clusterEvalQ(cl, source("xgboost_one_step.R"))
  parallel::clusterEvalQ(cl, source("neuralnet_one_step.R"))
  parallel::clusterEvalQ(cl, source("neuralnet_prediction_scheme.R"))
  parallel::clusterEvalQ(cl, source("regression_one_step.R"))
  parallel::clusterEvalQ(cl, source("regression_prediction_scheme.R"))
  parallel::clusterEvalQ(cl, source("regression_lsm_optimization.R"))
  #parallel::clusterEvalQ(cl, library("mrf"))
  lst_forecast = parallel::parLapply(cl, 1:window_size, rolling_window_single,
                                     data = data, ccps = ccps,
                                     agg_per_lvl = agg_per_lvl,
                                     horizon = horizon,
                                     window_size = window_size,
                                     method = method)
  array_forecast = array(unlist(lst_forecast))
  matrix_forecast = matrix(array_forecast, ncol = horizon, byrow = TRUE)
  parallel::stopCluster(cl)
  result = list("Rolling_Window" = matrix_forecast)    # Return
  return(result)
}

