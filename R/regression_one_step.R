regression_one_step <- function(data, ccps, agg_per_lvl){
  # INPUT
  # data[1:n]             Vector with n time series values.
  # ccps                  Vector with numbers which are associated with wavelet levels.
  #                       The last number is associated with the smooth level.
  #                       Each number determines the number of coefficient used per level.
  #                       The selection follows a specific scheme.
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  #
  # OUTPUT
  # forecast              Double one step forecast.
  # Author: QS, 02/2021
  if(!is.vector(data)){
    message("Data must be of type vector!")
    return()
  }
  if(!is.vector(ccps)){
    message("ccps must be of type vector")
    return()
  }
  if(!is.vector(agg_per_lvl)){
    message("agg_per_lvl must be of type vector")
    return()
  }
  
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

