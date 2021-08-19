regression_one_step <- function(UnivariateData, CoefficientCombination, Aggregation){
  # INPUT
  # UnivariateData[1:n]                   Numerical vector with n time series
  #                                       values.
  # CoefficientCombination[1:Scales+1]    Numerical vector with numbers which
  #                                       are associated with wavelet levels.
  #                                       The last number is associated with the
  #                                       smooth level. Each number determines
  #                                       the number of coefficient used per
  #                                       level. The selection follows a
  #                                       specific scheme.
  # Aggregation[1:Scales]    Numerical vector carrying numbers whose index is
  #                          associated with the wavelet level. The numbers
  #                          indicate the number of time in points used for
  #                          aggregation from the original time series.
  #
  # OUTPUT
  # forecast    Numerical value with one step forecast.
  # Author: QS, 02/2021
  if(!is.vector(UnivariateData)){
    message("Data must be of type vector!")
    return()
  }
  if(!is.vector(CoefficientCombination)){
    message("ccps must be of type vector")
    return()
  }
  if(!is.vector(Aggregation)){
    message("Aggregation must be of type vector")
    return()
  }

  if(is.null(Aggregation)){
    stop("Parameter Aggregation is not defined")
  }
  if((length(CoefficientCombination)-1) != length(Aggregation)){
    stop("Length of ccps must be longer than that of Aggregation by one.
    Parameter ccps needs one coefficient per wavelet level and one
         coefficient for the final smooth level.
         Parameter Aggregation defines the number of levels of the decomposition.")
  }
  scales = length(CoefficientCombination)-1
  # Decomposition
  dec_res <- decomposition(UnivariateData, Aggregation)
  # Training
  trs_res <- training(dec_res$UnivariateData,
                      dec_res$WaveletCoefficients,
                      dec_res$SmoothCoefficients,
                      dec_res$Scales,
                      CoefficientCombination,
                      Aggregation)
  arr_future_points = trs_res$points_in_future
  matrix = trs_res$lsmatrix
  # Optimization method
  Weights = regression_lsm_optimization(arr_future_points, matrix)
  # Forecast
  forecast = prediction_scheme(dec_res$WaveletCoefficients,
                               dec_res$SmoothCoefficients,
                               CoefficientCombination, Aggregation)
  forecast = Weights %*% forecast
  return(forecast)
}

