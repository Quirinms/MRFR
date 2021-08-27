mrf_one_step_forecast <- function(UnivariateData, CoefficientCombination,
                                  Aggregation, Method="r"){
  # DESCRIPTION
  # This function creates a one step forecast using the multiresolution
  # forecasting framework.
  #
  # INPUT
  # UnivariateData[1:n]    Vector with n time series values.
  #
  # OPTIONAL
  # CoefficientCombination    Vector with numbers which are associated with wavelet levels.
  #                           The last number is associated with the smooth level.
  #                           Each number determines the number of coefficient used per level.
  #                           The selection follows a specific scheme.
  # Aggregation               Vector carrying numbers whose index is associated with the
  #                           wavelet level. The numbers indicate the number of time in
  #                           points used for aggregation from the original time series.
  # Method           String indicating which method to use
  #                  Available methods: 'r'  = Regression
  #                                     'nn' = Neural Network
  # OUTPUT
  # forecast    Numerical value with one step forecast
  #
  # Author: QS, 02/2021
  if(!is.vector(UnivariateData)){
    message("Data must be of type vector")
    return()
  }
  if(!is.vector(CoefficientCombination)){
    message("ccps must be of type vector")
    return()
  }
  if(!is.vector(Aggregation)){
    message("agg_per_lvl must be of type vector")
    return()
  }
  if(Method == "r"){
    Forecast = mrf_regression_one_step_forecast(UnivariateData,
                                                CoefficientCombination,
                                                Aggregation)
  }else if(Method == "nn"){
    Forecast = mrf_neuralnet_one_step_forecast(UnivariateData,
                                               CoefficientCombination,
                                               Aggregation)
  }else{
    print("No valid methodname given => Returning.")
    Forecast = 0
  }
  return(Forecast)
}



#
#
#
#
#
