multi_step <- function(UnivariateData, Horizon, CoefficientCombination, Aggregation, Method = "r"){
  # DESCRIPTION
  # Computes multi-step forecasts with the multiresolution method.
  #
  # INPUT
  # UnivariateData[1:n]                   Numerical vector with n time series
  #                                       values.
  # Horizon          Number indicating horizon for forecast from 1 to horizon.
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
  # OPTIONAL
  # Method           String indicating which method to use
  #                  Available methods: 'r'  = Autoregression
  #                                     'nn' = Neural Network
  #
  # OUTPUT
  # multistep[1:Horizon]    Numerical vector with forecast of horizon according
  #                         to its index.
  #
  # Author: QS, 02/2021
  if(!is.vector(UnivariateData)){
    message("Data must be of type vector")
    return()
  }
  if(!is.double(Horizon)){
    message("steps must be of type double")
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
  if(!is.character(Method)){
    message("method must be of type character")
    return()
  }
  multistep = c()
  for(i in 1:Horizon){
    forecast = 0
    forecast = onestep(UnivariateData, CoefficientCombination, Aggregation, Method)
    multistep = c(multistep, forecast)
    UnivariateData = as.vector(unlist(c(UnivariateData, forecast)))
  }
  multistep = array(unlist(multistep))
  # Cap forecasts exceeding certain limits
  max_val = max(UnivariateData)
  upper_limit = 0
  if(max_val > 0){
    upper_limit = 1.3*max(UnivariateData)
  }else{
    upper_limit = 0.7*max(UnivariateData)
  }
  min_val = min(UnivariateData)
  lower_limit = 0
  if(min_val > 0){
    lower_limit = 0.7*min(UnivariateData)
  }else{
    lower_limit = 1.3*min(UnivariateData)
  }
  for(i in 1:Horizon){
    forecast = multistep[i]
    if(forecast > upper_limit){
      forecast = upper_limit
    }
    if(forecast < lower_limit){
      forecast = lower_limit
    }
    multistep[i] = forecast
  }
  return(multistep)
}


#
#
#
#
#
