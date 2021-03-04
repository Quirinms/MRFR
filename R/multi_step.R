multi_step <- function(data, steps, ccps, agg_per_lvl, method = "r"){
  # INPUT
  # data[1:n]             Vector with n time series values.
  #
  # OPTIONAL
  # steps                 Number indicating horizon for forecast from 1 to horizon
  # ccps                  Vector with numbers which are associated with wavelet levels.
  #                       The last number is associated with the smooth level.
  #                       Each number determines the number of coefficient used per level.
  #                       The selection follows a specific scheme.
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  # method                String indicating which method to use (r = Autoregression, nn = Neural Network)
  #
  #
  # OUTPUT
  # forecast              Array carrying the forecasts, where the index of the entry
  #                       is associated with the horizon of the forecast.
  # Author: QS, 02/2021
  if(!is.vector(data)){
    message("Data must be of type vector")
    return()
  }
  if(!is.double(steps)){
    message("steps must be of type double")
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
  if(!is.character(method)){
    message("method must be of type character")
    return()
  }
  multistep = c()
  for(i in 1:steps){
    forecast = 0
    switch(method,
           "r"={
             forecast = regression_one_step(data, ccps, agg_per_lvl)
             forecast = forecast$forecast
           },
           "nn"={
             forecast = neuralnet_one_step(data, ccps, agg_per_lvl)
             forecast = forecast$forecast
           },
           {
             message("Use 'r' for autoregression and 'nn' for neural network
              (multilayer perceptron)")
             return()
           }
    )
    #if(method == "xgboost"){
    #  forecast = xgboost_one_step(data, ccps, agg_per_lvl)
    #  forecast = forecast$forecast
    #}
    #if(method == "svm"){
    #  forecast = svm_one_step(data, ccps, agg_per_lvl)
    #  forecast = forecast$forecast
    #}

    multistep = c(multistep, forecast)
    data = as.vector(unlist(c(data, forecast)))
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
