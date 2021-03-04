rolling_window <- function(data, ccps, agg_per_lvl, horizon = 14,
                           window_size = 365, method = "r", numClusters = 1){
  # INPUT
  # data[1:n]             Vector with n time series values.
  # ccps[]                Vector with numbers which are associated with wavelet levels.
  #                       The last number is associated with the smooth level.
  #                       Each number determines the number of coefficient used per level.
  #                       The selection follows a specific scheme.
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  # 
  # OPTIONAL
  # horizon               Number indicating horizon for forecast from 1 to horizon.
  # window_size           Number indicating how many points are used to create cross validation.
  # method                String indicating which method to use (r = Autoregression, nn = Neural Network).
  # numClusters           Number of clusters used for parallel computing.
  # 
  #
  # OUTPUT
  # Rolling_Window[m,n]   Matrix with m many forecasts with horizon from 1 to n
  # 
  # Author: QS, 02/2021
  # Non parallel implementation
  if(!is.vector(data)){
    message("Data must be of type vector")
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
  if(!is.double(horizon)){
    message("horizon must be of type double")
    return()
  }
  if(!is.double(window_size)){
    message("window_size must be of type double")
    return()
  }
  if(!is.character(method)){
    message("method must be of type character")
    return()
  }
  #if(!is.double(numClusters)){
  #  message("numClusters must be of type double")
  #  return()
  #}
  if(numClusters == 1){
    int_total_length  = length(data)                        # Length time series
    mat_Error_Forecast = rbind()
    for(i in 1:(window_size)){

      int_CFCP = int_total_length - window_size - horizon + i # Current Forecast Position
      dfTrain  = data[1:int_CFCP]
      dfTest   = data[int_CFCP+1:horizon]
      forecast = multi_step(data = dfTrain, steps = horizon, ccps = ccps,
                            method = method, agg_per_lvl = agg_per_lvl)
      arr_Error = as.numeric(forecast$forecast) - dfTest
      #arr_Error = forecast - dfTest
      mat_Error_Forecast = rbind(mat_Error_Forecast, arr_Error)

    }
    matrix_forecast = matrix(mat_Error_Forecast, ncol = horizon, byrow = TRUE)
  }
  else{ # Parallel implementation

    if (!requireNamespace('parallel', quietly = TRUE)) {
      message(
        "Package parallel is missing in function rolling_window
      No computations are performed.
      Please install the packages which are defined in 'Suggests'"
      )
      return()
    }

    available_cores = parallel::detectCores()    # Number of cores available
    cores = available_cores[1]-1                 # Do not overuse => av_cores-1
    if(numClusters != "max"){
      if(numClusters > cores){
        message("There are not enough cores. Note that only maximum of detectCores() - 1 is allowed as maximum.")
        return()
      }
      if(numClusters < 1){
        message("Input smaller 1 is not allowed!")
        return()
      }
    }
    cl <- parallel::makeCluster(cores)
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
  }
  result = list("Rolling_Window" = matrix_forecast)    # Return
  return(result)
}

