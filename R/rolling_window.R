rolling_window <- function(UnivariateData, CoefficientCombination, Aggregation,
                           Horizon = 2, Window = 3, Method = "r",
                           NumClusters = 1){
  # DESCRIPTION
  # This function computes a rolling forecasting origin for one- or multi-step
  # forecasts with a specific method. Multi step forecasts are computed
  # recursively with the one step forecast method.
  #
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
  # OPTIONAL
  # Horizon          Number indicating horizon for forecast from 1 to horizon.
  # Window           Number indicating how many points are used for cross validation.
  # Method           String indicating which method to use
  #                  Available methods: 'r'  = Autoregression
  #                                     'nn' = Neural Network
  # NumClusters      Number of clusters used for parallel computing.
  #
  # OUTPUT
  # Error[1:Window,1:Horizon]       Numerical Matrix with 'Window' many rows
  #                                 entries indicating one time point with
  #                                 'Horizon' many forecast errors.
  # Forecast[1:Window,1:Horizon]    Numerical Matrix with 'Window' many rows
  #                                 entries indicating one time point with
  #                                 'Horizon' many forecasts.
  #
  # Author: QS, 02/2021
  # Non parallel implementation
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
  if(!is.double(Horizon)){
    message("horizon must be of type double")
    return()
  }
  if(!is.double(Window)){
    message("window_size must be of type double")
    return()
  }
  if(!is.character(Method)){
    message("method must be of type character")
    return()
  }
  #if(!is.double(numClusters)){
  #  message("numClusters must be of type double")
  #  return()
  #}
  if(NumClusters == 1){
    int_total_length  = length(UnivariateData)                        # Length time series
    matError = rbind()
    matForecast = rbind()
    for(i in 1:(Window)){

      int_CFCP = int_total_length - Window - Horizon + i # Current Forecast Position
      dfTrain  = UnivariateData[1:int_CFCP]
      dfTest   = UnivariateData[int_CFCP+1:Horizon]
      forecast = multi_step(UnivariateData = dfTrain, Horizon = Horizon,
                            CoefficientCombination = CoefficientCombination,
                            Aggregation = Aggregation, Method = Method)
      arr_Error = as.numeric(forecast) - dfTest
      #arr_Error = forecast - dfTest
      matError = rbind(matError, arr_Error)
      matForecast = rbind(matForecast, dfTest)
      #print("Training X")
      #print(dfTrain[(int_CFCP-5):int_CFCP],)
      #print("Training Y")
      #print(dfTest)
    }
    Error = matrix(matError, ncol = Horizon, byrow = TRUE)
    Forecast = matrix(matForecast, ncol = Horizon, byrow = TRUE)
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
    if(NumClusters != "max"){
      if(NumClusters > cores){
        message(paste0("There are only ", cores, " cores available. Falling back to maximum available number of cores."))
        #message("There are not enough cores. Note that only maximum of detectCores() - 1 is allowed as maximum.")
        #return()
      }
      if(NumClusters < 1){
        message("Input smaller 1 is not allowed!")
        return()
      }
    }
    cl <- parallel::makeCluster(cores)
    parallel::clusterEvalQ(cl, source("decomposition.R"))
    parallel::clusterEvalQ(cl, source("training.R"))
    parallel::clusterEvalQ(cl, source("multi_step.R"))
    parallel::clusterEvalQ(cl, source("onestep.R"))
    #parallel::clusterEvalQ(cl, source("svm_one_step.R"))
    #parallel::clusterEvalQ(cl, source("xgboost_one_step.R"))
    parallel::clusterEvalQ(cl, source("neuralnet_one_step.R"))
    parallel::clusterEvalQ(cl, source("prediction_scheme.R"))
    parallel::clusterEvalQ(cl, source("regression_one_step.R"))
    #parallel::clusterEvalQ(cl, source("regression_prediction_scheme.R"))
    parallel::clusterEvalQ(cl, source("regression_lsm_optimization.R"))
    #parallel::clusterEvalQ(cl, library("mrf"))
    lst_forecast = parallel::parLapply(cl, 1:Window, rolling_window_single,
                                       data = UnivariateData, ccps = CoefficientCombination,
                                       agg_per_lvl = Aggregation,
                                       horizon = Horizon,
                                       window_size = Window,
                                       method = Method)

    RawVector = unlist(lst_forecast)
    RawMatrix = matrix(RawVector, ncol = Horizon*2, byrow = TRUE)
    Error = RawMatrix[,1:Horizon]
    Forecast = RawMatrix[,(Horizon+1):(2*Horizon)]
    parallel::stopCluster(cl)
  }
  return(list("Error"=Error, "Forecast"=Forecast))
}

