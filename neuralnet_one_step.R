neuralnet_one_step <- function(data, ccps, agg_per_lvl){
  # INPUT
  # data[1:n]             Vector with n time series values.
  # 
  # OPTIONAL
  # ccps                  Vector with numbers which are associated with wavelet levels.
  #                       The last number is associated with the smooth level.
  #                       Each number determines the number of coefficient used per level.
  #                       The selection follows a specific scheme.
  # agg_per_lvl           Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  #
  # OUTPUT
  # forecast              Double one step forecast
  # Author: QS, 02/2021
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

  if((length(ccps)-1) != length(agg_per_lvl)){
    stop("Length of ccps must be longer than that of agg_per_level by one.
    Parameter ccps needs one coefficient per wavelet level and one
         coefficient for the final smooth level.
         Parameter agg_per_level defines the number of levels of the decomposition.")
  }
  dec_res <- decomposition(data, agg_per_lvl)                              # Decomposition
  trs_res <- training_scheme(dec_res$data,                                 # Create training scheme
                             dec_res$wmatrix,
                             dec_res$rhwtCoeff,
                             dec_res$scales,
                             ccps, agg_per_lvl)
  arr_future_points = trs_res$points_in_future
  lm_matrix = trs_res$lsmatrix
  num_feature = dim(lm_matrix)[2]
  num_trainingpoint = dim(lm_matrix)[1]
  forecast_scheme = neuralnet_prediction_scheme(dec_res$wmatrix,           # Create forecast scheme
                                                dec_res$rhwtCoeff,
                                                ccps, agg_per_lvl)
  arr_pred = array(unlist(forecast_scheme), dim=c(1,num_feature))
  matrix = rbind(lm_matrix, arr_pred)
  set.seed(8675309)
  if (!requireNamespace('monmlp', quietly = TRUE)) {
    message(
      "Package monmlp is missing in function neuralnet_one_step
      No computations are performed.
      Please install the packages which are defined in 'Suggests'"
    )
    return()
  }

  #model = monmlp::monmlp.fit(x = lm_matrix, y = as.matrix(arr_future_points),
  #                           hidden1 = 8, hidden2 = 0, iter.max = 500,
  #                           Th = my_tan, To = my_tan,
  #                           Th.prime = my_dtan, To.prime = my_dtan,
  #                           method ="BFGS")

  #model = monmlp::monmlp.fit(x = lm_matrix, y = as.matrix(arr_future_points),
  #                           hidden1 = 8, hidden2 = 0, iter.max = 500,
  #                           Th = my_relu, To = my_relu,
  #                           Th.prime = my_drelu, To.prime = my_drelu,
  #                           method ="BFGS")

  model = monmlp::monmlp.fit(x = lm_matrix, y = as.matrix(arr_future_points),
                             hidden1 = 8, hidden2 = 0, scale.y	= TRUE,
                             iter.max = 500, method ="BFGS")
  res <- monmlp::monmlp.predict(x = matrix, weights = model)
  len_res = length(res)
  forecast = res[len_res]
  result = list("forecast" = forecast)                                     # Return
  return(result)
}


my_relu <- function(x){
  return(max(0, x))
}
my_drelu <- function(x){
  if(x >= 0){
    return(1)
  }else{
    return(0)
  }
}

my_tan <- function(x){
  a = exp(x) - exp(-x)
  b = exp(x) + exp(-x)
  return(a/b)
}
my_dtan <- function(x){
  return(1-(tan(x)**2))
}


#
#
#
#
#
