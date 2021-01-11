#' One Step Forecast with Neural Network
#'
#' This function creates a one step forecast using a multi layer perceptron
#' with one hidden Layer. The number of input is the number of coefficients
#' chosen with the parameter ccps.
#'
#' @param data One dimensional array of signal values
#' @param ccps 1D array containing the number of coefficients chosen per each
#' scale. It accounts for all wavelet levels and the last smooth coefficient
#' level. The last number of the array is the number of coefficients for the
#' smooth level. Therefore the length of this array is larger by one than the
#' parameter "scales".
#' @return List of parameter with forecast.
#'
#' @examples
#' data(AirPassengers)
#' len_data = length(array(AirPassengers))
#' result = neuralnet_one_step(array(AirPassengers)[1:(len_data-1)], c(1,1,1))
#' forecast = result$forecast
#' true_value = array(AirPassengers)[len_data]
#' error = true_value - forecast
#'
#' @export
neuralnet_one_step <- function(data, ccps, agg_per_lvl){
  if(is.null(agg_per_lvl)){
    stop("Parameter agg_per_lvl is not defined")
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
  if(requireNamespace("monmlp", quietly = TRUE)){
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
  }
  else{
    stop("Package 'monmlp' is not installed")
  }
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
