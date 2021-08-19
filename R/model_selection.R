model_selection <- function(UnivariateData, Aggregation, Horizon = 14,
                            Window = 365, Method = "r", crit = "AIC",
                            itermax = 1, lower_limit = 1, upper_limit = 2,
                            NumClusters = 1){
  # INPUT
  # UnivariateData[1:n]      Numerical vector with n values
  # Aggregation[1:Scales]    Numerical vector carrying numbers whose index is
  #                          associated with the wavelet level. The numbers
  #                          indicate the number of time in points used for
  #                          aggregation from the original time series.
  #
  #
  # OPTIONAL
  # Horizon          Number indicating horizon for forecast from 1 to horizon.
  # Window           Number indicating how many points are used for cross validation.
  # Method           String indicating which method to use
  #                  Available methods: 'r'  = Autoregression
  #                                     'nn' = Neural Network
  # crit             String indicating which criterion to use.
  #                  Available criterion: AIC = Akaikes Information Criterion
  #                                       MRE = Mean Root Error
  # itermax          Number of iterations for evolutionary optimization method.
  # lower_limit      Lower limit for coefficients selected for each level.
  # upper_limit      Higher limit for coefficients selected for each level.
  # NumClusters      Number of clusters used for parallel computing.
  #
  #
  # OUTPUT
  # Error[1:Window,1:Horizon]       Numerical Matrix with 'Window' many rows
  #                                 entries indicating one time point with
  #                                 'Horizon' many forecast errors.
  # Best[1:Scales+1]                Numerical vector with integers associated
  #                                 with the best found number of coefficients
  #                                 per wavelet scale (1:Scales) and number of
  #                                 coefficients for the smooth approximation
  #                                 level in the last entry.
  #
  # Author: QS, 02/2021
  if(!is.vector(UnivariateData)){
    message("Data must be of type vector")
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
  if(!is.character(crit)){
    message("crit must be of type character")
    return()
  }
  if(!is.double(itermax)){
    message("itermax must be of type double")
    return()
  }
  if(!is.double(lower_limit)){
    message("lower_limit must be of type double")
    return()
  }
  if(!is.double(upper_limit)){
    message("upper_limit must be of type double")
    return()
  }
  #if(!is.double(numClusters)){
  #  message("numClusters must be of type double")
  #  return()
  #}
  len_data = length(UnivariateData)
  scales = length(Aggregation)
  len_ccps = scales + 1
  lower <- rep(lower_limit,len_ccps)
  upper <- rep(upper_limit,len_ccps)

  if (!requireNamespace('DEoptim', quietly = TRUE)) {
    message(
      "Package DEoptim is missing in function evolutionary_optim.
      No computations are performed.
      Please install the packages which are defined in 'Suggests'"
    )
    return()
  }
  res = DEoptim::DEoptim(crit_rolling_window, lower, upper, fnMap = round,
                         UnivariateData = UnivariateData, Aggregation = Aggregation,
                         Window = Window,
                         Horizon = Horizon, Method = Method, NumClusters = NumClusters,
                         crit = crit,
                         control = DEoptim::DEoptim.control(itermax = itermax))
  ccps = as.numeric(res$optim$bestmem)
  res = rolling_window(UnivariateData,
                       CoefficientCombination = ccps,
                       Aggregation = Aggregation,
                       Horizon = Horizon,
                       Window = Window,
                       Method = Method,
                       NumClusters = NumClusters)
  mat_error = res$Error
  return(list("Error"=mat_error, "Best"=ccps))
}

crit_rolling_window <- function(CoefficientCombination, Aggregation, UnivariateData, Window, Horizon,
                                Method, crit = "AIC", NumClusters = "max"){
  res = rolling_window(UnivariateData = UnivariateData,
                       CoefficientCombination = CoefficientCombination,
                       Aggregation = Aggregation,
                       Horizon = Horizon,
                       Window = Window,
                       Method = Method,
                       NumClusters = NumClusters)
  mat_error = res$Error
  MAE = sum(abs(mat_error))/(Window*Horizon)
  AIC = length(UnivariateData) * log(MAE^2) + 2*(sum(CoefficientCombination)+1)
  if(crit == "MAE"){
    return(MAE)
  }
  if(crit == "AIC"){
    return(AIC)
  }
  if(crit == "MRE"){
    SRE = (sum(sqrt(as.complex(((-1)*mat_error)))))
    MRE = SRE/(Window*Horizon)
    a = Re(MRE)
    b = Im(MRE)
    magnitude = abs(MRE)
    gamma = -1
    if(a > 0){
      gamma = atan(b/a)
    }else{
      if((a == 0) && (b > 0)){
        gamma = pi/2
      }
      else{
        if((a==0) && (b==0)){
          gamma = 0
        }
      }
    }
    kappa = 1 - ((4*gamma)/pi)
    qm_vec = c(magnitude, kappa)
    pareto = norm(as.matrix(qm_vec))
    return(pareto)
  }
}

write_rolling_window <- function(mat_table, str_name){
  utils::write.table(mat_table, str_name, row.names = FALSE)
}

get_savings_string <- function(method, horizon, data_name, dir_name, ccps, agg_per_lvl){
  date = substring(Sys.time(), 1, 10)
  str_dir_data_method = paste(dir_name, date, "_", data_name, "_", method, "_rw_", sep="")
  str_step = paste(as.character(horizon), "_step_", sep = "")
  str_ccps = paste(paste(ccps, collapse = "_"), "_ccps_", sep="")
  str_aggregation = paste(paste(agg_per_lvl, collapse = "_"), "_agg", sep="")
  str_final = paste(str_dir_data_method, str_step, str_ccps, str_aggregation, ".csv", sep="")
  return(str_final)
}

#
#
#
#
#
