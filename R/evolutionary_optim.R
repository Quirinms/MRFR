evolutionary_optim <- function(data,
                               agg_per_lvl,
                               data_name = "data",
                               dir_name = "",
                               horizon = 14,
                               window_size = 365,
                               method = "r",
                               crit = "AIC",
                               itermax = 1,
                               lower_limit = 1,
                               upper_limit = 2,
                               write = FALSE,
                               numClusters = 1){
  # INPUT
  # data[1:n]             Vector with n values
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  #
  # OPTIONAL
  # data_name             String with name for output file.
  # dir_name              String with directory for output file.
  # horizon               Number indicating horizon for forecast from 1 to horizon.
  # window_size           Number indicating how many points are used for cross validation.
  # method                String indicating which method to use (r = Autoregression, nn = Neural Network)
  # crit                  String indicating which criterion to use:
  #                       (AIC = Akaikes Information Criterion)
  #                       (MRE = Mean Root Error).
  # itermax               Number of iterations for evolutionary optimization method.
  # lower_limit           Lower limit for coefficients selected for each level.
  # upper_limit           Higher limit for coefficients selected for each level.
  # numClusters           Number of clusters used for parallel computing.
  #
  #
  # OUTPUT
  #
  # Author: QS, 02/2021
  if(!is.vector(data)){
    message("Data must be of type vector")
    return()
  }
  if(!is.vector(agg_per_lvl)){
    message("agg_per_lvl must be of type vector")
    return()
  }
  if(!is.character(data_name)){
    message("data_name must be of type character")
    return()
  }
  if(!is.character(dir_name)){
    message("dir_name must be of type character")
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
  if(!is.logical(write)){
    message("write must be of type logical")
    return()
  }
  #if(!is.double(numClusters)){
  #  message("numClusters must be of type double")
  #  return()
  #}
  len_data = length(data)
  scales = length(agg_per_lvl)
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
                         data = data, agg_per_lvl = agg_per_lvl,
                         window_size = window_size,
                         horizon = horizon, method = method,
                         crit = crit, numClusters = numClusters, control = DEoptim::DEoptim.control(itermax = itermax))
  ccps = as.numeric(res$optim$bestmem)
  mat_error = rolling_window(data,
                             ccps = ccps,
                             agg_per_lvl = agg_per_lvl,
                             horizon = horizon,
                             window_size = window_size,
                             method = method)
  mat_error = mat_error$Rolling_Window
  if(write){
    str_Error = get_savings_string(method = method,
                                   horizon = horizon,
                                   data_name = data_name,
                                   dir_name = dir_name,
                                   ccps = ccps,
                                   agg_per_lvl = agg_per_lvl)
    write_rolling_window(mat_error, str_Error)
  }
  return(mat_error)
}

crit_rolling_window <- function(ccps, agg_per_lvl, data, window_size, horizon,
                                method, crit = "AIC", numClusters = "max"){
  mat_error = rolling_window(data = data,
                             ccps = ccps,
                             agg_per_lvl = agg_per_lvl,
                             horizon = horizon,
                             window_size = window_size,
                             method = method,
                             numClusters = numClusters)
  mat_error = mat_error$Rolling_Window
  MAE = sum(abs(mat_error))/(365*14)
  AIC = length(data) * log(MAE^2) + 2*(sum(ccps)+1)
  if(crit == "MAE"){
    return(MAE)
  }
  if(crit == "AIC"){
    return(AIC)
  }
  if(crit == "MRE"){
    SRE = (sum(sqrt(as.complex(((-1)*mat_error)))))
    MRE = SRE/(window_size*horizon)
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
