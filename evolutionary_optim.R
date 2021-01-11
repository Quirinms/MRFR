#' Model Selection for Multiresolution Forecasts
#'
#' This function computes the errors of rolling windows for multiple settings
#' for a multi layer perceptron with one hidden Layer or a regression. The
#' procedures saves csv containing the different settings. The settings are
#' saved as String. The AIC of the rolling window result are computed.
#' Each setting is assigned one AIC value.
#'
#' @param data One dimensional array of signal values
#' @param data_name String indicating file name in order to assign a name for
#' the resulting csv file to save.
#' @param horizon Integer indicating the number of steps to forecast ahead
#' @param window_size Integer indicating the number of days, for which a forecast
#' should be evaluated. The procedure uses automatically the last part of the time
#' series
#' @param method Character indicating Regression ("r") or the Neural Network ("n")
#' @param level Integer indicating how many decomposition levels to produce. The
#' settings for the parameters is adjusted to exactly this level.
#' @return None
#' @export
evolutionary_optim <- function(data, data_name = "data", dir_name = "",
                               agg_per_lvl,
                               horizon = 14,
                               window_size = 365,
                               method = "r",
                               crit = "AIC", itermax = 1){
  len_data = length(data)
  if(len_data < (3*365)){
    warning("The length of time series is not long enough")
  }
  scales = length(agg_per_lvl)
  len_ccps = scales + 1
  lower <- rep(1,len_ccps)
  upper <- rep(7,len_ccps)
  if(requireNamespace("DEoptim", quietly = TRUE)){
    res = DEoptim::DEoptim(crit_rolling_window, lower, upper, fnMap = round,
                           data = data, agg_per_lvl = agg_per_lvl,
                           window_size = window_size,
                           horizon = horizon, method = method,
                           crit = crit, control = DEoptim::DEoptim.control(itermax = itermax))
  }
  else{
    warning("Package 'DEoptim' is not installed")
  }
  ccps = as.numeric(res$optim$bestmem)
  mat_error = rolling_window_parallel(data,
                                      ccps = ccps,
                                      agg_per_lvl = agg_per_lvl,
                                      horizon = horizon,
                                      window_size = window_size,
                                      method = method)
  mat_error = mat_error$Rolling_Window
  str_Error = get_savings_string(method = method,
                                 horizon = horizon,
                                 data_name = data_name,
                                 dir_name = dir_name,
                                 ccps = ccps,
                                 agg_per_lvl = agg_per_lvl)
  save_results(mat_error, str_Error)
}

crit_rolling_window <- function(ccps, agg_per_lvl, data, window_size, horizon,
                                method, crit = "AIC"){
  mat_error = rolling_window_parallel(data = data,
                                      ccps = ccps,
                                      agg_per_lvl = agg_per_lvl,
                                      horizon = horizon,
                                      window_size = window_size,
                                      method = method)
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

save_results <- function(df_table, str_name){
  if(requireNamespace("utils", quietly = TRUE)){
    utils::write.table(df_table, str_name, row.names = FALSE)
  }
  else{
    warning("Package 'utils' is not installed")
  }
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