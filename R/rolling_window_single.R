rolling_window_single <- function(i, data, ccps, agg_per_lvl, horizon = 14,
                                  window_size = 365, method = "r"){
  int_total_length  = length(data)                        # Length time series
  int_CFCP = int_total_length - window_size - horizon + i # Current Forecast Position
  dfTrain  = data[1:int_CFCP]
  dfTest   = data[int_CFCP+1:horizon]
  forecast = multi_step(data = dfTrain, steps = horizon, ccps = ccps,
                        method = method, agg_per_lvl = agg_per_lvl)
  arr_Error = as.numeric(forecast$forecast) - dfTest
  return(list("Error" = arr_Error))
}
