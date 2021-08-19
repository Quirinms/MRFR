rolling_window_single <- function(i, data, ccps, agg_per_lvl, horizon = 14,
                                  window_size = 365, method = "r"){
  int_total_length  = length(data)                        # Length time series
  int_CFCP = int_total_length - window_size - horizon + i # Current Forecast Position
  dfTrain  = data[1:int_CFCP]
  dfTest   = data[int_CFCP+1:horizon]
  forecast = multi_step(UnivariateData = dfTrain,
                        Horizon = horizon,
                        CoefficientCombination = ccps,
                        Aggregation = agg_per_lvl,
                        Method = method)
  arr_Error = as.numeric(forecast) - dfTest
  return(list("Error" = arr_Error, "Forecast"=dfTest))
}
