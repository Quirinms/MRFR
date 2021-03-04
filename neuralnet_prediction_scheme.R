neuralnet_prediction_scheme <- function(wmatrix, rhwtCoeff, ccps, agg_per_lvl){
  # INPUT
  # wmatrix[scales, n]    Matrix with wavelet coefficients
  # rhwtCoeff[scales, n]  Matrix with smooth coefficients
  # ccps                  Vector with numbers which are associated with wavelet levels.
  #                       The last number is associated with the smooth level.
  #                       Each number determines the number of coefficient used per level.
  #                       The selection follows a specific scheme.
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  #
  # OUTPUT
  # future_point          Double carrying one step forecast
  #
  # Author: QS, 02/2021
  
  #if(!is.matrix(wmatrix)){
  #  message("wmatrix must be of type matrix")
  #  return()
  #}
  #if(!is.matrix(rhwtCoeff)){
  #  message("rhwtCoeff must be of type matrix")
  #  return()
  #}
  #if(!is.vector(ccps)){
  #  message("ccps must be of type vector")
  #  return()
  #}
  #if(!is.vector(agg_per_lvl)){
  #  message("agg_per_lvl must be of type vector")
  #  return()
  #}
  
  scales = length(agg_per_lvl)
  time = dim(wmatrix)[2]
  future_point = c()
  counter = 1
  for(s in 0:scales){
    for(k in 0:(ccps[s+1]-1)){
      if(s != scales){
        index = time - k*agg_per_lvl[s+1]
        future_point = c(future_point, wmatrix[s+1, index])
        counter = counter + 1
      }
      else{
        index = time - k*agg_per_lvl[s]
        future_point = c(future_point, rhwtCoeff[s, index])
        counter = counter + 1
      }
    }
  }
  result = list("future_point" = future_point)
  return(result)
}

#
#
#
#
#
