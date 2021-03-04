decomposition <- function(data, agg_per_lvl = c(2,4,8,16,32)){
  # INPUT
  # data[1:n]             Vector with n time series values.
  #
  # OPTIONAL
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  #
  # OUTPUT
  # data                  Vector with n time series values.
  # wmatrix[scales, n]    Matrix with wavelet coefficients.
  # rhwtCoeff[scales, n]  Matrix with smooth coefficients.
  # scales                Number of wavelet levels.
  #
  # Author: QS, 02/2021
  if(!is.vector(data)){
    message("Data must be of type vector!")
    return()
  }
  if(!is.vector(agg_per_lvl)){
    message("agg_per_lvl must be of type vector")
    return()
  }
  intLenTS = length(data)
  if(intLenTS < max(agg_per_lvl)){
    warning("The length of the time series or data is not long enough for given aggregation")
  }
  scales = length(agg_per_lvl)
  sprintf("Decomposition is built for %i levels", scales)
  rhwtCoeff = dynamic_aggregation(data, agg_per_lvl)
  wmatrix   = compute_wavelet_matrix(data, intLenTS, scales, rhwtCoeff, wmatrix)
  result = list("data" = data, "wmatrix" = wmatrix, "rhwtCoeff" = rhwtCoeff, "scales" = scales)
  return(result)
}

first_level_aggregation <- function(data, time, number_aggregation_points){
  aggregation_sum = 0
  for(i in 0:(number_aggregation_points-1)){
    aggregation_sum = aggregation_sum + data[time-i]
  }
  averaged_aggregation = aggregation_sum/number_aggregation_points
  return(averaged_aggregation)
}

dynamic_aggregation <- function(data, agg_per_lvl){
  intLenTS = length(data)
  scales = length(agg_per_lvl)
  rhwtCoeff = matrix(0L, nrow = scales, ncol = intLenTS)    # Smooth coefficients
  for(level in 1:scales){
    start = agg_per_lvl[level]
    for(time in 1:intLenTS){
      if(time >= start){                                   # Construct only those time points, for which there is data (depends on scale)
        rhwtCoeff[level, time] = first_level_aggregation(data, time, agg_per_lvl[level])
      }
    }
  }
  return(rhwtCoeff)
}

compute_wavelet_matrix <- function(data, intLenTS, scales, rhwtCoeff, wmatrix){
  wmatrix   = matrix(0L, nrow = scales, ncol = intLenTS)    # Wavelet coefficients (Differences)
  for(time in 1:intLenTS){    # Consider all possible time points
    for(scale in 1:scales){   # Build all difference level
      if(scale == 1){
        wmatrix[scale, time] = data[time] - rhwtCoeff[scale, time]
      }
      else{
        wmatrix[scale, time] = rhwtCoeff[scale-1, time] - rhwtCoeff[scale, time]
      }
    }
  }
  return(wmatrix)
}

#
#
#
#
#
