wavelet_decomposition <- function(UnivariateData, Aggregation = c(2,4,8,16,32)){
  # DESCRIPTION
  # This function decomposes a time series in its wavelet and smooth
  # coefficients using the redundant Haar wavelet transform.
  #
  # INPUT
  # UnivariateData[1:n]           Numerical vector with n time series values.
  #
  # OPTIONAL
  # Aggregation[1:Scales]         Numerical vector carrying numbers whose index
  #                               is associated with the wavelet level. The
  #                               numbers indicate the number of values used for
  #                               aggregation from the original time series.
  #
  # OUTPUT
  # UnivariateData[1:n]               Numerical vector with n time series values.
  # WaveletCoefficients[Scales, n]    Numerical matrix with wavelet coefficients.
  # SmoothCoefficients[Scales, n]     Numerical matrix with smooth coefficients.
  # Scales                            Number of wavelet levels.
  #
  # DETAILS
  #
  #
  # Author: QS, 02/2021
  if(!is.vector(UnivariateData)){
    message("Data must be of type vector!")
    return()
  }
  if(!is.vector(Aggregation)){
    message("agg_per_lvl must be of type vector")
    return()
  }
  intLenTS = length(UnivariateData)
  if(intLenTS < max(Aggregation)){
    warning("The length of the time series or data is not long enough for given aggregation")
  }
  Scales = length(Aggregation)
  sprintf("Decomposition is built for %i levels", Scales)
  SmoothCoefficients = dynamic_aggregation(UnivariateData, Aggregation)
  WaveletCoefficients   = compute_wavelet_matrix(UnivariateData, intLenTS, Scales, SmoothCoefficients, WaveletCoefficients)
  return(list("UnivariateData" = UnivariateData,
              "WaveletCoefficients" = WaveletCoefficients,
              "SmoothCoefficients" = SmoothCoefficients,
              "Scales" = Scales))
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
