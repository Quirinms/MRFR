#' Decomposition
#'
#' This function decomposes a signal in its wavelet coefficients and its the
#' smooth coefficients. The number of scales/levels are controllable.
#'
#' @param data One dimensional array of signal values
#' @param scales Number of wavelet scales (smooth scale excluded!)
#' @return List of parameters with 1D Array with Original Values
#' 2D Matrix with Wavelet Coefficients, 2D Matrix with Smooth Coefficient
#' and number of scales
#'@examples
#'data(AirPassengers)
#'plot(AirPassengers, type = "l", col = "black")
#'result = decomposition(array(AirPassengers), scales = 2)
#'plot(result$rhwtCoeff[2,4:length(result$rhwtCoeff[2,])], type = "l", col = "blue")
#'lines(array(AirPassengers)[4:length(result$rhwtCoeff[2,])], col = "black")
#'
#'@export
decomposition <- function(data, agg_per_lvl = c(2,4,8,16,32)){
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



decomposition_deprecated <- function(data, scales = 2, dyadic = TRUE, agg_per_lvl = c(7,2,2,2,3,2)){
  intLenTS = length(data)
  if(dyadic == TRUE){
    sprintf("Decomposition is built for %i levels", scales)
    # Requirement for training
    numPoints   = scales-1                   # Number of points lost for building each level
    numMaxRange = 2**scales                  # Number of points needed for training biggest scale
    rangeScales = numMaxRange + numPoints    # Required number of swt coefficients on lowest scale
    if(intLenTS < rangeScales){
      warning("Length of Data is not sufficient for given level")
    }
    rhwtCoeff = dyadic_aggregation(data, scales, intLenTS)
    wmatrix   = compute_wavelet_matrix(data, intLenTS, scales, rhwtCoeff, wmatrix)
    # Following equation must be true
    # sum(wmatrix[:, idx]) + rhwtCoeff[scales-1, idx] = signal[idx] - idx = 0, ..., len(data)-1
    result = list("data" = data, "wmatrix" = wmatrix, "rhwtCoeff" = rhwtCoeff, "scales" = scales)
    return(result)
  }else{
    scales = length(agg_per_lvl)
    sprintf("Decomposition is built for %i levels", scales)
    construction_size = 1
    for(count in agg_per_lvl){
      construction_size = construction_size * count
    }
    # Requirement for training
    #numPoints   = scales-1                   # Number of points lost for building each level
    #numMaxRange = 2**scales                  # Number of points needed for training biggest scale
    #rangeScales = numMaxRange + numPoints    # Required number of swt coefficients on lowest scale
    if(intLenTS < construction_size){
      warning("Length of Data is not sufficient to construct custom aggregation")
    }
    rhwtCoeff = custom_aggregation(data, scales, intLenTS, agg_per_lvl)
    wmatrix   = compute_wavelet_matrix(data, intLenTS, scales, rhwtCoeff, wmatrix)
    # Following equation must be true
    # sum(wmatrix[:, idx]) + rhwtCoeff[scales-1, idx] = signal[idx] - idx = 0, ..., len(data)-1
    result = list("data" = data, "wmatrix" = wmatrix, "rhwtCoeff" = rhwtCoeff, "scales" = scales)
    return(result)
  }
  return()
}


dyadic_aggregation_deprecated <- function(data, scales, intLenTS){
  arrLB = numeric(scales)             # Array for indicating sizes needed for each level
  for(i in 0:(scales-1)){
    arrLB[i+1] = 2**i
  }
  rhwtCoeff = matrix(0L, nrow = scales, ncol = intLenTS)    # Smooth coefficients
  for(level in 1:scales){
    for(time in 1:intLenTS){
      if(time >= sum(arrLB[1:level])+1){                    # Use only the time points, for which there is data (depends on scale)
        idx = time - (2**(level-1))
        if(level == 1){                                     # Use the original time series for the first level of constructing
          rhwtCoeff[level, time] = 0.5*(data[idx] + data[time])
        }
        else{                                               # Use the coefficients from previous levels for further levels
          rhwtCoeff[level, time] = 0.5*(rhwtCoeff[level-1, idx] + rhwtCoeff[level-1, time])
        }
      }
    }
  }
  return(rhwtCoeff)
}

custom_aggregation_deprecated <- function(data, scales, intLenTS, agg_per_lvl){
  arrLB = numeric(scales)             # Array for indicating sizes needed for each level
  for(i in 0:(scales-1)){
    arrLB[i+1] = scales**i
  }
  rhwtCoeff = matrix(0L, nrow = scales, ncol = intLenTS)    # Smooth coefficients
  for(level in 1:scales){
    for(time in 1:intLenTS){
      limit = 1
      for(tc in agg_per_lvl[1:(level)]){                    # Example: agg_per_lvl = c(7, 2, 2, 2,  2)
        limit = limit * tc                                  #          Level index =   7,14,28,56,112
      }
      if(time >= limit){
        agg = agg_per_lvl[level]
        weight = 1/agg
        sum_points = 0
        if(level == 1){
          for(aggregation in 0:(agg-1)){
            sum_points = sum_points + data[time-aggregation]
          }
          rhwtCoeff[level, time] = weight*sum_points
        }
        else{
          translation = 1
          for(i in 1:(level-1)){
            translation = agg_per_lvl[i] * translation
          }
          for(aggregation in 0:(agg-1)){
            sum_points = sum_points + rhwtCoeff[(level-1), (time-(aggregation*translation))]
          }
          rhwtCoeff[level, time] = weight * sum_points
        }
      }
    }
  }
  return(rhwtCoeff)
}



#
#
#
#
#