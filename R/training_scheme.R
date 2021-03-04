training_scheme <- function(arrValues, wmatrix, rhwtCoeff, scales, ccps, agg_per_lvl){
  # INPUT
  # arrValues[1:n]        Time series with n values
  # wmatrix[scales, n]    Matrix with wavelet coefficients
  # rhwtCoeff[scales, n]  Matrix with smooth coefficients
  # scales                Number of wavelet levels
  # ccps                  Vector with numbers which are associated with wavelet levels.
  #                       The last number is associated with the smooth level.
  #                       Each number determines the number of coefficient used per level.
  #                       The selection follows a specific scheme.
  # agg_per_lvl[]         Vector carrying numbers whose index is associated with the
  #                       wavelet level. The numbers indicate the number of time in
  #                       points used for aggregation from the original time series.
  #
  #
  # OUTPUT
  # points_in_future[1:n] n many values of the time series, for which there
  #                       is an equation from a prediction scheme.
  # lsmatrix[m,n]         Matrix carrying predictive equations associated with a
  #                       specific value of the time series.
  #
  # Author: QS, 02/2021
  if(!is.vector(arrValues)){
    message("arrValues must be of type vector")
    return()
  }
  if(!is.matrix(wmatrix)){
    message("wmatrix must be of type matrix")
    return()
  }
  if(!is.matrix(rhwtCoeff)){
    message("rhwtCoeff must be of type matrix")
    return()
  }
  #if(!is.double(scales)){
  #  message("scales must be of type double")
  #  return()
  #}
  if(!is.vector(ccps)){
    message("ccps must be of type vector")
    return()
  }
  if(!is.vector(agg_per_lvl)){
    message("agg_per_lvl must be of type vector")
    return()
  }
  if(length(ccps) != (scales + 1)){
    warning("Length of ccps must be number of scales from decomposition + 1")
  }
  maxConLen = max(agg_per_lvl)                                          # Range needed for constructing decomposition
  maxReqLen = get_required_training_length(scales, ccps, agg_per_lvl)   # Range needed for constructing model
  startTraining = maxReqLen + maxConLen                                 # Offset needed at the start for choosing coefficients
  len_data = length(arrValues)                                          # Length of swt dec
  intNumEquations = round((len_data - startTraining-1), 0)              # Number of equations available for training
  numberWeights = sum(ccps)                                             # Number of equations needed
  if(numberWeights > intNumEquations){
    stop("There are not enough equations for training. Your time series is too short!")
  }
  lsmatrix = matrix(data = 0, nrow = intNumEquations, ncol = numberWeights)    # Matrix for equations
  for(i in 1:intNumEquations){
    time = startTraining + i
    counter = 1
    for(s in 0:(scales)){
      for(k in 0:(ccps[s+1]-1)){
        if(s != (scales)){
          index = time - k*agg_per_lvl[s+1]
          lsmatrix[i, counter] = wmatrix[s+1, index]
          counter = counter +  1
        }
        else{
          index = time - k*agg_per_lvl[s]
          lsmatrix[i, counter] = rhwtCoeff[s, index]
          counter = counter + 1
        }
      }
    }
  }
  #target_vector = arrValues[(length(arrValues)-intNumEquations+1):length(arrValues)]
  target_vector = as.vector(arrValues[(length(arrValues)-intNumEquations+1):length(arrValues)])
  result = list("points_in_future" = target_vector, "lsmatrix" = lsmatrix)
  return(result)
}

get_required_training_length <- function(scales, ccps, agg_per_lvl = NULL){
  minCut = numeric(length(ccps))
  for(i in 1:length(ccps)){
    if(i != length(ccps)){
      minCut[i] = ccps[i]*agg_per_lvl[i]
    }else{
      minCut[i] = ccps[i]*agg_per_lvl[i-1]
    }
  }
  return(max(minCut))
}


#
#
#
#
#
