#' Generic Training Scheme for wavelet framework
#' 
#' This function computes the input for the training phase required for
#' one step forecasts:
#' A 2D matrix containing wavelet and smooth coefficients and an 1D array
#' containing the target outcome.
#' This computational step is required for all one step forecast procedures
#' contained in this package. Meaning: auto regressive, neural network,
#' xgboost and support vector machine.
#' 
#' @param arrValues 1D array containing signal values
#' @param wmatrix 2D matrix containing wavelet coefficients
#' @param rhwtCoeff 2D matrix containing smooth coefficients
#' @param scales Number of wavelet scales (smooth scale excluded!)
#' @param ccps 1D array containing the number of coefficients chosen per each
#' scale. It accounts for all wavelet levels and the last smooth coefficient
#' level. The last number of the array is the number of coefficients for the
#' smooth level. Therefore the length of this array is larger by one than the
#' parameter "scales".
#' @return List of parameters with 2D matrix with wavelet and smooth
#' coefficients and a 1D array with target values for training.
#' @examples
#' data(AirPassengers)
#' len_data = length(array(AirPassengers))
#' ccps = c(1,1,1)
#' scales = length(ccps) - 1
#' # Decomposition
#' dec_res <- decomposition(array(AirPassengers), scales = 2)
#' # Training
#' trs_res <- training_scheme(dec_res$data, dec_res$wmatrix, dec_res$rhwtCoeff, dec_res$scales, ccps)
#' @export
training_scheme <- function(arrValues, wmatrix, rhwtCoeff, scales, ccps,
                            agg_per_lvl = NULL){
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
  target_vector = arrValues[(length(arrValues)-intNumEquations+1):length(arrValues)]
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





training_scheme_deprecated <- function(arrValues, wmatrix, rhwtCoeff, scales, ccps,
                            dyadic = TRUE, agg_per_lvl = NULL){
  # training = c("linear", "dyadic", "aggregation_per_lvl")
  # Error capturing
  if(length(ccps) != (scales + 1)){
    warning("Length of ccps must be number of scales from decomposition + 1")
  }
  # Highest range of coefficients for constructing model equations
  maxReqLen = get_required_training_length(scales, ccps, dyadic, agg_per_lvl)   # Range needed for constructing model
  maxConLen = get_required_construction_length(scales, dyadic, agg_per_lvl)
  startTraining = maxReqLen + maxConLen
  len_data = length(arrValues)                                                  # Length of swt dec
  # Consider zero coefficients at start
  #startTraining = round(2*maxReqLen, 0)                                        # Cut from start two times
  intNumEquations = round((len_data - startTraining-1), 0)                      # Number of equations available for training
  numberWeights = sum(ccps)                                                     # Number of equations needed
  if(numberWeights > intNumEquations){
    stop("There are not enough equations for training. Your time series is too short!")
  }
  lsmatrix = matrix(data = 0, nrow = intNumEquations, ncol = numberWeights)# Matrix for equations
  for(i in 1:intNumEquations){
    time = startTraining + i
    counter = 1
    for(s in 0:(scales)){
      for(k in 0:(ccps[s+1]-1)){
        if(s != (scales)){
          index = get_index(s, scales, time, k, dyadic, agg_per_lvl)
          lsmatrix[i, counter] = wmatrix[s+1, index]
          counter = counter +  1
        }
        else{
          index = get_index(s, scales, time, k, dyadic, agg_per_lvl)
          lsmatrix[i, counter] = rhwtCoeff[s, index]
          counter = counter + 1
        }
      }
    }
  }
  target_vector = get_target_vector(arrValues, startTraining, intNumEquations, dyadic)
  result = list("points_in_future" = target_vector, "lsmatrix" = lsmatrix)
  return(result)
}

get_target_vector_deprecated <- function(arrValues, startTraining, intNumEquations, dyadic){
  if(dyadic == TRUE){
    return(arrValues[(startTraining+2):length(arrValues)])
  }else{
    return(arrValues[(length(arrValues)-intNumEquations+1):length(arrValues)])
  }
}

get_required_construction_length_deprecated <- function(scales, dyadic, agg_per_lvl = NULL){
  req_cons_len = 1
  if(dyadic == TRUE){
    req_cons_len = ((scales * 2) - 1)
  }else{
    level = length(agg_per_lvl)
    translation = 1
    for(i in 1:level){
      translation = agg_per_lvl[i] * translation
    }
    req_cons_len = translation-1
  }
  return(req_cons_len)
}

get_required_training_length_deprecated <- function(scales, ccps, dyadic = TRUE, agg_per_lvl = NULL){
  start_time = 0
  if(dyadic == TRUE){
    minCut = rep(0, (scales+1))
    for(i in 1:(scales+1)){ # 1, ..., scales + 1
      if(i != (scales+1)){
        minCut[i] = ccps[i]*(2**(i))
      }
      else{
        minCut[i] = ccps[i]*(2**(i-1))
      }
    }
    start_time = max(minCut)
  }else{
    minCut = rep(0, (scales+1))
    for(i in 1:(scales+1)){ # 1, ..., scales + 1
      level = i
      if(i == (scales+1)){
        level = i-1
      }
      translation = 1
      for(j in 1:level){
        translation = agg_per_lvl[j] * translation
      }
      minCut[i] = (ccps[i]*translation)+1
    }
    start_time = max(minCut)
  }
  return(start_time)
}

get_index_deprecated <- function(s, scales, time, k, dyadic, agg_per_lvl = NULL){
  index = 0
  level = s
  if(s != scales){
    level = s + 1
  }
  if(dyadic == TRUE){
    index = time - ((k)*(2**level))
  }else{
    translation = 1
    for(i in 1:level){
      translation = agg_per_lvl[i] * translation
    }
    index = (time - (k * translation))
  }
  return(index)
}
