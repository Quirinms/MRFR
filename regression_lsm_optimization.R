#' Least Square Method for Regression
#'
#' This function computes the weights for the auto regression depending
#' on the given wavelet decomposition.
#'
#' @param points_in_future 1D array with the target values
#' @param lsmatrix 2D matrix containing the wavelet and smooth coefficients
#' chosen for training
#' @return List of parameter with 1D array carrying weights.
#' @examples
#' data(AirPassengers)
#' len_data = length(array(AirPassengers))
#' ccps = c(1,1,1)
#' scales = length(ccps) - 1
#' # Decomposition
#' dec_res <- decomposition(array(AirPassengers), scales = scales)
#' # Training
#' trs_res <- training_scheme(dec_res$data, dec_res$wmatrix, dec_res$rhwtCoeff, dec_res$scales, ccps)
#' arr_future_points = trs_res$points_in_future
#' matrix = trs_res$lsmatrix
#' # Optimization method
#' weights = regression_lsm_optimization(arr_future_points, matrix)
#'@export
regression_lsm_optimization <- function(points_in_future, lsmatrix){
  if(requireNamespace("limSolve", quietly = TRUE)){
    weights = limSolve::Solve(lsmatrix, points_in_future)
  }
  else{
    warning("Package 'limSolve' is not installed")
  }
  result = list("weights" = weights)
  return(result)
}



