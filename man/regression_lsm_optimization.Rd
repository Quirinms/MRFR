\name{regression_lsm_optimization}
\alias{regression_lsm_optimization}
\title{Least Square Method for Regression}
\description{
This function computes the weights for the autoregression depending
on the given wavelet decomposition.
It uses ordinary least square method for optimizing a linear equation system.
}
\usage{
regression_lsm_optimization(points_in_future, lsmatrix)
}
\arguments{
\item{points_in_future}{n many values of the time series, for which there is an
equation from a prediction scheme.}
\item{lsmatrix}{Matrix carrying predictive equations associated with a specific
value of the time series.}
}
\value{
List of
\item{weights}{Array of weights carrying the solution for a matrix problem,
which was solves with ordinary least squares.}
}
\references{
Aussem, A., Campbell, J., and Murtagh, F. Waveletbased Feature Extraction
and Decomposition Strategies for Financial Forecasting. International Journal of
Computational Intelligence in Finance, 6:5–12, 1998.

Renaud, O., Starck, J.-L., and Murtagh, F. Prediction based on a Multiscale De-
composition. International Journal of Wavelets, Multiresolution and Information
Processing, 1(2):217–232. doi:10.1142/S0219691303000153, 2003.

Murtagh, F., Starck, J.-L., and Renaud, O. On Neuro-Wavelet Modeling. Decision
Support Systems, 37(4):475–484. doi:10.1016/S0167-9236(03)00092-7, 2004.

Renaud, O., Starck, J.-L., and Murtagh, F. Wavelet-based combined Signal Filter-
ing and Prediction. IEEE Transactions on Systems, Man, and Cybernetics, Part
B (Cybernetics), 35(6):1241–1251. doi:10.1109/TSMCB.2005.850182, 2005.
}
\author{
Quirin Stier
}
\examples{
data(AirPassengers)
len_data = length(array(AirPassengers))
ccps = c(1,1,1)
agg_per_lvl = c(2,4)
# Decomposition
dec_res <- decomposition(as.vector(AirPassengers), Aggregation = agg_per_lvl)
# Training
trs_res <- training(dec_res$UnivariateData, dec_res$WaveletCoefficients,
dec_res$SmoothCoefficients,
dec_res$Scales, ccps, agg_per_lvl)
arr_future_points = trs_res$points_in_future
matrix = trs_res$lsmatrix
# Optimization method
weights = regression_lsm_optimization(arr_future_points, matrix)
# Forecast
scheme = prediction_scheme(dec_res$WaveletCoefficients,
dec_res$SmoothCoefficients, ccps, agg_per_lvl)
forecast = weights %*% scheme
}
\keyword{Ordinary least squares}
