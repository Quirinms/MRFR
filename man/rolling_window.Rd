\name{rolling_window}
\alias{rolling_window}
\title{Rolling forecasting origin for Multiresolution Forecasts}
\description{
This function creates a rolling forecasting origin for a specifc one step
forecast method. Multi step forecasts are computed recursively with the one step
forecast method.
Thus, h-step forecast for h = 1,..., horizon for window_size many steps can be
computed.
The forecasting method can be an autoregression or a neural network (multilayer
perceptron).
The ccps parameter controls the number of coefficients chosen for each wavelet
and smooth part level individually.
The numClusters parameter determines the number of cluster used for parallel
computation. numClusters = 1 performs a non parallel version.
numClusters is constrained by the maximum number of clusters available minus one
to prevent the machine to be overchallenged.
}
\usage{
rolling_window(
  data,
  ccps,
  agg_per_lvl,
  horizon = 14,
  window_size = 365,
  method = "r",
  numClusters = 1
)
}
\arguments{
\item{data}{Vector with n time series values.}
\item{ccps}{Vector with numbers which are associated with wavelet and smooth
levels from decomposition.
The last number is associated with the smooth part level.
The preceding numbers are associated with the wavelet levels which are ordered
increasingly.
Each number determines the number of coefficient used per level.
The coefficient selection follows a fixed scheme.}
\item{agg_per_lvl}{Vector carrying numbers whose index is associated with the
wavelet level. The numbers indicate the number of time in points used for
aggregation from the original time series.}
\item{horizon}{Number indicating horizon for forecast from 1 to horizon.}
\item{window_size}{Number indicating how many points are used to create cross
validation.}
\item{method}{String indicating which method to use (r = Autoregression, nn =
Neural Network).}
\item{numClusters}{Number of clusters used for parallel computing.}
}
\value{
List of
\item{Rolling_Window}{Matrix with m many forecasts with horizon from 1 to n}
}
\references{
Hyndman, R. and Athanasopoulos, G. Forecasting: principles and practice.
OTexts, 3 edition. 2018.
}
\author{
Quirin Stier
}
\examples{
data(AirPassengers)
plot(AirPassengers, type = "l", col = "black")
result = decomposition(as.vector(array(AirPassengers)), agg_per_lvl = c(2,4))
plot(result$rhwtCoeff[2,4:length(result$rhwtCoeff[2,])], type = "l", col = "blue")
lines(array(AirPassengers)[4:length(result$rhwtCoeff[2,])], col = "black")
}
\keyword{Rolling forecasting origin}
