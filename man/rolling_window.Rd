\name{rolling_window}
\alias{rolling_window}
\title{Rolling forecasting origin for Multiresolution Forecasts}
\description{
This function computes a rolling forecasting origin for one- or multi-step
forecasts with a specific method. Multi step forecasts are computed recursively
with the one step forecast method.
}
\usage{
rolling_window(UnivariateData, CoefficientCombination, Aggregation, Horizon = 2,
Window = 3, Method = "r", NumClusters = 1)
}
\arguments{
\item{UnivariateData}{[1:n] Numerical vector with n values.}
\item{CoefficientCombination}{[1:Scales+1] Numerical vector with numbers which
are associated with wavelet levels. The last number is associated with the
smooth level. Each number determines the number of coefficient used per level.
The selection follows a specific scheme.}
\item{Aggregation}{[1:Scales] Numerical vector carrying numbers whose index is
associated with the wavelet level. The numbers indicate the number of time in
points used for aggregation from the original time series.}
\item{Horizon}{Number indicating horizon for forecast from 1 to horizon.}
\item{Window}{Number indicating how many points are used for cross validation.}
\item{Method}{String indicating which method to use. Available methods: 'r'  =
Autoregression. 'nn' = Neural Network.}
\item{NumClusters}{Number of clusters used for parallel computing.}
}
\value{
List of
\item{Error}{[1:Window,1:Horizon] Numerical Matrix with 'Window' many row
entries indicating one time point with 'Horizon' many forecast errors.}
\item{Forecast}{[1:Window,1:Horizon] Numerical Matrix with 'Window' many row
entries indicating one time point with 'Horizon' many forecasts.}
}
\details{
Thus, h-step forecast for h = 1,..., horizon for window_size many steps can be
computed.
The forecasting method can be an autoregression or a neural network (multilayer
perceptron).
The CoefficientCombination parameter controls the number of coefficients chosen for each wavelet
and smooth part level individually.
The NumClusters parameter determines the number of cluster used for parallel
computation. NumClusters = 1 performs a non parallel version.
NumClusters is constrained by the maximum number of clusters available minus one
to prevent the machine to be overchallenged.
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
res = rolling_window(as.vector(array(AirPassengers)), c(10,10,10), c(2,4),
                           Horizon = 2, Window = 3, Method = "r",
                           NumClusters = 1)
Error = res$Error
Forecast = res$Forecast
}
\keyword{Rolling forecasting origin}
