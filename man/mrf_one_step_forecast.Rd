\name{mrf_one_step_forecast}
\alias{mrf_one_step_forecast}
\title{mrf_one_step_forecast Step Forecast}
\description{
This function creates a one step forecast using the multiresolution forecasting
framework based on the manuscript [Stier/Thrun, 2021] which is currently under
review.
}
\usage{
mrf_one_step_forecast(UnivariateData, CoefficientCombination, Aggregation,
Method="r")
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
\item{Method}{String indicating which method to use. Available methods: 'r'  =
Autoregression. 'nn' = Neural Network.}
}
\value{
\item{forecast}{Numerical value with one step forecast}
}
\references{
[Stier/Thrun, 2021] Stier, Q., and Thrun, M. C.: Multiresolution Forecasting for
Industrial Applications, Technical Report, 2021.

}
\author{
Quirin Stier
}
\examples{
data(AirPassengers)
len_data = length(array(AirPassengers))
UnivariateData = as.vector(AirPassengers)[1:(len_data-1)]
forecast = mrf_one_step_forecast(UnivariateData, c(1,1,1), c(2,4))
true_value = array(AirPassengers)[len_data]
error = true_value - forecast
}
\concept{One Step Forecasts}
\concept{Neural Networks}
\concept{Multilayer Perceptron}
\keyword{Regression}
