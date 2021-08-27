\name{mrf_multi_step_forecast}
\alias{mrf_multi_step_forecast}
\title{Multiresolution Forecast}
\description{
This function creates a multi step forecast for all horizons from 1 to steps
based on the manuscript [Stier/Thrun, 2021] which is currently under review.
The deployed forecast method can be an autoregression or a neural network
(multilayer perceptron with one hidden layer). Multi step forecasts are computed
recursively.
}
\usage{
mrf_multi_step_forecast(UnivariateData, Horizon, CoefficientCombination,
Aggregation, Method = "r")
}
\arguments{
\item{UnivariateData}{[1:n] Numerical vector with n values.}
\item{Horizon}{Number indicating horizon for forecast from 1 to horizon.}
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
List of
\item{multistep}{[1:Horizon] Numerical vector with forecast of horizon according
to its index.}
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
# One-step forecast (Multiresolution Forecast)
one_step = mrf_multi_step_forecast(UnivariateData = UnivariateData,
                                    Horizon = 2,
                                    CoefficientCombination = c(1,1,1),
                                    Aggregation = c(2,4),
                                    Method="r")
# Multi-step forecast (Multiresolution Forecast)
# Horizon = 2 => Forecast with Horizon 1 and 2 as vector
multi_step = mrf_multi_step_forecast(UnivariateData = UnivariateData,
                                      Horizon = 2,
                                      CoefficientCombination = c(1,1,1),
                                      Aggregation = c(2,4),
                                      Method="r")
}
\concept{Multi step forecasts}

