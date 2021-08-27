\name{mrf_model_selection}
\alias{mrf_model_selection}
\title{Model selection for Multiresolution Forecasts}
\description{
Evaluates the best coefficient combination for a given aggregation scheme based
on a rolling forecasting origin  based on the manuscript [Stier/Thrun, 2021]
which is currently under review.
}
\usage{
mrf_model_selection(UnivariateData, Aggregation, Horizon = 1, Window = 2,
Method = "r", crit = "AIC", itermax = 1, lower_limit = 1, upper_limit = 2,
NumClusters = 1)
}
\arguments{
\item{UnivariateData}{[1:n] Numerical vector with n values.}
\item{Aggregation}{[1:Scales] Numerical vector carrying numbers whose index is
associated with the wavelet level. The numbers indicate the number of time in
points used for aggregation from the original time series.}
\item{Horizon}{Number indicating horizon for forecast from 1 to horizon.}
\item{Window}{Number indicating how many points are used for cross validation.}
\item{Method}{String indicating which method to use. Available methods: 'r'  =
Autoregression. 'nn' = Neural Network.}
\item{crit}{String with criterion. Available criterions: "AIC" = Akaikes Info.
Crit. "MAE" = Mean Abs. Error. "MRE" = Mean Root Error. Default: crit = "AIC".}
\item{itermax}{Number of iterations used in the differential evolutionary
optimization algorithm. Default: itermax = 1.}
\item{lower_limit}{[1:Scales+1] Numeric vector: Lower limit for coefficients
selected for each level.}
\item{upper_limit}{[1:Scales+1] Numeric vector: Higher limit for coefficients
selected for each level.}
\item{NumClusters}{Number of clusters used for parallel computing. Default:
NumClusters = 1.}
}
\value{
\item{CoefficientCombination}{[1:Scales+1] Numerical vector with numbers which
are associated with wavelet levels. The last number is associated with the
smooth level. Each number determines the number of coefficient used per level.
The selection follows a specific scheme. Best combination of coefficients
found by the model selection algorithm.}
\item{Aggregation}{[1:Scales] Numerical vector carrying numbers whose index is
associated with the wavelet level. The numbers indicate the number of time in
points used for aggregation from the original time series. Best Aggregation
scheme found by the model selection algorithm.}
}
\details{
The evaluation function (optimization function) is built with a rolling
forecasting origin (rolling_window function), which computes a h-step ahead
forecast (for h = 1, ..., horizon) for 'Window' many steps.
The input space is searched with an evolutionary optimization method.
The search is restricted to one fixed aggregation scheme (parameter:
'Aggregation').
The deployed forecast method can be an autoregression or a neural network
(multilayer perceptron with one hidden layer).
}
\author{
Quirin Stier
}
\references{
[Stier/Thrun, 2021] Stier, Q., and Thrun, M. C.: Multiresolution Forecasting for
Industrial Applications, Technical Report, 2021.
}
\examples{
\donttest{
data(entsoe)
UnivariateData = entsoe$value
Aggregation = c(2,4)
res = mrf_model_selection(UnivariateData, Aggregation, Horizon = 1, Window = 2,
Method = "r", crit = "AIC", itermax = 1, lower_limit = 1, upper_limit = 2,
NumClusters = 1)
BestCoefficientCombination = res$CoefficientCombination
}
}
