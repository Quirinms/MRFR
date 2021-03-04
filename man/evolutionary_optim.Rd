\name{evolutionary_optim}
\alias{evolutionary_optim}
\title{Model Selection for Multiresolution Forecasts}
\description{
This function computes a model selection using a criterion (AIC, MRE).
The evaluation function (optimization function) is built with a rolling
forecasting origin (rolling_window function), which computes a h-step ahead
forecast (for h = 1, ..., horizon) for window_size many steps.
The input space is searched with an evolutionary optimization method.
The deployed forecast method can be an autoregression or a neural network
(multilayer perceptron with one hidden layer).
}
\usage{
evolutionary_optim(
  data,
  agg_per_lvl,
  data_name = "data",
  dir_name = "",
  horizon = 14,
  window_size = 365,
  method = "r",
  crit = "AIC",
  itermax = 1,
  lower_limit = 1,
  upper_limit = 2,
  write = FALSE,
  numClusters = 1)
}
\arguments{
\item{data}{Vector with n time series values.}
\item{agg_per_lvl}{Vector carrying numbers whose index is associated with the
wavelet level. The numbers indicate the number of time in points used for
aggregation from the original time series.}
\item{data_name}{String with name for output file.}
\item{dir_name}{String with directory for output file.}
\item{horizon}{Number indicating horizon for forecast from 1 to horizon.}
\item{window_size}{Number indicating how many points are used for cross
validation.}
\item{method}{String indicating which method to use (r = Autoregression, nn =
Neural Network).}
\item{crit}{String indicating which criterion to use:
(AIC = Akaikes Information Criterion)
(MRE = Mean Root Error)}
\item{itermax}{Number of iterations for evolutionary optimization method.}
\item{lower_limit}{Lower limit for coefficients selected for each level}
\item{upper_limit}{Higher limit for coefficients selected for each level}
\item{write}{Boolean indicating if the results are written to file (TRUE) or not
(FALSE)}
\item{numClusters}{Number of clusters used for parallel computing.}
}
\value{
None
}
\author{
Quirin Stier
}
\examples{
data(entsoe)
evolutionary_optim(data = entsoe$value, agg_per_lvl = c(2,4), horizon = 1,
window_size = 1, method = "r", crit = "AIC", itermax = 1, upper_limit = 1,
write = FALSE, numClusters = 1)
}
