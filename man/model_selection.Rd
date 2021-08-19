\name{model_selection}
\alias{model_selection}
\title{Model Selection for Multiresolution Forecasts}
\description{
This function computes a model selection using a criterion (AIC, MRE).
}
\usage{
model_selection(UnivariateData, Aggregation, Horizon = 14, Window = 365,
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
\item{crit}{String indicating which criterion to use. Available criterion: AIC =
Akaikes Information Criterion. MRE = Mean Root Error.}
\item{itermax}{Number of iterations for evolutionary optimization method.}
\item{lower_limit}{Lower limit for coefficients selected for each level}
\item{upper_limit}{Higher limit for coefficients selected for each level}
\item{NumClusters}{Number of clusters used for parallel computing.}
}
\value{
\item{Error}{[1:Window, 1:Horizon] Numerical Matrix with 'Window' many rows
entries indicating one time point with 'Horizon' many forecast errors.}
\item{Best}{[1:Scales+1] Numerical vector with integers associated with the
best found number of coefficients per wavelet scale (1:Scales) and number of
coefficients for the smooth approximation level in the last entry.}
}
\details{
The evaluation function (optimization function) is built with a rolling
forecasting origin (rolling_window function), which computes a h-step ahead
forecast (for h = 1, ..., horizon) for window_size many steps.
The input space is searched with an evolutionary optimization method.
The deployed forecast method can be an autoregression or a neural network
(multilayer perceptron with one hidden layer).
}
\author{
Quirin Stier
}
\references{
Hyndman, R. and Athanasopoulos, G. Forecasting: principles and practice.
OTexts, 3 edition. 2018.
}
\examples{
\donttest{
data(entsoe)
model_selection(UnivariateData = entsoe$value, Aggregation = c(2,4), Horizon = 1,
Window = 1, Method = "r", crit = "AIC", itermax = 1, upper_limit = 1,
NumClusters = 1)
}
}
