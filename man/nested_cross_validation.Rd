\name{nested_cross_validation}
\alias{nested_cross_validation}
\title{Nested cross validation for Multiresolution Forecasts}
\description{
This function computes a nested cross validation
(with the rolling forecasting origin).
The data is split into 3 datasets: training, test and evaluation dataset.
The best model is selected on the test and its performance is measured on the
evaluation dataset.
}
\usage{
nested_cross_validation(UnivariateData, Horizon=14, EvaluationLength=2, TestLength=2,
Method = "r", MultivariateData=NULL, NumMV=1, NumClusters = 1)
}
\arguments{
\item{UnivariateData}{[1:n] Numerical vector with n values.}
\item{Horizon}{Number indicating horizon for forecast from 1 to horizon.}
\item{EvaluationLength}{Number indicating how many points are used for cross
validation for the evaluation dataset.}
\item{TestLength}{Number indicating how many points are used for cross
validation for the test dataset.}
\item{Method}{String indicating which method to use. Available methods: 'r'  =
Autoregression. 'nn' = Neural Network.}
\item{MultivariateData}{Not implemented yet.}
\item{NumMV}{Not implemented yet.}
\item{NumClusters}{Number of clusters used for parallel computing.}
}
\value{
\item{Best}{[1:Scales+1] Numerical vector with integers associated with the
best found number of coefficients per wavelet scale (1:Scales) and number of
coefficients for the smooth approximation level in the last entry.}
\item{Error}{[1:Window, 1:Horizon] Numerical Matrix with 'Window' many rows
entries indicating one time point with 'Horizon' many forecast errors.}
\item{Forecast}{[1:Window, 1:Horizon]  Numerical Matrix with 'Window' many rows
entries indicating one time point with 'Horizon' many forecasts.}
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
res = nested_cross_validation(entsoe$value, Horizon = 2, EvaluationLength=2,
TestLength=2, Method="r", MultivariateData=NULL, NumMV=1, NumClusters=1)
BestCoefficientCombination = res$Best
Error = res$Error
Forecast = res$Forecast
}
}
