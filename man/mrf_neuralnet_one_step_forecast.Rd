\name{mrf_neuralnet_one_step_forecast}
\alias{mrf_neuralnet_one_step_forecast}
\title{One Step Forecast with Neural Network}
\description{
This function creates a one step forecast using a multi layer perceptron
with one hidden Layer. The number of input is the sum of all coefficients
chosen with the parameter CoefficientCombination. The CoefficientCombination
parameter controls the number of coefficients chosen for each wavelet and smooth
part level individually.
}
\usage{
mrf_neuralnet_one_step_forecast(UnivariateData, CoefficientCombination,
Aggregation)
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
}
\value{
\item{forecast}{Numerical value with one step forecast}
}
\references{
Aussem, A., Campbell, J., and Murtagh, F. Waveletbased Feature Extraction
and Decomposition Strategies for Financial Forecasting. International Journal of
Computational Intelligence in Finance, 6,5-12, 1998.

Renaud, O., Starck, J.-L., and Murtagh, F. Prediction based on a Multiscale De-
composition. International Journal of Wavelets, Multiresolution and Information
Processing, 1(2):217-232. doi:10.1142/S0219691303000153, 2003.

Murtagh, F., Starck, J.-L., and Renaud, O. On Neuro-Wavelet Modeling. Decision
Support Systems, 37(4):475-484. doi:10.1016/S0167-9236(03)00092-7, 2004.

Renaud, O., Starck, J.-L., and Murtagh, F. Wavelet-based combined Signal Filter-
ing and Prediction. IEEE Transactions on Systems, Man, and Cybernetics, Part
B (Cybernetics), 35(6):1241-1251. doi:10.1109/TSMCB.2005.850182, 2005.
}
\author{
Quirin Stier
}
\examples{
data(AirPassengers)
len_data = length(as.vector(array(AirPassengers)))
UnivariateData = as.vector(AirPassengers)[1:(len_data-1)]
CoefficientCombination = c(1,1,1)
Aggregation = c(2,4)
forecast = mrf_neuralnet_one_step_forecast(UnivariateData,
                                           CoefficientCombination,
                                           Aggregation)
true_value = array(AirPassengers)[len_data]
error = true_value - forecast
}
\keyword{One Step Forecasts, Neural Networks, Multilayer Perceptron}
