\name{multi_step}
\alias{multi_step}
\title{Multi Step Forecast}
\description{
This function creates a multi step forecast for all horizons from 1 to steps.
The deployed forecast method can be an autoregression or a neural network
(multilayer perceptron with one hidden layer). Multi step forecasts are computed
recursively.
}
\usage{
multi_step(data, steps, ccps, agg_per_lvl, method = "r")
}
\arguments{
\item{data}{Vector with n time series values.}
\item{steps}{Number indicating horizon for forecast from 1 to horizon.}
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
\item{method}{String indicating which method to use (r = Autoregression, nn =
Neural Network).}
}
\value{
List of
\item{forecast}{Vector carrying the forecasts, where the index of the
entry is associated with the horizon of the forecast.}
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
result = multi_step(array(AirPassengers)[1:(len_data-1)], 2, c(1,1,1), c(2,4), method="r")
forecast = result$forecast
}
\keyword{Multi step forecasts}

