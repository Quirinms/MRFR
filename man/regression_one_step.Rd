\name{regression_one_step}
\alias{regression_one_step}
\title{One Step Forecast with Regression}
\description{
This function creates a one step forecast using an autoregression method.
The ccps parameter controls the number of coefficients chosen for each wavelet
and smooth part level individually.
}
\usage{
regression_one_step(data, ccps, agg_per_lvl)
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
}
\value{
List of
\item{forecast}{Double one step forecast}
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
len_data = length(as.vector(array(AirPassengers)))
result = regression_one_step(array(AirPassengers)[1:(len_data-1)], c(1,1,1), c(2,4))
forecast = result$forecast
}
\keyword{Autoregression, One step forecast}
