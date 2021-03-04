\name{neuralnet_prediction_scheme}
\alias{neuralnet_prediction_scheme}
\title{One Step Prediction Scheme for Neural Network}
\description{
This function computes the 1D array input for the prediction phase specific
for the Machine Learning one step forecast methods (Neural Network, SVM,
XGBoost).
}
\usage{
neuralnet_prediction_scheme(wmatrix, rhwtCoeff, ccps, agg_per_lvl)
}
\arguments{
\item{wmatrix}{Matrix with wavelet coefficients.}
\item{rhwtCoeff}{Matrix with smooth coefficients.}
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
\item{future_point}{Double carrying one step forecast.}
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
ccps = c(1,1,1)
scales = length(ccps) - 1
agg_per_lvl = c(2,4)
dec_res <- decomposition(as.vector(array(AirPassengers)), agg_per_lvl = agg_per_lvl)
trs_res <- training_scheme(dec_res$data,dec_res$wmatrix,dec_res$rhwtCoeff,
dec_res$scales, ccps,agg_per_lvl)
arr_future_points = trs_res$points_in_future
lm_matrix = trs_res$lsmatrix
num_feature = dim(lm_matrix)[2]
num_trainingpoint = dim(lm_matrix)[1]
forecast_scheme = neuralnet_prediction_scheme(dec_res$wmatrix,
dec_res$rhwtCoeff,ccps,agg_per_lvl)
}
\keyword{Wavelets}
