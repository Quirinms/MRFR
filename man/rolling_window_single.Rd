\name{rolling_window_single}
\alias{rolling_window_single}
\title{Rolling Window for Multiresolution Forecasts}
\description{
This function creates a single step for a rolling forecasting origin for a
specifc one step forecast method. Multi step forecasts are computed recursively
with the one step forecast method.
Thus, h-step forecast for h = 1,..., horizon for window_size many steps can be
computed.
The forecasting method can be an autoregression or a neural network (multilayer
perceptron).
The ccps parameter controls the number of coefficients chosen for each wavelet
and smooth part level individually.
}
\usage{
rolling_window_single(
  i,
  data,
  ccps,
  agg_per_lvl,
  horizon = 14,
  window_size = 365,
  method = "r"
)
}
\arguments{
\item{i}{number indicating index for parallel computation.}
\item{data}{Time series with n values.}
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
\item{horizon}{Number indicating horizon for forecast from 1 to horizon.}
\item{window_size}{Number indicating how many points are used to create cross
validation.}
\item{method}{String indicating which method to use (r = Autoregression, nn =
Neural Network).}
}
\value{
List of parameter with a 2D matrix of the forecast error.
}
\author{
Quirin Stier
}

