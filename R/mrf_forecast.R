mrf_forecast <- function(Model, Horizon=1){
  # DESCRIPTION
  # Computes multi-step forecasts with a given multiresolution model.
  #
  # INPUT
  # Model        List containing model specifications from mrf_train().
  #
  # OPTIONAL
  # Horizon      Number indicating forecast horizon. Horizon = 1 means
  #              one-step forecast and Horizon > 1 means a one-step forecast
  #              and all multi-step forecasts from horizon 2 to 'Horizon'.
  #
  # OUTPUT
  # List of 2 elements:
  # Forecast[1:Horizon]    Numerical vector with forecast with 'Horizon' many
  #                        entries.
  # Model                  List containing model specifications from mrf_train().
  #
  # Author: QS, 08/2021
  Data         = Model$Data
  Coefficients = Model$CoefficientCombination
  Aggregation  = Model$Aggregation
  Method       = Model$Method
  if(!is.vector(Data)){
    message("Data must be of type vector")
    return()
  }
  if(!is.double(Horizon)){
    message("Horizon must be of type double")
    return()
  }
  if(!is.vector(Coefficients)){
    message("Coefficients must be of type vector")
    return()
  }
  if(!is.vector(Aggregation)){
    message("Aggregation must be of type vector")
    return()
  }
  if(!is.character(Method)){
    message("Method must be of type character")
    return()
  }
  Forecast = mrf_multi_step_forecast(Data, Horizon, Coefficients,
                                     Aggregation, Method)
  return(list("Forecast"=Forecast,
              "Model"=Model))
}


#
#
#
#
#
