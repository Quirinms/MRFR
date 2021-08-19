nested_cross_validation <- function(UnivariateData, Horizon = 14,
                                    EvaluationLength=2, TestLength=2,
                                    Method="r", MultivariateData=NULL, NumMV=1,
                                    NumClusters=1){
  # INPUT
  # UnivariateData[1:n]      Numerical vector with n values
  #
  # OPTIONAL
  # Horizon             Number indicating horizon for forecast from 1 to horizon.
  # EvaluationLength    Number indicating how many points are used for cross
  #                     validation for the evaluation dataset.
  # TestLength          Number indicating how many points are used for cross
  #                     validation for the test dataset.
  # Method              String indicating which method to use
  #                     Available methods: 'r'  = Autoregression
  #                                        'nn' = Neural Network
  # MultivariateData    Not implemented yet.
  # NumMV               Not implemented yet.
  # NumClusters         Number of clusters used for parallel computing.
  #
  #
  # OUTPUT
  # Best[1:Scales+1]                Numerical vector with integers associated
  #                                 with the best found number of coefficients
  #                                 per wavelet scale (1:Scales) and number of
  #                                 coefficients for the smooth approximation
  #                                 level in the last entry.
  # Error[1:Window,1:Horizon]       Numerical Matrix with 'Window' many rows
  #                                 entries indicating one time point with
  #                                 'Horizon' many forecast errors.
  # Forecast[1:Window,1:Horizon]    Numerical Matrix with 'Window' many rows
  #                                 entries indicating one time point with
  #                                 'Horizon' many forecasts.
  #
  # Author: QS, 02/2021
  DataLength = length(UnivariateData)
  # Create Test data
  # => There is no need for explicitly defining an evaluation dataset in case of a split in test and evaluation data.
  TestUnivariateData = UnivariateData[0:(DataLength - EvaluationLength)]


  if (length(TestUnivariateData) != (DataLength - EvaluationLength)){
    print("Something went wrong when splitting Data in Test and Evaluation part in modelSelection.py!")
  }
  if(length(TestUnivariateData) == 0){
    print("modelSelection.py: There is no Testdata for computing a model!")
  }
  #if(typeof(MultivariateData) == np.ndarray){
  #  TestMultivariateData = MultivariateData[0:(DataLength - EvaluationLength)]
  #}
  #if(TestMultivariateData.shape[0] != (DataLength - EvaluationLength)){
  #  print("Something went wrong when splitting MultivariateData in Test and Evaluation part in nested_cross_validation.R!")
  #}

  #if(length(TestUnivariateData) == 0){
  #  print("nested_cross_validation.R: There is no Testdata (Multivariate data part) for computing a model!")
  #}else{
  #  TestMultivariateData = NULL
  #}

  AllAggregations = list(c(2,4), c(2,4,8), c(2,4,8,16), c(2,4,8,16,32))

  lst_Results = list()
  selectMAE = c()
  counter = 1
  for(i in 1:length(AllAggregations)){
    #NumCoeffs = length(Aggregation) + 1
    lower_limit = 1
    if(Method=="r"){
      upper_limit = 15
    }else{
      upper_limit = 8
    }
    Aggregation = AllAggregations[[i]]
    res = model_selection(UnivariateData = TestUnivariateData,
                          Aggregation=Aggregation,
                          Horizon = Horizon,
                          Window = TestLength,
                          Method = Method,
                          crit = "AIC",
                          itermax = 1, lower_limit = lower_limit,
                          upper_limit = upper_limit,
                          NumClusters = NumClusters)
    MAE = sum(abs(res$Error))/(TestLength*Horizon)
    lst_Results[[i]] = res$Best
    selectMAE = c(selectMAE, MAE)
  }
  IdxMin = which(selectMAE==min(selectMAE))
  Best = lst_Results[[IdxMin]]

  res = rolling_window(UnivariateData = UnivariateData,
                       CoefficientCombination = Best,
                       Aggregation = Aggregation,
                       Horizon = Horizon,
                       Window = EvaluationLength,
                       Method = Method,
                       NumClusters = NumClusters)
  return(list("Best"=Best, "Error"=res$Error, "Forecast"=res$Forecast))
}



#
#
#
#
#
