# borrowed from https://github.com/Larsvanderlaan/CDML/blob/main/R/isotonic_regression.R
isoreg_with_xgboost <- function(x, y, max_depth = 15, min_child_weight = 20, weights = NULL) {
  if(is.null(weights)) {
    weights <- rep(1, length(y))
  }
  # Create an XGBoost DMatrix object from the data, including weights
  data <- xgboost::xgb.DMatrix(data = as.matrix(x), label = as.vector(y), weight = weights)
  
  # Set parameters for the monotonic XGBoost model
  params <- list(
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    monotone_constraints = 1,  # Enforce monotonic increase
    eta = 1,
    gamma = 0,
    lambda = 0
  )
  
  # Train the model with one boosting round
  iso_fit <- xgboost::xgb.train(params = params, data = data, nrounds = 1)
  
  # Prediction function for new data
  fun <- function(x) {
    data_pred <- xgboost::xgb.DMatrix(data = as.matrix(x))
    pred <- predict(iso_fit, data_pred)
    return(pred)
  }
  
  return(fun)
}