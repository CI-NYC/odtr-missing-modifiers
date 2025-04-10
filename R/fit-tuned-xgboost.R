library(mlr3verse)

fit_xgboost_with_tuning <- function(data, target, n_evals = 100) {
  task <- as_task_regr(data, target)
  search_space <- lts("regr.xgboost.default")
  tuner_grid_search = tnr("grid_search", resolution = 5, batch_size = 5)
  xgboost_auto_tuner <- auto_tuner(
    tuner = tuner_grid_search,
    learner = lrn("regr.xgboost"), 
    resampling = rsmp("cv", folds = 10), 
    measure = msr("regr.mse"),
    terminator = trm("evals", n_evals = n_evals), 
    search_space = search_space
  )
  xgboost_auto_tuner$train(task)
  return(xgboost_auto_tuner)
}

predict_tuned_xgboost <- function(object, newdata) {
  object$predict_newdata(newdata)$response
}
