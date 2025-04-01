cv_glmnet_formula <- function(formula, data, family) {
  y <- data[[as.character(formula[[2]])]]
  formula[[2]] <- NULL
  X <- model.matrix(formula, data = data)
  
  out <- list(fit = cv.glmnet(X, y, family = family), 
              formula = formula)
  class(out) <- "cv_glmnet_formula"
  out
}

predict.cv_glmnet_formula <- function(object, newdata, ...) {
  X <- model.matrix(object$formula, data = newdata)
  predict(object$fit, newx = X, s = "lambda.1se", type = "response", ...) |> 
    as.vector()
}