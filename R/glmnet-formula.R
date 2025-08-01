cv_glmnet_formula <- function(formula, data, family, subset, ...) {
  y <- data[[as.character(formula[[2]])]]
  
  if (is.null(y)) {
    y <- get(as.character(formula[[2]]))
  }
  
  formula[[2]] <- NULL
  X <- model.matrix(formula, data = data)
  
  if (missing(family)) family <- "gaussian"
  
  if (!missing(subset)) {
    X <- X[subset, ]
    y <- y[subset]
  }
  
  out <- list(fit = cv.glmnet(X, y, family = family), 
              formula = formula)
  class(out) <- "cv_glmnet_formula"
  out
}

predict.cv_glmnet_formula <- function(object, newdata, s = c("lambda.1se", "lambda.min"), ...) {
  X <- model.matrix(object$formula, data = newdata)
  predict(object$fit, newx = X, s = match.arg(s), type = "response", ...) |> 
    as.vector()
}