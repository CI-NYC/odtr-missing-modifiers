library(mlr3superlearner)
library(mlr3extralearners)
library(lmtp)
library(dplyr)
library(kernelshap)
library(shapviz)
library(ggplot2)
library(future)
library(patchwork)
library(purrr)

source("R/load-data.R")
source("R/fit-tuned-xgboost.R")

# Load CTN data
data <- load_analysis_data()
vars <- data$vars
data <- data$data

# Load the model for the CATE based on doubly-robust unbiased transformation
model <- readRDS("data/application/drv/psuedo-regression-tuned-xgboost-cate-model.rds")

lookup_cate <-
  select(data, who, all_of(vars$V1)) |>
  select(-who) |> 
  unique()

lookup_cate$cate <- predict_tuned_xgboost(model, newdata = lookup_cate)

data <- left_join(data, lookup_cate)
shifted <- mutate(data, bupnx = as.numeric(cate < 0))

# Superlearner library
learners <- c("mean", "glm", "earth", "ranger", "bart", "cv_glmnet")

set.seed(634)

lmtp_no_homeless_odtr <- lmtp_tmle(
  data,
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shifted = shifted, 
  learners_trt = learners, 
  learners_outcome = learners, 
  folds = 20
)

saveRDS(lmtp_no_homeless_odtr, "data/application/drv/lmtp-no-homeless.rds")
