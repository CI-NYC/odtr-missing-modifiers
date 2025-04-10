library(mlr3superlearner)
library(mlr3extralearners)
library(lmtp)
library(dplyr)
library(kernelshap)
library(shapviz)
library(ggplot2)
library(future)

source("R/load-data.R")
source("R/fit-tuned-xgboost.R")

# load CTN data
data <- load_analysis_data()
vars <- data$vars
data <- data$data

# load pseudo regression models
models <- readRDS("data/application/drv/psuedo-regression-tuned-xgboost-models.rds")

cate <- function(object, newdata) {
  hbar_1_homeless <- predict_tuned_xgboost(object$kappa_model_1_1, newdata = newdata)
  hbar_0_homeless <- predict_tuned_xgboost(object$kappa_model_1_0, newdata = newdata)
  qbar_homeless <- predict_tuned_xgboost(object$lambda_model_1, newdata = newdata)
  
  hbar_1_not_homeless <- predict_tuned_xgboost(object$kappa_model_0_1, newdata = newdata)
  hbar_0_not_homeless <- predict_tuned_xgboost(object$kappa_model_0_0, newdata = newdata)
  qbar_not_homeless <- predict_tuned_xgboost(object$lambda_model_0, newdata = newdata)
  
  psi_homeless <- (hbar_1_homeless - hbar_0_homeless) / qbar_homeless
  psi_not_homeless <- (hbar_1_not_homeless - hbar_0_not_homeless) / qbar_not_homeless
  
  newdata$homeless*psi_homeless + 
    (1 - newdata$homeless)*psi_not_homeless
}

ctn30_homeless <- 
  ctn30_no_homeless <-
  filter(data, project == 0) |>
  select(who, all_of(c(vars$V1, vars$V2)))

ctn30_homeless$homeless <- 1
ctn30_no_homeless$homeless <- 0

lookup_cate <-
  filter(data, project == 1) |>
  select(who, all_of(c(vars$V1, vars$V2))) |>
  rbind(ctn30_homeless, ctn30_no_homeless) |>
  select(-who) |> 
  unique()

# shap values -------------------------------------------------------------

off_row <- which(rowSums(lookup_cate == 0) == ncol(lookup_cate))[1]
shap <- kernelshap(object = models, pred_fun = cate, 
                   X = lookup_cate[-off_row, ], 
                   bg_X = lookup_cate[off_row, ])

shapviz(shap) |> 
  sv_importance(kind = "bee") + 
  scale_y_discrete(labels = c(
    "homeless" = "Homeless",
    "ivdrug" = "IV drug use", 
    "hwithdraw_4" = "Severe withdrawal", 
    "hwithdraw_3" = "Medium withdrawal",
    "hwithdraw_2" = "Mild withdrawal",
    "bbenzo30_base" = "Benzo. use", 
    "cocdisorder" = "Cocaine use disorder", 
    "xrace_3" = "Hispanic", 
    "hasEpilepsy" = "Epilepsy", 
    "xrace_4" = "Other race", 
    "xrace_2" = "Non-Hispanic, Black"
  )) + 
  theme_bw()

shapviz(shap) |> 
  sv_importance() + 
  scale_y_discrete(labels = c(
    "homeless" = "Homeless",
    "ivdrug" = "IV drug use", 
    "hwithdraw_4" = "Severe withdrawal", 
    "hwithdraw_3" = "Medium withdrawal",
    "hwithdraw_2" = "Mild withdrawal",
    "bbenzo30_base" = "Benzo. use", 
    "cocdisorder" = "Cocaine use disorder", 
    "xrace_3" = "Hispanic", 
    "hasEpilepsy" = "Epilepsy", 
    "xrace_4" = "Other race", 
    "xrace_2" = "Non-Hispanic, Black"
  )) + 
  theme_bw()

negative_homeless_shap_rows <- which(shap$S[, "homeless"] < 0)
positive_homeless_shap_rows <- which(shap$S[, "homeless"] > 0)
summary(lookup_cate[-off_row, ][negative_homeless_shap_rows, ])
summary(lookup_cate[-off_row, ][positive_homeless_shap_rows, ])

# relapse under the rule --------------------------------------------------

lookup_cate$cate <- cate(models, lookup_cate)

ctn30 <- left_join(ctn30_homeless, lookup_cate) |>
  rename(cate_homeless = cate) |>
  select(who, cate_homeless) |>
  cbind({
    left_join(ctn30_no_homeless, lookup_cate) |>
      rename(cate_no_homeless = cate) |>
      select(cate_no_homeless)
  })

ctn30$cate_min <- pmin(ctn30$cate_no_homeless, ctn30$cate_homeless)
ctn30$cate_max <- pmax(ctn30$cate_no_homeless, ctn30$cate_homeless)

ctn51 <- filter(data, project == 1) |> 
  left_join(lookup_cate) |> 
  mutate(bupnx = as.numeric(cate < 0))

ctn30_decisive <- 
  filter(ctn30, sign(cate_min) == sign(cate_max)) |> 
  select(who, cate = cate_max) |> 
  inner_join(data, by = "who") |> 
  mutate(bupnx = as.numeric(cate < 0))

ctn30_ambiguous_bupnx <- 
  filter(ctn30, sign(cate_min) != sign(cate_max)) |> 
  select(who) |> 
  inner_join(data, by = "who") |> 
  mutate(bupnx = 1)

ctn30_ambiguous_no_bupnx <- 
  filter(ctn30, sign(cate_min) != sign(cate_max)) |> 
  select(who) |> 
  inner_join(data, by = "who") |> 
  mutate(bupnx = 0)

shifted_ambiguous_bupnx <- bind_rows(ctn51, ctn30_decisive, ctn30_ambiguous_bupnx)
shifted_ambiguous_no_bupnx <- bind_rows(ctn51, ctn30_decisive, ctn30_ambiguous_no_bupnx)

# superlearner library
learners <- c("mean", "glm", "earth", "ranger", "bart", "cv_glmnet")

lmtp_all_bupnx <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shift = static_binary_on, 
  learners_trt = learners, 
  learners_outcome = learners
)

lmtp_all_naltrexone <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shift = static_binary_off, 
  learners_trt = learners, 
  learners_outcome = learners
)

lmtp_ambiguous_bupnx <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shifted = arrange(shifted_ambiguous_bupnx, who), 
  learners_trt = learners, 
  learners_outcome = learners
)

lmtp_ambiguous_no_bupnx <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shifted = arrange(shifted_ambiguous_no_bupnx, who), 
  learners_trt = learners, 
  learners_outcome = learners
)
