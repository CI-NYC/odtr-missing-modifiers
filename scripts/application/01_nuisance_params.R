library(origami)
library(mlr3superlearner)
library(mlr3extralearners)
library(dplyr)

source("R/isotonic-calibration.R")
source("R/fit-tuned-xgboost.R")
source("R/load-data.R")

# load CTN data
data <- load_analysis_data()
vars <- data$vars
data <- data$data

# number of crossfit folds
folds <- 20

# superlearner library
learners <- c("mean", "glm", "earth", "ranger", "bart", "cv_glmnet")

set.seed(4534)

crossfit_folds <- make_folds(data, V = folds)
n <- nrow(data)

b <- m <- matrix(nrow = n, ncol = 2)
colnames(b) <- colnames(m) <- c("A=0", "A=1")

gs <- g <- q <- th <- matrix(nrow = n, ncol = 1)

# only used for rapid prototyping
if (folds == 1) {
  crossfit_folds[[1]]$training_set <- crossfit_folds[[1]]$validation_set
}

set.seed(431)

for (v in 1:folds) {
  train <- training(data, crossfit_folds[[v]])
  valid <- validation(data, crossfit_folds[[v]])
  valid_1 <- mutate(valid, bupnx = 1)
  valid_0 <- mutate(valid, bupnx = 0)
  
  valid_indexes <- crossfit_folds[[v]]$validation_set
  
  m_model <- mlr3superlearner(
    data = train[, c(vars$Y, vars$A, vars$V1, vars$W)], 
    target = vars$Y, 
    outcome_type = "binomial",
    library = learners, 
    discrete = FALSE
  )
  
  m[valid_indexes, "A=0"] <- predict(m_model, newdata = valid_0)
  m[valid_indexes, "A=1"] <- predict(m_model, newdata = valid_1)
  
  b_model <- mlr3superlearner(
    data = filter(train, week_12_relapse == 1 & project == 1) |> 
      select(all_of(c(vars$V2, vars$A, vars$V1, vars$W))),
    target = vars$V2, 
    outcome_type = "binomial", 
    library = learners, 
    discrete = FALSE
  )
  
  b[valid_indexes, "A=0"] <- predict(b_model, newdata = valid_0)
  b[valid_indexes, "A=1"] <- predict(b_model, newdata = valid_1)
  
  q_model <- mlr3superlearner(
    data = filter(train, project == 1) |> 
      select(all_of(c(vars$V2, vars$V1, vars$W))), 
    target = vars$V2,
    outcome_type = "binomial",
    library = learners,
    discrete = FALSE
  )
  
  q[valid_indexes, 1] <- predict(q_model, newdata = valid)
  
  g_model <- mlr3superlearner(
    data = select(train, all_of(c(vars$A, vars$W, vars$V1))), 
    target = vars$A,
    outcome_type = "binomial",
    library = learners,
    discrete = FALSE
  )

  g[valid_indexes, 1] <- predict(g_model, newdata = valid)
  
  gs_model <- mlr3superlearner(
    data = filter(train, project == 1) |> 
      select(all_of(c(vars$A, vars$W, vars$V1))), 
    target = vars$A,
    outcome_type = "binomial",
    library = learners,
    discrete = FALSE
  )
  
  gs[valid_indexes, 1] <- predict(gs_model, newdata = valid)
  
  t_model <- mlr3superlearner(
    data = select(train, all_of(c(vars$S, vars$V1, vars$W))), 
    target = vars$S, 
    outcome_type = "binomial",
    library = learners,
    discrete = FALSE
  )
  
  th[valid_indexes, 1] <- predict(t_model, newdata = valid)
}

Ia1 <- as.numeric(data$bupnx == 1)
Ia0 <- as.numeric(data$bupnx == 0)
Is1 <- as.numeric(data$project == 1)
Iv2_1 <- as.numeric(data$homeless == 1)
Iv2_1[is.na(Iv2_1)] <- -999
Iv2_0 <- as.numeric(data$homeless == 0)
Iv2_0[is.na(Iv2_0)] <- -999
Iy1 <- as.numeric(data$week_12_relapse == 1)

# isotonic calibrations
m[, "A=0"] <- isoreg_with_xgboost(m[Ia1 == 0, "A=0"], data$week_12_relapse[Ia1 == 0])(m[, "A=0"])
m[, "A=1"] <- isoreg_with_xgboost(m[Ia1 == 1, "A=1"], data$week_12_relapse[Ia1 == 1])(m[, "A=1"])
g[, 1] <- isoreg_with_xgboost(g, data$bupnx)(g)
gs[, 1] <- isoreg_with_xgboost(gs[Is1 == 1], data$bupnx[Is1 == 1])(gs)
q[, 1] <- isoreg_with_xgboost(q[Is1 == 1], Iv2_1[Is1 == 1])(q)
th[, 1] <- isoreg_with_xgboost(th, Is1)(th)
b[, "A=0"] <- isoreg_with_xgboost(b[Ia0*Iy1*Is1 == 1, "A=0"], Iv2_1[Ia0*Iy1*Is1 == 1])(b[, "A=0"])
b[, "A=1"] <- isoreg_with_xgboost(b[Ia1*Iy1*Is1 == 1, "A=1"], Iv2_1[Ia1*Iy1*Is1 == 1])(b[, "A=1"])

# abusing the fact that the DAG implies s-admissibility
ms <- m

# construct r
r1 <- m[, "A=1"]*gs*th 
r0 <- m[, "A=0"]*(1 - gs)*th 

# truncating weights
g <- pmin(pmax(g, 0.025), 0.975)
th <- pmax(th, 0.025)
r1 <- pmax(r1, 0.025)
r0 <- pmax(r0, 0.025)

# save nuisance parameters
list(
  m = m,
  b = b,
  g = g,
  gs = gs,
  q = q,
  t = th,
  r1 = r1,
  r0 = r0
) |> 
  saveRDS("data/application/drv/nuisance-parameters-calibrated-updated.rds")

# riesz representers
alpha_g1 <- Ia1 / g
alpha_g0 <- Ia0 / (1 - g)
alpha_t <- Is1 / th
alpha_r1 <- (Is1*Iy1*Ia1) / r1
alpha_r0 <- (Is1*Iy1*Ia0) / r0

# uncentered eifs
D_k1_v21_uc <- alpha_r1*(m[, "A=1"]*Iv2_1 - m[, "A=1"]*b[, "A=1"]) + 
  alpha_g1*(b[, "A=1"]*data$week_12_relapse - b[, "A=1"]*m[, "A=1"]) +
  b[, "A=1"]*m[, "A=1"]

D_k0_v21_uc <- alpha_r0*(m[, "A=0"]*Iv2_1 - m[, "A=0"]*b[, "A=0"]) + 
  alpha_g0*(b[, "A=0"]*data$week_12_relapse - b[, "A=0"]*m[, "A=0"]) +
  b[, "A=0"]*m[, "A=0"]

D_k1_v20_uc <- alpha_r1*(m[, "A=1"]*Iv2_0 - m[, "A=1"]*(1 - b[, "A=1"])) + 
  alpha_g1*((1 - b[, "A=1"])*data$week_12_relapse - (1 - b[, "A=1"])*m[, "A=1"]) +
  (1 - b[, "A=1"])*m[, "A=1"]

D_k0_v20_uc <- alpha_r0*(m[, "A=0"]*Iv2_0 - m[, "A=0"]*(1 - b[, "A=0"])) + 
  alpha_g0*((1 - b[, "A=0"])*data$week_12_relapse - (1 - b[, "A=0"])*m[, "A=0"]) +
  (1 - b[, "A=0"])*m[, "A=0"]

D_q1_uc <- alpha_t*(Iv2_1 - q) + q
D_q0_uc <- alpha_t*(Iv2_0 - (1 - q)) + (1 - q)

set.seed(7563)

# fit pseudo regressions
kappa_model_1_1 <- fit_xgboost_with_tuning(
  mutate(data, D_k1_v21_uc = D_k1_v21_uc) |> 
    select(D_k1_v21_uc, all_of(vars$V1)), 
  "D_k1_v21_uc", 100
)

kappa_model_1_0 <- fit_xgboost_with_tuning(
  mutate(data, D_k0_v21_uc = D_k0_v21_uc) |> 
    select(D_k0_v21_uc, all_of(vars$V1)), 
  "D_k0_v21_uc", 100
)

kappa_model_0_1 <- fit_xgboost_with_tuning(
  mutate(data, D_k1_v20_uc = D_k1_v20_uc) |> 
    select(D_k1_v20_uc, all_of(vars$V1)), 
  "D_k1_v20_uc", 100
)

kappa_model_0_0 <- fit_xgboost_with_tuning(
  mutate(data, D_k0_v20_uc = D_k0_v20_uc) |> 
    select(D_k0_v20_uc, all_of(vars$V1)), 
  "D_k0_v20_uc", 100
)

lambda_model_1 <- fit_xgboost_with_tuning(
  mutate(data, D_q1_uc = D_q1_uc) |> 
    select(D_q1_uc, all_of(vars$V1)), 
  "D_q1_uc", 100
)

lambda_model_0 <- fit_xgboost_with_tuning(
  mutate(data, D_q0_uc = D_q0_uc) |> 
    select(D_q0_uc, all_of(vars$V1)), 
  "D_q0_uc", 100
)

# save pseudo regression models
list(
  kappa_model_1_1 = kappa_model_1_1,
  kappa_model_1_0 = kappa_model_1_0,
  kappa_model_0_1 = kappa_model_0_1,
  kappa_model_0_0 = kappa_model_0_0,
  lambda_model_1 = lambda_model_1,
  lambda_model_0 = lambda_model_0
) |> 
  saveRDS("data/application/drv/psuedo-regression-tuned-xgboost-models-updated.rds")

# ignoring homeless status ------------------------------------------------

alpha <- (2*data$bupnx - 1) / (data$bupnx*g + (1 - data$bupnx)*(1 - g))
D_ate <- alpha * (data$week_12_relapse - (data$bupnx*m[, "A=1"] + (1 - data$bupnx)*m[, "A=0"])) + 
  m[, "A=1"] - m[, "A=0"]

set.seed(53443)
model_cate <- fit_xgboost_with_tuning(
  mutate(data, D_ate = D_ate) |> 
    select(D_ate, all_of(vars$V1)), 
  "D_ate", 100
)

saveRDS(model_cate, "data/application/drv/psuedo-regression-tuned-xgboost-cate-model.rds")
