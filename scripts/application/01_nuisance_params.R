library(origami)
library(mlr3superlearner)
library(mlr3extralearners)
library(dplyr)
library(fastDummies)

source("R/isotonic-calibration.R")

# load CTN data
data <- readRDS("data/application/drv/ctn_data.rds")

# TODO - figure out what to do with these variables
data <- filter(data, !is.na(hasBipolar), !is.na(hasSchiz), !is.na(hwithdraw))

# load yaml files of variables to use
vars <- yaml::read_yaml("scripts/application/vars.yaml")

# dummy coding categorical variables
V1 <- dummy_cols(select(data, all_of(vars$V1)), 
                 remove_first_dummy = TRUE, 
                 remove_selected_columns = TRUE)

data <- bind_cols(select(data, -all_of(vars$V1)), V1)

vars$V1 <- names(V1)

# number of crossfit folds
folds <- 1

# superlearner library
learners <- c("mean", "glm", "earth", "ranger", "bart", "cv_glmnet")

crossfit_folds <- make_folds(data, V = folds)
n <- nrow(data)

m <- matrix(nrow = n, ncol = 2)
colnames(m) <- c("A=0", "A=1")
mc <- m

g <- q <- e <- b <- th <- matrix(nrow = n, ncol = 1)
gc <- qc <- ec <- bc <- thc <- g

crossfit_folds[[1]]$training_set <- crossfit_folds[[1]]$validation_set

for (v in 1:folds) {
  train <- training(data, crossfit_folds[[v]])
  valid <- validation(data, crossfit_folds[[v]])
  valid_1 <- mutate(valid, bupnx = 1)
  valid_0 <- mutate(valid, bupnx = 0)
  
  valid_indexes <- crossfit_folds[[v]]$validation_set
  
  m_model <- mlr3superlearner(
    data = train[, c(vars$Y, vars$A, vars$W)], 
    target = vars$Y, 
    outcome_type = "binomial",
    library = learners, 
    discrete = FALSE
  )
  
  m[valid_indexes, "A=0"] <- predict(m_model, newdata = valid_0)
  m[valid_indexes, "A=1"] <- predict(m_model, newdata = valid_1)
  
  b_model <- mlr3superlearner(
    data = filter(train, week_12_relapse == 1 & project == 1) |> 
      select(all_of(c(vars$V2, vars$A, vars$W))),
    target = vars$V2, 
    outcome_type = "binomial", 
    library = learners, 
    discrete = FALSE
  )
  
  b[valid_indexes, 1] <- predict(b_model, newdata = valid)
  
  q_model <- mlr3superlearner(
    data = filter(train, project == 1) |> 
      select(all_of(c(vars$V2, vars$W))), 
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
  
  t_model <- mlr3superlearner(
    data = select(train, all_of(c(vars$S, vars$W))), 
    target = vars$S, 
    outcome_type = "binomial",
    library = learners,
    discrete = FALSE
  )
  
  th[valid_indexes, 1] <- predict(t_model, newdata = valid)
  
  e_model <- mlr3superlearner(
    data = filter(train, week_12_relapse == 1) |> 
      select(all_of(c(vars$S, vars$W, vars$A))), 
    target = vars$S, 
    outcome_type = "binomial",
    library = learners,
    discrete = FALSE
  )

  e[valid_indexes, 1] <- predict(e_model, newdata = valid)
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
mc[, "A=0"] <- isoreg_with_xgboost(m[Ia1 == 0, "A=0"], data$week_12_relapse[Ia1 == 0])(m[, "A=0"])
mc[, "A=1"] <- isoreg_with_xgboost(m[Ia1 == 1, "A=1"], data$week_12_relapse[Ia1 == 1])(m[, "A=1"])
gc[, 1] <- isoreg_with_xgboost(g, data$bupnx)(g)
qc[, 1] <- isoreg_with_xgboost(q[Is1 == 1], Iv2_1[Is1 == 1])(q)
thc[, 1] <- isoreg_with_xgboost(th, Is1)(th)
bc[, 1] <- isoreg_with_xgboost(b[Iy1*Is1 == 1], Iv2_1[Iy1*Is1 == 1])(b)
ec[, 1] <- isoreg_with_xgboost(e[Iy1 == 1], Is1[Iy1 == 1])(e)

# truncating weights
gc <- pmin(pmax(gc, 0.025), 0.975)
thc <- pmax(thc, 0.025)
ec <- pmax(ec, 0.025)

# riesz representers
alpha_m1 <- Ia1 / g
alpha_m0 <- Ia0 / (1 - g)
alpha_q <- Is1 / th
alpha_b <- Is1 / e

# uncentered eifs
D_m1_uc <- alpha_m1*(data$week_12_relapse - m[, "A=1"]) + m[, "A=1"]
D_m0_uc <- alpha_m0*(data$week_12_relapse - m[, "A=0"]) + m[, "A=0"]
D_q1_uc <- alpha_q*(Iv2_1 - q) + q
D_q0_uc <- alpha_q*(Iv2_0 - (1 - q)) + (1 - q)
D_b1_uc <- alpha_b*(Iv2_1 - b) + b
D_b0_uc <- alpha_b*(Iv2_0 - (1 - b)) + (1 - b)

fit_softbart_mlr3 <- function(data, target) {
  task <- as_task_regr(data, target)
  learner <- LearnerRegrSoftBart$new()
  learner$train(task)
  return(learner)
}

predict_softbart_mlr3 <- function(object, newdata) {
  object$predict_newdata(newdata)$response
}

mbar_model_1 <- fit_softbart_mlr3(
  mutate(data, D_m1_uc = D_m1_uc) |> 
    select(D_m1_uc, all_of(vars$V1)), 
  "D_m1_uc"
)

# mbar_model_1 <- mlr3superlearner(
#   data = mutate(data, D_m1_uc = D_m1_uc) |> 
#     select(D_m1_uc, all_of(vars$V1)), 
#   target = "D_m1_uc", 
#   outcome_type = "continuous", 
#   library = c("glm", "earth", "ranger", "bart"),
#   discrete = FALSE
# )

mbar_model_0 <- fit_softbart_mlr3(
  mutate(data, D_m0_uc = D_m0_uc) |> 
    select(D_m0_uc, all_of(vars$V1)), 
  "D_m0_uc"
)

# mbar_model_0 <- mlr3superlearner(
#   data = mutate(data, D_m0_uc = D_m0_uc) |> 
#     select(D_m0_uc, all_of(vars$V1)), 
#   target = "D_m0_uc", 
#   outcome_type = "continuous", 
#   library = c("glm", "earth", "ranger", "bart"),
#   discrete = FALSE
# )

qbar_model_1 <- fit_softbart_mlr3(
  mutate(data, D_q1_uc = D_q1_uc) |> 
    select(D_q1_uc, all_of(vars$V1)), 
  "D_q1_uc"
)

# qbar_model_1 <- mlr3superlearner(
#   data = mutate(data, D_q1_uc = D_q1_uc) |> 
#     select(D_q1_uc, all_of(vars$V1)), 
#   target = "D_q1_uc", 
#   outcome_type = "continuous", 
#   library = c("glm", "earth", "ranger", "bart"),
#   discrete = FALSE
# )

qbar_model_0 <- fit_softbart_mlr3(
  data = mutate(data, D_q0_uc = D_q0_uc) |> 
    select(D_q0_uc, all_of(vars$V1)), 
  target = "D_q0_uc"
)

# qbar_model_0 <- mlr3superlearner(
#   data = mutate(data, D_q0_uc = D_q0_uc) |> 
#     select(D_q0_uc, all_of(vars$V1)), 
#   target = "D_q0_uc", 
#   outcome_type = "continuous", 
#   library = c("glm", "earth", "ranger", "bart"),
#   discrete = FALSE
# )

bbar_model_1 <- fit_softbart_mlr3(
  data = mutate(data, D_b1_uc = D_b1_uc) |> 
    filter(week_12_relapse == 1) |> 
    select(D_b1_uc, all_of(c(vars$A, vars$V1))), 
  target = "D_b1_uc"
)

# bbar_model_1 <- mlr3superlearner(
#   data = mutate(data, D_b1_uc = D_b1_uc) |> 
#     filter(week_12_relapse == 1) |> 
#     select(D_b1_uc, all_of(c(vars$A, vars$V1))), 
#   target = "D_b1_uc", 
#   outcome_type = "continuous", 
#   library = c("glm", "earth", "ranger", "bart"),
#   discrete = FALSE
# )

bbar_model_0 <- fit_softbart_mlr3(
  data = mutate(data, D_b0_uc = D_b0_uc) |>
    filter(week_12_relapse == 1) |> 
    select(D_b0_uc, all_of(c(vars$A, vars$V1))), 
  target = "D_b0_uc"
)

# bbar_model_0 <- mlr3superlearner(
#   data = mutate(data, D_b0_uc = D_b0_uc) |>
#     filter(week_12_relapse == 1) |> 
#     select(D_b0_uc, all_of(c(vars$A, vars$V1))), 
#   target = "D_b0_uc", 
#   outcome_type = "continuous", 
#   library = c("glm", "earth", "ranger", "bart"),
#   discrete = FALSE
# )

cate <- function(newdata) {
  mbar_1 <- predict_softbart_mlr3(mbar_model_1, newdata = newdata)
  mbar_0 <- predict_softbart_mlr3(mbar_model_0, newdata = newdata)
  
  if (newdata$homeless == 1) {
    bbar_1 <- predict_softbart_mlr3(bbar_model_1, newdata = mutate(newdata, bupnx = 1))
    bbar_0 <- predict_softbart_mlr3(bbar_model_1, newdata = mutate(newdata, bupnx = 0))
    qbar <- predict_softbart_mlr3(qbar_model_1, newdata = newdata)
  }
  
  if (newdata$homeless == 0) {
    bbar_1 <- predict_softbart_mlr3(bbar_model_0, newdata = mutate(newdata, bupnx = 1))
    bbar_0 <- predict_softbart_mlr3(bbar_model_0, newdata = mutate(newdata, bupnx = 0))
    qbar <- predict_softbart_mlr3(qbar_model_0, newdata = newdata)
  }
  
  (bbar_1*mbar_1 - bbar_0*mbar_0) / qbar
}

s0_v2_1 <- s0_v2_0 <- 
  filter(data, project == 0) |> 
  select(who, all_of(c(vars$V1, vars$V2)))

s0_v2_1$homeless <- 1
s0_v2_0$homeless <- 0

lookup <- 
  filter(data, project == 1) |> 
  select(who, all_of(c(vars$V1, vars$V2))) |> 
  rbind(s0_v2_1, s0_v2_0) |> 
  unique()

cates <- vector("numeric", nrow(lookup))
for (i in 1:nrow(lookup)) {
  cates[i] <- cate(lookup[i, ])
}

lookup$cate <- cates

s0 <- left_join(s0_v2_1, lookup) |> 
  rename(cate_1 = cate) |> 
  select(who, cate_1) |> 
  cbind({
    left_join(s0_v2_0, lookup) |> 
      rename(cate_0 = cate) |> 
      select(cate_0)
  })

s0$cate_min <- pmin(s0$cate_0, s0$cate_1)
s0$cate_max <- pmax(s0$cate_0, s0$cate_1)

shifted_d <- 
  left_join(filter(data, project == 1), lookup) |> 
  mutate(bupnx = as.numeric(cate < 0))

lmtp_tmle(
  data,
  vars$A, 
  vars$Y, 
  vars$W,
  shift = static_binary_on
)

