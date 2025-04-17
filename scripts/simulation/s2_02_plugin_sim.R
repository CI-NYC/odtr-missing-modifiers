library(origami)
library(glmnet)
library(dplyr)
library(glue)

source("../../R/dgp-02.R")
source("../../R/glmnet-formula.R")
source("../../R/isotonic-calibration.R")

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (id == "undefined" || id == "") id <- 1

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- list(1e5, 2)
}

n <- as.numeric(args[[1]])
folds <- as.numeric(args[[2]])

set.seed(id)
data <- generate(n)

crossfit_folds <- make_folds(data, V = folds)

b <- m <- matrix(nrow = n, ncol = 2)
colnames(b) <- colnames(m) <- c("A=0", "A=1")

q <- matrix(nrow = n, ncol = 1)

for (v in 1:folds) {
  train <- training(data, crossfit_folds[[v]])
  valid <- validation(data, crossfit_folds[[v]])
  valid_1 <- mutate(valid, A = 1)
  valid_0 <- mutate(valid, A = 0)
  
  valid_indexes <- crossfit_folds[[v]]$validation_set
  
  m_model <- cv_glmnet_formula(Y ~ -1 + W1*W2*A*V1_1*V1_2*V1_3, 
                               data = train, family = "binomial")
  
  m[valid_indexes, "A=0"] <- predict(m_model, newdata = valid_0)
  m[valid_indexes, "A=1"] <- predict(m_model, newdata = valid_1)
  
  b_model <- cv_glmnet_formula(V2 ~ -1 + W1*W2*A*V1_1*V1_2*V1_3, 
                               data = filter(train, Y == 1 & S == 1), 
                               family = "binomial")
  
  b[valid_indexes, "A=0"] <- predict(b_model, newdata = valid_0)
  b[valid_indexes, "A=1"] <- predict(b_model, newdata = valid_1)
  
  q_model <- cv_glmnet_formula(V2 ~ W1*W2*V1_1*V1_2*V1_3, 
                               data = filter(train, S == 1), 
                               family = "binomial")
  
  q[valid_indexes, 1] <- predict(q_model, newdata = valid)
}

Ia1 <- as.numeric(data$A == 1)
Ia0 <- as.numeric(data$A == 0)
Is1 <- as.numeric(data$S == 1)
Iv2_1 <- as.numeric(data$V2 == 1)
Iv2_1[is.na(Iv2_1)] <- -999
Iv2_0 <- as.numeric(data$V2 == 0)
Iv2_0[is.na(Iv2_0)] <- -999
Iy1 <- as.numeric(data$Y == 1)

# isotonic calibrations
m[, "A=0"] <- isoreg_with_xgboost(m[Ia1 == 0, "A=0"], data$Y[Ia1 == 0])(m[, "A=0"])
m[, "A=1"] <- isoreg_with_xgboost(m[Ia1 == 1, "A=1"], data$Y[Ia1 == 1])(m[, "A=1"])
q[, 1] <- isoreg_with_xgboost(q[Is1 == 1], Iv2_1[Is1 == 1])(q)
b[, "A=0"] <- isoreg_with_xgboost(b[Ia0*Iy1*Is1 == 1, "A=0"], Iv2_1[Ia0*Iy1*Is1 == 1])(b[, "A=0"])
b[, "A=1"] <- isoreg_with_xgboost(b[Ia1*Iy1*Is1 == 1, "A=1"], Iv2_1[Ia1*Iy1*Is1 == 1])(b[, "A=1"])

# Fit pseudo regressions
h_1_1 <- lm(m[, "A=1"]*b[, "A=1"] ~ V1_1*V1_2*V1_3, data = data)
h_1_0 <- lm(m[, "A=0"]*b[, "A=0"] ~ V1_1*V1_2*V1_3, data = data)
h_0_1 <- lm(m[, "A=1"]*(1 - b[, "A=1"]) ~ V1_1*V1_2*V1_3, data = data)
h_0_0 <- lm(m[, "A=0"]*(1 - b[, "A=0"]) ~ V1_1*V1_2*V1_3, data = data)
qbar <- lm(q ~ V1_1*V1_2*V1_3, data = data)

cate <- function(newdata) {
  hbar_1_v21 <- predict(h_1_1, newdata = newdata)
  hbar_0_v21 <- predict(h_1_0, newdata = newdata)
  qbar_v21 <- predict(qbar, newdata = newdata)
  
  hbar_1_v20 <- predict(h_0_1, newdata = newdata)
  hbar_0_v20 <- predict(h_0_0, newdata = newdata)
  qbar_v20 <- 1 - qbar_v21
  
  psi_v21 <- (hbar_1_v21 - hbar_0_v21) / qbar_v21
  psi_v20 <- (hbar_1_v20 - hbar_0_v20) / qbar_v20
  
  newdata$V2*psi_v21 + (1 - newdata$V2)*psi_v20
}

res <- expand.grid(V1_1 = 0:1, V1_2 = 0:1, V1_3 = 0:1, V2 = 0:1)
res$cate <- cate(res)

saveRDS(res, glue("../../data/sim/sim_plugin_dag2_{n}_{id}.rds"))
