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

b <- ms1 <- m <- matrix(nrow = n, ncol = 2)
colnames(ms1) <- colnames(b) <- colnames(m) <- c("A=0", "A=1")

gs <- g <- q <- th <- matrix(nrow = n, ncol = 1)

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
  
  ms_model <- cv_glmnet_formula(Y ~ -1 + W1*W2*A*V1_1*V1_2*V1_3, 
                               data = filter(train, S == 1), family = "binomial")
  
  ms1[valid_indexes, "A=0"] <- predict(ms_model, newdata = valid_0)
  ms1[valid_indexes, "A=1"] <- predict(ms_model, newdata = valid_1)
  
  b_model <- cv_glmnet_formula(V2 ~ -1 + W1*W2*A*V1_1*V1_2*V1_3, 
                               data = filter(train, Y == 1 & S == 1), 
                               family = "binomial")
  
  b[valid_indexes, "A=0"] <- predict(b_model, newdata = valid_0)
  b[valid_indexes, "A=1"] <- predict(b_model, newdata = valid_1)
  
  q_model <- cv_glmnet_formula(V2 ~ W1*W2*V1_1*V1_2*V1_3, 
                               data = filter(train, S == 1), 
                               family = "binomial")
  
  q[valid_indexes, 1] <- predict(q_model, newdata = valid)
  
  g_model <- cv_glmnet_formula(A ~ -1 + W1*W2*V1_1*V1_2*V1_3, family = "binomial", data = train)
  g[valid_indexes, 1] <- predict(g_model, newdata = valid)
  
  gs_model <- cv_glmnet_formula(A ~ -1 + W1*W2*V1_1*V1_2*V1_3, family = "binomial", 
                                data = filter(train, S == 1))
  gs[valid_indexes, 1] <- predict(gs_model, newdata = valid)
  
  t_model <- cv_glmnet_formula(S ~ -1 + W1*W2*V1_1*V1_2*V1_3, family = "binomial", data = train)
  th[valid_indexes, 1] <- predict(t_model, newdata = valid)
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
ms1[, "A=0"] <- isoreg_with_xgboost(ms1[Ia0*Is1 == 1, "A=0"], data$Y[Ia0*Is1 == 1])(ms1[, "A=0"])
ms1[, "A=1"] <- isoreg_with_xgboost(ms1[Ia1*Is1 == 1, "A=1"], data$Y[Ia1*Is1 == 1])(ms1[, "A=1"])
g[, 1] <- isoreg_with_xgboost(g, data$A)(g)
gs[, 1] <- isoreg_with_xgboost(gs[Is1 == 1], data$A[Is1 == 1])(gs)
q[, 1] <- isoreg_with_xgboost(q[Is1 == 1], Iv2_1[Is1 == 1])(q)
th[, 1] <- isoreg_with_xgboost(th, Is1)(th)
b[, "A=0"] <- isoreg_with_xgboost(b[Ia0*Iy1*Is1 == 1, "A=0"], Iv2_1[Ia0*Iy1*Is1 == 1])(b[, "A=0"])
b[, "A=1"] <- isoreg_with_xgboost(b[Ia1*Iy1*Is1 == 1, "A=1"], Iv2_1[Ia1*Iy1*Is1 == 1])(b[, "A=1"])

# riesz representers
alpha_m1 <- Ia1 / g
alpha_m0 <- Ia0 / (1 - g)
alpha_ms1 <- Ia1 / gs
alpha_ms0 <- Ia0 / (1 - gs)
alpha_q <- Is1 / th
alpha_b <- alpha_q
alpha_foo1 <- Iy1 / ms1[, "A=1"]
alpha_foo0 <- Iy1 / ms1[, "A=0"]

# uncentered eifs
D_k1_v21_uc <- alpha_foo1*alpha_ms1*alpha_b*(m[, "A=1"]*Iv2_1 - m[, "A=1"]*b[, "A=1"]) + 
  alpha_m1*(b[, "A=1"]*data$Y - b[, "A=1"]*m[, "A=1"]) +
  b[, "A=1"]*m[, "A=1"]

D_k0_v21_uc <- alpha_foo0*alpha_ms0*alpha_b*(m[, "A=0"]*Iv2_1 - m[, "A=0"]*b[, "A=0"]) + 
  alpha_m0*(b[, "A=0"]*data$Y - b[, "A=0"]*m[, "A=0"]) +
  b[, "A=0"]*m[, "A=0"]

D_k1_v20_uc <- alpha_foo1*alpha_ms1*alpha_b*(m[, "A=1"]*Iv2_0 - m[, "A=1"]*(1 - b[, "A=1"])) + 
  alpha_m1*((1 - b[, "A=1"])*data$Y - (1 - b[, "A=1"])*m[, "A=1"]) +
  (1 - b[, "A=1"])*m[, "A=1"]

D_k0_v20_uc <- alpha_foo0*alpha_ms0*alpha_b*(m[, "A=0"]*Iv2_0 - m[, "A=0"]*(1 - b[, "A=0"])) + 
  alpha_m0*((1 - b[, "A=0"])*data$Y - (1 - b[, "A=0"])*m[, "A=0"]) +
  (1 - b[, "A=0"])*m[, "A=0"]

# D_q1_uc <- alpha_q*(Iv2_1 - q) + q
# D_q0_uc <- alpha_q*(Iv2_0 - (1 - q)) + (1 - q)

# Fit pseudo regressions
kappa_model_1_1 <- lm(D_k1_v21_uc ~ V1_1*V1_2*V1_3, data = data)
kappa_model_1_0 <- lm(D_k0_v21_uc ~ V1_1*V1_2*V1_3, data = data)
kappa_model_0_1 <- lm(D_k1_v20_uc ~ V1_1*V1_2*V1_3, data = data)
kappa_model_0_0 <- lm(D_k0_v20_uc ~ V1_1*V1_2*V1_3, data = data)
# lambda_model_1 <- lm(D_q1_uc ~ V1_1*V1_2*V1_3, data = data)
# lambda_model_0 <- lm(D_q0_uc ~ V1_1*V1_2*V1_3, data = data)

pseudo_cate <- function(newdata) {
  if (newdata$V2 == 1) {
    hbar_1 <- predict(kappa_model_1_1, newdata = newdata)
    hbar_0 <- predict(kappa_model_1_0, newdata = newdata)
    # qbar <- predict(lambda_model_1, newdata = newdata)
  }
  
  if (newdata$V2 == 0) {
    hbar_1 <- predict(kappa_model_0_1, newdata = newdata)
    hbar_0 <- predict(kappa_model_0_0, newdata = newdata)
    # qbar <- predict(lambda_model_0, newdata = newdata)
  }
  
  hbar_1 - hbar_0
}

res <- expand.grid(
  V1_1 = 0:1, V1_2 = 0:1, V1_3 = 0:1, V2 = 0:1
)

pseudo_cates <- vector("numeric", nrow(res))
for (i in 1:nrow(res)) {
  pseudo_cates[i] <- pseudo_cate(res[i, ])
}

res$psuedo_cate <- pseudo_cates

saveRDS(res, glue("../../data/sim/sim_numeratorLearner_dag2_{n}_{id}.rds"))
