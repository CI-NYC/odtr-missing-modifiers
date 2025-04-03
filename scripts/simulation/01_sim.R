library(origami)
library(glmnet)
library(dplyr)
library(glue)

source("../../R/dgp.R")
source("../../R/glmnet-formula.R")
source("../../R/isotonic-calibration.R")

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (id == "undefined" || id == "") id <- 1

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- list(1500, 10)
}

n <- as.numeric(args[[1]])
folds <- as.numeric(args[[2]])
alpha <- 0.75

data <- generate(n, alpha)

crossfit_folds <- make_folds(data, V = folds)

m <- matrix(nrow = n, ncol = 2)
colnames(m) <- c("A=0", "A=1")

g <- q <- e <- b <- th <- matrix(nrow = n, ncol = 1)

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
  
  b[valid_indexes, 1] <- predict(b_model, newdata = valid)
  
  q_model <- cv_glmnet_formula(V2 ~ W1*W2*V1_1*V1_2*V1_3, 
                               data = filter(train, S == 1), 
                               family = "binomial")
  
  q[valid_indexes, 1] <- predict(q_model, newdata = valid)
  
  g_model <- cv_glmnet_formula(A ~ -1 + W1*W2*V1_1*V1_2*V1_3, family = "binomial", data = train)
  g[valid_indexes, 1] <- predict(g_model, newdata = valid)
  
  t_model <- cv_glmnet_formula(S ~ -1 + W1*W2*V1_1*V1_2*V1_3, family = "binomial", data = train)
  th[valid_indexes, 1] <- predict(t_model, newdata = valid)
  
  e_model <- cv_glmnet_formula(S ~  -1 + W1*W2*A*V1_1*V1_2*V1_3, family = "binomial", 
                 data = filter(train, Y == 1))
  e[valid_indexes, 1] <- predict(e_model, newdata = valid)
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
g[, 1] <- isoreg_with_xgboost(g, data$A)(g)
q[, 1] <- isoreg_with_xgboost(q[Is1 == 1], Iv2_1[Is1 == 1])(q)
th[, 1] <- isoreg_with_xgboost(th, Is1)(th)
b[, 1] <- isoreg_with_xgboost(b[Iy1*Is1 == 1], Iv2_1[Iy1*Is1 == 1])(b)
e[, 1] <- isoreg_with_xgboost(e[Iy1 == 1], Is1[Iy1 == 1])(e)

# riesz representers
alpha_m1 <- Ia1 / g
alpha_m0 <- Ia0 / (1 - g)
alpha_q <- Is1 / th
alpha_b <- Is1 / e

# uncentered eifs
D_m1_uc <- alpha_m1*(data$Y - m[, "A=1"]) + m[, "A=1"]
D_m0_uc <- alpha_m0*(data$Y - m[, "A=0"]) + m[, "A=0"]
D_q1_uc <- alpha_q*(Iv2_1 - q) + q
D_q0_uc <- alpha_q*(Iv2_0 - (1 - q)) + (1 - q)
D_b1_uc <- alpha_b*(Iv2_1 - b) + b
D_b0_uc <- alpha_b*(Iv2_0 - (1 - b)) + (1 - b)

mbar_model_1 <- glm(D_m1_uc ~ V1_1*V1_2*V1_3, data = data)
mbar_model_0 <- glm(D_m0_uc ~ V1_1*V1_2*V1_3, data = data)
qbar_model_1 <- glm(D_q1_uc ~ V1_1*V1_2*V1_3, data = data)
qbar_model_0 <- glm(D_q0_uc ~ V1_1*V1_2*V1_3, data = data)
bbar_model_1 <- glm(D_b1_uc ~ A*V1_1*V1_2*V1_3, data = data, subset = Iy1 == 1)
bbar_model_0 <- glm(D_b0_uc ~ A*V1_1*V1_2*V1_3, data = data, subset = Iy1 == 1)

cate <- function(newdata) {
  mbar_1 <- predict(mbar_model_1, newdata = newdata)
  mbar_0 <- predict(mbar_model_0, newdata = newdata)
  
  if (newdata$V2 == 1) {
    bbar_1 <- predict(bbar_model_1, newdata = mutate(newdata, A = 1))
    bbar_0 <- predict(bbar_model_1, newdata = mutate(newdata, A = 0))
    qbar <- predict(qbar_model_1, newdata = newdata)
  }
  
  if (newdata$V2 == 0) {
    bbar_1 <- predict(bbar_model_0, newdata = mutate(newdata, A = 1))
    bbar_0 <- predict(bbar_model_0, newdata = mutate(newdata, A = 0))
    qbar <- predict(qbar_model_0, newdata = newdata)
  }
  
  (bbar_1*mbar_1 - bbar_0*mbar_0) / qbar
}

res <- expand.grid(
  V1_1 = 0:1, V1_2 = 0:1, V1_3 = 0:1, V2 = 0:1
)

cates <- vector("numeric", nrow(res))
for (i in 1:nrow(res)) {
  cates[i] <- cate(res[i, ])
}

res$cate <- cates

saveRDS(res, glue("../../data/sim/sim_drcLearner_{alpha}_{n}_{id}.rds"))
