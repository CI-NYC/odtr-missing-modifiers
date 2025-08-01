library(data.table)

source("R/dgp-01.R")

set.seed(54334)

tmp1 <- generate(1e7, FALSE, 1)
tmp0 <- generate(1e7, FALSE, 0)

setDT(tmp1)
setDT(tmp0)

mean_of_y <- function(a, v1_1, v1_2, v1_3, v2) {
  if (a == 1) {
    tmp1[V1_1 == v1_1 & V1_2 == v1_2 & V1_3 == v1_3 & V2 == v2, mean(Y)]
  } else {
    tmp0[V1_1 == v1_1 & V1_2 == v1_2 & V1_3 == v1_3 & V2 == v2, mean(Y)]
  }
}

prob_v2 <- function(v2, v1_1, v1_2, v1_3) {
  tmp1[V1_1 == v1_1 & V1_2 == v1_2 & V1_3 == v1_3, mean(V2 == v2)]
}

pseudo_CATE <- function(v1_1, v1_2, v1_3, v2) {
  cate <- mean_of_y(1, v1_1, v1_2, v1_3, v2) - mean_of_y(0, v1_1, v1_2, v1_3, v2)
  round(cate * prob_v2(v2, v1_1, v1_2, v1_3), 3)
}

vals <- expand.grid(
  v1_1 = c(0, 1),
  v1_2 = c(0, 1),
  v1_3 = c(0, 1),
  v2 = c(0, 1)
)

vals$pseudo_cate <- unlist(.mapply(pseudo_CATE, dots = vals, MoreArgs = NULL))
