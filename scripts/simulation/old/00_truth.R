library(tidyverse)

source("R/dgp.R")

params <- expand.grid(A = c(0, 1), V1 = c(0, 1), V2 = c(0, 1))
paramsS1 <- params 

for (i in 1:nrow(params)) {
  tmp <- generate(
    1e7,
    0.5,
    FALSE,
    A = params$A[i],
    V1 = params$V1[i],
    V2 = params$V2[i]
  )
  
  params$psi[i] <- mean(tmp$Y)
  paramsS1$psi[i] <- mean(tmp$Y[tmp$S == 1])
}

truth <- 
  pivot_wider(params, names_from = "A", values_from = "psi") |> 
  mutate(cate = `1` - `0`) |> 
  arrange(V1, V2)
