library(tidyverse)
library(glue)
library(data.table)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

# load true values
source("scripts/simulation/00_truth.R")

dag <- 1

summarize_sim <- function(n) {
  res <- 
    read_zip_rds(glue("data/sim/sim_drcLearner_dag{dag}_{n}.zip")) |> 
    data.table::rbindlist()
  
  group_by(res, V1_1, V1_2, V1_3, V2) |>
    rename(v1_1 = V1_1, v1_2 = V1_2, v1_3 = V1_3, v2 = V2) |>
    left_join(rename(vals, truth = cate)) |> 
    summarize(abs_bias = abs(mean(cate) - mean(truth)),
              rmse = sqrt(mean((cate - truth)^2))) |> 
    ungroup() |> 
    mutate(n = n)
}

res <- map(c(500, 1000, 2500, 10000), summarize_sim) |> 
  list_rbind()

tmp <- generate(1e7, FALSE)
setDT(tmp)

weights <- tmp[, .N, by = .(V1_1, V1_2, V1_3, V2)][, prop := N / sum(N)][]
setnames(weights, 
         c("V1_1", "V1_2", "V1_3", "V2"), 
         c("v1_1", "v1_2", "v1_3", "v2"))

weighted_res <- 
  left_join(res, weights) |> 
  group_by(n) |> 
  summarize(abs_bias = weighted.mean(abs_bias, prop),
            rmse = weighted.mean(rmse, prop))
