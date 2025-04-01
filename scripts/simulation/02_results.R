library(tidyverse)
library(glue)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

summarize_sim <- function(n) {
  res <- 
    read_zip_rds(glue("data/sim/sim_drLearner_components_0.75_{n}.zip")) |> 
    data.table::rbindlist()
  
  group_by(res, V1_1, V1_2, V1_3, V2) |>
    rename(v1_1 = V1_1, v1_2 = V1_2, v1_3 = V1_3, v2 = V2) |>
    left_join(rename(vals, truth = cate)) |> 
    summarize(abs_bias = abs(mean(cate) - mean(truth)),
              mse = mean((cate - truth)^2)) |> 
    ungroup() |> 
    mutate(n = n)
}

map(c(500, 1000, 2500, 10000), summarize_sim) |> 
  list_rbind() |> 
  group_by(n) |> 
  summarize(median_abs_bias = median(abs_bias),
            median_mse = median(mse))

            