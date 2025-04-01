library(mingle)
# library(mlr3extralearners)
# library(mlr3pipelines)
library(glue)
library(ife)

source("../R/dgp.R")

# options(future.globals.maxSize = 2 * 1024^3) # Example: 2 GB (adjust as needed)

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (id == "undefined" || id == "") id <- 1

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- list(1500, 1)
}

# [500, 1000, 5000]
n <- as.numeric(args[[1]])
folds <- as.numeric(args[[2]])
# [-2, 0]
alpha <- 0.5 

# learners <- list("mean", "glm",
#                  list("glm", filter = po("modelmatrix", formula = ~ 0 + .^2)),
#                  list("earth", degree = 2), 
#                  list("earth", degree = 3))

# set.seed(id)
tmp <- generate(n, alpha)

set.seed(id)
res <- print(mingle_cate(
    tmp, 
    "A", "Y", "W", "S", "V1", "V2", 
    folds = folds, 
    learners_b = "glm",
    learners_m = "glm",
    learners_q = "glm",
    learners_t = "glm",
    learners_e = "glm",
    learners_g = "glm",
  )
)

saveRDS(res, glue("../data/sim_{estimator}_{alpha}_{n}_{id}.rds"))
