library(tidyverse)
library(mice)
library(yaml)

# load yaml files of variables to use
vars <- read_yaml("scripts/application/vars.yaml")

# Extract homeless variable for CTN0051
homeless <- 
  read_csv("data/application/src/ctn51.csv") |> 
  select(PATID, homeless) |> 
  mutate(homeless = ifelse(homeless == "Yes", 1, 0), 
         PATID = factor(str_c("0", PATID))) |> 
  rename(who = PATID)

# clean data per AJE paper
load("data/application/src/clean_combined_imputed_data9-8-21.Rdata", verbose = TRUE)
raw_data <- ALT_patients_imputed_03[[1]]

weeks_to_relapse <- function(rand, relapse) {
  time <- lubridate::interval(rand, relapse)
  days_to_relapse <- lubridate::as.period(time, unit = "day")
  days_to_relapse / lubridate::as.period(lubridate::dweeks(1))
}

ctn_data <- 
  mutate(raw_data,
         across(starts_with("has"), ~ as.numeric(.x == "yes")),
         across(ends_with("disorder"), ~ as.numeric(.x == "yes")),
         across(c("bamphetamine30_base", "bcannabis30_base", "bbenzo30_base", "ivdrug"),
                ~ as.numeric(.x == "yes")),
         across(c("switched_meds", "never_initiated"), as.numeric),
         weeks_to_relapse = weeks_to_relapse(rand_dt, relapse_date),
         week_12_relapse = as.numeric(weeks_to_relapse <= 12),
         sex = if_else(sex == "female", 2, 1)) |> 
  select(all_of(c(vars$key, vars$S, vars$A, vars$Y, vars$W))) |> 
  filter(project %in% c("30", "51"), 
         !is.na(medicine))

# join homeless indicator 
ctn_data <- left_join(ctn_data, homeless)

ctn_data <- mutate(ctn_data, 
                   medicine = ifelse(medicine == "bup", 1, 0), 
                   sex = sex - 1)

saveRDS(ctn_data, "data/application/drv/ctn_data.rds")
