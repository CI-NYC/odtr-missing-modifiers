load_analysis_data <- function() {
  # load CTN data
  data <- readRDS("data/application/drv/ctn_data.rds")
  # load yaml files of variables to use
  vars <- yaml::read_yaml("scripts/application/vars.yaml")
  
  # Removes 5 observations (<1%)
  data <- dplyr::filter(data, !is.na(hasBipolar), !is.na(hasSchiz), !is.na(hwithdraw))
  
  # dummy coding categorical variables
  V1 <- fastDummies::dummy_cols(dplyr::select(data, dplyr::all_of(vars$V1)), 
                                remove_first_dummy = TRUE, 
                                remove_selected_columns = TRUE)
  
  data <- dplyr::bind_cols(dplyr::select(data, -dplyr::all_of(vars$V1)), V1)
  
  vars$V1 <- names(V1)
  list(data = data, vars = vars)
}