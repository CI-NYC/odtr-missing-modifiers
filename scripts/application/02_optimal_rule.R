library(mlr3superlearner)
library(mlr3extralearners)
library(lmtp)
library(dplyr)
library(kernelshap)
library(shapviz)
library(ggplot2)
library(future)
library(patchwork)
library(purrr)

source("R/load-data.R")
source("R/fit-tuned-xgboost.R")

pt_to_mm <- function(pt) {
  pt / ggplot2::.pt
}

# load CTN data
data <- load_analysis_data()
vars <- data$vars
data <- data$data

# load pseudo regression models
models <- readRDS("data/application/drv/psuedo-regression-tuned-xgboost-models.rds")

cate <- function(object, newdata) {
  hbar_1_homeless <- predict_tuned_xgboost(object$kappa_model_1_1, newdata = newdata)
  hbar_0_homeless <- predict_tuned_xgboost(object$kappa_model_1_0, newdata = newdata)
  qbar_homeless <- predict_tuned_xgboost(object$lambda_model_1, newdata = newdata)
  
  hbar_1_not_homeless <- predict_tuned_xgboost(object$kappa_model_0_1, newdata = newdata)
  hbar_0_not_homeless <- predict_tuned_xgboost(object$kappa_model_0_0, newdata = newdata)
  qbar_not_homeless <- predict_tuned_xgboost(object$lambda_model_0, newdata = newdata)
  
  psi_homeless <- (hbar_1_homeless - hbar_0_homeless) / qbar_homeless
  psi_not_homeless <- (hbar_1_not_homeless - hbar_0_not_homeless) / qbar_not_homeless
  
  newdata$homeless*psi_homeless + 
    (1 - newdata$homeless)*psi_not_homeless
}

ctn30_homeless <- 
  ctn30_no_homeless <-
  filter(data, project == 0) |>
  select(who, all_of(c(vars$V1, vars$V2)))

ctn30_homeless$homeless <- 1
ctn30_no_homeless$homeless <- 0

lookup_cate <-
  filter(data, project == 1) |>
  select(who, all_of(c(vars$V1, vars$V2))) |>
  rbind(ctn30_homeless, ctn30_no_homeless) |>
  select(-who) |> 
  unique()

# shap values -------------------------------------------------------------

off_row <- which(rowSums(lookup_cate == 0) == ncol(lookup_cate))[1]
shap <- kernelshap(object = models, pred_fun = cate, 
                   X = lookup_cate[-off_row, ], 
                   bg_X = lookup_cate[off_row, ])

shap_beeswarm <- shapviz(shap) |> 
  sv_importance(kind = "bee") + 
  scale_y_discrete(labels = c(
    "homeless" = "Homeless",
    "ivdrug" = "IV drug use", 
    "hwithdraw_4" = "Severe withdrawal", 
    "hwithdraw_3" = "Medium withdrawal",
    "hwithdraw_2" = "Mild withdrawal",
    "bbenzo30_base" = "Benzo. use", 
    "cocdisorder" = "Cocaine use disorder", 
    "xrace_3" = "Hispanic", 
    "hasEpilepsy" = "Epilepsy", 
    "xrace_4" = "Other race", 
    "xrace_2" = "Non-Hispanic, Black"
  )) + 
  theme_classic(base_size = 8,
                base_line_size = 0.4,
                base_rect_size = 0.4) +
  theme(
    legend.position = c(0.9, 0.5), 
    legend.background = element_rect(fill = alpha('white', 0.33)), 
    text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.title = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.line = element_line(size = pt_to_mm(1)),
    axis.ticks = element_line(size = pt_to_mm(1))
  )

shap_importance <- shapviz(shap) |> 
  sv_importance(show_numbers = FALSE, 
                number_size = 1.5) + 
  scale_y_discrete(labels = c(
    "homeless" = "Homeless",
    "ivdrug" = "IV drug use", 
    "hwithdraw_4" = "Severe withdrawal", 
    "hwithdraw_3" = "Medium withdrawal",
    "hwithdraw_2" = "Mild withdrawal",
    "bbenzo30_base" = "Benzo. use", 
    "cocdisorder" = "Cocaine use disorder", 
    "xrace_3" = "Hispanic", 
    "hasEpilepsy" = "Epilepsy", 
    "xrace_4" = "Other race", 
    "xrace_2" = "Non-Hispanic, Black"
  )) + 
  theme_classic() + 
  theme(
    axis.title.y = element_blank(),        # Remove y-axis title
    axis.text.y = element_blank(),         # Remove y-axis text labels
    axis.ticks.y = element_blank(),        # Remove y-axis ticks
    axis.line.y = element_blank(),         # Remove y-axis line
    text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.text.x = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.title.x = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.line.x = element_line(size = pt_to_mm(1)),
    axis.ticks.x = element_line(size = pt_to_mm(1))
  ) + 
  scale_x_continuous(expand = expansion(mult = c(0, .15)))

png("figures/shap-values.png", width = 6, height = 3, units = "in", res = 600)
shap_beeswarm + shap_importance + plot_layout(widths = c(2/3, 1/3))
dev.off()

lookup_cate[-off_row, ] |> 
  mutate(row = row_number()) |> 
  arrange(desc(homeless), desc(ivdrug), desc(hwithdraw_4), desc(hwithdraw_3), desc(cocdisorder))

lookup_cate[-off_row, ] |> 
  mutate(row = row_number()) |> 
  filter(homeless == 1, ivdrug == 1, hwithdraw_4 == 1, cocdisorder == 1, bbenzo30_base == 1) |> 
  (\(x) x$row)()

  arrange(desc(homeless), desc(ivdrug), desc(hwithdraw_4), desc(hwithdraw_3), desc(cocdisorder))


sv_waterfall(shapviz(shap), row_id = .Last.value, max_display = Inf) +
  theme(axis.text = element_text(size = 11))

negative_homeless_shap_rows <- which(shap$S[, "homeless"] < 0)
positive_homeless_shap_rows <- which(shap$S[, "homeless"] > 0)
summary(lookup_cate[-off_row, ][negative_homeless_shap_rows, ])
summary(lookup_cate[-off_row, ][positive_homeless_shap_rows, ])

# relapse under the rule --------------------------------------------------

lookup_cate$cate <- cate(models, lookup_cate)

ctn30 <- left_join(ctn30_homeless, lookup_cate) |>
  rename(cate_homeless = cate) |>
  select(who, cate_homeless) |>
  cbind({
    left_join(ctn30_no_homeless, lookup_cate) |>
      rename(cate_no_homeless = cate) |>
      select(cate_no_homeless)
  })

ctn30$cate_min <- pmin(ctn30$cate_no_homeless, ctn30$cate_homeless)
ctn30$cate_max <- pmax(ctn30$cate_no_homeless, ctn30$cate_homeless)

ctn51 <- filter(data, project == 1) |> 
  left_join(lookup_cate) |> 
  mutate(bupnx = as.numeric(cate < 0))

ctn30_decisive <- 
  filter(ctn30, sign(cate_min) == sign(cate_max)) |> 
  select(who, cate = cate_max) |> 
  inner_join(data, by = "who") |> 
  mutate(bupnx = as.numeric(cate < 0))

ctn30_ambiguous_bupnx <- 
  filter(ctn30, sign(cate_min) != sign(cate_max)) |> 
  select(who) |> 
  inner_join(data, by = "who") |> 
  mutate(bupnx = 1)

ctn30_ambiguous_no_bupnx <- 
  filter(ctn30, sign(cate_min) != sign(cate_max)) |> 
  select(who) |> 
  inner_join(data, by = "who") |> 
  mutate(bupnx = 0)

shifted_ambiguous_bupnx <- bind_rows(ctn51, ctn30_decisive, ctn30_ambiguous_bupnx)
shifted_ambiguous_no_bupnx <- bind_rows(ctn51, ctn30_decisive, ctn30_ambiguous_no_bupnx)

# superlearner library
learners <- c("mean", "glm", "earth", "ranger", "bart", "cv_glmnet")

lmtp_all_bupnx <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shift = static_binary_on, 
  learners_trt = learners, 
  learners_outcome = learners
)

saveRDS(lmtp_all_bupnx, "data/application/drv/lmtp-all-bupnx.rds")

lmtp_all_naltrexone <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shift = static_binary_off, 
  learners_trt = learners, 
  learners_outcome = learners
)

saveRDS(lmtp_all_naltrexone, "data/application/drv/lmtp-all-naltrexone.rds")

lmtp_ambiguous_bupnx <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shifted = arrange(shifted_ambiguous_bupnx, who), 
  learners_trt = learners, 
  learners_outcome = learners
)

saveRDS(lmtp_ambiguous_bupnx, "data/application/drv/lmtp-ambiguous-bupnx.rds")

lmtp_ambiguous_no_bupnx <- lmtp_tmle(
  arrange(data, who),
  vars$A,
  vars$Y,
  c(vars$W, vars$V1),
  shifted = arrange(shifted_ambiguous_no_bupnx, who), 
  learners_trt = learners, 
  learners_outcome = learners
)

saveRDS(lmtp_ambiguous_no_bupnx, "data/application/drv/lmtp-ambiguous-naltrexone.rds")

# paper results -----------------------------------------------------------

# % decreases
1 - lmtp_ambiguous_bupnx$estimate / lmtp_randomized$estimate
1 - lmtp_ambiguous_bupnx$estimate / lmtp_all_bupnx$estimate
1 - lmtp_ambiguous_bupnx$estimate / lmtp_all_naltrexone$estimate

1 - lmtp_ambiguous_no_bupnx$estimate / lmtp_randomized$estimate
1 - lmtp_ambiguous_no_bupnx$estimate / lmtp_all_bupnx$estimate
1 - lmtp_ambiguous_no_bupnx$estimate / lmtp_all_naltrexone$estimate

all_point_estimates <- 
  map(list("BUP-NX" = lmtp_all_bupnx, 
           "XR-NTX" = lmtp_all_naltrexone, 
           "d1" = lmtp_ambiguous_bupnx, 
           "d0" = lmtp_ambiguous_no_bupnx), 
      tidy) |> 
  list_rbind(names_to = "rule")

point_estimate_plot <- ggplot(all_point_estimates, aes(
  x = factor(rule, levels = c("XR-NTX", "BUP-NX", "d1", "d0")),
  y = estimate
)) +
  geom_point(position = position_dodge(.75), size = pt_to_mm(3)) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0.15,
    position = position_dodge(.75),
    linewidth = pt_to_mm(1)
  ) +
  coord_cartesian(ylim = c(0.2, 1)) +
  labs(
    x = "Treatment rule",
    y = "Expected Risk of Relapse by 12 Weeks"
  ) +
  theme_classic(base_size = 8,
                base_line_size = 0.4,
                base_rect_size = 0.4) +
  theme(legend.position = c(0.15, 0.775),
        legend.background = element_rect(fill = "white",
                                         color = "black",
                                         size = pt_to_mm(1)),
        text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        axis.text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        axis.title = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        axis.line = element_line(size = pt_to_mm(1)),
        axis.ticks = element_line(size = pt_to_mm(1)),
        legend.text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        legend.text.align = 0,
        legend.title.align = 0.5)

all_risk_difference_estimates <- 
  map(list("d1 - BUP-NX" = lmtp_ambiguous_bupnx$estimate - lmtp_all_bupnx$estimate, 
           "d1 - XR-NTX" = lmtp_ambiguous_bupnx$estimate - lmtp_all_naltrexone$estimate, 
           "d0 - BUP-NX" = lmtp_ambiguous_no_bupnx$estimate - lmtp_all_bupnx$estimate, 
           "d0 - XR-NTX" = lmtp_ambiguous_no_bupnx$estimate - lmtp_all_naltrexone$estimate), 
      ife::tidy) |> 
  list_rbind(names_to = "contrast")

risk_difference_plot <- 
  ggplot(all_risk_difference_estimates, aes(
    x = factor(contrast, levels = c("d1 - XR-NTX", "d0 − XR-NTX", "d1 - BUP-NX", "d0 - BUP-NX")),
    y = estimate
  )) +
  geom_point(position = position_dodge(0.75), size = pt_to_mm(3)) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
    ),
    width = 0.15,
    position = position_dodge(0.75),
    linewidth = pt_to_mm(1)
  ) +
  geom_hline(yintercept = 0, size = pt_to_mm(1), linetype = "dashed") +
  scale_y_continuous(breaks = c(0.2, 0, -0.2, -0.4),
                     labels = c("0.2", "0", "-0.2", "-0.4")) +
  #scale_x_discrete(guide = guide_axis(angle = -45)) +
  labs(
    x = "Contrast",
    y = "Expected Difference in Risk \nof Relapse by 12 Weeks",
    linetype = NULL
  ) +
  coord_cartesian(ylim = c(-0.5, 0.1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        axis.text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        axis.title = element_text(size = 7, color = "black", face = "plain", family = "sans"),
        axis.line = element_line(size = pt_to_mm(1)),
        axis.ticks = element_line(size = pt_to_mm(1)))

png("figures/point-estimates.png", width = 6, height = 3, units = "in", res = 600)
print({
  point_estimate_plot + risk_difference_plot + 
    plot_layout(widths = c(0.4, 0.6))
})
dev.off()
