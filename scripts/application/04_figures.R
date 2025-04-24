library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(lmtp)

pt_to_mm <- function(pt) {
  pt / ggplot2::.pt
}

# Load all lmtp point estimates under the given rule
estimate_bupnx <- readRDS("data/application/drv/lmtp-all-bupnx.rds")
estimate_xrntx <- readRDS("data/application/drv/lmtp-all-naltrexone.rds")
estimate_odtr_no_homeless <- readRDS("data/application/drv/lmtp-no-homeless.rds")
estimate_odtr_ambiguous_bupnx <- readRDS("data/application/drv/lmtp-ambiguous-bupnx-updated.rds")
estimate_odtr_ambiguous_xrntx <- readRDS("data/application/drv/lmtp-ambiguous-naltrexone-updated.rds")

# # % decreases
# 1 - lmtp_ambiguous_bupnx$estimate / lmtp_all_bupnx$estimate
# 1 - lmtp_ambiguous_bupnx$estimate / lmtp_all_naltrexone$estimate
# 
# 1 - lmtp_ambiguous_no_bupnx$estimate / lmtp_all_bupnx$estimate
# 1 - lmtp_ambiguous_no_bupnx$estimate / lmtp_all_naltrexone$estimate

# Extract point estimates into a data.frame
all_point_estimates <- 
  map(list("BUP-NX" = estimate_bupnx, 
           "XR-NTX" = estimate_xrntx, 
           "d*" = estimate_odtr_no_homeless,
           "d1" = estimate_odtr_ambiguous_bupnx, 
           "d0" = estimate_odtr_ambiguous_xrntx), 
      tidy) |> 
  list_rbind(names_to = "rule")

percent_decrease_estimates <- 
  map(list("d* / BUP-NX" = 1 - estimate_odtr_no_homeless$estimate / estimate_bupnx$estimate, 
           "d* / XR-NTX" = 1 - estimate_odtr_no_homeless$estimate / estimate_xrntx$estimate, 
           "d1 / BUP-NX" = 1 - estimate_odtr_ambiguous_bupnx$estimate / estimate_bupnx$estimate, 
           "d1 / XR-NTX" = 1 - estimate_odtr_ambiguous_bupnx$estimate / estimate_xrntx$estimate, 
           "d0 / BUP-NX" = 1 - estimate_odtr_ambiguous_xrntx$estimate / estimate_bupnx$estimate, 
           "d0 / XR-NTX" = 1 - estimate_odtr_ambiguous_xrntx$estimate / estimate_xrntx$estimate), 
      ife::tidy) |> 
  list_rbind(names_to = "contrast")

# Plot of point estimates with confidence intervals
point_estimate_plot <- 
  ggplot(
    all_point_estimates, 
    aes(x = factor(rule, levels = c("XR-NTX", "BUP-NX", "d*", "d1", "d0")), y = estimate)
  ) +
  geom_point(position = position_dodge(.75), size = pt_to_mm(3)) +
  geom_errorbar(
    aes(ymin = conf.low,ymax = conf.high),
    width = 0.15,
    position = position_dodge(.75),
    linewidth = pt_to_mm(1)
  )

point_estimate_plot <- 
  point_estimate_plot + 
  coord_cartesian(ylim = c(0.2, 1)) +
  labs(
    x = NULL,
    y = "Expected Risk of Relapse by 12 Weeks"
  ) +
  theme_classic(
    base_size = 8,
    base_line_size = 0.4,
    base_rect_size = 0.4
  ) +
  theme(
    legend.position = c(0.15, 0.775),
    legend.background = element_rect(fill = "white", color = "black", size = pt_to_mm(1)),
    text = element_text(size = 7, color = "black", face = "plain", family = "sans"),
    axis.text = element_text(size = 6, color = "black", face = "plain", family = "sans"),
    axis.title = element_text(size = 6, color = "black", face = "plain", family = "sans"),
    axis.line = element_line(size = pt_to_mm(1)),
    axis.ticks = element_line(size = pt_to_mm(1)),
    legend.text = element_text(size = 6, color = "black", face = "plain", family = "sans"),
    legend.text.align = 0,
    legend.title.align = 0.5
  )

percent_decrease_plot <- 
  ggplot(
    percent_decrease_estimates, 
    aes(
      x = factor(
        contrast, 
        levels = c("d* / XR-NTX", "d1 / XR-NTX", "d0 / XR-NTX", "d* / BUP-NX", "d1 / BUP-NX", "d0 / BUP-NX")),
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
  scale_y_continuous(labels = scales::percent) + 
  labs(
    x = NULL,
    y = "Percent decrease in Expected Risk \nof Relapse by 12 Weeks",
    linetype = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 6, color = "black", face = "plain", family = "sans"),
    axis.text = element_text(size = 6, color = "black", face = "plain", family = "sans"),
    axis.title = element_text(size = 6, color = "black", face = "plain", family = "sans"),
    axis.line = element_line(size = pt_to_mm(1)),
    axis.ticks = element_line(size = pt_to_mm(1)), 
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

png("figures/point-estimates-updated.png", width = 6, height = 3, units = "in", res = 600)
print({
  point_estimate_plot + percent_decrease_plot + 
    plot_layout(widths = c(0.4, 0.6))
})
dev.off()
