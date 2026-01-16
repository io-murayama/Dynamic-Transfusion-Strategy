### visualization (forest plots) ###
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

data_dir   <- './data/'
output_dir <- './output/'
date       <- '20251011'
time_window_width <- 24
cov_order_label <- "fwd"
sg <- "all"
xlim <- 0.04
x_breaks_by <- 0.01

sg_to_plot <- c(
  all                   = "All patients",
  age_70_or_older       = "Age ≥ 70",
  age_under_70          = "Age < 70",
  high_base_mbp         = "Initial MAP ≥ 65",
  low_base_mbp          = "Initial MAP < 65",
  acute_coronary_syndrome = "ACS",
  septic_shock          = "Septic shock"
)

load(paste0(output_dir, date, "_gformula_ci_", time_window_width, "hr_", cov_order_label, "_", sg, ".RData"))

key_combined <- paste0("combined_", sg)
dfc          <- results[[key_combined]]

x_lim    <- c(-xlim, xlim)
x_breaks <- seq(-xlim, xlim, by = x_breaks_by)
dpi_out  <- 600

df_fp <- dfc %>% 
  filter(time_points %in% c(7, 28)) %>% 
  mutate(
    tp = factor(
      time_points,
      levels = c(28, 7),
      labels = c("Day 28", "Day 7")
    )
  )

theme_forest <- 
  theme_classic() +
  theme(
    text            = element_text(family = "Arial"),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 16),
    plot.title      = element_text(size = 18)
  )

p_top <- ggplot(df_fp, aes(y = tp)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ll_risk_diff_9, xmax = ul_risk_diff_9),
                 height = 0.25, linewidth = 0.9) +
  geom_point(aes(x = risk_diff_9_mean), size = 3.2, shape = 16, color = "black") +
  scale_x_continuous(limits = x_lim, breaks = x_breaks, labels = function(x) x * 100) +
  labs(
    title = paste0(sg_to_plot[[sg]], " — 9 g/dL vs 7 g/dL"),
    x = "Risk difference (percentage points)", 
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_forest

p_bot <- ggplot(df_fp, aes(y = tp)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ll_risk_diff_8, xmax = ul_risk_diff_8),
                 height = 0.25, linewidth = 0.9) +
  geom_point(aes(x = risk_diff_8_mean), size = 3.2, shape = 16, color = "black") +
  scale_x_continuous(limits = x_lim, breaks = x_breaks, labels = function(x) x * 100) +
  labs(
    title = paste0(sg_to_plot[[sg]], " — 8 g/dL vs 7 g/dL"),
    x = "Risk difference (percentage points)", 
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_forest

p_forest <- p_top / p_bot
print(p_forest)

f_forest <- file.path(
  output_dir,
  sprintf("%s_g_forest_rd_%shr_%s_%s.png", date, time_window_width, cov_order_label, sg)
)

ggsave(
  filename = f_forest,
  plot     = p_forest,
  width    = 8, 
  height   = 6, 
  dpi      = dpi_out
)
