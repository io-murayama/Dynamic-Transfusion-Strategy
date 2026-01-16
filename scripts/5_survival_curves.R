### visualization (survival curves) ###
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

#========================================#
# Configurations                         #
#========================================#
data_dir   <- './data/'
output_dir <- './output/'
date       <- '20251011'
time_window_width <- 6
cov_order_label <- "fwd"
sg <- "all"
ylim <- 0.00
y_breaks_by <- 0.20

sg_to_plot <- c(
  all                   = "All patients",
  age_70_or_older       = "Age ≥ 70",
  age_under_70          = "Age < 70",
  high_base_mbp         = "Initial MAP ≥ 65",
  low_base_mbp          = "Initial MAP < 65",
  acute_coronary_syndrome = "ACS",
  septic_shock          = "Septic shock"
)

# 1. Load Output
load(paste0(output_dir, date, "_gformula_ci_", time_window_width, "hr_", cov_order_label, "_", sg, ".RData"))

# 2. Visualization
key_combined <- paste0("combined_", sg)
dfc          <- results[[key_combined]]

pal_col   <- c("7 g/dL" = "#A6CEE3",
               "8 g/dL" = "#1F78B4",
               "9 g/dL" = "#08306B")
dpi_out   <- 600

p_combined <- ggplot(dfc, aes(x = time_points)) +
  geom_line(aes(y = surv7_mean, colour = "7 g/dL"), linewidth = 0.8) +
  geom_line(aes(y = surv8_mean, colour = "8 g/dL"), linewidth = 0.8) +
  geom_line(aes(y = surv9_mean, colour = "9 g/dL"), linewidth = 0.8) +
  scale_x_continuous(
    limits = c(0, followup_length),
    breaks = seq(0, followup_length, by = 7)
  ) +
  scale_y_continuous(
    limits = c(ylim, 1.00),
    breaks = seq(ylim, 1.00, by = y_breaks_by)
  ) +
  scale_colour_manual(values = pal_col) +
  labs(
    x      = "Days Since First Hb < 9 g/dL in ICU",
    y      = "Survival Proportion",
    colour = "Hb Threshold",
    title  = sg_to_plot[[sg]]
  ) +
  theme_classic() +
  theme(
    text            = element_text(family = "Arial"),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 16),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 18),
    legend.position = "right",
    plot.title      = element_text(size = 18)
  )

print(p_combined)


f_combined <- file.path(
  output_dir,
  sprintf(
    "%s_g_surv_combined_%shr_%s_%s.png",
    date, time_window_width, cov_order_label, sg
  )
)


ggsave(
  filename = f_combined,
  plot     = p_combined,
  width    = 8,
  height   = 6,
  dpi      = dpi_out
)