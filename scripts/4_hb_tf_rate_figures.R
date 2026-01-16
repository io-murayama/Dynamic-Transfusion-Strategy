### visualization (hb, transfusion rate) ###
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

# 1. Load Output
load(paste0(output_dir, date, "_sim_subfigures_", time_window_width, "hr_", cov_order_label, ".RData"))

# 2. Visualization
pal_col   <- c("7 g/dL" = "#A6CEE3",
               "8 g/dL" = "#1F78B4",
               "9 g/dL" = "#08306B")
dpi_out   <- 600
pd_width  <- 0.2
pd        <- position_dodge(width = pd_width)

x_max_full    <- followup_length
x_max_full_pd <- x_max_full + pd_width

p1 <- ggplot(tp_all, aes(x = time_points,
                         y = hb_median,
                         colour = strategy)) +
  geom_line(position = pd) +
  geom_point(size = 1, position = pd) +
  geom_errorbar(aes(ymin = hb_p25, ymax = hb_p75),
                width = 0,
                alpha = 0.9,
                show.legend = FALSE,
                position = pd) +
  scale_x_continuous(
    limits = c(-pd_width, x_max_full_pd),
    breaks = seq(0, floor(x_max_full_pd), by = 7)
  ) +
  scale_y_continuous(
    limits = c(7.5, 11),
    breaks = seq(7.5, 11, by = 0.5)
  ) +
  scale_colour_manual(values = pal_col) +
  labs(
    x = "Days Since First Hb < 9 g/dL in ICU",
    y = "Blood Hemoglobin (g/dL)",
    colour = NULL
  ) +
  theme_classic() +
  theme(
    text            = element_text(family = "Arial"),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 16),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 18),
    legend.position = "right"
  )

print(p1)

f1 <- file.path(
  output_dir,
  paste0(date, "_hb_timepoints_",
         time_window_width, "hr_", cov_order_label, ".png")
)
ggsave(filename = f1,
       plot     = p1,
       width    = 12,
       height   = 5,
       dpi      = dpi_out)


p2 <- ggplot(
  tp_all,
  aes(
    x      = time_points,
    y      = transfusion_rate,
    colour = strategy
  )
) +
  geom_line() +
  geom_point(size = 1) +
  scale_x_continuous(
    limits = c(0, x_max_full),
    breaks = seq(0, floor(x_max_full), 7)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  scale_colour_manual(values = pal_col) +
  labs(
    x      = "Days Since First Hb < 9 g/dL in ICU",
    y      = "Proportion transfused",
    colour = NULL
  ) +
  theme_classic() +
  theme(
    text         = element_text(family = "Arial"),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.position = "right"
  )

print(p2)


f2 <- file.path(
  output_dir,
  paste0(date, "_transfusion_rate_timepoints_",
         time_window_width, "hr_", cov_order_label, ".png")
)

ggsave(
  filename = f2,
  plot     = p2,
  width    = 12,
  height   = 5,
  dpi      = dpi_out
)


p_combined <- p1 / p2 +
  plot_layout(heights = c(1, 1))

print(p_combined)


f_combined <- file.path(
  output_dir,
  paste0(date, "_combined_hb_txrate_timepoints_",
         time_window_width, "hr_", cov_order_label, ".png")
)

ggsave(
  filename = f_combined,
  plot     = p_combined,
  width    = 12,
  height   = 10,
  dpi      = dpi_out
)
