### visualization (forest plots) ###
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
library(tidyverse)
if (!requireNamespace("forestplot", quietly = TRUE)) install.packages("forestplot")
library(forestplot)
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
library(grid)

data_dir   <- "./data/"
output_dir <- "./output/"
date       <- "20251011"
time_window_width <- 6
cov_order_label   <- "fwd"

sg_to_plot <- c(
  all             = "All patients",
  age_70_or_older = "Age ≥ 70",
  age_under_70    = "Age < 70",
  high_base_mbp   = "Initial MAP ≥ 65",
  low_base_mbp    = "Initial MAP < 65",
  acute_coronary_syndrome = "ACS",
  septic_shock    = "Septic shock"
)
days_to_keep <- c(7, 28)
sg_keys      <- names(sg_to_plot)

rows_spec <- tibble(
  sg_key = c(
    "all",
    "age_header",
    "age_70_or_older",
    "age_under_70",
    "map_header",
    "high_base_mbp",
    "low_base_mbp",
    "acute_coronary_syndrome",
    "septic_shock"
  ),
  label = c(
    "All patients",
    "Age",
    "  ≥70 years",
    "  <70 years",
    "MAP at time zero",
    "  ≥65 mmHg",
    "  <65 mmHg",
    "ACS",
    "Septic shock"
  ),
  row_type = c(
    "data",
    "header",
    "data",
    "data",
    "header",
    "data",
    "data",
    "data",
    "data"
  ),
  row_order = seq_len(9)
)

extract_one_sg <- function(sg, output_dir, date, time_window_width, cov_order_label) {
  load_path <- paste0(
    output_dir, date, "_gformula_ci_",
    time_window_width, "hr_", cov_order_label, "_", sg, ".RData"
  )
  load(load_path)
  combined_name <- paste0("combined_", sg)
  combined_df <- results[[combined_name]]
  
  combined_df %>%
    mutate(time_points = as.numeric(time_points)) %>%
    filter(round(time_points, 2) %in% days_to_keep) %>%
    transmute(
      sg_key = sg,
      day = as.integer(round(time_points, 0)),
      mort_7 = 1 - surv7_mean,
      mort_8 = 1 - surv8_mean,
      mort_9 = 1 - surv9_mean,
      rd_8vs7 = risk_diff_8_mean,
      rd_lcl_8vs7 = ll_risk_diff_8,
      rd_ucl_8vs7 = ul_risk_diff_8,
      rd_9vs7 = risk_diff_9_mean,
      rd_lcl_9vs7 = ll_risk_diff_9,
      rd_ucl_9vs7 = ul_risk_diff_9
    )
}

sg_data_long <- map_dfr(
  sg_keys,
  ~ extract_one_sg(.x, output_dir, date, time_window_width, cov_order_label)
)

load_path_df <- paste0(data_dir, "df_", date, "_", time_window_width, "hr.RData")
load(load_path_df)

df_base <- df %>%
  filter(time_window_index == 0)

n_table <- tibble(
  sg_key = c(
    "all",
    "age_70_or_older",
    "age_under_70",
    "high_base_mbp",
    "low_base_mbp",
    "acute_coronary_syndrome",
    "septic_shock"
  ),
  n_patients = c(
    nrow(df_base),
    nrow(filter(df_base, age >= 70)),
    nrow(filter(df_base, age < 70)),
    nrow(filter(df_base, mbp >= 65)),
    nrow(filter(df_base, mbp < 65)),
    nrow(filter(df_base, acute_coronary_syndrome == 1)),
    nrow(filter(df_base, septic_shock == 1))
  )
)

forestplot_input <- rows_spec %>%
  crossing(day = days_to_keep) %>%
  left_join(n_table, by = "sg_key") %>%                  # ← これを追加
  left_join(sg_data_long, by = c("sg_key", "day")) %>%
  mutate(
    mort_7 = ifelse(row_type == "data", mort_7, NA_real_),
    mort_8 = ifelse(row_type == "data", mort_8, NA_real_),
    mort_9 = ifelse(row_type == "data", mort_9, NA_real_),
    rd_8vs7 = ifelse(row_type == "data", rd_8vs7, NA_real_),
    rd_lcl_8vs7 = ifelse(row_type == "data", rd_lcl_8vs7, NA_real_),
    rd_ucl_8vs7 = ifelse(row_type == "data", rd_ucl_8vs7, NA_real_),
    rd_9vs7 = ifelse(row_type == "data", rd_9vs7, NA_real_),
    rd_lcl_9vs7 = ifelse(row_type == "data", rd_lcl_9vs7, NA_real_),
    rd_ucl_9vs7 = ifelse(row_type == "data", rd_ucl_9vs7, NA_real_)
  ) %>%
  arrange(day, row_order) %>%
  select(
    day, row_order, label, row_type, n_patients,
    mort_7, mort_8, mort_9,
    rd_9vs7, rd_lcl_9vs7, rd_ucl_9vs7,
    rd_8vs7, rd_lcl_8vs7, rd_ucl_8vs7
  )

cfg <- expand_grid(
  day = c(7, 28),
  strategy = c("9vs7", "8vs7")
)

out_paths <- character(nrow(cfg))

for (i in seq_len(nrow(cfg))) {
  
  day <- cfg$day[[i]]
  strategy <- cfg$strategy[[i]]
  
  df <- forestplot_input %>%
    filter(day == !!day) %>%
    arrange(row_order)

  xlim   <- c(-8, 8)
  xticks <- seq(-8, 8, by = 2)
  
  if (strategy == "9vs7") {
    col_a <- "mort_9"
    col_b <- "mort_7"
    label_a <- "9 g/dL (%)"
    label_b <- "7 g/dL (%)"
    title_strategy <- "9 vs 7 g/dL"
    rd_col <- "rd_9vs7"
    lcl_col <- "rd_lcl_9vs7"
    ucl_col <- "rd_ucl_9vs7"
  } else {
    col_a <- "mort_8"
    col_b <- "mort_7"
    label_a <- "8 g/dL (%)"
    label_b <- "7 g/dL (%)"
    title_strategy <- "8 vs 7 g/dL"
    rd_col <- "rd_8vs7"
    lcl_col <- "rd_lcl_8vs7"
    ucl_col <- "rd_ucl_8vs7"
  }
  
  df2 <- df %>%
    mutate(
      n_txt      = ifelse(is.na(n_patients), "", as.character(n_patients)),
      mort_a_txt = ifelse(is.na(.data[[col_a]]), "", sprintf("%.1f", 100 * .data[[col_a]])),
      mort_b_txt = ifelse(is.na(.data[[col_b]]), "", sprintf("%.1f", 100 * .data[[col_b]])),
      rd_ci_txt  = ifelse(
        is.na(.data[[rd_col]]),
        "",
        sprintf(
          "%.1f (%.1f to %.1f)",
          100 * .data[[rd_col]],
          100 * .data[[lcl_col]],
          100 * .data[[ucl_col]]
        )
      )
    )
  
  mort_header <- paste0(day, "-day mortality")
  
  tabletext_header <- rbind(
    c("", "", mort_header, mort_header, "Risk difference"),
    c("Subgroup","No. of patients", label_a, label_b, "(pp; 95% CI)")
  )
  
  tabletext_body <- df2 %>%
    select(label, n_txt, mort_a_txt, mort_b_txt, rd_ci_txt) %>%
    as.matrix()
  
  tabletext <- rbind(tabletext_header, tabletext_body)
  
  mean_vec_data  <- ifelse(df2$row_type == "data", 100 * df2[[rd_col]],  NA_real_)
  lower_vec_data <- ifelse(df2$row_type == "data", 100 * df2[[lcl_col]], NA_real_)
  upper_vec_data <- ifelse(df2$row_type == "data", 100 * df2[[ucl_col]], NA_real_)
  
  mean_vec_data_plot  <- pmin(pmax(mean_vec_data,  xlim[1], na.rm = FALSE), xlim[2], na.rm = FALSE)
  lower_vec_data_plot <- pmax(lower_vec_data, xlim[1], na.rm = FALSE)
  upper_vec_data_plot <- pmin(upper_vec_data, xlim[2], na.rm = FALSE)
  
  mean_vec  <- c(NA_real_, NA_real_, mean_vec_data_plot)
  lower_vec <- c(NA_real_, NA_real_, lower_vec_data_plot)
  upper_vec <- c(NA_real_, NA_real_, upper_vec_data_plot)
  
  is_summary_vec <- c(TRUE, TRUE, rep(FALSE, nrow(df2)))
  
  out_path <- file.path(output_dir, paste0(date, "_forestplot_day", day, "_", strategy, ".png"))
  out_paths[[i]] <- out_path
  
  png(out_path, width = 2300, height = 820, res = 200)
  
  fp <- forestplot(
    labeltext  = tabletext,
    mean       = mean_vec,
    lower      = lower_vec,
    upper      = upper_vec,
    is.summary = is_summary_vec,
    zero       = 0,
    xlog       = FALSE,
    xticks     = xticks,
    xlim       = xlim,
    graphwidth = unit(95, "mm"),
    colgap     = unit(8, "mm"),
    align      = c("l", "r", "r", "r", "r"),
    boxsize    = 0.16,
    line.margin = 0.02,
    xlab       = "",
    title = paste0(
      "Subgroup analysis of ",
      day,
      "-day mortality: risk difference (",
      title_strategy,
      ")"
    ),
    col        = fpColors(box = "black", line = "black", summary = "black"),
    fn.ci_norm = fpDrawNormalCI,
    hrzl_lines = list("3" = gpar(lwd = 1)),
    txt_gp     = fpTxtGp(
         title = gpar(fontface = "bold", cex = 1.10),
         label = gpar(cex = 0.86),
         ticks = gpar(cex = 0.88)
    )
  )
  fp <- fp_set_zebra_style(fp, "#F5F5F5", ignore_subheaders = TRUE)
  
  print(fp)
  dev.off()
  
  message("Saved: ", out_path)
}

print(out_paths)