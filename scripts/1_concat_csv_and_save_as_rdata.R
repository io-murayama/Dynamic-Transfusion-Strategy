### concat csv and save as rdata ###
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# Config
date <- "20251011"
time_window_width <- 6
data_dir <- paste0("./data/", date, "_", time_window_width, "hr/")
export_dir <- "./data/" # RDataの保存先

# ファイルの読み込みと結合
file_list <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
sorted_file_list <- sort(file_list) # 名前順に結合する

df_list <- lapply(
  sorted_file_list,
  # icu_stay_idは桁割れ防止のためcharacterで読み込む
  function(f) read.csv(f, colClasses = c(icu_stay_id = "character"))
)

df <- bind_rows(df_list)

# データ型の変更
sapply(df, class)

df$hospital_id             <- as.character(df$hospital_id)
df$hospital_id             <- factor(df$hospital_id, levels = c("1","2","3","5","6","8","9","11"))
levels(df$hospital_id)

df$female                  <- as.integer(df$female)
df$noradrenaline           <- as.integer(df$noradrenaline)
df$adrenaline              <- as.integer(df$adrenaline)
df$dopamine                <- as.integer(df$dopamine)
df$vasopressin             <- as.integer(df$vasopressin)
df$label_mechanical_ventilation <- as.integer(df$label_mechanical_ventilation)

df$trauma                  <- as.integer(df$trauma)
df$circulatory             <- as.integer(df$circulatory)
df$infection               <- as.integer(df$infection)
df$digestive               <- as.integer(df$digestive)
df$neoplasms               <- as.integer(df$neoplasms)
df$respiratory             <- as.integer(df$respiratory)
df$neurological            <- as.integer(df$neurological)
df$external_causes         <- as.integer(df$external_causes)
df$others                  <- as.integer(df$others)
df$acute_coronary_syndrome <- as.integer(df$acute_coronary_syndrome)
df$septic_shock            <- as.integer(df$septic_shock)

df$label_blood_transfusion <- as.integer(df$label_blood_transfusion)
df$icu_death               <- as.integer(df$icu_death)
df$icu_discharge_alive     <- as.integer(df$icu_discharge_alive)
df$survival_after_icu_discharge <- as.integer(df$survival_after_icu_discharge)

sapply(df, class)

# 本解析に必要な列だけ残す
columns_to_use <- c(
  "icu_stay_id", "time_window_index", "age", "female", "hospital_id",
  "infection", "trauma", "circulatory", "digestive", "neoplasms",
  "neurological", "external_causes", "respiratory", "others",
  "acute_coronary_syndrome", "septic_shock", "bt", "hr", "rr", "mbp",
  "spo2", "hemoglobin", "ph", "lactate", "noradrenaline",
  "adrenaline", "vasopressin", "dopamine", "label_mechanical_ventilation", "label_blood_transfusion",
  "icu_death", "icu_discharge_alive", "survival_after_icu_discharge"
)

df <- df[, columns_to_use]

# survival_after_icu_discharge < 0 の場合、ICU内死亡とみなす
# ICUの情報と退院後の情報のデータ不整合であり、退院後死亡とみなすより
# ICU内死亡とみなす方が適切である
print(paste0("Number of unique icu_stay_id: ", length(unique(df$icu_stay_id))))
print(paste0("Number of unique icu_stay_id with survival_after_icu_discharge < 0: ", 
             length(unique(df$icu_stay_id[df$survival_after_icu_discharge < 0]))))

df <- df %>%
  group_by(icu_stay_id) %>%
  mutate(
    icu_discharge_alive = ifelse(
      survival_after_icu_discharge < 0 & time_window_index == max(time_window_index), 0L, icu_discharge_alive
    ),
    icu_death = ifelse(
      survival_after_icu_discharge < 0 & time_window_index == max(time_window_index), 1L, icu_death
    ),
    survival_after_icu_discharge = ifelse(
      survival_after_icu_discharge < 0, 0L, survival_after_icu_discharge
    )
  ) %>% 
  ungroup()

# RDataで保存
save(df, file = paste0(export_dir,"df_",date,"_",time_window_width,"hr.RData"))

# EDA用にCSVでも保存
write.csv(df, file = paste0(export_dir,"df_",date,"_",time_window_width,"hr.csv"), row.names = FALSE)
