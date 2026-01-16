### parametric gformula simulations ###
if (!require("Amelia")) install.packages("Amelia")
library(Amelia)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("gfoRmula")) install.packages("gfoRmula")
library(gfoRmula)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("dtplyr")) install.packages("dtplyr")
library(dtplyr)
if (!require("future.apply")) install.packages("future.apply")
library(future.apply)
if (!require("parallel")) install.packages("parallel")
library(parallel) 
if (!require("progressr")) install.packages("progressr")
library(progressr)
if (!require("RhpcBLASctl")) install.packages("RhpcBLASctl")
library(RhpcBLASctl)


#========================================#
# Configurations                         #
#========================================#
# 進捗表示の有効化
handlers(handler_txtprogressbar)
handlers(global = TRUE)
options(
  progressr.enable = TRUE,
  progressr.clear = FALSE
)

data_dir   <- './data/'
output_dir <- './output/'
date       <- '20251011'

set.seed(813)
followup_length   <- 28        # followup期間（日）
time_window_width <- 6         # time windowの幅（時間）
thresholds        <- c(7,8,9)  # 動的治療戦略のHb閾値（g/dL）
n_iter            <- 1000      # bootstrapの反復回数
size              <- NULL      # bootstrapサンプルサイズ（NULLなら元の症例数を使う）

# memory管理
total_ram     <- 90 * 1024^3  # 90GB, 環境に合わせて適切に設定する
max_workers   <- detectCores() - 1 
print(sprintf("Detected %d CPU cores, using up to %d workers", detectCores(), max_workers))
ram_use_frac  <- 0.7          # 利用するRAMの割合
expand_factor <- 80           # 処理中に膨らむ分（中間オブジェクトなど）を見込む係数

subgroup_filters  <- rlang::quos(
  all                     = (TRUE),
  age_70_or_older         = (age >= 70),
  age_under_70            = (age < 70),
  high_base_mbp           = (mbp >= 65),
  low_base_mbp            = (mbp < 65),
  acute_coronary_syndrome = (acute_coronary_syndrome == 1),
  septic_shock            = (septic_shock == 1)
)

# flags
args          <- commandArgs(trailingOnly = TRUE)
simulate_main <- !("--no-sim-main" %in% args)     # main（bootstrap）はデフォルトで実行
simulate_sub  <- "--sim-sub" %in% args            # sub（resamplingなし1回simulation）はオプションで実行
use_cov_inv   <- "--cov-inv" %in% args            # sensitivity analysis用
if (use_cov_inv) {
  cov_order_label <- "inv"
} else {
  cov_order_label <- "fwd"
}

sg_idx        <- which(args == "--sg")

if (length(sg_idx) == 0) {
  chosen_sg   <- "all" # 指定がない場合は all（全症例）を選択
} else if (length(sg_idx) > 1) {
  stop("Error: --sg must be specified only once.")
} else if (length(args) < sg_idx + 1) {
  stop("Error: --sg must be followed by a subgroup name.")
} else {
  chosen_sg <- args[sg_idx + 1]
}

if ("--time-window-width" %in% args) {            # sensitivity analysis用(時間幅の変更)
  idx <- which(args == "--time-window-width")
  if (length(idx) > 0 && idx < length(args)) {
    time_window_width <- as.numeric(args[idx + 1])
  } else {
    stop("Error: --time-window-width requires a numeric value")
  }
}

if (24 %% time_window_width != 0) {
  stop(sprintf("Error: time_window_width (%d) must be a divisor of 24", time_window_width))
}

message(sprintf(
  "[CLI] simulate_main=%s, simulate_sub=%s, subgroup=%s, cov=%s, time_window_width=%d",
  simulate_main,
  simulate_sub,
  chosen_sg,
  if (use_cov_inv) "inv" else "fwd",
  time_window_width
))


#========================================#
# 0. parallel planning                   #
#    subgroupごとに適切な並列数を推定    #
#========================================#
estimate_bytes_per_row <- function(df_filtered, sample_n = 50000L) {
  n <- nrow(df_filtered)
  # データから最大 sample_n 行（デフォルト5万行）をランダムにサンプリング
  # (全体が sample_n 行未満なら全行を使う)
  s <- df_filtered[sample.int(n, min(sample_n, n)), ]
  return(as.numeric(object.size(s)) / nrow(s))
}

workers_for_sg <- function(n_rows_sg,
                           bytes_per_row,
                           n_iter,
                           total_ram,
                           max_workers,
                           ram_use_frac,
                           expand_factor) {
  budget         <- total_ram * ram_use_frac                   # 予算メモリ
  bytes_per_iter <- n_rows_sg * bytes_per_row * expand_factor  # 1回のsimulationに必要なメモリ
  w              <- floor(budget / bytes_per_iter)             # 同時に動かせるワーカー数
  return(max(1L, min(max_workers, w, n_iter)))
}


#========================================#
# 1. Bootstrap, Imputation               #
#    icu_stay_id単位でリサンプリング     #
#    単一代入法により欠損を処理          #
#========================================#
bootstrap_by_icu <- function(df) {
  orig_ids <- df %>% distinct(icu_stay_id)
  id_map   <- orig_ids %>%
    slice_sample(n = nrow(orig_ids), replace = TRUE) %>%
    mutate(new_id = row_number())
  df_boot  <- id_map %>%
    left_join(df, by = "icu_stay_id", relationship = "many-to-many") %>% 
    mutate(icu_stay_id = as.integer(new_id)) %>%
    select(-new_id) %>%
    arrange(icu_stay_id, time_window_index)

  return(df_boot)
}

single_imputation <- function(df) {
  # 補完対象列
  target_cols <- c(
    'bt', 'hr', 'rr', 'mbp', 'spo2', 'hemoglobin', 'ph', 'lactate'
  )
  # 以下のbinary variableは定義上欠損し得ない
  # 'noradrenaline', 'adrenaline', 'vasopressin', 'dopamine',
  # 'label_mechanical_ventilation', 'label_blood_transfusion'
  
  other_cols <- setdiff(names(df), target_cols)
  if (any(sapply(df[other_cols], function(x) any(is.na(x))))) {
    stop("Error in single_imputation: target_cols 以外に欠損値が含まれています。")
  }

  sg_idvars <- if (!is.null(chosen_sg) && chosen_sg %in% names(df)) {
    if (chosen_sg == "acute_coronary_syndrome" && "circulatory" %in% names(df)) {
      c(chosen_sg, "circulatory")
    } else {
      chosen_sg
    }
  } else {
    NULL
  }
  
  amel <- amelia(
    x = df,
    m = 1,
    cs = 'icu_stay_id',
    ts = 'time_window_index',
    polytime = 2,
    noms = 'hospital_id',
    idvars = sg_idvars,
    empri = 0.005 * nrow(df),
    p2s = 2,
    parallel = 'no'
  )
  
  df_imputed <- amel$imputations[[1]]
  
  return(df_imputed)
}

#========================================#
# 2. G-formula survival model            #
#    ICU退室後のモデルを作成             #
#========================================#
create_cumulatives <- function(df, df_last) {
  counts <- df %>%
    group_by(icu_stay_id) %>%
    summarise(
      num_blood_transfusion = sum(label_blood_transfusion == 1, na.rm = TRUE),
      pct_noradrenaline     = mean(noradrenaline == 1,          na.rm = TRUE),
      pct_adrenaline        = mean(adrenaline == 1,             na.rm = TRUE),
      pct_dopamine          = mean(dopamine == 1,               na.rm = TRUE),
      pct_vasopressin       = mean(vasopressin == 1,            na.rm = TRUE),
      pct_mv                = mean(label_mechanical_ventilation == 1, na.rm = TRUE),
      # count_mv, count_catecholamineはモデルの変数としては使用せず、secondary endpoint解析のみに使用
      count_mv              = sum(label_mechanical_ventilation == 1,  na.rm = TRUE),
      count_catecholamine   = sum(noradrenaline == 1 | adrenaline == 1 | dopamine == 1 | vasopressin == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_last_with_cum <- df_last %>%
    left_join(counts, by = "icu_stay_id")
  
  return(df_last_with_cum)
}

extract_last_time_window <- function(df, cum = FALSE) {
  df_last <- df %>%
    group_by(icu_stay_id) %>%
    slice_max(time_window_index, n=1) %>%
    ungroup()
  if (cum) df_last <- create_cumulatives(df, df_last)
  
  return(df_last)
}

get_gformula_pooled_model <- function(df, followup_length, time_window_width){
  df_last_with_cum <- extract_last_time_window(df, cum = TRUE)
  
  df_surv_g <- df_last_with_cum %>% 
    # ICUを退室した症例のみを抽出
    # = 28日目までICUに残り続けた人と途中でICU内死亡した人を除外
    filter(icu_discharge_alive == 1) %>% 
    # 退室後モデルの作成に使用するデータの範囲をICU内でHb<9となった時刻から
    # 28日後までに限定するための変数を作成
    mutate(
      icu_length_of_stay = floor((time_window_index + 1) * time_window_width / 24),
      followup_after_icu_discharge = pmin(
        # 退室時~1日目までの死亡はsurvival_after_icu_discharge = 0 と記録されている
        survival_after_icu_discharge + 1, 
        pmax(followup_length - icu_length_of_stay, 0))
    ) %>%
    uncount(weights = followup_after_icu_discharge, .remove = FALSE) %>% 
    group_by(icu_stay_id) %>%
    mutate(
      time = row_number() - 1,
      event = as.integer(time == survival_after_icu_discharge)
    ) %>%
    ungroup()
  
  # 説明変数のリスト
  static_vars       <- c("age", "I(age^2)", "female", "hospital_id",
                         "infection","trauma","circulatory","digestive","neoplasms",
                         "neurological","external_causes","respiratory","others",
                         "acute_coronary_syndrome","septic_shock")
  time_varying_vars <- c("bt",         "I(bt^2)",
                         "hr",         "I(hr^2)",
                         "rr",         "I(rr^2)",
                         "mbp",        "I(mbp^2)",
                         "spo2",       "I(spo2^2)",
                         "ph",         "I(ph^2)",
                         "hemoglobin", "I(hemoglobin^2)",
                         "lactate",    "I(lactate^2)")
  cum_vars          <- c("num_blood_transfusion","I(num_blood_transfusion^2)",
                         "pct_noradrenaline","pct_adrenaline","pct_dopamine",
                         "pct_vasopressin","pct_mv")
  time_vars         <- c("time","I(time^2)",
                         "time_window_index","I(time_window_index^2)","I(time_window_index^3)",
                         "time_window_index:time","time_window_index:I(time^2)","I(time_window_index^2):time")
  interaction_vars  <- c("num_blood_transfusion:time",
                         "num_blood_transfusion:I(time^2)",
                         "num_blood_transfusion:mbp",
                         "num_blood_transfusion:hemoglobin",
                         "num_blood_transfusion:time_window_index",
                         "hemoglobin:time_window_index",
                         "num_blood_transfusion:hemoglobin:time_window_index")
  
  all_vars <- c(static_vars, time_varying_vars, cum_vars, time_vars, interaction_vars)
  
  # サブグループごとに説明変数を調整
  # 退室後予測モデルは病名カテゴリーと別の変数のproduct termを含まないので
  # 病名カテゴリ項だけ処理すれば良い
  if (chosen_sg %in% all_vars) {
    all_vars <- setdiff(all_vars, chosen_sg)
    message(sprintf("[get_gformula_pooled_model] subgroup=%s → 説明変数 '%s' を除外", chosen_sg, chosen_sg))
    # acute_coronary_syndromeならcirculatoryも除外する
    if (chosen_sg == "acute_coronary_syndrome") {
      all_vars <- setdiff(all_vars, "circulatory")
      message("[get_gformula_pooled_model] subgroup=acute_coronary_syndrome → 説明変数 'circulatory' も除外")
    }
  } else {
    message(sprintf("[get_gformula_pooled_model] chosen_sg='%s' は 説明変数 に存在しません。除外せずに実行", chosen_sg))
  }
  
  formula      <- as.formula(paste("event ~", paste(all_vars, collapse = " + ")))
  g_pool_model <- glm(
    formula,
    data = df_surv_g,
    family = binomial()
  )
  
  message("[get_gformula_pooled_model] g-formula pooled logistic modelを構築")
  return(g_pool_model)
}


#========================================#
# 3. g‐formula simulation                #
#    ICU内シミュレーション               #
#========================================#
make_dynamic_interventions <- function(thresholds) {
  dynamic_strategy <- function(newdf, pool, intvar, intvals, time_name, t) {
    threshold      <- intvals
    newdf[hemoglobin >= threshold,   (intvar) := 0]
    newdf[hemoglobin <  threshold,   (intvar) := 1]
  }
  int_list <- list()
  for (th in thresholds) {
    int_list[[paste0("below_", th)]] <- list(c(dynamic_strategy, th))
  }
  message("[make_dynamic_interventions] 動的戦略を定義: ", paste(thresholds, collapse = ", "))
  return(int_list)
}

run_gformula_simulation <- function(df_boot, interventions, 
                                    followup_length, time_window_width, keep_model_fits = FALSE) {
  # gfoRmula packageはdata.tableを採用しているので、data.tableに変換しておく
  dt     <- as.data.table(df_boot)
  rm(df_boot); gc()

  # gfoRmula package 内で実装される id の map 対応表を作成
  # https://github.com/CausalInference/gfoRmula/blob/7b0fa3d68bc6637124ac2ae9ea23dd32b486db1c/R/gformula.R#L2707C1-L2712C57
  # 後でremapした時の妥当性検証のため、ageを保持しておく
  id_map <- dt %>%
    lazy_dt() %>%
    select(icu_stay_id, age) %>%
    group_by(icu_stay_id) %>%
    summarise(age_ref = first(age), .groups = "drop") %>%
    arrange(icu_stay_id) %>%
    mutate(id = row_number()) %>%
    select(id, icu_stay_id, age_ref) %>%
    as.data.table()   
  setkey(id_map, id)
  
  n_unique <- length(unique(dt$icu_stay_id))
  
  # gfoRmulaで必要な列だけ残す（処理の効率化）
  columns_to_use <- c(
    "icu_stay_id", "time_window_index", "icu_death", "icu_discharge_alive",
    "age","female","infection", "trauma", "circulatory", "digestive", "neoplasms",
    "neurological", "external_causes", "respiratory", "others",
    "hospital_id","acute_coronary_syndrome","septic_shock",
    "bt","hr","rr","mbp","spo2","ph","hemoglobin","lactate",
    "noradrenaline","adrenaline","dopamine","vasopressin","label_mechanical_ventilation","label_blood_transfusion"
  )
  
  # サブグループごとに使用する列を調整
  if (chosen_sg %in% columns_to_use) {
    columns_to_use <- setdiff(columns_to_use, chosen_sg)
    message(sprintf("[run_gformula_simulation] subgroup=%s → 変数 '%s' を除外", chosen_sg, chosen_sg))
    # acute_coronary_syndromeならcirculatoryも除外する
    if (chosen_sg == "acute_coronary_syndrome") {
      columns_to_use <- setdiff(columns_to_use, "circulatory")
      message("[run_gformula_simulation] subgroup=acute_coronary_syndrome → 変数 'circulatory' も除外")
    }
  } else {
    message(sprintf("[run_gformula_simulation] chosen_sg='%s' は 使用する列 に存在しません。除外せずに実行", chosen_sg))
  }
  
  dt <- dt[, ..columns_to_use]
  
  # basecovs
  basecovs <- c("age", "female","hospital_id",
                "infection", "trauma", "circulatory", "digestive", "neoplasms",
                "neurological", "external_causes", "respiratory", "others",
                "acute_coronary_syndrome", "septic_shock")
  
  if (chosen_sg %in% basecovs) {
    basecovs <- setdiff(basecovs, chosen_sg)
    message(sprintf("[run_gformula_simulation] subgroup=%s → basecovs から変数 '%s' を除外", chosen_sg, chosen_sg))
    if (chosen_sg == "acute_coronary_syndrome") {
      basecovs <- setdiff(basecovs, "circulatory")
      message("[run_gformula_simulation] subgroup=acute_coronary_syndrome → basecovs から変数 'circulatory' も除外")
    }
  } else {
    message(sprintf("[run_gformula_simulation] chosen_sg='%s' は basecovs に存在しません。除外せずに実行。", chosen_sg))
  }

  covnames_fwd  <- c("bt", "hr", "rr", "mbp", "spo2", "ph", 
                    "hemoglobin", "lactate", "noradrenaline", "adrenaline", 
                    "dopamine", "vasopressin", "label_mechanical_ventilation", "label_blood_transfusion")
  covtypes_fwd  <- c("normal", "normal", "normal","normal", "normal", "normal","normal", "normal",
                     "binary", "binary", "binary", "binary", "binary", "binary")

  covnames <- if (use_cov_inv) {
    rev(covnames_fwd)
  } else {
    covnames_fwd
  }
  covtypes <- if (use_cov_inv) {
    rev(covtypes_fwd)
  } else {
    covtypes_fwd
  }
  
  # model for Y
  y_static_vars       <- c("age", "I(age^2)", "female", "hospital_id",
                          "infection", "trauma", "circulatory", "digestive", "neoplasms",
                          "neurological", "external_causes", "respiratory", "others",
                          "acute_coronary_syndrome", "septic_shock")
  y_time_varying_vars <- c("bt",
                           "hr", "I(hr^2)",
                           "rr", "I(rr^2)",
                           "mbp", "I(mbp^2)",
                           "spo2", "I(spo2^2)",
                           "ph", "I(ph^2)",
                           "hemoglobin", "I(hemoglobin^2)",
                           "lactate", "I(lactate^2)",
                           "noradrenaline", "adrenaline", "dopamine", "vasopressin",
                           "label_mechanical_ventilation", "label_blood_transfusion",                   
                           "lag1_label_blood_transfusion")
  y_time_vars         <- c("time_window_index", "I(time_window_index^2)", "I(time_window_index^3)")
  y_interaction_vars  <- c("infection:time_window_index",
                           "trauma:time_window_index",
                           "circulatory:time_window_index",
                           "digestive:time_window_index",
                           "neoplasms:time_window_index",
                           "neurological:time_window_index",
                           "external_causes:time_window_index",
                           "respiratory:time_window_index",
                           "others:time_window_index",
                           "label_blood_transfusion:time_window_index",
                           "label_blood_transfusion:I(time_window_index^2)",
                           "label_blood_transfusion:I(time_window_index^3)")

  if (chosen_sg %in% y_static_vars) {
    y_static_vars <- setdiff(y_static_vars, chosen_sg)
    # 対応する product term も除外
    y_interaction_vars <- y_interaction_vars[!grepl(paste0("^", chosen_sg, ":"), y_interaction_vars)]
    message(sprintf("[run_gformula_simulation(y_model)] subgroup=%s → 変数 '%s' を除外", chosen_sg, chosen_sg))
    if (chosen_sg == "acute_coronary_syndrome") {
      y_static_vars <- setdiff(y_static_vars, "circulatory")
      y_interaction_vars <- y_interaction_vars[!grepl("^circulatory:", y_interaction_vars)]
      message("[run_gformula_simulation(y_model)] subgroup=acute_coronary_syndrome → 変数 'circulatory' も除外")
    }
  } else {
    message(sprintf("[run_gformula_simulation(y_model)] chosen_sg='%s' は 説明変数 に存在しません。除外せずに実行", chosen_sg))
  }
    
  y_all_vars <- c(y_static_vars, y_time_varying_vars, y_time_vars, y_interaction_vars)
  y_formula  <- as.formula(paste("icu_death ~", paste(y_all_vars, collapse = " + ")))
  
  # model for L, A
  l_static_vars       <- c("age", "I(age^2)", "female", "hospital_id",
                           "infection", "trauma", "circulatory", "digestive", "neoplasms",
                           "neurological", "external_causes", "respiratory", "others",
                           "acute_coronary_syndrome", "septic_shock")
  l_lagged_vars       <- c("lag1_bt", 
                           "lag1_hr", "I(lag1_hr^2)", 
                           "lag1_rr", "I(lag1_rr^2)", 
                           "lag1_mbp","I(lag1_mbp^2)",
                           "lag1_spo2","I(lag1_spo2^2)",
                           "lag1_ph","I(lag1_ph^2)",
                           "lag1_hemoglobin","I(lag1_hemoglobin^2)",
                           "lag1_lactate","I(lag1_lactate^2)",
                           "lag1_noradrenaline","lag1_adrenaline","lag1_dopamine","lag1_vasopressin",
                           "lag1_label_mechanical_ventilation","lag1_label_blood_transfusion")
  l_time_vars         <- c("time_window_index", "I(time_window_index^2)")
  
  a_static_vars       <- l_static_vars
  a_time_varying_vars <- c("bt",
                           "hr", "I(hr^2)",
                           "rr", "I(rr^2)",
                           "mbp", "I(mbp^2)",
                           "spo2", "I(spo2^2)",
                           "ph", "I(ph^2)",
                           "hemoglobin", "I(hemoglobin^2)",
                           "lactate", "I(lactate^2)",
                           "noradrenaline", "adrenaline", "dopamine", "vasopressin",
                           "label_mechanical_ventilation")
  a_hist_vars         <- c("lag1_label_blood_transfusion", "lag2_label_blood_transfusion")
  a_time_vars         <- l_time_vars

  if (chosen_sg %in% l_static_vars) {
    l_static_vars <- setdiff(l_static_vars, chosen_sg)
    a_static_vars <- setdiff(a_static_vars, chosen_sg)
    message(sprintf("[run_gformula_simulation(l_model, a_model)] subgroup=%s → 変数 '%s' を除外", chosen_sg, chosen_sg))
    if (chosen_sg == "acute_coronary_syndrome") {
      l_static_vars <- setdiff(l_static_vars, "circulatory")
      a_static_vars <- setdiff(a_static_vars, "circulatory")
      message("[run_gformula_simulation(l_model, a_model)] subgroup=acute_coronary_syndrome → 変数 'circulatory' も除外")
    }
  } else {
    message(sprintf("[run_gformula_simulation(l_model, a_model)] chosen_sg='%s' は 説明変数 に存在しません。除外せずに実行", chosen_sg))
  }
  
  l_common_vars <- c(l_static_vars, l_lagged_vars, l_time_vars)
  a_all_vars    <- c(a_static_vars, a_time_varying_vars, a_hist_vars, a_time_vars)

  cov_formula_vars_fwd <- list(
    bt                           =   l_common_vars,
    hr                           = c(l_common_vars, "bt"),
    rr                           = c(l_common_vars, "bt", "hr"),
    mbp                          = c(l_common_vars, "bt", "hr", "rr"),
    spo2                         = c(l_common_vars, "bt", "hr", "rr", "mbp"),
    ph                           = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2"),
    hemoglobin                   = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph"),
    lactate                      = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph", "hemoglobin"),
    noradrenaline                = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph", "hemoglobin", "lactate"),
    adrenaline                   = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph", "hemoglobin", "lactate", "noradrenaline"),
    dopamine                     = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph", "hemoglobin", "lactate", "noradrenaline", "adrenaline"),
    vasopressin                  = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph", "hemoglobin", "lactate", "noradrenaline", "adrenaline", "dopamine"),
    label_mechanical_ventilation = c(l_common_vars, "bt", "hr", "rr", "mbp", "spo2", "ph", "hemoglobin", "lactate",
                                     "noradrenaline", "adrenaline", "dopamine", "vasopressin"),
    label_blood_transfusion      =   a_all_vars
  )
  
  # sensitivity analysis用
  cov_formula_vars_inv <- list(
    label_blood_transfusion      = a_all_vars,
    label_mechanical_ventilation = l_common_vars,
    vasopressin                  = c(l_common_vars, "label_mechanical_ventilation"),
    dopamine                     = c(l_common_vars, "label_mechanical_ventilation", "vasopressin"),
    adrenaline                   = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine"),
    noradrenaline                = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline"),
    lactate                      = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline"),
    hemoglobin                   = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate"),
    ph                           = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate", "hemoglobin"),
    spo2                         = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate", "hemoglobin", "ph"),
    mbp                          = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate", "hemoglobin", "ph", "spo2"),
    rr                           = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate", "hemoglobin", "ph", "spo2", "mbp"),
    hr                           = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate", "hemoglobin", "ph", "spo2", "mbp", "rr"),
    bt                           = c(l_common_vars, "label_mechanical_ventilation", "vasopressin", "dopamine", "adrenaline", "noradrenaline", "lactate", "hemoglobin", "ph", "spo2", "mbp", "rr", "hr")
  )
  
  cov_formula_vars <- if (use_cov_inv) {
    cov_formula_vars_inv
  } else {
    cov_formula_vars_fwd
  }

  cov_formulas <- lapply(names(cov_formula_vars), function(var){
    rhs <- paste(cov_formula_vars[[var]], collapse = " + ")
    as.formula(paste(var, "~", rhs))
  })
  names(cov_formulas) <- names(cov_formula_vars)
  
  
  res <- gformula(
    seed          = 813,
    obs_data      = dt,
    id            = "icu_stay_id",
    time_name     = "time_window_index",
    time_points   = followup_length*(24/time_window_width),
    outcome_name  = "icu_death",
    outcome_type  = "survival",
    ymodel        = y_formula,
    basecovs      = basecovs,
    covnames      = covnames,
    covtypes      = covtypes,  
    histories     = c(lagged),
    histvars      = list(covnames),
    covparams     = list(covmodels = cov_formulas),
    intvars       = list("label_blood_transfusion", "label_blood_transfusion","label_blood_transfusion"),
    interventions = interventions,                    
    int_descript  = c("Initiate blood transfusions below 7 g/dL", "Initiate blood transfusions below 8 g/dL", 
                      "Initiate blood transfusions below 9 g/dL"),
    ref_int       = 1,
    nsamples      = 0,
    nsimul        = n_unique,
    model_fits    = keep_model_fits,
    sim_data_b    = TRUE,
    parallel      = FALSE,
    threads       = 1
  )
  message("[run_gformula_simulation] g-formula 実行完了")
  
  columns_to_keep <- c(columns_to_use, c("Py", "id", "lag1_label_blood_transfusion"))
  columns_to_keep <- setdiff(columns_to_keep, c("icu_stay_id", "icu_death", "icu_discharge_alive"))

  sims_merged <- lapply(res$sim_data[setdiff(names(res$sim_data), "Natural course")], function(x_dt) {
    x_dt <- x_dt[, ..columns_to_keep]
    x_dt <- id_map[x_dt, on = "id"]
    badN <- x_dt[is.na(age_ref) | is.na(age) | age != age_ref, .N]
    if (badN > 0L) {
      stop(sprintf("[run_gformula_simulation] age mismatch detected: %d rows", badN))
    }
    
    x_dt[, c("id", "age_ref") := NULL]  # 不要になった列を削除
    return(x_dt[])
  })
  names(sims_merged) <- setdiff(names(res$sim_data), "Natural course")
  
  return(list(
    sims = sims_merged,
    model_fits = if (keep_model_fits) res else NULL
  ))
}

get_discharge_model <- function(df) {
  df <- df %>%
    group_by(icu_stay_id) %>%
    arrange(icu_stay_id, time_window_index) %>% 
    mutate(lag1_label_blood_transfusion = lag(label_blood_transfusion, n = 1L, default = 0)) %>%
    ungroup()
  
  # model for ICU discharge
  c_static_vars       <- c("age", "I(age^2)", "female", "hospital_id",
                           "infection", "trauma", "circulatory", "digestive", "neoplasms",
                           "neurological", "external_causes", "respiratory", "others",
                           "acute_coronary_syndrome", "septic_shock")
  c_time_varying_vars <- c("bt",         "I(bt^2)",
                           "hr",         "I(hr^2)",
                           "rr",         "I(rr^2)",
                           "mbp",        "I(mbp^2)",
                           "spo2",       "I(spo2^2)",
                           "ph",         "I(ph^2)",
                           "hemoglobin", "I(hemoglobin^2)",
                           "lactate",    "I(lactate^2)",
                           "noradrenaline", "adrenaline", "dopamine", "vasopressin",
                           "label_mechanical_ventilation", "label_blood_transfusion",
                           "lag1_label_blood_transfusion")
  c_time_vars         <- c("time_window_index", "I(time_window_index^2)", "I(time_window_index^3)")
  c_interaction_vars  <- c("label_blood_transfusion:time_window_index",
                           "label_blood_transfusion:I(time_window_index^2)",
                           "label_blood_transfusion:I(time_window_index^3)")
  
  if (chosen_sg %in% c_static_vars) {
    c_static_vars <- setdiff(c_static_vars, chosen_sg)
    message(sprintf("[get_discharge_model] subgroup=%s → 変数 '%s' を除外", chosen_sg, chosen_sg))
    if (chosen_sg == "acute_coronary_syndrome") {
      c_static_vars <- setdiff(c_static_vars, "circulatory")
      message("[get_discharge_model] subgroup=acute_coronary_syndrome → 変数 'circulatory' も除外")
    }
  } else {
    message(sprintf("[get_discharge_model] chosen_sg='%s' は 説明変数 に存在しません。除外せずに実行", chosen_sg))
  }
  
  c_all_vars <- c(c_static_vars, c_time_varying_vars, c_time_vars, c_interaction_vars)
  discharge_formula <- as.formula(paste("icu_discharge_alive ~", paste(c_all_vars, collapse = " + ")))
  
  discharge_model <- glm(
    formula = discharge_formula,
    data = df,
    family = binomial()
  )
  
  return(discharge_model)
}

process_first_simulation <- function(g_results, 
                                     discharge_model, 
                                     return_all_time_points = FALSE) { # sub figure用のsimulationではTRUEにして全time pointを取り出す
  sim_icu_death_and_icu_discharge <- list()
  sim_icu_death_only              <- list() # ICU内simulationだけのsurvival curveを得るための処理
  
  for (nm in names(g_results$sims)) {
    dt_sim <- g_results$sims[[nm]] %>% 
      as_tibble() %>% # data.table → data.frame
      mutate(
        icu_discharge_alive_prob = predict(
          discharge_model,
          newdata = .,
          type    = "response"
        ),
        icu_discharge_alive = rbinom(nrow(.), 1, icu_discharge_alive_prob),
        icu_death           = rbinom(nrow(.), 1, Py)
      ) %>% 
      select(-starts_with("lag")) # simulation結果として帰ってくるこれ以降不要な列を削除
    
    # icu_discharge_alive, icu_death以降の行を削除
    dt_sim_icu_death_and_icu_discharge <- dt_sim %>% 
      group_by(icu_stay_id) %>% 
      arrange(time_window_index, .by_group = TRUE) %>% 
      mutate(event_flag = (icu_discharge_alive == 1 | icu_death == 1)) %>%
      slice(
        # 最初に event_flag が TRUE になる行までを取得
        seq_len(if (any(event_flag)) which(event_flag)[1] else n())
      ) %>%
      mutate(
        # 最初のイベント行で両方 TRUE の場合 icu_discharge_alive を 0 に変える
        icu_discharge_alive = if_else(row_number() == replace_na(which(event_flag)[1], -1) & icu_death == 1, 0, icu_discharge_alive)
      ) %>%
      ungroup() %>% 
      select(-event_flag) 
    
    # icu_death以降の行を削除
    dt_sim_icu_death_only <- dt_sim %>%
      group_by(icu_stay_id) %>% 
      arrange(time_window_index, .by_group = TRUE) %>% 
      mutate(event_flag = (icu_death == 1)) %>%
      slice(
        # 最初に event_flag が TRUE になる行までを取得
        seq_len(if (any(event_flag)) which(event_flag)[1] else n())
      ) %>%
      ungroup() %>% 
      select(-event_flag)
    
    rm(dt_sim); gc()
    
    if (return_all_time_points) {
      # sub figure用のsimulationでは全time pointを取り出す
      
      if (nm == "Initiate blood transfusions below 7 g/dL") {
        sim_icu_death_and_icu_discharge$sim7 <- dt_sim_icu_death_and_icu_discharge
        sim_icu_death_only$sim7              <- dt_sim_icu_death_only
      } else if (nm == "Initiate blood transfusions below 8 g/dL") {
        sim_icu_death_and_icu_discharge$sim8 <- dt_sim_icu_death_and_icu_discharge
        sim_icu_death_only$sim8              <- dt_sim_icu_death_only
      } else if (nm == "Initiate blood transfusions below 9 g/dL") {
        sim_icu_death_and_icu_discharge$sim9 <- dt_sim_icu_death_and_icu_discharge
        sim_icu_death_only$sim9              <- dt_sim_icu_death_only
      }
      next
    }
    
    # secondary endpoint解析に必要な累積変数を残して最後のtime windowだけ取り出す
    dt_sim_icu_death_and_icu_discharge <- extract_last_time_window(dt_sim_icu_death_and_icu_discharge, cum = TRUE)
    dt_sim_icu_death_only              <- extract_last_time_window(dt_sim_icu_death_only, cum = TRUE)
    
    if (nm == "Initiate blood transfusions below 7 g/dL") {
      sim_icu_death_and_icu_discharge$sim7 <- dt_sim_icu_death_and_icu_discharge
      sim_icu_death_only$sim7              <- dt_sim_icu_death_only
    } else if (nm == "Initiate blood transfusions below 8 g/dL") {
      sim_icu_death_and_icu_discharge$sim8 <- dt_sim_icu_death_and_icu_discharge
      sim_icu_death_only$sim8              <- dt_sim_icu_death_only
    } else if (nm == "Initiate blood transfusions below 9 g/dL") {
      sim_icu_death_and_icu_discharge$sim9 <- dt_sim_icu_death_and_icu_discharge
      sim_icu_death_only$sim9              <- dt_sim_icu_death_only
    }
  }
  
  message("[process_first_simulation] ICU内simulationを実行")
  return(list(
    sim_icu_death_and_icu_discharge = sim_icu_death_and_icu_discharge,
    sim_icu_death_only              = sim_icu_death_only
  ))
}


#======================================================#
# 4. g‐formula survival prediction with simulated data #
#    ICU退室後シミュレーション                         #
#======================================================#
build_second_stage_individuals <- function(sims1, g_pool_model, followup_length){
  
  process_second_simulation <- function(sim_data, g_pool_model, followup_length){
    surv_results <- sim_data %>%  
      # ICUを退室した症例のみを抽出
      # = 最後のtime windowでもICUに残り続けた人、ICU内死亡した人を除外
      filter(icu_discharge_alive == 1) %>% 
      # 実際に生存関数の描画に使うのはfollowup_length - icu_length_of_stay分だけだが
      # 多くestimateする分には結果に影響しないので一律followup_length分展開する
      uncount(weights = followup_length, .remove = F) %>% 
      group_by(icu_stay_id) %>%
      mutate(time = row_number() - 1) %>%
      ungroup()
    
    # g_pool_modelの予測値をbernoulli分布に基づいてhazardを0,1に変換（simulation）
    surv_results[["hazard"]] <- predict(g_pool_model, newdata = surv_results, type = "response")
    surv_results[["death_after_discharge"]] <- rbinom(nrow(surv_results), 1, surv_results[["hazard"]])
    
    # 必要最小限の列への圧縮
    out <- surv_results %>%
      mutate(time_updated = time + 1) %>% 
      select(icu_stay_id, time_window_index, time_updated, 
             hazard, death_after_discharge, count_mv, count_catecholamine) %>%
      arrange(icu_stay_id, time_window_index)
    
    return(out)
  }
  
  list(
    sim7 = process_second_simulation(sims1$sim7, g_pool_model, followup_length),
    sim8 = process_second_simulation(sims1$sim8, g_pool_model, followup_length),
    sim9 = process_second_simulation(sims1$sim9, g_pool_model, followup_length)
  )
}


#======================================================#
# 5. summarize simulation results                      #
#    生存曲線/secondary endpointの計算                 #
#======================================================#
# 5-0. 共通関数
calc_survival <- function(df, time_points) {
  sapply(time_points, function(t) {
    1 - mean(df$time_since_baseline <= t, na.rm = TRUE)
  })
}

convert_to_surv_data <- function(sim_results, followup_length, time_window_width) {
  time_points <- seq(0, followup_length, by = time_window_width / 24)
  
  surv7  <- calc_survival(sim_results$sim7, time_points)
  surv8  <- calc_survival(sim_results$sim8, time_points)
  surv9  <- calc_survival(sim_results$sim9, time_points)
  
  result <- tibble(
    time_points = time_points,
    surv7 = surv7,
    surv8 = surv8,
    surv9 = surv9
  )
  return(result)
}

# 5-1. ICU内シミュレーションのみの結果（生存曲線）を得るための処理（combinedには不要）
get_first_stage_surv_data <- function(sim_icu_death_only, time_window_width) {
  out_list <- list()
  
  for (nm in names(sim_icu_death_only)) {
    # ICU内で死亡した症例
    sim_1_death <- sim_icu_death_only[[nm]] %>% 
      filter(icu_death == 1) %>% 
      mutate(time_since_baseline = (time_window_index + 1) * time_window_width / 24) %>% 
      select(icu_stay_id, time_since_baseline)
    
    # ICU内で生存したままだった症例
    sim_1_survive <- sim_icu_death_only[[nm]] %>% 
      filter(icu_death == 0) %>% 
      mutate(time_since_baseline = 100 * 365) %>% # 解析に影響しない十分長い時間（100年）を代入
      select(icu_stay_id, time_since_baseline)
    
    # sim_1_death, sim_1_surviveを縦方向に結合
    sim_all <- bind_rows(sim_1_death, sim_1_survive)
    if (!all(unique(sim_all$icu_stay_id) %in% unique(sim_icu_death_only[[nm]]$icu_stay_id))) {
      stop(sprintf("Error: icu_stay_id mismatch in %s", nm))
    }
    
    if (nm == "sim7") {
      out_list$sim7 <- sim_all
    } else if (nm == "sim8") {
      out_list$sim8 <- sim_all
    } else if (nm == "sim9") {
      out_list$sim9 <- sim_all
    }
  }
  
  return(out_list)
}

# 5-2. ICU退室後シミュレーションのみの結果（生存曲線）を得るための処理
get_second_stage_surv_data <- function(sim_list_2){
  
  calculate_second_stage_surv_func <- function(sim_data, surv_col_name) {
    surv_data <- sim_data %>% 
      arrange(icu_stay_id, time_updated) %>%
      group_by(icu_stay_id) %>% 
      mutate(surv_indv = cumprod(1-hazard)) %>% 
      ungroup() %>% 
      group_by(time_updated) %>%
      summarise("{surv_col_name}" := mean(surv_indv))
    
    return(surv_data)
  }
  
  s7 <- calculate_second_stage_surv_func(sim_list_2$sim7, "surv7")
  s8 <- calculate_second_stage_surv_func(sim_list_2$sim8, "surv8")
  s9 <- calculate_second_stage_surv_func(sim_list_2$sim9, "surv9")
  
  g_surv_data <- s7 %>%
    select(time_updated, surv7) %>%
    inner_join(s8 %>% select(time_updated, surv8), by = "time_updated") %>%
    inner_join(s9 %>% select(time_updated, surv9), by = "time_updated") %>% 
    add_row(time_updated=0, surv7=1, surv8=1, surv9=1, .before = 1)
  
  return(g_surv_data)
}

# 5-3. ICU内・退室後シミュレーションを組み合わせた結果（生存曲線）を得るための処理
combine_simulation_results <- function(sim_list_1, sim_list_2, time_window_width, followup_length) {
  out_list <- list()
  
  for (nm in names(sim_list_2)) {
    # ICU内で死亡した症例
    sim_1_death <- sim_list_1[[nm]] %>% 
      filter(icu_death == 1) %>% 
      mutate(time_since_baseline = (time_window_index + 1) * time_window_width / 24,
             icu_free_days = 0,
             catecholamine_free_days = 0,
             mv_free_days = 0) %>% 
      select(icu_stay_id, time_since_baseline, icu_free_days, catecholamine_free_days, mv_free_days)
    
    # ICU内で生存したまま退室もしなかった症例
    sim_1_survive <- sim_list_1[[nm]] %>% 
      filter(icu_death == 0 & icu_discharge_alive == 0) %>% 
      mutate(time_since_baseline = 100 * 365, # 解析に影響しない十分長い時間（100年）を代入
             icu_free_days = 0,
             catecholamine_free_days = icu_free_days + ((time_window_index + 1) - count_catecholamine) * time_window_width / 24,
             mv_free_days = icu_free_days + ((time_window_index + 1) - count_mv) * time_window_width / 24) %>% 
      select(icu_stay_id, time_since_baseline, icu_free_days, catecholamine_free_days, mv_free_days)
    
    # ICUを退室した症例
    sim_2 <- sim_list_2[[nm]] %>% 
      group_by(icu_stay_id) %>% 
      arrange(time_updated, .by_group = TRUE) %>% 
      slice(
        seq_len(
          if (any(death_after_discharge == 1)) 
            which(death_after_discharge == 1)[1] 
          else n())
      ) %>%
      filter(time_updated == max(time_updated)) %>% 
      ungroup() %>% 
      mutate(
        # death_after_discharge = 1なら ICU在室時間 + 退院後時間, 
        # death_after_discharge = 0なら解析に影響しない十分長い時間（100年）を代入
        time_since_baseline = if_else(
          death_after_discharge == 1, 
          (time_window_index + 1) * time_window_width / 24 + time_updated, 
          100 * 365 # 解析に影響しない十分長い時間（100年）を代入
          ),
        icu_free_days           = if_else(
          time_since_baseline <= followup_length, 
          0,  # followup期間内に死亡した場合は0日
          pmax(followup_length - (time_window_index + 1) * time_window_width / 24, 0)  # followup期間 - ICU滞在期間
          ), 
        catecholamine_free_days = if_else(
          time_since_baseline <= followup_length, 
          0,  # followup期間内に死亡した場合は0日
          pmax(icu_free_days + ((time_window_index + 1) - count_catecholamine) * time_window_width / 24, 0)
          ),
        mv_free_days            = if_else(
          time_since_baseline <= followup_length, 
          0,  # followup期間内に死亡した場合は0日
          pmax(icu_free_days + ((time_window_index + 1) - count_mv) * time_window_width / 24, 0)
        )
      ) %>%
      select(icu_stay_id, time_since_baseline, icu_free_days, catecholamine_free_days, mv_free_days)
    
    # sim_1_death, sim_1_survive, sim_2を縦方向に結合
    sim_all <- bind_rows(sim_1_death, sim_1_survive, sim_2)
    if (!all(unique(sim_all$icu_stay_id) %in% unique(sim_list_1[[nm]]$icu_stay_id))) {
      stop(sprintf("Error: icu_stay_id mismatch in %s", nm))
    }
    
    if (nm == "sim7") {
      out_list$sim7 <- sim_all
    } else if (nm == "sim8") {
      out_list$sim8 <- sim_all
    } else if (nm == "sim9") {
      out_list$sim9 <- sim_all
    }
  }
  
  return(out_list)
}

# 5-4. secondary endpointの集計
# sim_results$sim7, sim8, sim9 に対して、icu_free_days, catecholamine_free_days, mv_free_daysの平均値を計算
summarize_secondary_endpoints <- function(sim_results) {
  summary_df <- tibble(
    strategy = c("7 g/dL", "8 g/dL", "9 g/dL"),
    icu_free_days = c(
      mean(sim_results$sim7$icu_free_days, na.rm = TRUE),
      mean(sim_results$sim8$icu_free_days, na.rm = TRUE),
      mean(sim_results$sim9$icu_free_days, na.rm = TRUE)
    ),
    catecholamine_free_days = c(
      mean(sim_results$sim7$catecholamine_free_days, na.rm = TRUE),
      mean(sim_results$sim8$catecholamine_free_days, na.rm = TRUE),
      mean(sim_results$sim9$catecholamine_free_days, na.rm = TRUE)
    ),
    mv_free_days = c(
      mean(sim_results$sim7$mv_free_days, na.rm = TRUE),
      mean(sim_results$sim8$mv_free_days, na.rm = TRUE),
      mean(sim_results$sim9$mv_free_days, na.rm = TRUE)
    )
  )
  return(summary_df)
}


#==================#
# 6. bootstrap     #
#==================#
run_one_iter_for_one_sg <- function(df_filtered, sg_ids, i,
                                    interventions, followup_length, time_window_width){
  # worker内での処理を1スレッドに制限
  try({
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    data.table::setDTthreads(1)
  }, silent = TRUE)

  # 1. bootstrap集団を作成し、single imputationを実行
  df_boot      <- bootstrap_by_icu(df_filtered)
  df_boot      <- single_imputation(df_boot)
  
  # 2. ICU退室後死亡のpooled logistic modelを作成
  g_pool_model <- get_gformula_pooled_model(df_boot, followup_length, time_window_width)
  
  # 3. ICU内シミュレーション
  g_results       <- run_gformula_simulation(df_boot, interventions, 
                                             followup_length, time_window_width)
  discharge_model <- get_discharge_model(df_boot)
  sims1           <- process_first_simulation(g_results, discharge_model)
  
  # 4. ICU退室後シミュレーション
  sims2           <- build_second_stage_individuals(sims1$sim_icu_death_and_icu_discharge, 
                                                    g_pool_model, followup_length)
  
  # 5. 生存曲線/secondary endpointの計算
  # 5-1. ICU内シミュレーションのみの結果（生存曲線）を得るための処理（combinedには不要）
  first_sg   <- get_first_stage_surv_data(sims1$sim_icu_death_only, time_window_width)
  surv_first <- convert_to_surv_data(
    sim_results       = first_sg, 
    followup_length   = followup_length, 
    time_window_width = time_window_width)
  
  # 5-2. ICU退室後シミュレーションのみの結果（生存曲線）を得るための処理（combinedには不要）
  surv_second <- get_second_stage_surv_data(sims2) %>% rename(time_points = time_updated)
  
  # 5-3. ICU内・退室後シミュレーションを組み合わせた結果（生存曲線）を得るための処理
  sim_comb <- combine_simulation_results(
    sim_list_1        = sims1$sim_icu_death_and_icu_discharge,
    sim_list_2        = sims2,
    time_window_width = time_window_width,
    followup_length   = followup_length
  )
  surv_comb    <- convert_to_surv_data(
    sim_results       = sim_comb, 
    followup_length   = followup_length, 
    time_window_width = time_window_width)
  
  # 5-4. secondary endpointの計算
  secondary_endpoints <- summarize_secondary_endpoints(sim_comb)
  
  # 中間オブジェクトの掃除（ピークメモリ削減）
  rm(g_pool_model, g_results, sims1, sims2); gc()
  
  list(
    first = surv_first,
    second = surv_second,
    combined = surv_comb,
    secondary_endpoints = secondary_endpoints
  )
}

calculate_confidence_intervals_survival <- function(surv_df){
  surv_df <- surv_df %>% 
    group_by(time_points) %>%
    summarise(
      surv7_mean = mean(surv7),
      surv7_sd   = sd(surv7),
      surv8_mean = mean(surv8),
      surv8_sd   = sd(surv8),
      surv9_mean = mean(surv9),
      surv9_sd   = sd(surv9),
      risk_diff_8_mean = mean((1-surv8) - (1-surv7)),
      risk_diff_8_sd   = sd((1-surv8) - (1-surv7)),
      risk_diff_9_mean = mean((1-surv9) - (1-surv7)),
      risk_diff_9_sd   = sd((1-surv9) - (1-surv7))
    ) %>% 
    mutate(
      ul7 = surv7_mean + qnorm(0.975) * surv7_sd,
      ll7 = surv7_mean - qnorm(0.975) * surv7_sd,
      ul8 = surv8_mean + qnorm(0.975) * surv8_sd,
      ll8 = surv8_mean - qnorm(0.975) * surv8_sd,
      ul9 = surv9_mean + qnorm(0.975) * surv9_sd,
      ll9 = surv9_mean - qnorm(0.975) * surv9_sd,
      ul_risk_diff_8 = risk_diff_8_mean + qnorm(0.975) * risk_diff_8_sd,
      ll_risk_diff_8 = risk_diff_8_mean - qnorm(0.975) * risk_diff_8_sd,
      ul_risk_diff_9 = risk_diff_9_mean + qnorm(0.975) * risk_diff_9_sd,
      ll_risk_diff_9 = risk_diff_9_mean - qnorm(0.975) * risk_diff_9_sd
    )
  
  return(surv_df)
}

calculate_confidence_intervals_secondary_endpoints <- function(sec_df){
  sec_summary <- sec_df %>%
    group_by(strategy) %>%
    summarise(
      icu_free_days_mean           = mean(icu_free_days),
      icu_free_days_sd             = sd(icu_free_days),
      catecholamine_free_days_mean = mean(catecholamine_free_days),
      catecholamine_free_days_sd   = sd(catecholamine_free_days),
      mv_free_days_mean            = mean(mv_free_days),
      mv_free_days_sd              = sd(mv_free_days)
    ) %>%
    mutate(
      ul_icu_free_days           = icu_free_days_mean + qnorm(0.975) * icu_free_days_sd,
      ll_icu_free_days           = icu_free_days_mean - qnorm(0.975) * icu_free_days_sd,
      ul_catecholamine_free_days = catecholamine_free_days_mean + qnorm(0.975) * catecholamine_free_days_sd,
      ll_catecholamine_free_days = catecholamine_free_days_mean - qnorm(0.975) * catecholamine_free_days_sd,
      ul_mv_free_days            = mv_free_days_mean + qnorm(0.975) * mv_free_days_sd,
      ll_mv_free_days            = mv_free_days_mean - qnorm(0.975) * mv_free_days_sd
    )
  
  return(sec_summary)
}

get_gformula_ci_single_sg <- function(
    df_filtered, interventions,
    followup_length, time_window_width,
    n_iter, size,
    total_ram, max_workers, ram_use_frac, expand_factor
) {
  
  # 指定された症例数にランダムサンプリング
  # 実験用に使うだけで、本解析ではsize = NULLとする
  if (!is.null(size)) {
    sampled_ids <- df_filtered %>%
      distinct(icu_stay_id) %>%
      slice_sample(n = size, replace = FALSE) %>%
      pull(icu_stay_id)
    df_filtered <- df_filtered %>%
      filter(icu_stay_id %in% sampled_ids)
  }
  
  # 最適なworker数を決定
  bytes_per_row <- estimate_bytes_per_row(df_filtered)
  ids           <- df_filtered %>%                                    
    distinct(icu_stay_id) %>%                               
    pull(icu_stay_id)                                       
  n_rows_sg     <- nrow(df_filtered)
  n_workers     <- workers_for_sg(
    n_rows_sg     = n_rows_sg, 
    bytes_per_row = bytes_per_row, 
    n_iter        = n_iter, 
    total_ram     = total_ram, 
    max_workers   = max_workers, 
    ram_use_frac  = ram_use_frac, 
    expand_factor = expand_factor
  )
  future::plan(future::multisession, workers = n_workers) 
  message(sprintf(                                          
    "[scheduler] sg=%s, ids=%d, rows≈%d, data≈%.3f GB -> workers=%d",
    chosen_sg, length(ids), n_rows_sg, (bytes_per_row * n_rows_sg) / 1024^3, n_workers
  ))
  
  # bootstrap処理を並列実行
  pieces <- with_progress({
    p <- progressor(steps = n_iter)
    future_lapply(seq_len(n_iter), function(i) {
      p(sprintf("sg=%s iter=%d/%d", chosen_sg, i, n_iter))
      out <- tryCatch(
        {
          run_one_iter_for_one_sg(
            df_filtered = df_filtered, 
            sg_ids = ids, 
            i = i,
            interventions = interventions, 
            followup_length = followup_length,
            time_window_width = time_window_width
          )
        },
        error = function(e) {
          message(sprintf("[bootstrap] sg=%s iter=%d FAILED: %s",
                          chosen_sg, i, e$message))
          return(NULL)
        }
      )
      return(out)
    }, 
    future.seed = TRUE
    )
  })

  pieces_valid <- Filter(Negate(is.null), pieces)
  n_success    <- length(pieces_valid)
  if (n_success == 0L) {
    stop("[get_gformula_ci_single_sg] All bootstrap iterations failed (n_success = 0).")
  }
  message(sprintf("[get_gformula_ci_single_sg] Successful iterations: %d / %d",
                  n_success, n_iter))
  
  # 結果を集計
  surv_results_1    <- bind_rows(lapply(pieces_valid, `[[`, "first"))
  surv_results_2    <- bind_rows(lapply(pieces_valid, `[[`, "second"))
  surv_results_c    <- bind_rows(lapply(pieces_valid, `[[`, "combined"))
  surv_results_se   <- bind_rows(lapply(pieces_valid, `[[`, "secondary_endpoints"))
  rm(pieces, pieces_valid); gc()
  
  g_surv_data_1     <- calculate_confidence_intervals_survival(surv_results_1)
  g_surv_data_2     <- calculate_confidence_intervals_survival(surv_results_2)
  g_surv_data_c     <- calculate_confidence_intervals_survival(surv_results_c)
  g_surv_data_se    <- calculate_confidence_intervals_secondary_endpoints(surv_results_se)
  
  results           <- list()
  results[[paste0("first_", chosen_sg)]]               <- tibble(g_surv_data_1)
  results[[paste0("second_", chosen_sg)]]              <- tibble(g_surv_data_2)
  results[[paste0("combined_", chosen_sg)]]            <- tibble(g_surv_data_c)
  results[[paste0("secondary_endpoints_", chosen_sg)]] <- tibble(g_surv_data_se)
  results[["n_success"]]                               <- n_success   
  
  # 結果をRDataで保存
  rdata_file <- file.path(output_dir, paste0(date, "_gformula_ci_", time_window_width, "hr_", cov_order_label, "_", chosen_sg, ".RData"))
  save(
    list     = c("results", "thresholds", "followup_length",
               "time_window_width", "n_iter", "n_success", "size", "chosen_sg", "cov_order_label"),
    file     = rdata_file
  )
  message("[get_gformula_ci_single_sg] Results saved for subgroup: ",
          chosen_sg, " -> ", rdata_file)
  
  return(results)
}


#==========================#
# 7. sim for sub figures   #
#==========================#
sim_without_bootstrap <- function(df_filtered, interventions, 
                              followup_length, time_window_width){
  # single imputationを実行
  df_filtered <- single_imputation(df_filtered)
  
  # ICU内シミュレーション
  g_results   <- run_gformula_simulation(df_filtered, interventions,
                                       followup_length, time_window_width)
  discharge_model <- get_discharge_model(df_filtered)
  sims1           <- process_first_simulation(g_results, 
                                              discharge_model, 
                                              return_all_time_points = TRUE)$sim_icu_death_and_icu_discharge
  
  # time window ごとに subplot に必要な情報を集計
  summarise_per_tw <- function(sim_data, strategy_label) {
    sim_data <- sim_data %>%
      # ・Day0に値が表示された方がグラフが綺麗に見える
      # ・Hemoglobinは各time windowにおいてfirst valueを取っている
      # ので、各time windowにおけるHbの値, 輸血割合は左寄せでstart_time上に表示している
      mutate(
        time_points = time_window_index * time_window_width / 24
      ) %>%
      group_by(time_points) %>%
      summarise(
        hb_median          = median(hemoglobin, na.rm = TRUE),
        hb_p25             = quantile(hemoglobin, 0.25, na.rm = TRUE),
        hb_p75             = quantile(hemoglobin, 0.75, na.rm = TRUE),
        transfusion_rate   = mean(label_blood_transfusion == 1, na.rm = TRUE),
        .groups            = "drop"
      ) %>%
      mutate(strategy = strategy_label)
  }
  
  tp7    <- summarise_per_tw(sims1$sim7, "7 g/dL")
  tp8    <- summarise_per_tw(sims1$sim8, "8 g/dL")
  tp9    <- summarise_per_tw(sims1$sim9, "9 g/dL")
  tp_all <- bind_rows(tp7, tp8, tp9)
  

  rdata_file <- file.path(output_dir, paste0(date, "_sim_subfigures_", time_window_width, "hr_", cov_order_label, ".RData"))
  save(list = c("tp_all", "followup_length", "time_window_width", "interventions", "cov_order_label"),
       file = rdata_file)
  message("[sim_without_bootstrap] Results saved: ", rdata_file)
  
  return(tp_all)
}


#==================#
# 8. main flow     #
#==================#
# データの読み込み
load(paste0(data_dir, "df_", date, "_", time_window_width, "hr.RData"))
if (any(is.na(df$survival_after_icu_discharge))) {
  stop("Error: survival_after_icu_discharge contains NA values.")
}
sg_ids <- df %>%
  filter(time_window_index == 0) %>%
  filter(!!subgroup_filters[[chosen_sg]]) %>%
  distinct(icu_stay_id) %>%
  pull(icu_stay_id)
df_filtered   <- df %>% 
  filter(icu_stay_id %in% sg_ids) %>%
  filter(time_window_index <= followup_length*(24/time_window_width) - 1) # フォローアップ期間に合わせてデータを抽出
rm(df); gc()

# 動的介入の定義
interventions <- make_dynamic_interventions(thresholds)

# main simulation 
if (simulate_main) {
  get_gformula_ci_single_sg(
    df_filtered       = df_filtered, 
    interventions     = interventions, 
    followup_length   = followup_length, 
    time_window_width = time_window_width, 
    n_iter            = n_iter, 
    size              = size,
    total_ram         = total_ram,
    max_workers       = max_workers,
    ram_use_frac      = ram_use_frac,
    expand_factor     = expand_factor
  )
} else {
  message("[SKIP] main simulation")
}

# simulation for sub figures
if (simulate_sub) {
  sim_without_bootstrap(
    df_filtered, 
    interventions, 
    followup_length, 
    time_window_width
  )
} else {
  message("[SKIP] sub simulation")
}


#=======================================#
# 9. 検証用（ブートストラップなし）     #
#=======================================#
# # データの読み込み
# load(paste0(data_dir, "df_", date, "_", time_window_width, "hr.RData"))
# if (any(is.na(df$survival_after_icu_discharge))) {
#   stop("Error: survival_after_icu_discharge contains NA values.")
# }
# df_filtered   <- df %>%
#   filter(!!subgroup_filters[[chosen_sg]]) %>%
#   filter(time_window_index <= followup_length*(24/time_window_width) - 1) # フォローアップ期間に合わせてデータを抽出
# rm(df); gc()
# 
# # 動的介入の定義
# interventions <- make_dynamic_interventions(thresholds)
# 
# if (!is.null(size)) {
#   sampled_ids <- df_filtered %>%
#     distinct(icu_stay_id) %>%
#     slice_sample(n = size, replace = FALSE) %>%
#     pull(icu_stay_id)
#   df_filtered <- df_filtered %>%
#     filter(icu_stay_id %in% sampled_ids)
# }
# 
# # 1. bootstrap集団を作成し、single imputationを実行
# df_boot      <- bootstrap_by_icu(df_filtered)
# df_boot      <- single_imputation(df_boot)
# 
# # 2. ICU退室後死亡のpooled logistic modelを作成
# g_pool_model <- get_gformula_pooled_model(df_boot, followup_length, time_window_width)
# 
# # 3. ICU内シミュレーション
# g_results       <- run_gformula_simulation(df_boot, interventions,
#                                            followup_length, time_window_width, keep_model_fits = TRUE)
# discharge_model <- get_discharge_model(df_boot)
# sims1           <- process_first_simulation(g_results, discharge_model)
# 
# # 4. ICU退室後シミュレーション
# sims2           <- build_second_stage_individuals(sims1$sim_icu_death_and_icu_discharge,
#                                                   g_pool_model, followup_length)
# 
# # 5. 生存曲線/secondary endpointの計算
# # 5-1. ICU内シミュレーションのみの結果（生存曲線）を得るための処理（combinedには不要）
# first_sg   <- get_first_stage_surv_data(sims1$sim_icu_death_only, time_window_width)
# surv_first <- convert_to_surv_data(
#   sim_results       = first_sg,
#   followup_length   = followup_length,
#   time_window_width = time_window_width)
# 
# # 5-2. ICU退室後シミュレーションのみの結果（生存曲線）を得るための処理（combinedには不要）
# surv_second <- get_second_stage_surv_data(sims2) %>% rename(time_points = time_updated)
# 
# # 5-3. ICU内・退室後シミュレーションを組み合わせた結果（生存曲線）を得るための処理
# sim_comb <- combine_simulation_results(
#   sim_list_1        = sims1$sim_icu_death_and_icu_discharge,
#   sim_list_2        = sims2,
#   time_window_width = time_window_width,
#   followup_length   = followup_length
# )
# surv_comb    <- convert_to_surv_data(
#   sim_results       = sim_comb,
#   followup_length   = followup_length,
#   time_window_width = time_window_width)
# 
# # 5-4. secondary endpointの計算
# secondary_endpoints <- summarize_secondary_endpoints(sim_comb)
