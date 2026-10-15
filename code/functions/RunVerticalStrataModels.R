# ==============================================================================
# RunVerticalStrataModels.R
# Runs Random Forest models across Shannon VSD and the 5 vertical vegetation
# strata (av1-av5) for both S1+S2 (M6) and S2-only (M7), and computes the
# side-by-side comparison and Delta metrics (S1+S2 minus S2-only).
# ==============================================================================

RunVerticalStrataModels <- function(
  df_input,
  s1s2_predictors,
  s2_predictors,
  subset_name,
  ntree = 1000,
  seed = 123
) {
  shannon_and_strata <- c("shannon", "av1", "av2", "av3", "av4", "av5")
  strata_labels <- c(
    shannon = "Shannon VSD (Overall)",
    av1     = "Stratum 1 (0–0.5 m)",
    av2     = "Stratum 2 (0.5–1 m)",
    av3     = "Stratum 3 (1–2 m)",
    av4     = "Stratum 4 (2–4 m)",
    av5     = "Stratum 5 (>4 m)"
  )

  make_empty_summary <- function(var, model_id, model_desc, n_preds) {
    tibble::tibble(
      model_id = model_id,
      model_description = model_desc,
      subset_name = subset_name,
      response = var,
      n_obs = nrow(df_input),
      n_predictors = n_preds,
      ntree = ntree,
      mtry = NA_real_,
      n_cv_folds = 0L,
      Rsquared_mean = NA_real_,
      Rsquared_sd = NA_real_,
      Rsquared_median = NA_real_,
      Rsquared_q25 = NA_real_,
      Rsquared_q75 = NA_real_,
      Rsquared_iqr = NA_real_,
      Rsquared_min = NA_real_,
      Rsquared_max = NA_real_,
      RMSE_mean = NA_real_,
      RMSE_sd = NA_real_,
      RMSE_median = NA_real_,
      RMSE_q25 = NA_real_,
      RMSE_q75 = NA_real_,
      RMSE_iqr = NA_real_,
      RMSE_min = NA_real_,
      RMSE_max = NA_real_,
      MAE_mean = NA_real_,
      MAE_sd = NA_real_,
      MAE_median = NA_real_,
      MAE_q25 = NA_real_,
      MAE_q75 = NA_real_,
      MAE_iqr = NA_real_,
      MAE_min = NA_real_,
      MAE_max = NA_real_,
      Rep20_R2_cor_mean = NA_real_,
      Rep20_R2_cor_sd = NA_real_,
      Rep20_R2_trad_mean = NA_real_,
      Rep20_R2_trad_sd = NA_real_,
      Rep20_RMSE_mean = NA_real_,
      Rep20_RMSE_sd = NA_real_,
      Rep20_MAE_mean = NA_real_,
      Rep20_MAE_sd = NA_real_,
      Rsquared_mean_sd_str = NA_character_,
      Rsquared_med_iqr_str = NA_character_,
      RMSE_mean_sd_str = NA_character_,
      MAE_mean_sd_str = NA_character_,
      Rep20_R2_cor_str = NA_character_,
      Rep20_R2_trad_str = NA_character_
    )
  }

  s1s2_list <- list()
  s2_list   <- list()

  for (var in shannon_and_strata) {
    n_nonzero <- sum(df_input[[var]] > 0, na.rm = TRUE)

    if (n_nonzero >= 3) {
      s1s2_list[[var]] <- FitRepeatedCVRandomForest(
        data = df_input,
        response_col = var,
        predictor_cols = s1s2_predictors,
        model_id = paste0("M6_S1S2_", var),
        model_description = "S1 + S2 Vertical Strata (M6)",
        subset_name = subset_name,
        ntree = ntree,
        seed = seed
      )

      s2_list[[var]] <- FitRepeatedCVRandomForest(
        data = df_input,
        response_col = var,
        predictor_cols = s2_predictors,
        model_id = paste0("M7_S2_", var),
        model_description = "S2 Only Vertical Strata (M7)",
        subset_name = subset_name,
        ntree = ntree,
        seed = seed
      )
    } else {
      s1s2_list[[var]] <- list(
        model = NULL,
        summary = make_empty_summary(var, paste0("M6_S1S2_", var), "S1 + S2 Vertical Strata (M6)", length(s1s2_predictors)),
        folds_raw = tibble::tibble(),
        importance = tibble::tibble()
      )
      s2_list[[var]] <- list(
        model = NULL,
        summary = make_empty_summary(var, paste0("M7_S2_", var), "S2 Only Vertical Strata (M7)", length(s2_predictors)),
        folds_raw = tibble::tibble(),
        importance = tibble::tibble()
      )
    }
  }

  s1s2_summary <- dplyr::bind_rows(lapply(s1s2_list, `[[`, "summary")) %>%
    dplyr::mutate(sensor_set = "S1+S2 (M6)", stratum_label = dplyr::recode(response, !!!strata_labels))
  s2_summary   <- dplyr::bind_rows(lapply(s2_list, `[[`, "summary")) %>%
    dplyr::mutate(sensor_set = "S2 Only (M7)", stratum_label = dplyr::recode(response, !!!strata_labels))

  long_summary <- dplyr::bind_rows(s1s2_summary, s2_summary) %>%
    dplyr::select(subset_name, sensor_set, response, stratum_label, dplyr::everything())

  # Side-by-side comparison table with Delta R2, Delta RMSE, Delta MAE
  comparison_df <- s1s2_summary %>%
    dplyr::select(
      subset_name,
      response,
      stratum_label,
      n_obs,
      S1S2_mtry = mtry,
      S1S2_Rsquared_mean = Rsquared_mean,
      S1S2_Rsquared_sd = Rsquared_sd,
      S1S2_RMSE_mean = RMSE_mean,
      S1S2_RMSE_sd = RMSE_sd,
      S1S2_MAE_mean = MAE_mean,
      S1S2_MAE_sd = MAE_sd,
      S1S2_Rep20_R2_cor_mean = Rep20_R2_cor_mean,
      S1S2_Rep20_R2_trad_mean = Rep20_R2_trad_mean
    ) %>%
    dplyr::left_join(
      s2_summary %>%
        dplyr::select(
          subset_name,
          response,
          S2_mtry = mtry,
          S2_Rsquared_mean = Rsquared_mean,
          S2_Rsquared_sd = Rsquared_sd,
          S2_RMSE_mean = RMSE_mean,
          S2_RMSE_sd = RMSE_sd,
          S2_MAE_mean = MAE_mean,
          S2_MAE_sd = MAE_sd,
          S2_Rep20_R2_cor_mean = Rep20_R2_cor_mean,
          S2_Rep20_R2_trad_mean = Rep20_R2_trad_mean
        ),
        by = c("subset_name", "response")
    ) %>%
    dplyr::mutate(
      Delta_Rsquared_mean = S1S2_Rsquared_mean - S2_Rsquared_mean,
      Delta_RMSE_mean = S1S2_RMSE_mean - S2_RMSE_mean,
      Delta_MAE_mean = S1S2_MAE_mean - S2_MAE_mean,
      Delta_Rep20_R2_cor = S1S2_Rep20_R2_cor_mean - S2_Rep20_R2_cor_mean,
      Delta_Rep20_R2_trad = S1S2_Rep20_R2_trad_mean - S2_Rep20_R2_trad_mean
    )

  folds_raw <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(s1s2_list, `[[`, "folds_raw")) %>% dplyr::mutate(sensor_set = "S1+S2 (M6)"),
    dplyr::bind_rows(lapply(s2_list, `[[`, "folds_raw")) %>% dplyr::mutate(sensor_set = "S2 Only (M7)")
  )

  importance_all <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(s1s2_list, `[[`, "importance")),
    dplyr::bind_rows(lapply(s2_list, `[[`, "importance"))
  )

  return(list(
    s1s2_models = lapply(s1s2_list, `[[`, "model"),
    s2_models = lapply(s2_list, `[[`, "model"),
    long_summary = long_summary,
    comparison_df = comparison_df,
    folds_raw = folds_raw,
    importance = importance_all
  ))
}
