# ==============================================================================
# FitRepeatedCVRandomForest.R
# Trains and evaluates Random Forest regression models using 5-fold x 20-repeat
# cross-validation (100 folds), returning both 100-fold distribution statistics
# and 20-repeat pooled out-of-fold metrics.
# ==============================================================================

FitRepeatedCVRandomForest <- function(
  data,
  response_col = "shannon",
  predictor_cols,
  model_id,
  model_description,
  subset_name = "All_Plots",
  ntree = 1000,
  n_folds = 5,
  n_repeats = 20,
  seed = 123
) {
  mod_data <- data %>%
    dplyr::select(dplyr::all_of(c(response_col, predictor_cols))) %>%
    stats::na.omit()

  form <- stats::as.formula(paste(response_col, "~ ."))

  set.seed(seed)
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = n_folds,
    repeats = n_repeats,
    savePredictions = "final"
  )

  set.seed(seed)
  rf_fit <- caret::train(
    form,
    data = mod_data,
    method = "rf",
    trControl = ctrl,
    ntree = ntree,
    importance = TRUE
  )

  best_mtry <- rf_fit$bestTune$mtry
  best_row  <- rf_fit$results[rf_fit$results$mtry == best_mtry, ]
  res       <- rf_fit$resample %>%
    dplyr::mutate(
      model_id = model_id,
      model_description = model_description,
      subset_name = subset_name,
      response = response_col,
      n_obs = nrow(mod_data),
      mtry = best_mtry,
      ntree = ntree
    )

  # Compute 20-repeat pooled out-of-fold metrics (pools all n_obs predictions per repeat)
  preds <- rf_fit$pred
  preds$Rep <- sub("^Fold\\d+\\.", "", preds$Resample)
  rep_mat <- do.call(rbind, lapply(split(preds, preds$Rep), function(df) {
    c(
      R2_cor  = stats::cor(df$obs, df$pred)^2,
      R2_trad = 1 - sum((df$obs - df$pred)^2) / sum((df$obs - mean(df$obs))^2),
      RMSE    = sqrt(mean((df$obs - df$pred)^2)),
      MAE     = mean(abs(df$obs - df$pred))
    )
  }))

  summary_row <- tibble::tibble(
    model_id = model_id,
    model_description = model_description,
    subset_name = subset_name,
    response = response_col,
    n_obs = nrow(mod_data),
    n_predictors = length(predictor_cols),
    ntree = ntree,
    mtry = best_mtry,
    n_cv_folds = nrow(res),
    # 100-fold CV metrics
    Rsquared_mean = best_row$Rsquared,
    Rsquared_sd = best_row$RsquaredSD,
    Rsquared_median = stats::median(res$Rsquared, na.rm = TRUE),
    Rsquared_q25 = unname(stats::quantile(res$Rsquared, 0.25, na.rm = TRUE)),
    Rsquared_q75 = unname(stats::quantile(res$Rsquared, 0.75, na.rm = TRUE)),
    Rsquared_iqr = stats::IQR(res$Rsquared, na.rm = TRUE),
    Rsquared_min = min(res$Rsquared, na.rm = TRUE),
    Rsquared_max = max(res$Rsquared, na.rm = TRUE),
    RMSE_mean = best_row$RMSE,
    RMSE_sd = best_row$RMSESD,
    RMSE_median = stats::median(res$RMSE, na.rm = TRUE),
    RMSE_q25 = unname(stats::quantile(res$RMSE, 0.25, na.rm = TRUE)),
    RMSE_q75 = unname(stats::quantile(res$RMSE, 0.75, na.rm = TRUE)),
    RMSE_iqr = stats::IQR(res$RMSE, na.rm = TRUE),
    RMSE_min = min(res$RMSE, na.rm = TRUE),
    RMSE_max = max(res$RMSE, na.rm = TRUE),
    MAE_mean = best_row$MAE,
    MAE_sd = best_row$MAESD,
    MAE_median = stats::median(res$MAE, na.rm = TRUE),
    MAE_q25 = unname(stats::quantile(res$MAE, 0.25, na.rm = TRUE)),
    MAE_q75 = unname(stats::quantile(res$MAE, 0.75, na.rm = TRUE)),
    MAE_iqr = stats::IQR(res$MAE, na.rm = TRUE),
    MAE_min = min(res$MAE, na.rm = TRUE),
    MAE_max = max(res$MAE, na.rm = TRUE),
    # 20-repeat pooled out-of-fold metrics
    Rep20_R2_cor_mean = mean(rep_mat[, "R2_cor"], na.rm = TRUE),
    Rep20_R2_cor_sd = stats::sd(rep_mat[, "R2_cor"], na.rm = TRUE),
    Rep20_R2_trad_mean = mean(rep_mat[, "R2_trad"], na.rm = TRUE),
    Rep20_R2_trad_sd = stats::sd(rep_mat[, "R2_trad"], na.rm = TRUE),
    Rep20_RMSE_mean = mean(rep_mat[, "RMSE"], na.rm = TRUE),
    Rep20_RMSE_sd = stats::sd(rep_mat[, "RMSE"], na.rm = TRUE),
    Rep20_MAE_mean = mean(rep_mat[, "MAE"], na.rm = TRUE),
    Rep20_MAE_sd = stats::sd(rep_mat[, "MAE"], na.rm = TRUE),
    # Formatted strings for manuscript tables
    Rsquared_mean_sd_str = sprintf("%.3f ± %.3f", best_row$Rsquared, best_row$RsquaredSD),
    Rsquared_med_iqr_str = sprintf(
      "%.3f (%.3f–%.3f)",
      stats::median(res$Rsquared, na.rm = TRUE),
      stats::quantile(res$Rsquared, 0.25, na.rm = TRUE),
      stats::quantile(res$Rsquared, 0.75, na.rm = TRUE)
    ),
    RMSE_mean_sd_str = sprintf("%.2f ± %.2f", best_row$RMSE, best_row$RMSESD),
    MAE_mean_sd_str = sprintf("%.2f ± %.2f", best_row$MAE, best_row$MAESD),
    Rep20_R2_cor_str = sprintf(
      "%.3f ± %.3f",
      mean(rep_mat[, "R2_cor"], na.rm = TRUE),
      stats::sd(rep_mat[, "R2_cor"], na.rm = TRUE)
    ),
    Rep20_R2_trad_str = sprintf(
      "%.3f ± %.3f",
      mean(rep_mat[, "R2_trad"], na.rm = TRUE),
      stats::sd(rep_mat[, "R2_trad"], na.rm = TRUE)
    )
  )

  # Extract Variable Importance (both unscaled %IncMSE and scaled 0-100)
  var_names_map <- c(
    vh  = "S1-VH",
    vv  = "S1-VV",
    b01 = "S2-B1",
    b02 = "S2-B2",
    b03 = "S2-B3",
    b04 = "S2-B4",
    b05 = "S2-B5",
    b06 = "S2-B6",
    b07 = "S2-B7",
    b08 = "S2-B8",
    b8a = "S2-B8A",
    b09 = "S2-B9",
    b11 = "S2-B11",
    b12 = "S2-B12"
  )

  vip_unscaled <- caret::varImp(rf_fit, scale = FALSE)$importance %>%
    tibble::rownames_to_column("predictor") %>%
    dplyr::rename(IncMSE_unscaled = Overall)

  vip_scaled <- caret::varImp(rf_fit, scale = TRUE)$importance %>%
    tibble::rownames_to_column("predictor") %>%
    dplyr::rename(Importance_scaled_0_100 = Overall)

  vip_df <- vip_unscaled %>%
    dplyr::left_join(vip_scaled, by = "predictor") %>%
    dplyr::mutate(
      model_id = model_id,
      model_description = model_description,
      subset_name = subset_name,
      response = response_col,
      predictor_label = dplyr::recode(predictor, !!!var_names_map)
    ) %>%
    dplyr::select(
      model_id, model_description, subset_name, response,
      predictor, predictor_label, IncMSE_unscaled, Importance_scaled_0_100
    ) %>%
    dplyr::arrange(dplyr::desc(IncMSE_unscaled))

  return(list(
    model = rf_fit,
    summary = summary_row,
    folds_raw = res,
    importance = vip_df
  ))
}
