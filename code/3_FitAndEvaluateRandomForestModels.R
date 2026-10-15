###########################################################################
# 3. Fit and Evaluate Random Forest Models (M1-M7) & Test Residual Moran's I
###########################################################################

# 1. Dependencies
    # Packages
        library(tidyverse)
        library(terra)
        library(randomForest)
        library(caret)

    # Load project functions
        function_files <- list.files(file.path("code", "functions"), pattern = "\\.[rR]$")
        sapply(file.path("code", "functions", function_files), source)

    # Global Modelling Parameters (matching Section 2.4.1 of manuscript)
        NTREE     <- 1000
        N_FOLDS   <- 5
        N_REPEATS <- 20
        SEED      <- 123

# 2. Load Extracted Plot Data and Plot Polygons from Step 2
    model_inputs <- readRDS(file.path("outputs", "model_input_data.rds"))
    rf_data      <- model_inputs$extracted_df
    s1_preds     <- model_inputs$s1_predictors
    s2_preds     <- model_inputs$s2_predictors
    s1s2_preds   <- model_inputs$s1s2_predictors
    plots_vect   <- terra::vect(file.path("outputs", "field_plots_30m_squares.gpkg"))

# 3. Fit and Evaluate Overall and Stratified VSD Models (M1 - M5)
    # M3: Full Combined Model (S1 + S2, All Plots, n = 74)
    m3_full <- FitRepeatedCVRandomForest(
      data              = rf_data,
      response_col      = "shannon",
      predictor_cols    = s1s2_preds,
      model_id          = "M3_Full_S1S2",
      model_description = "Combined S1 + S2 (M3)",
      subset_name       = "All_Plots (n=74)",
      ntree             = NTREE,
      n_folds           = N_FOLDS,
      n_repeats         = N_REPEATS,
      seed              = SEED
    )

    # M1: Sentinel-1 Only (All Plots, n = 74)
    m1_s1 <- FitRepeatedCVRandomForest(
      data              = rf_data,
      response_col      = "shannon",
      predictor_cols    = s1_preds,
      model_id          = "M1_S1_Only",
      model_description = "Sentinel-1 Only (M1)",
      subset_name       = "All_Plots (n=74)",
      ntree             = NTREE,
      n_folds           = N_FOLDS,
      n_repeats         = N_REPEATS,
      seed              = SEED
    )

    # M2: Sentinel-2 Only (All Plots, n = 74)
    m2_s2 <- FitRepeatedCVRandomForest(
      data              = rf_data,
      response_col      = "shannon",
      predictor_cols    = s2_preds,
      model_id          = "M2_S2_Only",
      model_description = "Sentinel-2 Only (M2)",
      subset_name       = "All_Plots (n=74)",
      ntree             = NTREE,
      n_folds           = N_FOLDS,
      n_repeats         = N_REPEATS,
      seed              = SEED
    )

    # M4 & M5 (UKCEH LCM Forest Cover > 0.2 [n = 56] and <= 0.2 [n = 18])
    m4_forest_lcm <- FitRepeatedCVRandomForest(
      data              = rf_data %>% dplyr::filter(lcm_is_forest),
      response_col      = "shannon",
      predictor_cols    = s1s2_preds,
      model_id          = "M4_Forest_LCM",
      model_description = "Forested Plots - UKCEH LCM >20% (M4)",
      subset_name       = "LCM_Forest (n=56)",
      ntree             = NTREE,
      n_folds           = N_FOLDS,
      n_repeats         = N_REPEATS,
      seed              = SEED
    )

    m5_nonforest_lcm <- FitRepeatedCVRandomForest(
      data              = rf_data %>% dplyr::filter(!lcm_is_forest),
      response_col      = "shannon",
      predictor_cols    = s1s2_preds,
      model_id          = "M5_NonForest_LCM",
      model_description = "Non-Forested Plots - UKCEH LCM <=20% (M5)",
      subset_name       = "LCM_NonForest (n=18)",
      ntree             = NTREE,
      n_folds           = N_FOLDS,
      n_repeats         = N_REPEATS,
      seed              = SEED
    )

    models_m1_m5_list <- list(
      M3_Full_S1S2     = m3_full,
      M1_S1_Only       = m1_s1,
      M2_S2_Only       = m2_s2,
      M4_Forest_LCM    = m4_forest_lcm,
      M5_NonForest_LCM = m5_nonforest_lcm
    )

    summary_m1_m5 <- dplyr::bind_rows(lapply(models_m1_m5_list, `[[`, "summary"))
    folds_m1_m5   <- dplyr::bind_rows(lapply(models_m1_m5_list, `[[`, "folds_raw"))
    vip_m1_m5     <- dplyr::bind_rows(lapply(models_m1_m5_list, `[[`, "importance"))

    print(summary_m1_m5 %>% dplyr::select(model_id, n_obs, mtry, Rsquared_mean_sd_str, RMSE_mean_sd_str, Rep20_R2_cor_str))

# 4. Test Spatial Autocorrelation (Moran's I) on Model Residuals (M1, M2, M3)
    moran_m3 <- TestResidualMoransI(plots_vect, m3_full$model, rf_data, model_id = "M3_Full_S1S2", seed = SEED)
    moran_m1 <- TestResidualMoransI(plots_vect, m1_s1$model, rf_data, model_id = "M1_S1_Only", seed = SEED)
    moran_m2 <- TestResidualMoransI(plots_vect, m2_s2$model, rf_data, model_id = "M2_S2_Only", seed = SEED)

    morans_i_summary <- dplyr::bind_rows(
      moran_m3$moran_summary,
      moran_m1$moran_summary,
      moran_m2$moran_summary
    )
    residuals_all <- dplyr::bind_rows(
      moran_m3$residuals_df,
      moran_m1$residuals_df,
      moran_m2$residuals_df
    )
    print(morans_i_summary)

# 5. Fit and Evaluate Vertical Strata Models (M6: S1+S2 vs M7: S2 Only across av1-av5)
    # (A) All Plots (n = 74)
    strata_all_plots <- RunVerticalStrataModels(
      df_input        = rf_data,
      s1s2_predictors = s1s2_preds,
      s2_predictors   = s2_preds,
      subset_name     = "All Plots (n=74)",
      ntree           = NTREE,
      seed            = SEED
    )

    # (B) UKCEH LCM Forested Plots (n = 56)
    strata_lcm_forest <- RunVerticalStrataModels(
      df_input        = rf_data %>% dplyr::filter(lcm_is_forest),
      s1s2_predictors = s1s2_preds,
      s2_predictors   = s2_preds,
      subset_name     = "Forest (n=56)",
      ntree           = NTREE,
      seed            = SEED
    )

    # (C) UKCEH LCM Non-Forested Plots (n = 18)
    strata_lcm_nonforest <- RunVerticalStrataModels(
      df_input        = rf_data %>% dplyr::filter(!lcm_is_forest),
      s1s2_predictors = s1s2_preds,
      s2_predictors   = s2_preds,
      subset_name     = "Non-forest (n=18)",
      ntree           = NTREE,
      seed            = SEED
    )

    strata_long_summary <- dplyr::bind_rows(
      strata_all_plots$long_summary,
      strata_lcm_forest$long_summary,
      strata_lcm_nonforest$long_summary
    )

    strata_comparison_delta <- dplyr::bind_rows(
      strata_all_plots$comparison_df,
      strata_lcm_forest$comparison_df,
      strata_lcm_nonforest$comparison_df
    )

    strata_vip_all <- dplyr::bind_rows(
      strata_all_plots$importance,
      strata_lcm_forest$importance,
      strata_lcm_nonforest$importance
    )

    print(strata_comparison_delta %>% dplyr::select(subset_name, response, S1S2_Rsquared_mean, S2_Rsquared_mean, Delta_Rsquared_mean))

# 6. Save Fitted Models and All Statistical Output CSVs
    saveRDS(
      list(
        models_m1_m5         = models_m1_m5_list,
        strata_all_plots     = strata_all_plots,
        strata_lcm_forest    = strata_lcm_forest,
        strata_lcm_nonforest = strata_lcm_nonforest
      ),
      file.path("outputs", "rf_models_all.rds")
    )

    write_csv(
      summary_m1_m5,
      file.path("outputs", "stats_model_performance_M1_M5.csv")
    )
    write_csv(
      folds_m1_m5,
      file.path("outputs", "stats_model_performance_M1_M5_100folds_raw.csv")
    )
    write_csv(
      m3_full$importance,
      file.path("outputs", "stats_variable_importance_M3.csv")
    )
    write_csv(
      vip_m1_m5,
      file.path("outputs", "stats_variable_importance_M1_M5.csv")
    )
    write_csv(
      morans_i_summary,
      file.path("outputs", "stats_morans_i_residuals.csv")
    )
    write_csv(
      residuals_all,
      file.path("outputs", "stats_plot_predictions_and_residuals.csv")
    )
    write_csv(
      strata_long_summary,
      file.path("outputs", "stats_strata_performance_M6_M7_long.csv")
    )
    write_csv(
      strata_comparison_delta,
      file.path("outputs", "stats_strata_comparison_M6_vs_M7_delta.csv")
    )
    write_csv(
      strata_vip_all,
      file.path("outputs", "stats_strata_variable_importance_M6_M7.csv")
    )
