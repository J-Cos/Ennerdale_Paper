#####################################################
# 5. Generate Publication Figures (300 DPI PNG)
#####################################################

# 1. Dependencies
    # Packages
        library(tidyverse)
        library(terra)
        library(tidyterra)
        library(sf)
        library(ggspatial)
        library(ggdist)
        library(patchwork)

    # Load project functions
        function_files <- list.files(file.path("code", "functions"), pattern = "\\.[rR]$")
        sapply(file.path("code", "functions", function_files), source)

# 2. Load Saved Outputs from Previous Steps
    vip_m3_df      <- read_csv(file.path("outputs", "stats_variable_importance_M3.csv"), show_col_types = FALSE)
    folds_m1_m5_df <- read_csv(file.path("outputs", "stats_model_performance_M1_M5_100folds_raw.csv"), show_col_types = FALSE)
    strata_long_df <- read_csv(file.path("outputs", "stats_strata_performance_M6_M7_long.csv"), show_col_types = FALSE)
    area_mask      <- terra::vect(file.path("outputs", "area_mask.gpkg"))
    out_rast       <- terra::rast(file.path("outputs", "predicted_shannon.tif"))
    out_sd_rast    <- terra::rast(file.path("outputs", "predicted_shannon_sd.tif"))
    plot_df        <- readRDS(file.path("outputs", "pixel_predictions_by_lcm.rds"))

# 3. Construct Figures
    fig_vip_m3 <- MakeVariableImportancePlot(
      vip_df = vip_m3_df
    )

    fig_prediction_map <- MakePredictionMap(
      out_rast  = out_rast,
      area_mask = area_mask
    )

    fig_uncertainty_map <- MakeUncertaintyMap(
      out_sd_rast = out_sd_rast,
      area_mask   = area_mask
    )

    fig_vsd_violin <- MakeVegRaincloudPlot(
      plot_df = plot_df,
      y_col   = "prediction",
      y_label = "Predicted VSD",
      x_angle = 50
    )

    fig_uncertainty_violin <- MakeVegRaincloudPlot(
      plot_df = plot_df,
      y_col   = "uncertainty",
      y_label = "Standard Deviation",
      x_angle = 60
    )

    fig_combined_4panel <- MakeCombinedLandscapeAndRaincloudFigure(
      out_rast    = out_rast,
      out_sd_rast = out_sd_rast,
      area_mask   = area_mask,
      plot_df     = plot_df
    )

    fig_cv_distributions <- MakeModelCVComparisonPlot(
      folds_raw_df = folds_m1_m5_df
    )

    fig_combined_cv_vip <- MakeCombinedModelCVAndVIPFigure(
      folds_raw_df = folds_m1_m5_df,
      vip_df       = vip_m3_df
    )

    fig_strata_comparison <- MakeStrataComparisonPlot(
      strata_long_df = strata_long_df
    )

# 4. Save All Figures to figures/ as 300 DPI PNG
    ggsave(
      filename = file.path("figures", "Variable_importance.png"),
      plot     = fig_vip_m3,
      width    = 240,
      height   = 180,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "prediction_map.png"),
      plot     = fig_prediction_map,
      width    = 25,
      height   = 20,
      units    = "cm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "uncertainty_map.png"),
      plot     = fig_uncertainty_map,
      width    = 25,
      height   = 20,
      units    = "cm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "VSD_violin_plot.png"),
      plot     = fig_vsd_violin,
      width    = 240,
      height   = 180,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "uncertainty_violin_plot.png"),
      plot     = fig_uncertainty_violin,
      width    = 240,
      height   = 180,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "combined_predictions_uncertainty_4panel.png"),
      plot     = fig_combined_4panel,
      width    = 340,
      height   = 230,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "model_cv_r2_distributions.png"),
      plot     = fig_cv_distributions,
      width    = 220,
      height   = 145,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "combined_model_cv_and_vip_2panel.png"),
      plot     = fig_combined_cv_vip,
      width    = 330,
      height   = 155,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )

    ggsave(
      filename = file.path("figures", "vertical_strata_comparison.png"),
      plot     = fig_strata_comparison,
      width    = 270,
      height   = 135,
      units    = "mm",
      dpi      = 300,
      bg       = "white"
    )
