##############################################################################
# 4. Predict Wall-to-Wall Landscape VSD & Uncertainty and Run Spatial Sampling
##############################################################################

# 1. Dependencies
    # Packages
        library(tidyverse)
        library(terra)
        library(sf)
        library(randomForest)
        library(exactextractr)
        library(effectsize)

    # Load project functions
        function_files <- list.files(file.path("code", "functions"), pattern = "\\.[rR]$")
        sapply(file.path("code", "functions", function_files), source)

# 2. Load Spatial Layers and Fitted Full Random Forest Model (M3)
    area_mask         <- terra::vect(file.path("outputs", "area_mask.gpkg"))
    sentinel_1_masked <- terra::rast(file.path("outputs", "sentinel_1_masked.tif"))
    sentinel_2_masked <- terra::rast(file.path("outputs", "sentinel_2_masked.tif"))
    rf_models_all     <- readRDS(file.path("outputs", "rf_models_all.rds"))
    full_model_m3     <- rf_models_all$models_m1_m5$M3_Full_S1S2$model

# 3. Predict Wall-to-Wall Shannon VSD and Tree-Level Standard Deviation
    landscape_preds <- PredictLandscapeVSDAndUncertainty(
      s1_masked = sentinel_1_masked,
      s2_masked = sentinel_2_masked,
      rf_model  = full_model_m3,
      agg_fact  = 3
    )

    out_rast    <- landscape_preds$out_rast
    out_sd_rast <- landscape_preds$out_sd_rast

# 4. Summarise Predictions & Uncertainty by UKCEH LCM Vegetation Class
    lcm_path <- file.path("data", "UKCEH_LCM", "LCM.tif")

    lcm_results <- SummariseLandscapeByLCM(
      lcm_path    = lcm_path,
      out_rast    = out_rast,
      out_sd_rast = out_sd_rast,
      area_mask   = area_mask
    )

    print(lcm_results$lcm_area_summary)
    print(lcm_results$forest_stats)
    print(lcm_results$whole_raster_tests)

# 5. Spatially Separated Random Sampling (200 Forest vs 200 Non-Forest 30m Plots)
    sampling_results <- SampleLandscapeAndTest(
      lcm_path    = lcm_path,
      out_rast    = out_rast,
      out_sd_rast = out_sd_rast,
      area_mask   = area_mask,
      target_n    = 200,
      min_dist    = 100,
      plot_side   = 30,
      seed        = 123
    )

    print(sampling_results$sampled_tests_summary)

# 6. Save Rasters, Vectors, and Statistical CSVs
    terra::writeRaster(
      out_rast,
      file.path("outputs", "predicted_shannon.tif"),
      overwrite = TRUE
    )
    terra::writeRaster(
      out_sd_rast,
      file.path("outputs", "predicted_shannon_sd.tif"),
      overwrite = TRUE
    )
    terra::writeRaster(
      lcm_results$LCM_aligned,
      file.path("outputs", "lcm_aligned.tif"),
      overwrite = TRUE
    )
    terra::writeVector(
      sampling_results$sample_plots_vect,
      file.path("outputs", "sampled_400_plots.gpkg"),
      overwrite = TRUE
    )
    saveRDS(
      lcm_results$plot_df,
      file.path("outputs", "pixel_predictions_by_lcm.rds")
    )

    # Save structured statistical CSVs
    write_csv(
      lcm_results$lcm_area_summary,
      file.path("outputs", "stats_lcm_area_and_zonal_summary.csv")
    )
    write_csv(
      lcm_results$veg_stats,
      file.path("outputs", "stats_landscape_vegetation_class_summary.csv")
    )
    write_csv(
      lcm_results$forest_stats,
      file.path("outputs", "stats_landscape_forest_vs_nonforest_summary.csv")
    )
    write_csv(
      lcm_results$whole_raster_tests,
      file.path("outputs", "stats_landscape_whole_raster_tests.csv")
    )
    write_csv(
      sampling_results$plot_stats,
      file.path("outputs", "stats_sampled_400_plots_values.csv")
    )
    write_csv(
      sampling_results$sampled_tests_summary,
      file.path("outputs", "stats_sampled_400_plots_tests.csv")
    )
