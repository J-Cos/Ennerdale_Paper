#######################################################################
# 2. Buffer Field Plots and Extract Native-Resolution Raster Predictors
#######################################################################

# 1. Dependencies
    # Packages
        library(tidyverse)
        library(terra)
        library(tidyterra)

    # Load project functions
        function_files <- list.files(file.path("code", "functions"), pattern = "\\.[rR]$")
        sapply(file.path("code", "functions", function_files), source)

# 2. Load Cleaned Field Points and Masked Rasters from Previous Steps
    field_points_vect <- terra::vect(file.path("outputs", "field_plots_points.gpkg"))
    area_mask         <- terra::vect(file.path("outputs", "area_mask.gpkg"))
    sentinel_1_masked <- terra::rast(file.path("outputs", "sentinel_1_masked.tif"))
    sentinel_2_masked <- terra::rast(file.path("outputs", "sentinel_2_masked.tif"))
    forest_mask       <- terra::rast(file.path("outputs", "forest_mask.tif"))

# 3. Create 30m x 30m Square Plot Polygons and Verify Centroid Alignment
    plots_vect <- SquareBufferLatLon(
      pts        = field_points_vect,
      side       = 30,
      crs_metric = terra::crs(area_mask)
    )

    max_centroid_offset_m <- max(abs(terra::distance(
      terra::centroids(plots_vect),
      terra::project(field_points_vect, terra::crs(plots_vect)),
      pairwise = TRUE
    )))
    cat("Max centroid alignment offset (m):", max_centroid_offset_m, "\n")

# 4. Extract Area-Weighted Pixel Means at Native Resolution
    extracted_df <- ExtractPlotPredictors(
      plots_vect  = plots_vect,
      s1_masked   = sentinel_1_masked,
      s2_masked   = sentinel_2_masked,
      forest_mask = forest_mask
    )

    # Define predictor sets for downstream Random Forest modelling
    s1_predictors   <- c("vh", "vv")
    s2_predictors   <- c("b01", "b02", "b03", "b04", "b05", "b06", "b07", "b08", "b8a", "b09", "b11", "b12")
    s1s2_predictors <- c(s1_predictors, s2_predictors)

    # Compute descriptive summary and correlation matrix of extracted predictors
    predictor_summary <- extracted_df %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(c(s1s2_predictors, "msavi2", "ndvi", "nbr", "forest_cover")),
        names_to = "predictor",
        values_to = "value"
      ) %>%
      dplyr::group_by(predictor) %>%
      dplyr::summarise(
        n_plots = sum(!is.na(value)),
        mean = mean(value, na.rm = TRUE),
        sd = stats::sd(value, na.rm = TRUE),
        median = stats::median(value, na.rm = TRUE),
        iqr = stats::IQR(value, na.rm = TRUE),
        min = min(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        .groups = "drop"
      )

    cor_mat <- stats::cor(
      extracted_df %>% dplyr::select(dplyr::all_of(c("shannon", s1s2_predictors, "msavi2", "ndvi", "nbr"))),
      use = "complete.obs"
    ) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("variable")

# 5. Save Outputs
    terra::writeVector(
      plots_vect,
      file.path("outputs", "field_plots_30m_squares.gpkg"),
      overwrite = TRUE
    )
    write_csv(
      extracted_df,
      file.path("outputs", "plot_extracted_all_predictors.csv")
    )
    saveRDS(
      list(
        extracted_df    = extracted_df,
        s1_predictors   = s1_predictors,
        s2_predictors   = s2_predictors,
        s1s2_predictors = s1s2_predictors
      ),
      file.path("outputs", "model_input_data.rds")
    )
    write_csv(
      predictor_summary,
      file.path("outputs", "stats_extracted_predictors_summary.csv")
    )
    write_csv(
      cor_mat,
      file.path("outputs", "stats_predictor_correlation_matrix.csv")
    )
