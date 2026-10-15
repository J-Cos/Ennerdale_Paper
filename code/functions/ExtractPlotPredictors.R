# ==============================================================================
# ExtractPlotPredictors.R
# Creates 30m x 30m square plot polygons in metric CRS and extracts area-weighted
# pixel values from native-resolution Sentinel-1, Sentinel-2, and LCM rasters
# without resampling the rasters.
# ==============================================================================

SquareBufferLatLon <- function(pts, side, crs_metric) {
  d <- side / 2
  pts_proj <- terra::project(pts, crs_metric)
  xy <- terra::crds(pts_proj)
  n  <- nrow(xy)

  geom <- cbind(
    geom = rep(seq_len(n), each = 5),
    part = 1,
    x    = as.vector(t(cbind(
      xy[, 1] - d, xy[, 1] + d, xy[, 1] + d, xy[, 1] - d, xy[, 1] - d
    ))),
    y    = as.vector(t(cbind(
      xy[, 2] - d, xy[, 2] - d, xy[, 2] + d, xy[, 2] + d, xy[, 2] - d
    ))),
    hole = 0
  )

  v <- terra::vect(geom, type = "polygons", crs = crs_metric)
  terra::values(v) <- terra::values(pts)
  return(v)
}

ExtractPlotPredictors <- function(plots_vect, s1_masked, s2_masked, forest_mask) {
  # Extract S1 at native resolution (EPSG:32630)
  s1_extract <- terra::extract(
    s1_masked,
    plots_vect %>% tidyterra::select(plot_id, shannon),
    fun = mean,
    exact = TRUE,
    na.rm = TRUE,
    bind = TRUE
  )

  # Extract S2 at native resolution (EPSG:32630)
  s2_extract <- terra::extract(
    s2_masked,
    plots_vect %>% tidyterra::select(plot_id, shannon),
    fun = mean,
    exact = TRUE,
    na.rm = TRUE,
    bind = TRUE
  )

  # Project vector plots to LCM CRS (EPSG:27700) to extract exact fraction of binary forest_mask
  plots_27700 <- terra::project(
    plots_vect %>% tidyterra::select(plot_id, shannon),
    terra::crs(forest_mask)
  )
  forest_extract <- terra::extract(
    forest_mask,
    plots_27700,
    fun = mean,
    exact = TRUE,
    na.rm = TRUE,
    bind = TRUE
  )

  shannon_and_strata <- c("shannon", "av1", "av2", "av3", "av4", "av5")

  combined_df <- plots_vect %>%
    tibble::as_tibble() %>%
    dplyr::select(
      plot_id,
      dplyr::all_of(shannon_and_strata)
    ) %>%
    dplyr::left_join(
      s1_extract %>% tibble::as_tibble() %>% dplyr::select(-shannon),
      by = "plot_id"
    ) %>%
    dplyr::left_join(
      s2_extract %>% tibble::as_tibble() %>% dplyr::select(-shannon),
      by = "plot_id"
    ) %>%
    dplyr::left_join(
      forest_extract %>% tibble::as_tibble() %>% dplyr::select(plot_id, forest_cover = LCM_1),
      by = "plot_id"
    ) %>%
    dplyr::rename_with(.fn = tolower) %>%
    dplyr::mutate(
      lcm_is_forest = forest_cover > 0.2
    )

  return(combined_df)
}
