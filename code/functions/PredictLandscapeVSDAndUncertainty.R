# ==============================================================================
# PredictLandscapeVSDAndUncertainty.R
# Predicts Vegetation Structural Diversity (Shannon index) and tree-level
# standard deviation (prediction uncertainty) across the Wild Ennerdale study
# area using the combined S1+S2 Random Forest model (M3).
# ==============================================================================

PredictLandscapeVSDAndUncertainty <- function(
  s1_masked,
  s2_masked,
  rf_model,
  agg_fact = 3
) {
  base_frame <- terra::as.polygons(
    terra::aggregate(s1_masked, fact = agg_fact),
    aggregate = FALSE
  )
  base_frame$id <- seq_len(nrow(base_frame))

  # Extract S2 values directly from native EPSG:32630 raster
  s2_points <- exactextractr::exact_extract(
    s2_masked,
    sf::st_as_sf(base_frame),
    fun = "mean",
    append_cols = "id"
  )

  x_names <- rf_model$finalModel$xNames

  fit_frame <- base_frame %>%
    as.data.frame() %>%
    dplyr::mutate(
      VH = dplyr::if_else(is.infinite(VH), NA_real_, VH),
      VV = dplyr::if_else(is.infinite(VV), NA_real_, VV)
    ) %>%
    dplyr::left_join(s2_points, by = "id") %>%
    dplyr::rename_with(
      ~ stringr::str_remove(.x, "^mean\\."),
       dplyr::starts_with("mean.")
    ) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::select(id, dplyr::all_of(x_names))

  pred_data <- fit_frame %>%
    dplyr::filter(dplyr::if_all(-id, ~ !is.na(.)))

  # 1. Mean prediction across all trees
  rf_pred_all <- stats::predict(
    rf_model$finalModel,
    newdata = pred_data %>% dplyr::select(-id),
    predict.all = TRUE
  )

  pred_data$predicted_shannon <- rf_pred_all$aggregate
  pred_data$shannon_sd        <- apply(rf_pred_all$individual, 1, stats::sd)

  base_frame$predicted_shannon <- NA_real_
  base_frame$shannon_sd        <- NA_real_

  match_idx <- match(pred_data$id, base_frame$id)
  base_frame$predicted_shannon[match_idx] <- pred_data$predicted_shannon
  base_frame$shannon_sd[match_idx]        <- pred_data$shannon_sd

  template_rast <- terra::aggregate(s1_masked$VH, fact = agg_fact)

  out_rast <- terra::rasterize(
    base_frame,
    template_rast,
    field = "predicted_shannon"
  )
  names(out_rast) <- "predicted_shannon"

  out_sd_rast <- terra::rasterize(
    base_frame,
    template_rast,
    field = "shannon_sd"
  )
  names(out_sd_rast) <- "shannon_sd"

  return(list(
    out_rast = out_rast,
    out_sd_rast = out_sd_rast
  ))
}
