# ==============================================================================
# SummariseLandscapeByLCM.R
# Aligns UKCEH Land Cover Map (LCM) to the prediction grid and calculates area
# coverage, zonal statistics, and whole-landscape Forest vs Non-forest tests.
# ==============================================================================

GetLCMLabels <- function() {
  data.frame(
    class_value = 1:21,
    class_name = c(
      "Broadleaved woodland",
      "Coniferous woodland",
      "Arable and horticulture",
      "Improved grassland",
      "Neutral grassland",
      "Calcareous grassland",
      "Acid grassland",
      "Fen, marsh and swamp",
      "Heather",
      "Heather grassland",
      "Bog",
      "Inland rock",
      "Saltwater",
      "Freshwater",
      "Supra-littoral rock",
      "Supra-littoral sediment",
      "Littoral rock",
      "Littoral sediment",
      "Saltmarsh",
      "Urban",
      "Suburban"
    ),
    stringsAsFactors = FALSE
  )
}

GetVegClassesOrdered <- function() {
  c(
    "Broadleaved woodland",
    "Coniferous woodland",
    "Heather",
    "Heather grassland",
    "Acid grassland",
    "Calcareous grassland",
    "Neutral grassland",
    "Improved grassland",
    "Arable and horticulture",
    "Bog",
    "Fen, marsh and swamp"
  )
}

SummariseLandscapeByLCM <- function(lcm_path, out_rast, out_sd_rast, area_mask) {
  LCM_labels <- GetLCMLabels()
  veg_classes_ordered <- GetVegClassesOrdered()

  LCM_aligned <- terra::rast(lcm_path)$LCM_1 %>%
    terra::as.factor() %>%
    terra::project(out_sd_rast, method = "near") %>%
    terra::mask(area_mask)

  # 1. Zonal means and area coverage across all LCM classes present
  pred_means_lcm <- terra::zonal(out_rast, LCM_aligned, fun = "mean", na.rm = TRUE) %>%
    dplyr::rename(class_value = "LCM_1") %>%
    dplyr::mutate(class_value = as.numeric(class_value))

  uncert_means_lcm <- terra::zonal(out_sd_rast, LCM_aligned, fun = "mean", na.rm = TRUE) %>%
    dplyr::rename(class_value = "LCM_1") %>%
    dplyr::mutate(class_value = as.numeric(class_value))

  lcm_area_summary <- terra::freq(LCM_aligned) %>%
    as.data.frame() %>%
    dplyr::rename(
      class_value = value,
      n_cells = count
    ) %>%
    dplyr::mutate(
      class_value = as.numeric(class_value),
      area_m2 = n_cells * prod(terra::res(LCM_aligned)),
      area_ha = area_m2 / 10000
    ) %>%
    dplyr::left_join(LCM_labels, by = "class_value") %>%
    dplyr::mutate(
      percent_cover = 100 * area_ha / sum(area_ha)
    ) %>%
    dplyr::left_join(pred_means_lcm, by = "class_value") %>%
    dplyr::left_join(uncert_means_lcm, by = "class_value") %>%
    dplyr::arrange(dplyr::desc(area_ha))

  # 2. Pixel-level table for vegetated classes (1:11)
  plot_df <- tibble::tibble(
    prediction  = terra::values(out_rast)[, 1],
    uncertainty = terra::values(out_sd_rast)[, 1],
    class_value = terra::values(LCM_aligned)[, 1]
  ) %>%
    dplyr::filter(
      !is.na(prediction),
      !is.na(uncertainty),
      !is.na(class_value)
    ) %>%
    dplyr::mutate(class_value = as.numeric(class_value)) %>%
    dplyr::left_join(LCM_labels, by = "class_value") %>%
    dplyr::filter(class_name %in% veg_classes_ordered) %>%
    dplyr::mutate(
      forest_group = dplyr::if_else(
        class_name %in% c("Broadleaved woodland", "Coniferous woodland"),
        "Forested",
        "Non-forested"
      )
    )

  veg_stats <- plot_df %>%
    dplyr::group_by(forest_group, class_value, class_name) %>%
    dplyr::summarise(
      n_pixels = dplyr::n(),
      mean_prediction = mean(prediction),
      sd_prediction = stats::sd(prediction),
      median_prediction = stats::median(prediction),
      iqr_prediction = stats::IQR(prediction),
      mean_uncertainty = mean(uncertainty),
      sd_uncertainty = stats::sd(uncertainty),
      median_uncertainty = stats::median(uncertainty),
      iqr_uncertainty = stats::IQR(uncertainty),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_prediction))

  forest_stats <- plot_df %>%
    dplyr::group_by(forest_group) %>%
    dplyr::summarise(
      n_pixels = dplyr::n(),
      mean_prediction = mean(prediction),
      sd_prediction = stats::sd(prediction),
      median_prediction = stats::median(prediction),
      iqr_prediction = stats::IQR(prediction),
      mean_uncertainty = mean(uncertainty),
      sd_uncertainty = stats::sd(uncertainty),
      median_uncertainty = stats::median(uncertainty),
      iqr_uncertainty = stats::IQR(uncertainty),
      .groups = "drop"
    )

  # Whole-raster tests (Prediction and Uncertainty: Forested vs Non-forested)
  t_pred <- stats::t.test(prediction ~ forest_group, data = plot_df)
  w_pred <- stats::wilcox.test(prediction ~ forest_group, data = plot_df, exact = FALSE)
  d_pred <- effectsize::cohens_d(prediction ~ forest_group, data = plot_df)

  t_unc <- stats::t.test(uncertainty ~ forest_group, data = plot_df)
  w_unc <- stats::wilcox.test(uncertainty ~ forest_group, data = plot_df, exact = FALSE)
  d_unc <- effectsize::cohens_d(uncertainty ~ forest_group, data = plot_df)

  whole_raster_tests <- tibble::tibble(
    metric = c("Predicted_Shannon_VSD", "Prediction_Uncertainty_SD"),
    comparison = "Whole-raster pixels: Forested vs Non-forested",
    n_forested = sum(plot_df$forest_group == "Forested"),
    n_non_forested = sum(plot_df$forest_group == "Non-forested"),
    mean_forested = c(
      mean(plot_df$prediction[plot_df$forest_group == "Forested"]),
      mean(plot_df$uncertainty[plot_df$forest_group == "Forested"])
    ),
    sd_forested = c(
      stats::sd(plot_df$prediction[plot_df$forest_group == "Forested"]),
      stats::sd(plot_df$uncertainty[plot_df$forest_group == "Forested"])
    ),
    mean_non_forested = c(
      mean(plot_df$prediction[plot_df$forest_group == "Non-forested"]),
      mean(plot_df$uncertainty[plot_df$forest_group == "Non-forested"])
    ),
    sd_non_forested = c(
      stats::sd(plot_df$prediction[plot_df$forest_group == "Non-forested"]),
      stats::sd(plot_df$uncertainty[plot_df$forest_group == "Non-forested"])
    ),
    t_statistic = c(unname(t_pred$statistic), unname(t_unc$statistic)),
    t_df = c(unname(t_pred$parameter), unname(t_unc$parameter)),
    t_p_value = c(t_pred$p.value, t_unc$p.value),
    wilcox_W = c(unname(w_pred$statistic), unname(w_unc$statistic)),
    wilcox_p_value = c(w_pred$p.value, w_unc$p.value),
    cohens_d = c(d_pred$Cohens_d, d_unc$Cohens_d),
    cohens_d_ci_low = c(d_pred$CI_low, d_unc$CI_low),
    cohens_d_ci_high = c(d_pred$CI_high, d_unc$CI_high)
  )

  return(list(
    LCM_aligned = LCM_aligned,
    lcm_area_summary = lcm_area_summary,
    plot_df = plot_df,
    veg_stats = veg_stats,
    forest_stats = forest_stats,
    whole_raster_tests = whole_raster_tests
  ))
}
