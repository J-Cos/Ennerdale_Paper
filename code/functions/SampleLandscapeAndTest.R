# ==============================================================================
# SampleLandscapeAndTest.R
# Performs spatially separated random sampling (200 Forest + 200 Non-forest
# 30m x 30m plots with >= 100m minimum distance) across valid vegetated pixels
# and tests differences in predicted VSD and prediction uncertainty.
# ==============================================================================

SampleLandscapeAndTest <- function(
  lcm_path,
  out_rast,
  out_sd_rast,
  area_mask,
  target_n = 200,
  min_dist = 100,
  plot_side = 30,
  seed = 123
) {
  set.seed(seed)
  half_side <- plot_side / 2

  LCM <- terra::rast(lcm_path)
  lcm_32630 <- terra::project(LCM$LCM_1, terra::crs(out_rast), method = "near") %>%
    terra::mask(area_mask)

  valid_pred_lcm_grid <- !is.na(terra::resample(out_rast, lcm_32630, method = "near"))

  # 1 = Forest (LCM 1-2), 0 = Non-forest vegetation (LCM 3-11), NA = Rock/Water/Urban/Outside
  veg_class_mask <- terra::ifel(
    !valid_pred_lcm_grid, NA,
    terra::ifel(
      lcm_32630 %in% c(1, 2), 1,
      terra::ifel(lcm_32630 %in% 3:11, 0, NA)
    )
  )

  sample_df <- data.frame()

  for (cl in c(0, 1)) {
    candidate_cells <- which(terra::values(veg_class_mask) == cl)
    if (length(candidate_cells) == 0) next
    candidate_cells <- sample(candidate_cells)
    selected <- numeric(0)

    for (cell in candidate_cells) {
      xy_new <- terra::xyFromCell(veg_class_mask, cell)
      if (length(selected) == 0) {
        selected <- cell
      } else {
        xy_existing <- terra::xyFromCell(veg_class_mask, selected)
        d <- sqrt(
          (xy_existing[, 1] - xy_new[1, 1])^2 +
            (xy_existing[, 2] - xy_new[1, 2])^2
        )
        if (min(d) >= min_dist) {
          selected <- c(selected, cell)
        }
      }
      if (length(selected) >= target_n) break
    }

    xy <- terra::xyFromCell(veg_class_mask, selected)
    sample_df <- rbind(
      sample_df,
      data.frame(
        x     = xy[, 1],
        y     = xy[, 2],
        class = cl
      )
    )
  }

  samples <- terra::vect(sample_df, geom = c("x", "y"), crs = terra::crs(out_rast))
  sample_plots <- terra::buffer(samples, width = half_side, capstyle = "square")

  extracted_means <- terra::extract(
    c(out_rast, out_sd_rast),
    sample_plots,
    fun   = mean,
    exact = TRUE,
    na.rm = TRUE
  )

  plot_stats <- tibble::tibble(
    sample_id = extracted_means$ID,
    x = sample_df$x,
    y = sample_df$y,
    class_code = sample_plots$class,
    class_label = factor(
      sample_plots$class,
      levels = c(0, 1),
      labels = c("Non-forest", "Forest")
    ),
    mean_predicted_shannon = extracted_means$predicted_shannon,
    mean_shannon_sd = extracted_means$shannon_sd
  )

  sample_plots$mean_predicted_shannon <- plot_stats$mean_predicted_shannon
  sample_plots$mean_shannon_sd        <- plot_stats$mean_shannon_sd
  sample_plots$class_label            <- as.character(plot_stats$class_label)

  # Statistical tests on the 400 spatially separated sample plots
  t_vsd <- stats::t.test(mean_predicted_shannon ~ class_label, data = plot_stats)
  w_vsd <- stats::wilcox.test(mean_predicted_shannon ~ class_label, data = plot_stats, exact = FALSE)
  d_vsd <- effectsize::cohens_d(mean_predicted_shannon ~ class_label, data = plot_stats)

  t_sd  <- stats::t.test(mean_shannon_sd ~ class_label, data = plot_stats)
  w_sd  <- stats::wilcox.test(mean_shannon_sd ~ class_label, data = plot_stats, exact = FALSE)
  d_sd  <- effectsize::cohens_d(mean_shannon_sd ~ class_label, data = plot_stats)

  sampled_tests_summary <- tibble::tibble(
    metric = c("Sampled_Predicted_Shannon_VSD", "Sampled_Prediction_Uncertainty_SD"),
    comparison = "Spatially separated 30m x 30m plots (>=100m apart): Non-forest (n=200) vs Forest (n=200)",
    n_non_forest = sum(!is.na(plot_stats$mean_predicted_shannon[plot_stats$class_label == "Non-forest"])),
    n_forest = sum(!is.na(plot_stats$mean_predicted_shannon[plot_stats$class_label == "Forest"])),
    mean_non_forest = c(
      mean(plot_stats$mean_predicted_shannon[plot_stats$class_label == "Non-forest"], na.rm = TRUE),
      mean(plot_stats$mean_shannon_sd[plot_stats$class_label == "Non-forest"], na.rm = TRUE)
    ),
    sd_non_forest = c(
      stats::sd(plot_stats$mean_predicted_shannon[plot_stats$class_label == "Non-forest"], na.rm = TRUE),
      stats::sd(plot_stats$mean_shannon_sd[plot_stats$class_label == "Non-forest"], na.rm = TRUE)
    ),
    mean_forest = c(
      mean(plot_stats$mean_predicted_shannon[plot_stats$class_label == "Forest"], na.rm = TRUE),
      mean(plot_stats$mean_shannon_sd[plot_stats$class_label == "Forest"], na.rm = TRUE)
    ),
    sd_forest = c(
      stats::sd(plot_stats$mean_predicted_shannon[plot_stats$class_label == "Forest"], na.rm = TRUE),
      stats::sd(plot_stats$mean_shannon_sd[plot_stats$class_label == "Forest"], na.rm = TRUE)
    ),
    t_statistic = c(unname(t_vsd$statistic), unname(t_sd$statistic)),
    t_df = c(unname(t_vsd$parameter), unname(t_sd$parameter)),
    t_p_value = c(t_vsd$p.value, t_sd$p.value),
    wilcox_W = c(unname(w_vsd$statistic), unname(w_sd$statistic)),
    wilcox_p_value = c(w_vsd$p.value, w_sd$p.value),
    cohens_d = c(d_vsd$Cohens_d, d_sd$Cohens_d),
    cohens_d_ci_low = c(d_vsd$CI_low, d_sd$CI_low),
    cohens_d_ci_high = c(d_vsd$CI_high, d_sd$CI_high)
  )

  return(list(
    sample_plots_vect = sample_plots,
    plot_stats = plot_stats,
    sampled_tests_summary = sampled_tests_summary
  ))
}
