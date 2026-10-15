# ==============================================================================
# PlotHelpers.R
# Modular ggplot2 figure constructors used by 5_MakeFigures.R.
# ==============================================================================

GetVegPalette <- function() {
  c(
    "Broadleaved woodland"    = "#8E0C3A",
    "Coniferous woodland"     = "#D53E4F",
    "Heather"                 = "#F46D43",
    "Heather grassland"       = "#FDAE61",
    "Acid grassland"          = "#FEE08B",
    "Calcareous grassland"    = "#D9EF8B",
    "Improved grassland"      = "#A6D96A",
    "Neutral grassland"       = "#66C2A5",
    "Arable and horticulture" = "#3288BD",
    "Bog"                     = "#5E4FA2",
    "Fen, marsh and swamp"    = "#7B3294"
  )
}

GetRobustLegend <- function(plot) {
  g <- ggplot2::ggplotGrob(plot)
  grob_names <- sapply(g$grobs, function(x) x$name)
  idx <- grep("guide-box", grob_names)
  if (length(idx) > 0) {
    for (i in idx) {
      if (!inherits(g$grobs[[i]], "zeroGrob")) {
        return(g$grobs[[i]])
      }
    }
  }
  return(NULL)
}

MakeVariableImportancePlot <- function(vip_df, title_text = NULL) {
  ggplot2::ggplot(
    vip_df,
    ggplot2::aes(
      x = stats::reorder(predictor_label, IncMSE_unscaled),
      y = IncMSE_unscaled
    )
  ) +
    ggplot2::geom_col(
      fill = "#3A7D44",
      colour = "black",
      linewidth = 0.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = "Variable Importance"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12.5, hjust = 0.5),
      axis.text  = ggplot2::element_text(size = 10.5, colour = "black"),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      axis.line  = ggplot2::element_line(colour = "black")
    )
}

MakePredictionMap <- function(out_rast, area_mask) {
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = out_rast, na.rm = TRUE) +
    ggplot2::scale_fill_viridis_c(
      name = "Predicted VSD",
      na.value = "transparent"
    ) +
    ggplot2::geom_sf(
      data = area_mask,
      fill = NA,
      colour = "black",
      linewidth = 0.2
    ) +
    ggspatial::annotation_scale(location = "bl") +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    ggplot2::coord_sf() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 10),
      legend.title     = ggplot2::element_text(size = 11, face = "bold"),
      legend.text      = ggplot2::element_text(size = 9)
    )
}

MakeUncertaintyMap <- function(out_sd_rast, area_mask) {
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = out_sd_rast, na.rm = TRUE) +
    ggplot2::geom_sf(
      data = area_mask,
      fill = NA,
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_gradientn(
      name = "Prediction Uncertainty (SD)",
      colours = c("white", "red", "black"),
      na.value = "transparent"
    ) +
    ggspatial::annotation_scale(location = "bl") +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    ggplot2::coord_sf() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 10),
      legend.title     = ggplot2::element_text(size = 11, face = "bold"),
      legend.text      = ggplot2::element_text(size = 9)
    )
}

MakeVegRaincloudPlot <- function(plot_df, y_col, y_label, x_angle = 50, min_pixels = 10) {
  veg_classes_ordered <- GetVegClassesOrdered()
  veg_palette <- GetVegPalette()

  df_sub <- plot_df %>%
    dplyr::add_count(class_name, name = "n_class_pixels") %>%
    dplyr::filter(n_class_pixels >= min_pixels) %>%
    dplyr::mutate(
      class_name = factor(class_name, levels = intersect(veg_classes_ordered, unique(class_name))),
      y_val = .data[[y_col]]
    )

  ggplot2::ggplot(
    df_sub,
    ggplot2::aes(x = class_name, y = y_val, fill = class_name)
  ) +
    ggdist::stat_halfeye(
      adjust = 0.5,
      width = 0.7,
      .width = 0,
      justification = -0.2,
      point_colour = NA,
      alpha = 0.8
    ) +
    ggplot2::geom_boxplot(
      width = 0.12,
      fill = "white",
      colour = "black",
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      width = 0.05,
      alpha = 0.03,
      size = 0.2,
      colour = "grey40"
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2.5,
      fill = "gold",
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = veg_palette) +
    ggplot2::labs(x = "Vegetation Class", y = y_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1, size = 10, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      axis.title  = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

MakeCombinedLandscapeAndRaincloudFigure <- function(out_rast, out_sd_rast, area_mask, plot_df) {
  # 1. Build horizontal-legend donor plots and extract standalone legend grobs
  p_vsd_with_leg <- MakePredictionMap(out_rast, area_mask) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(6.5, "cm"),
        barheight = grid::unit(0.35, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )

  p_sd_with_leg <- MakeUncertaintyMap(out_sd_rast, area_mask) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(6.5, "cm"),
        barheight = grid::unit(0.35, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )

  leg_vsd <- GetRobustLegend(p_vsd_with_leg)
  leg_sd  <- GetRobustLegend(p_sd_with_leg)

  # 2. Hide inline legends on main map panels and add panel tags A and B
  tag_theme <- ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 15, colour = "black"),
    plot.tag.position = c(0.02, 0.98),
    plot.margin = ggplot2::margin(6, 14, 6, 6)
  )

  p_map_a <- MakePredictionMap(out_rast, area_mask) +
    ggplot2::labs(tag = "A") +
    ggplot2::theme(legend.position = "none") +
    tag_theme

  p_map_b <- MakeUncertaintyMap(out_sd_rast, area_mask) +
    ggplot2::labs(tag = "B") +
    ggplot2::theme(legend.position = "none") +
    tag_theme

  # 3. Constrain aspect.ratio on raincloud plots to match map proportions and add tags C and D
  p_rain_c <- MakeVegRaincloudPlot(plot_df, y_col = "prediction", y_label = "Predicted VSD", x_angle = 45) +
    ggplot2::labs(tag = "C") +
    ggplot2::theme(aspect.ratio = 0.52) +
    tag_theme

  p_rain_d <- MakeVegRaincloudPlot(plot_df, y_col = "uncertainty", y_label = "Prediction Uncertainty (SD)", x_angle = 45) +
    ggplot2::labs(tag = "D") +
    ggplot2::theme(aspect.ratio = 0.52) +
    tag_theme

  # 4. Assemble using a flat grid layout (no column nesting)
  fig_combined <- p_map_a + p_map_b +
    patchwork::wrap_elements(full = leg_vsd) + patchwork::wrap_elements(full = leg_sd) +
    p_rain_c + p_rain_d +
    patchwork::plot_layout(
      ncol = 2,
      heights = c(1, 0.16, 1.05)
    )

  return(fig_combined)
}

MakeModelCVComparisonPlot <- function(folds_raw_df, title_text = NULL) {
  model_order <- c(
    "M1_S1_Only",
    "M2_S2_Only",
    "M3_Full_S1S2",
    "M4_Forest_LCM",
    "M5_NonForest_LCM"
  )
  model_labels <- c(
    "M1_S1_Only"       = "M1: S1 Only (All, n=74)",
    "M2_S2_Only"       = "M2: S2 Only (All, n=74)",
    "M3_Full_S1S2"     = "M3: S1+S2 (All, n=74)",
    "M4_Forest_LCM"    = "M3 Forest (n=56)",
    "M5_NonForest_LCM" = "M3 Non-forest (n=18)"
  )
  spec_groups <- c(
    "M1_S1_Only"       = "S1 Only (M1)",
    "M2_S2_Only"       = "S2 Only (M2)",
    "M3_Full_S1S2"     = "S1+S2 Full Dataset (M3)",
    "M4_Forest_LCM"    = "S1+S2 Data Subsets of M3",
    "M5_NonForest_LCM" = "S1+S2 Data Subsets of M3"
  )
  spec_palette <- c(
    "S1 Only (M1)"             = "#4A7C9B",
    "S2 Only (M2)"             = "#8C6BB1",
    "S1+S2 Full Dataset (M3)"  = "#3A7D44",
    "S1+S2 Data Subsets of M3" = "#88C492"
  )

  df_plot <- folds_raw_df %>%
    dplyr::filter(model_id %in% model_order) %>%
    dplyr::mutate(
      model_label = factor(
        dplyr::recode(model_id, !!!model_labels),
        levels = rev(unname(model_labels))
      ),
      spec_group = factor(
        dplyr::recode(model_id, !!!spec_groups),
        levels = names(spec_palette)
      )
    )

  ggplot2::ggplot(df_plot, ggplot2::aes(x = model_label, y = Rsquared, fill = spec_group)) +
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = 3.48,
      ymin = -Inf, ymax = Inf,
      fill = "#E8F5E9", alpha = 0.65
    ) +
    ggplot2::geom_boxplot(
      alpha = 0.9,
      outlier.alpha = 0.4,
      width = 0.52,
      colour = "black",
      linewidth = 0.35
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2.8,
      fill = "gold",
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = spec_palette,
      name = "Model Specification"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = "Cross-Validation R²"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12.5, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 10.5, colour = "black"),
      axis.text.x = ggplot2::element_text(size = 10.5, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

MakeCombinedModelCVAndVIPFigure <- function(folds_raw_df, vip_df) {
  tag_theme <- ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 15, colour = "black"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = ggplot2::margin(20, 12, 8, 8)
  )

  p_a <- MakeModelCVComparisonPlot(folds_raw_df, title_text = NULL) +
    ggplot2::labs(tag = "A") +
    tag_theme

  p_b <- MakeVariableImportancePlot(vip_df, title_text = NULL) +
    ggplot2::labs(tag = "B") +
    tag_theme

  fig_combined <- p_a + p_b +
    patchwork::plot_layout(
      ncol = 2,
      widths = c(1.15, 1)
    )

  return(fig_combined)
}

MakeStrataComparisonPlot <- function(strata_long_df) {
  strata_order <- c("av1", "av2", "av3", "av4", "av5")
  strata_names <- c(
    "av1" = "1: 0–0.5 m",
    "av2" = "2: 0.5–1 m",
    "av3" = "3: 1–2 m",
    "av4" = "4: 2–4 m",
    "av5" = "5: >4 m"
  )
  subset_map <- c(
    "All_Plots (n=74)"     = "All Plots (n=74)",
    "All Plots (n=74)"     = "All Plots (n=74)",
    "LCM_Forest (n=56)"    = "Forest (n=56)",
    "Forest (n=56)"        = "Forest (n=56)",
    "LCM_NonForest (n=18)" = "Non-forest (n=18)",
    "Non-forest (n=18)"    = "Non-forest (n=18)"
  )
  subset_levels <- c("All Plots (n=74)", "Forest (n=56)", "Non-forest (n=18)")

  df_plot <- strata_long_df %>%
    dplyr::filter(response %in% strata_order) %>%
    dplyr::mutate(
      stratum = factor(dplyr::recode(response, !!!strata_names), levels = unname(strata_names)),
      subset_clean = factor(dplyr::recode(subset_name, !!!subset_map), levels = subset_levels),
      sensor_clean = dplyr::recode(
        sensor_set,
        "S1+S2 (M6)"   = "S1+S2 (M6)",
        "S2_Only (M7)" = "S2 Only (M7)",
        "S2 Only (M7)" = "S2 Only (M7)"
      )
    )

  pd <- ggplot2::position_dodge(width = 0.2)

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = stratum, y = Rsquared_mean, colour = sensor_clean, group = sensor_clean)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, Rsquared_mean - Rsquared_sd),
        ymax = pmin(1, Rsquared_mean + Rsquared_sd)
      ),
      width = 0.18,
      linewidth = 0.55,
      alpha = 0.75,
      position = pd,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(linewidth = 0.9, position = pd, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.8, position = pd, na.rm = TRUE) +
    ggplot2::facet_wrap(~ subset_clean, nrow = 1) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_colour_manual(
      values = c("S1+S2 (M6)" = "#3A7D44", "S2 Only (M7)" = "#8C6BB1"),
      name = "Model"
    ) +
    ggplot2::labs(
      title = NULL,
      x = "Vertical Vegetation Stratum",
      y = "Cross-Validation R²"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey94", colour = "black", linewidth = 0.5),
      strip.text       = ggplot2::element_text(face = "bold", size = 11, colour = "black"),
      axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1, size = 10, colour = "black"),
      axis.text.y      = ggplot2::element_text(size = 10, colour = "black"),
      axis.title       = ggplot2::element_text(size = 11, face = "bold"),
      legend.title     = ggplot2::element_text(size = 10.5, face = "bold"),
      legend.text      = ggplot2::element_text(size = 10),
      legend.position  = "bottom",
      panel.spacing    = grid::unit(1, "lines")
    )
}
