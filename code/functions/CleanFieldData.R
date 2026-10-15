# ==============================================================================
# CleanFieldData.R
# Standardizes column names in the Ennerdale field survey dataset, creates a
# spatial point vector (EPSG:4326), and summarises observed VSD and strata cover.
# ==============================================================================

NamesFix <- function(x) {
  x %>%
    tolower() %>%
    gsub(pattern = " ", replacement = "_", x = .)
}

CleanFieldData <- function(csv_path, crs_epsg = "EPSG:4326") {
  raw_df <- readr::read_csv(
    csv_path,
    name_repair = NamesFix,
    show_col_types = FALSE
  )

  cleaned_df <- raw_df %>%
    dplyr::select(
      plot_id,
      lat,
      lon,
      elev,
      dplyr::matches("^[abcd][1-5]$"),
      av1,
      av2,
      av3,
      av4,
      av5,
      shannon
    )

  pts_vect <- terra::vect(
    cleaned_df,
    geom = c("lon", "lat"),
    crs = crs_epsg
  )

  return(list(
    data_df = cleaned_df,
    points_vect = pts_vect
  ))
}

SummariseFieldData <- function(cleaned_df) {
  response_vars <- c("shannon", "av1", "av2", "av3", "av4", "av5")

  summary_all <- cleaned_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(response_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      n_plots = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      q25 = stats::quantile(value, 0.25, na.rm = TRUE),
      q75 = stats::quantile(value, 0.75, na.rm = TRUE),
      iqr = stats::IQR(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      .groups = "drop"
    )

  return(summary_all)
}
