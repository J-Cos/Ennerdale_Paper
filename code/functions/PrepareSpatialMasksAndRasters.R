# ==============================================================================
# PrepareSpatialMasksAndRasters.R
# Loads study boundary, erases buildings, loads and names Sentinel-1 and
# Sentinel-2 rasters, computes spectral indices (with corrected MSAVI2), and
# prepares UKCEH LCM forest and non-forest masks.
# ==============================================================================

CreateStudyAreaMask <- function(boundary_shp_path, buildings_shp_path) {
  enn_subLake   <- terra::vect(boundary_shp_path)
  enn_buildings <- terra::vect(buildings_shp_path)
  area_mask     <- terra::erase(enn_subLake, enn_buildings)
  if ("fid" %in% names(area_mask)) {
    area_mask$fid <- NULL
  }
  return(area_mask)
}

LoadAndMaskSentinel1 <- function(raster_dir, area_mask) {
  s1_files <- list.files(path = raster_dir, pattern = "Sentinel-1", full.names = TRUE)
  s1_rast  <- terra::rast(s1_files)
  names(s1_rast) <- ifelse(
    grepl("VH", basename(s1_files), ignore.case = TRUE), "VH",
    ifelse(grepl("VV", basename(s1_files), ignore.case = TRUE), "VV", basename(s1_files))
  )
  s1_masked <- terra::mask(s1_rast, area_mask, inverse = FALSE)
  return(s1_masked)
}

LoadAndMaskSentinel2 <- function(raster_dir, area_mask) {
  s2_files <- list.files(path = raster_dir, pattern = "Sentinel-2", full.names = TRUE)
  s2_rast  <- terra::rast(s2_files)
  names(s2_rast) <- stringr::str_extract(basename(s2_files), "B(\\d{2}|8A)")

  # Compute vegetation indices with analytically corrected MSAVI2 on [0, 1] reflectance scale
  s2_rast <- s2_rast %>%
    tidyterra::mutate(
      MSAVI2 = (2 * (B08 / 10000) + 1 - sqrt(pmax(0, (2 * (B08 / 10000) + 1)^2 - 8 * ((B08 - B04) / 10000)))) / 2,
      NDVI   = (B08 - B04) / (B08 + B04),
      NBR    = (B08 - B12) / (B08 + B12)
    )

  s2_masked <- terra::mask(s2_rast, area_mask, inverse = FALSE)
  return(s2_masked)
}

CreateLCMMasks <- function(lcm_path, area_mask) {
  LCM <- terra::rast(lcm_path)
  area_mask_lcm <- terra::project(area_mask, terra::crs(LCM))

  forest_mask <- (LCM$LCM_1 %in% c(1, 2)) %>%
    terra::mask(area_mask_lcm, inverse = FALSE)
  names(forest_mask) <- "LCM_1"

  non_forest_mask <- (LCM$LCM_1 %in% c(3, 4, 5, 6, 7, 8, 9, 10, 11)) %>%
    terra::mask(area_mask_lcm, inverse = FALSE)
  names(non_forest_mask) <- "LCM_1"

  return(list(
    forest_mask = forest_mask,
    non_forest_mask = non_forest_mask
  ))
}
