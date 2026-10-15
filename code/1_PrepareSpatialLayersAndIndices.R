##############################################################
# 1. Prepare Spatial Masks, Sentinel Rasters & Spectral Indices
##############################################################

# 1. Dependencies
    # Packages
        library(tidyverse)
        library(terra)
        library(tidyterra)

    # Load project functions
        function_files <- list.files(file.path("code", "functions"), pattern = "\\.[rR]$")
        sapply(file.path("code", "functions", function_files), source)

# 2. Create Study Area Mask (Wild Ennerdale minus Lake and Buildings)
    area_mask <- CreateStudyAreaMask(
      boundary_shp_path  = file.path("data", "Masks", "WE_boundary_subLake.shp"),
      buildings_shp_path = file.path("data", "Masks", "Ennerdale_buildings.shp")
    )

# 3. Load, Compute Spectral Indices, and Mask Sentinel-1 & Sentinel-2 Rasters
    raster_dir <- file.path("data", "rasters_S1_S2")

    sentinel_1_masked <- LoadAndMaskSentinel1(
      raster_dir = raster_dir,
      area_mask  = area_mask
    )

    sentinel_2_masked <- LoadAndMaskSentinel2(
      raster_dir = raster_dir,
      area_mask  = area_mask
    )

# 4. Prepare UKCEH Land Cover Map (LCM) Forest and Non-Forest Masks
    lcm_masks <- CreateLCMMasks(
      lcm_path  = file.path("data", "UKCEH_LCM", "LCM.tif"),
      area_mask = area_mask
    )

# 5. Save Spatial Outputs for Next Steps
    terra::writeVector(
      area_mask,
      file.path("outputs", "area_mask.gpkg"),
      overwrite = TRUE
    )
    terra::writeRaster(
      sentinel_1_masked,
      file.path("outputs", "sentinel_1_masked.tif"),
      overwrite = TRUE
    )
    terra::writeRaster(
      sentinel_2_masked,
      file.path("outputs", "sentinel_2_masked.tif"),
      overwrite = TRUE
    )
    terra::writeRaster(
      lcm_masks$forest_mask,
      file.path("outputs", "forest_mask.tif"),
      overwrite = TRUE
    )
    terra::writeRaster(
      lcm_masks$non_forest_mask,
      file.path("outputs", "non_forest_mask.tif"),
      overwrite = TRUE
    )
