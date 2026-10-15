#########################################################
# 0. Clean Field Survey Data and Compute Field Statistics
#########################################################

# 1. Dependencies
    # Packages
        library(tidyverse)
        library(terra)

    # Load project functions
        function_files <- list.files(file.path("code", "functions"), pattern = "\\.[rR]$")
        sapply(file.path("code", "functions", function_files), source)

# 2. Load and Clean Raw Field Data
    raw_csv_path <- file.path("data", "Clean_field_data_Ennerdale_22042026.csv")
    cleaned_field <- CleanFieldData(csv_path = raw_csv_path, crs_epsg = "EPSG:4326")

# 3. Compute Descriptive Statistics on Field Observations
    field_stats <- SummariseFieldData(cleaned_df = cleaned_field$data_df)
    print(field_stats)

# 4. Save Outputs
    # Save cleaned tabular and spatial field data for downstream scripts
        write_csv(cleaned_field$data_df, file.path("outputs", "field_data_cleaned.csv"))
        terra::writeVector(
          cleaned_field$points_vect,
          file.path("outputs", "field_plots_points.gpkg"),
          overwrite = TRUE
        )

    # Save structured statistical CSV
        write_csv(
          field_stats,
          file.path("outputs", "stats_field_summary.csv")
        )
