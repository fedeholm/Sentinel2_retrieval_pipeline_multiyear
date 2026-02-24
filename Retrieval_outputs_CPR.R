# OPTIONAL: set RETICULATE_PYTHON before loading reticulate (adjust path to your python)
Sys.setenv(RETICULATE_PYTHON = "YOUR_PYTHON_PATH")
library(reticulate)
reticulate::py_config()   # confirm points to the desired python

# INSTALL ALL NECESSARY PACKAGES

# 1) Initialize rgee & authenticate if needed (interactive)

library(rgee)
# If you haven't authenticated or want to re-authenticate:
rgee::ee_Authenticate()      # interactive browser auth, one-time
# Initialize EE (drive = TRUE if you intend to use Drive)
rgee::ee_Initialize(drive = TRUE)


# Full production example: export high-resolution (10 m) 5 km Sentinel-2 tiles for
# "now" and one year ago for a given center using rgee -> Google Drive export.
# Then download the Drive-exported GeoTIFFs robustly using rclone (example commands below).
#
# Key steps:
#  1) Initialize rgee (interactive authentication if needed).
#  2) Run export_county_tile(..., download_via = "drive", drive_folder = "rgee_backup").
#  3) Use the returned drive_name / drive_folder to locate the file in Drive.
#  4) Use rclone (CLI) to list and copy the file from Drive to local disk (reliable for large files).
#
# Notes:
#  - For very large exports use download_via = "drive" (the export happens to Google Drive).
#  - rgee will create an Earth Engine export task and move the resulting file to your Drive folder.
#  - After the task completes, rgee's wrapper typically attempts to move it to a local path;
#    if rgee doesn't download it automatically or you prefer rclone, you can download it from Drive.
#


# 2) Configure these values before running:

county_fips <- "#####"
center_coords <- c() # Coordinate pair
size_km <- 6           # production tile side length (km)
scale <- 10            # 10 m (Sentinel-2 native resolution)
window_days <- 30      # Time window to allow for low cloud coverage days
png_max_dim <- 1600   
output_dir <- "production_outputs"    # local folder for previews/metadata
drive_folder <- "rgee_backup"         # Drive folder where rgee will place exported files
prev_months <- c(12)                  # Time lag between images, include 1 year ago
bands <- c("B4","B3","B2")            # RGB
name_tag <- "Compass_Statesville "   # <- set your tag here (or NULL to omit)
filename_template <- "{fips}_{center}_{date}_{tag}"  # <- set to NULL to use default naming

# name_tag: a short label to append to output filenames (e.g. "warehouse", "sunglasses")
# filename_template: optional template string. Use tokens {fips}, {center}, {date}, {tag}.
# Example: "{fips}_{center}_{date}_{tag}" -> "37145_c_-78.98973_36.47930_2025-10-31_warehouse.tif"



# 3) Functions: Source your pipeline implementation (this defines export_county_tile and helpers).
#    Make sure these paths point to your actual files (the Core functions must include generate_output_filename).
source("YOUR_LOCAL_PATH/Core_functions.R")
source("YOUR_LOCAL_PATH/Single_file_pipeline.R")

# Ensure your get_s2_composite_for_date is set to harmonized collection (optional override)
get_s2_composite_for_date <- (function() {
  function(target_date, aoi_ee, window_days = 90, max_cloud_pct = 60) {
    target_date <- as.Date(target_date)
    start_date <- as.character(target_date - window_days + 1)
    end_date <- as.character(target_date)
    message(sprintf("Building composite for %s (window %s to %s)", target_date, start_date, end_date))
    col <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
      filterDate(start_date, end_date)$
      filterBounds(aoi_ee)$
      filter(ee$Filter$lte("CLOUDY_PIXEL_PERCENTAGE", max_cloud_pct))$
      map(mask_s2_sr)
    median_img <- col$median()$clip(aoi_ee)
    median_img
  }
})()
assign("get_s2_composite_for_date", get_s2_composite_for_date, envir = .GlobalEnv)


# 4) Run the high-resolution export (Drive). This will create EE tasks and upload TIFF(s) to Drive_folder.
#    Pass filename_template and name_tag so exported filenames include the provided tag.
res_prod <- export_county_tile(
  county_fips = county_fips,         # county by FIPS code
  center = center_coords,            # lon/lat center
  size_km = size_km,                 # tile size (km)
  target_date = Sys.Date(),          # "now"
  window_days = window_days,         # composite window
  prev_months = prev_months,         # include 1-year-ago period
  include_now = TRUE,
  bands = bands,
  scale = scale,                     # 10 m => large files
  max_cloud_pct = 50,
  output_dir = output_dir,
  download_via = "drive",            # export via Google Drive
  drive_folder = drive_folder,       # folder in Drive where rgee will place the exports
  png_max_dim = png_max_dim,
  snap_to = "nearest",
  filename_template = filename_template,  # optional template controlling filename
  name_tag = name_tag                       # custom tag included in filename
)

# 5) Inspect results
print(res_prod)

# End of script.