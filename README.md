# Sentinel2_retrieval_pipeline_multiyear
This repository contains R scripts to build Sentinel-2 composites (now and historical periods) for county tiles and export GeoTIFFs/PNGs via Earth Engine (rgee). It includes helpers to create square tiles around a coordinate or county centroid, optional snapping behavior when a requested center falls outside the county, and Drive/GCS download helpers.

Contents
- `Core_functions.R` — Core pipeline: utilities, EE compositing, raster export, PNG preview creation, square geometry helpers, and `export_county_tile()` (with `snap_to` option).
- `Single_file_pipeline_CPR.R` — Drive-aware helpers and an alternate `export_county_tile()` that supports multiple previous periods (`prev_months`) and robust retrieval helpers (googledrive + rclone).
- `Retrieval_outputs_CPR.R` — Wrapper that runs `export_county_tile()` for now + previous period.


Quick overview
- Build a small square tile (size_km × size_km) inside a county (selected by FIPS).
- Composite Sentinel-2 SR images over a time window ending on a target date (median composite).
- Export GeoTIFFs via Earth Engine; optionally use Google Drive or Google Cloud Storage as the transfer medium.
- Produce downsampled PNG previews and an optional comparison montage.
- Flexible options:
  - `scale` — export pixel size (meters). Use `scale = 10` for Sentinel-2 native resolution.
  - `window_days` — number of days in the composite window (default 90).
  - `prev_months` (drive-aware) — vector of month offsets for previous comparators (e.g., `c(18, 24)`).
  - `snap_to` — how to handle a provided `center` that lies outside the county: `"nearest"` (default), `"centroid"`, or `"none"`.

Prerequisites

1. R >= 4.x (tested on modern R)
2. Required R packages:
   - sf, tigris, lubridate, dplyr, terra, png
   - rgee (Earth Engine wrapper) and reticulate (Python interop)
   - (optional) googledrive for in-R Drive downloads
   - (optional) magick for higher-quality montages
3. Python environment with Earth Engine API and required dependencies (commonly installed in a conda env used by reticulate).
4. (Optional) rclone (external CLI) for robust Drive downloads.
5. Google account with Earth Engine access. For Drive exports, the same account should be available to rgee and (if used) googledrive or rclone.

Setup — R / rgee / Python

1. Create or identify a Python environment that has Earth Engine API, numpy, etc. Example (conda):
   - `conda create -n rgee_py311 python=3.11`
   - `conda activate rgee_py311`
   - `pip install earthengine-api numpy`

2. In R, point `reticulate` to that Python env (set before loading reticulate):
```r
Sys.setenv(RETICULATE_PYTHON = "###YOUR_PATH/python.exe")
library(reticulate)
reticulate::py_config()
```

3. Authenticate & initialize rgee (interactive one-time):
```r
library(rgee)
rgee::ee_Authenticate()       # if required (browser flow)
rgee::ee_Initialize(drive = TRUE)  # enable Drive exports
```

4. Install other R packages if needed:
```r
install.packages(c("sf","tigris","lubridate","dplyr","terra","png"))
# optional:
install.packages("googledrive")
install.packages("magick")
```

Load and source the pipeline

1. Start a fresh R session (recommended).
2. Set `RETICULATE_PYTHON` as above and initialize rgee.
3. Source the core functions first, then the Drive-aware file:
```r
source("Core_functions.R", local = globalenv(), echo = TRUE)
source("Single_file_pipeline_CPR.R", local = globalenv(), echo = TRUE)
```


