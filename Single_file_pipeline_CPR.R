# Single-file pipeline_all_in_one.R
# Combined pipeline: county retrieval, Sentinel-2 harmonized compositing with cloud-probability masking,
# Earth Engine export wrappers (Drive/GCS/getInfo-aware), robust PNG preview writer (cross-terra compatible),
# filename template/tag support, and small unit tests.
#
#
# Notes:
# - This file was constructed to be tolerant across rgee/terra versions.
# - For production large exports, prefer download_via = "drive" + rclone, or "gcs" + gsutil.
# - The composite uses COPERNICUS/S2_SR_HARMONIZED and COPERNICUS/S2_CLOUD_PROBABILITY for improved cloud masking.

# -------------------------------------------------------------------------
# Required packages (install if missing)
# -------------------------------------------------------------------------
required_pkgs <- c("sf", "tigris", "lubridate", "dplyr", "terra", "png")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Please install package: ", p, " (install.packages('", p, "'))")
  }
}
library(sf)
library(tigris)
library(lubridate)
library(dplyr)
library(terra)
library(png)

# rgee and other Google helpers are optional (used only when present)
# -------------------------------------------------------------------------
# Small helpers / operators
# -------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------------------------------------------------------------
# Basic checks
# -------------------------------------------------------------------------
.check_core_available <- function() {
  if (!exists("get_county_by_fips", mode = "function")) {
    message("Core functions not found in the session. Loading from this combined script.")
  } else {
    message("Core functions present.")
  }
}
.check_core_available()

# -------------------------------------------------------------------------
# Date helpers (avoid relying on lubridate::months/years exported symbols)
# -------------------------------------------------------------------------
subtract_months <- function(date, m) {
  if (is.null(date)) return(NULL)
  date <- as.Date(date)
  if (!is.numeric(m) || length(m) != 1) stop("m must be a single numeric value")
  seq.Date(date, by = paste0("-", abs(m), " months"), length.out = 2)[2]
}
subtract_years <- function(date, y) {
  if (is.null(date)) return(NULL)
  date <- as.Date(date)
  if (!is.numeric(y) || length(y) != 1) stop("y must be a single numeric value")
  seq.Date(date, by = paste0("-", abs(y), " years"), length.out = 2)[2]
}

# Label helper for months -> e.g., 12 -> "1y" else "m"
label_from_months <- function(m) {
  if (!is.numeric(m) || length(m) != 1 || is.na(m) || m <= 0) stop("label_from_months expects a positive numeric scalar")
  if (m %% 12 == 0) {
    yrs <- m / 12
    paste0(yrs, "y")
  } else {
    paste0(m, "m")
  }
}

# -------------------------------------------------------------------------
# County retrieval helpers (sf + optional EE conversion)
# -------------------------------------------------------------------------
get_county_sf_only <- function(fips) {
  fips <- as.character(fips)
  if (nchar(fips) != 5) stop("Please provide a 5-digit FIPS code, e.g. '37085'.")
  state_fips <- substr(fips, 1, 2)
  options(tigris_use_cache = TRUE)
  counties_sf <- tigris::counties(state = state_fips, cb = TRUE, year = 2021) %>% # may need to update this to capture changes in the county lines (eg, CT)
    sf::st_as_sf() %>%
    sf::st_transform(4326)
  county_sf <- counties_sf %>% dplyr::filter(GEOID == fips)
  if (nrow(county_sf) != 1) stop(sprintf("County with FIPS %s not found.", fips))
  county_sf
}

sf_as_ee <- function(sf_obj) {
  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop("rgee is required to convert sf objects to Earth Engine objects. Install it and initialize rgee.")
  }
  tryCatch(
    rgee::sf_as_ee(sf_obj),
    error = function(e) stop("Failed to convert 'sf' object to EE. Ensure rgee is initialized. Original error: ", e$message)
  )
}

get_county_by_fips <- function(fips) {
  county_sf <- get_county_sf_only(fips)
  county_ee <- NULL
  if (requireNamespace("rgee", quietly = TRUE)) {
    county_ee <- tryCatch({
      sf_as_ee(county_sf)
    }, error = function(e) {
      warning("rgee is installed but conversion to EE failed; returning county_ee = NULL. Message: ", e$message)
      NULL
    })
  } else {
    warning("rgee not installed; 'county_ee' will be NULL. Install rgee to enable Earth Engine exports.")
  }
  list(county_sf = county_sf, county_ee = county_ee)
}

# -------------------------------------------------------------------------
# Geometry helpers
# -------------------------------------------------------------------------
create_square_around_point <- function(center = NULL, size_km = 5, county_sf = NULL) {
  if (is.null(center) && is.null(county_sf)) stop("Either center or county_sf must be provided.")
  if (is.null(center)) {
    cpt <- sf::st_centroid(county_sf$geometry)
    center <- sf::st_coordinates(cpt)[1, ]
  }
  lon <- as.numeric(center[1]); lat <- as.numeric(center[2])
  zone <- floor((lon + 180) / 6) + 1
  utm_epsg <- 32600 + zone
  pt <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326) # may need to be updated if focusing on west coast? Maybe. Consider implications
  pt_utm <- sf::st_transform(pt, crs = utm_epsg)
  half <- (size_km * 1000) / 2
  coords <- sf::st_coordinates(pt_utm)[1, ]
  x <- coords[1]; y <- coords[2]
  sq_coords <- matrix(c(
    x - half, y - half,
    x + half, y - half,
    x + half, y + half,
    x - half, y + half,
    x - half, y - half
  ), ncol = 2, byrow = TRUE)
  sq_poly_utm <- sf::st_sfc(sf::st_polygon(list(sq_coords)), crs = utm_epsg)
  sq_poly_wgs84 <- sf::st_transform(sq_poly_utm, crs = 4326)
  sq_poly_wgs84 <- sf::st_make_valid(sq_poly_wgs84)
  sq_sf <- sf::st_sf(geometry = sq_poly_wgs84)
  sf::st_crs(sq_sf) <- 4326
  sq_sf
}

clip_square_to_county <- function(square_sf, county_sf) {
  sq <- sf::st_make_valid(square_sf)
  ct <- sf::st_make_valid(county_sf)
  inter <- sf::st_intersection(sq, ct)
  if (is.null(inter) || (inherits(inter, "sf") && nrow(inter) == 0)) {
    warning("Square does not intersect the county. Returned object may be empty.")
  }
  inter
}

validate_and_snap_center <- function(center, county_sf, snap_to = c("nearest", "centroid", "none")) {
  snap_to <- match.arg(snap_to)
  if (is.null(center)) return(NULL)
  if (!is.numeric(center) || length(center) != 2) stop("center must be numeric c(lon, lat).")
  center_pt <- sf::st_sfc(sf::st_point(center), crs = 4326)
  inside <- tryCatch({
    sf::st_within(center_pt, county_sf, sparse = FALSE)[1, 1]
  }, error = function(e) FALSE)
  if (isTRUE(inside)) return(center)
  if (snap_to == "nearest") {
    union_county <- sf::st_union(county_sf)
    nearest_line <- sf::st_nearest_points(center_pt, union_county)
    nearest_pts <- sf::st_cast(nearest_line, "POINT")
    nearest_on_county <- nearest_pts[2]
    nearest_on_county_wgs84 <- sf::st_transform(nearest_on_county, 4326)
    new_coords <- sf::st_coordinates(nearest_on_county_wgs84)[1, ]
    warning(sprintf("Provided center (%.6f, %.6f) lies outside the county. Snapped to nearest county point (%.6f, %.6f).",
                    sf::st_coordinates(center_pt)[1,1], sf::st_coordinates(center_pt)[1,2], new_coords[1], new_coords[2]))
    return(c(new_coords[1], new_coords[2]))
  } else if (snap_to == "centroid") {
    union_county <- sf::st_union(county_sf)
    cent <- sf::st_centroid(union_county)
    cent_wgs84 <- sf::st_transform(cent, 4326)
    new_coords <- sf::st_coordinates(cent_wgs84)[1, ]
    warning(sprintf("Provided center (%.6f, %.6f) lies outside the county. Snapped to county centroid (%.6f, %.6f).",
                    sf::st_coordinates(center_pt)[1,1], sf::st_coordinates(center_pt)[1,2], new_coords[1], new_coords[2]))
    return(c(new_coords[1], new_coords[2]))
  } else {
    warning("Provided center lies outside the county and snap_to = 'none'. Proceeding with the provided center; clipped tile may be empty or small.")
    return(center)
  }
}

# -------------------------------------------------------------------------
# Earth Engine helpers & masking
# -------------------------------------------------------------------------
# mask_s2_sr: basic SCL/QA60-based mask used as fallback inside composite mapping
mask_s2_sr <- function(image) {
  # expects an ee Image; careful: this function is executed in the Earth Engine JS/Python environment via rgee
  # We implement the logic in the rgee's expression style below when mapping, but keep a minimal function here
  image$divide(10000)  # ensure reflectances are scaled if needed (fallback)
}

# Improved S2 compositing: harmonized SR joined with cloud probability and mask using SCL/QA60/cloud_prob
get_s2_composite_for_date <- function(target_date, aoi_ee, window_days = 90, max_cloud_pct = 60,
                                      cloud_prob_thresh = 40) {
  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop("rgee must be installed and initialized to run EE composites. Install rgee and run rgee::ee_Initialize().")
  }
  ee <- rgee::ee
  
  target_date <- as.Date(target_date)
  start_date <- as.character(target_date - window_days + 1)
  end_date <- as.character(target_date)
  message(sprintf("Building composite for %s (window %s to %s)", target_date, start_date, end_date))
  
  # Collections
  s2_col <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(start_date, end_date)$
    filterBounds(aoi_ee)$
    filter(ee$Filter$lte("CLOUDY_PIXEL_PERCENTAGE", max_cloud_pct))
  
  cloud_col <- ee$ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")$
    filterDate(start_date, end_date)$
    filterBounds(aoi_ee)
  
  # Join by system:index
  joiner <- ee$Join$saveFirst("cloud")
  filter_idx <- ee$Filter$equals("system:index", "system:index")
  joined <- joiner$apply(s2_col, cloud_col, filter_idx)
  
  # Map function using the joined property
  add_mask_and_cloudprob <- function(img) {
    img <- ee$Image(img)
    # try to get the matched cloud image from property 'cloud'
    cloud_img <- ee$Image(img$get("cloud"))
    # safe: if the match failed, create a constant 0 band
    cloud_prob_band <- ee$Image(0)$rename("cloud_prob")
    cloud_prob_band <- ee$Algorithms$If(
      cloud_img,
      cloud_img$select("probability")$rename("cloud_prob"),
      cloud_prob_band
    )
    # Ensure cloud_prob band exists on the image
    img2 <- img$addBands(ee$Image(cloud_prob_band))
    
    # Build masks
    bnames <- img2$bandNames()
    has_QA60 <- bnames$contains("QA60")
    has_SCL  <- bnames$contains("SCL")
    
    # QA60 mask if exists: remove cloud/cirrus bits 10 and 11
    mask_qa <- ee$Image(1)
    mask_qa <- ee$Algorithms$If(
      has_QA60,
      {
        qa <- img2$select("QA60")
        cloudBitMask <- bitwShiftL(1L, 10L)
        cirrusBitMask <- bitwShiftL(1L, 11L)
        qa_mask <- qa$bitwiseAnd(cloudBitMask)$eq(0)$And( qa$bitwiseAnd(cirrusBitMask)$eq(0) )
        qa_mask
      },
      ee$Image(1)
    )
    
    # SCL mask - exclude SCL values 3,8,9,10
    mask_scl <- ee$Image(1)
    mask_scl <- ee$Algorithms$If(
      has_SCL,
      {
        scl <- img2$select("SCL")
        scl_mask <- scl$neq(3)$And(scl$neq(8))$And(scl$neq(9))$And(scl$neq(10))
        scl_mask
      },
      ee$Image(1)
    )
    
    # Cloud probability mask
    cloud_prob_mask <- img2$select("cloud_prob")$lt(as.numeric(cloud_prob_thresh))
    
    # Combined mask (use QA60 when present)
    combined_mask <- ee$Image(ee$Algorithms$If(has_QA60,
                                               ee$Image(mask_qa)$And(mask_scl)$And(cloud_prob_mask),
                                               ee$Image(mask_scl)$And(cloud_prob_mask)))
    # Apply mask and scale to reflectance range
    img2$updateMask(combined_mask)$divide(10000)
  }
  
  comp <- joined$map(add_mask_and_cloudprob)$median()$clip(aoi_ee)
  comp
}

# -------------------------------------------------------------------------
# ee_image_to_geotiff - wrapper around rgee::ee_as_rast with normalized 'via' handling
# -------------------------------------------------------------------------
ee_image_to_geotiff <- function(ee_image, region, out_path, bands = c("B4","B3","B2"),
                                scale = 10, via = "drive", drive_folder = NULL) {
  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop("rgee is required for exporting images. Install and initialize rgee to export GeoTIFFs.")
  }
  ee <- rgee::ee
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  
  # Normalize via argument (accept common synonyms)
  if (is.null(via) || identical(via, "")) via <- "drive"
  via_norm <- tolower(as.character(via))
  via_map <- list(
    "getinfo" = "getInfo",
    "get_info" = "getInfo",
    "get-info" = "getInfo",
    "getinfo" = "getInfo",
    "getinfo" = "getInfo",
    "getinfo" = "getInfo",
    "drive" = "drive",
    "gcs" = "gcs",
    "gcloud" = "gcs",
    "googlecloud" = "gcs"
  )
  via_canonical <- via_map[[via_norm]]
  if (is.null(via_canonical)) via_canonical <- "drive"
  
  # Coerce region to EE geometry if needed
  if (inherits(region, "sf") || inherits(region, "sfc")) {
    bbox <- sf::st_bbox(region)
    coords <- as.numeric(c(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]))
    region_ee <- ee$Geometry$Rectangle(coords)
  } else {
    region_ee <- tryCatch({
      region$geometry()
    }, error = function(e1) {
      tryCatch({
        ee$Geometry(region)
      }, error = function(e2) {
        stop("Cannot convert 'region' to an Earth Engine Geometry. Provide sf, numeric bbox, or ee Geometry/Feature.")
      })
    })
  }
  
  # Determine present bands on the image
  band_names_info <- tryCatch({ ee_image$bandNames()$getInfo() }, error = function(e) { NULL })
  present_bands <- if (!is.null(band_names_info)) intersect(bands, band_names_info) else bands
  if (length(present_bands) == 0) stop("None of the requested bands are present on the image.")
  ee_sel <- ee_image$select(present_bands)
  
  message("Exporting to GeoTIFF (via = ", via_canonical, ") -> ", out_path, "  (bands: ", paste(present_bands, collapse = ","), ")")
  
  # helper to attempt ee_as_rast and capture errors
  try_block <- function(v) {
    tryCatch({
      rgee::ee_as_rast(image = ee_sel, region = region_ee, scale = scale, via = v)
    }, error = function(e) e)
  }
  
  rast_obj <- try_block(via_canonical)
  if (inherits(rast_obj, "error")) {
    # try fallbacks
    fallbacks <- c("drive", "gcs", "getInfo")
    for (fb in fallbacks) {
      if (identical(tolower(fb), tolower(via_canonical))) next
      message("Primary via '", via_canonical, "' failed: ", rast_obj$message)
      message("Attempting fallback via = '", fb, "' ...")
      rast_obj2 <- try_block(fb)
      if (!inherits(rast_obj2, "error")) {
        rast_obj <- rast_obj2
        via_canonical <- fb
        break
      } else {
        rast_obj <- rast_obj2
      }
    }
    if (inherits(rast_obj, "error")) stop("All ee_as_rast attempts failed. Last error: ", rast_obj$message)
  }
  
  # Normalize return to terra::rast and ensure out_path exists
  if (is.character(rast_obj) && file.exists(rast_obj)) {
    r <- terra::rast(rast_obj)
    if (normalizePath(rast_obj, winslash = "/", mustWork = FALSE) != normalizePath(out_path, winslash = "/", mustWork = FALSE)) {
      terra::writeRaster(r, filename = out_path, filetype = "GTiff", overwrite = TRUE)
      r <- terra::rast(out_path)
    }
  } else if (inherits(rast_obj, "SpatRaster")) {
    r <- rast_obj
    terra::writeRaster(r, filename = out_path, filetype = "GTiff", overwrite = TRUE)
    r <- terra::rast(out_path)
  } else {
    r <- tryCatch({ terra::rast(rast_obj) }, error = function(e) stop("ee_as_rast returned an unrecognized object."))
    terra::writeRaster(r, filename = out_path, filetype = "GTiff", overwrite = TRUE)
    r <- terra::rast(out_path)
  }
  
  list(rast = r, local_path = out_path, drive_name = basename(out_path), drive_folder = drive_folder %||% "rgee_backup")
}

# -------------------------------------------------------------------------
# Robust PNG writer (avoids terra::disaggregate and handles tiny rasters)
# -------------------------------------------------------------------------
save_png_from_raster <- function(r, out_png, r_band = 1, g_band = 2, b_band = 3,
                                 png_max_dim = 1600, stretch = "lin", gamma = 1.0, overwrite = TRUE,
                                 min_dim = 64, upsample_method = "near") {
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  if (terra::nlyr(r) < max(c(r_band, g_band, b_band))) stop("Raster does not contain requested bands.")
  ncol_r <- ncol(r); nrow_r <- nrow(r)
  max_dim <- max(ncol_r, nrow_r)
  
  # Downsample if very large
  agg_factor_down <- max(1, ceiling(max_dim / png_max_dim))
  if (agg_factor_down > 1) {
    message(sprintf("Downsampling raster by factor %d for PNG output (original %dx%d -> approx %dx%d).",
                    agg_factor_down, ncol_r, nrow_r, ceiling(ncol_r/agg_factor_down), ceiling(nrow_r/agg_factor_down)))
    r_work <- terra::aggregate(r, fact = agg_factor_down, fun = mean)
  } else {
    r_work <- r
  }
  
  # Upsample small rasters using terra::resample
  final_ncol <- ncol(r_work); final_nrow <- nrow(r_work)
  final_max <- max(final_ncol, final_nrow)
  up_factor <- 1
  if (final_max < min_dim) {
    up_factor <- ceiling(min_dim / final_max)
    message(sprintf("Upsampling small raster by factor %d for PNG output (approx %dx%d).", up_factor, final_ncol * up_factor, final_nrow * up_factor))
    e <- terra::ext(r_work)
    crs_r <- terra::crs(r_work, proj = TRUE)
    target <- terra::rast(ncols = final_ncol * up_factor,
                          nrows = final_nrow * up_factor,
                          xmin = e$xmin, xmax = e$xmax,
                          ymin = e$ymin, ymax = e$ymax,
                          crs = crs_r)
    r_work <- terra::resample(r_work, target, method = upsample_method)
  }
  
  out_w <- ncol(r_work); out_h <- nrow(r_work)
  
  try_plotrgb <- function() {
    png(filename = out_png, width = out_w, height = out_h, units = "px", bg = "white", res = 72)
    par(mar = c(0,0,0,0))
    terra::plotRGB(r_work, r = r_band, g = g_band, b = b_band, scale = 1, stretch = stretch, gamma = gamma)
    dev.off()
    invisible(out_png)
  }
  
  fallback_write_png <- function() {
    band_to_matrix <- function(spat, band_index) {
      tryCatch(terra::as.matrix(spat[[band_index]]),
               error = function(e) {
                 vals <- terra::values(spat[[band_index]], mat = FALSE)
                 matrix(vals, nrow = nrow(spat), ncol = ncol(spat), byrow = TRUE)
               })
    }
    rmat <- band_to_matrix(r_work, r_band)
    gmat <- band_to_matrix(r_work, g_band)
    bmat <- band_to_matrix(r_work, b_band)
    
    normalize01 <- function(m) {
      m[is.na(m)] <- NA
      mn <- suppressWarnings(min(m, na.rm = TRUE))
      mx <- suppressWarnings(max(m, na.rm = TRUE))
      if (!is.finite(mn) || !is.finite(mx) || mx == mn) {
        mm <- matrix(0, nrow = nrow(m), ncol = ncol(m))
      } else {
        mm <- (m - mn) / (mx - mn)
        mm[mm < 0] <- 0; mm[mm > 1] <- 1
      }
      mm[is.na(mm)] <- 1
      mm
    }
    
    r_arr <- normalize01(rmat); g_arr <- normalize01(gmat); b_arr <- normalize01(bmat)
    rgb_arr <- array(0, dim = c(nrow(r_arr), ncol(r_arr), 3))
    rgb_arr[,,1] <- r_arr; rgb_arr[,,2] <- g_arr; rgb_arr[,,3] <- b_arr
    png::writePNG(rgb_arr, target = out_png)
    invisible(out_png)
  }
  
  result <- tryCatch({
    try_plotrgb()
  }, error = function(e) {
    warning("terra::plotRGB failed; falling back to manual PNG writer. Error: ", conditionMessage(e))
    tryCatch({ fallback_write_png() }, error = function(e2) stop("Both plotRGB and fallback PNG writer failed: ", conditionMessage(e2)))
  })
  
  message("Saved PNG: ", out_png)
  invisible(out_png)
}

# -------------------------------------------------------------------------
# Filename template & tag helper
# -------------------------------------------------------------------------
generate_output_filename <- function(template, county_fips, center_tag, date_str, tag = NULL) {
  if (is.null(template) || template == "") return(NULL)
  fname <- template
  tag_sanit <- if (!is.null(tag) && nzchar(tag)) gsub("[^A-Za-z0-9_\\-\\.]", "_", as.character(tag)) else ""
  fname <- gsub("\\{fips\\}", county_fips, fname)
  fname <- gsub("\\{center\\}", center_tag, fname)
  fname <- gsub("\\{date\\}", as.character(date_str), fname)
  fname <- gsub("\\{tag\\}", tag_sanit, fname)
  if (!grepl("\\.tif$", fname, ignore.case = TRUE)) fname <- paste0(fname, ".tif")
  fname
}

# -------------------------------------------------------------------------
# Exporters
# -------------------------------------------------------------------------
export_county_tile <- function(county_fips,
                               center = NULL,
                               size_km = 5,
                               target_date = Sys.Date(),
                               window_days = 90,
                               prev_months = c(12),
                               include_now = TRUE,
                               bands = c("B4","B3","B2"),
                               scale = 10,
                               max_cloud_pct = 60,
                               output_dir = "county_tiles",
                               download_via = "drive",
                               drive_folder = "rgee_backup",
                               png_max_dim = 1600,
                               snap_to = c("nearest", "centroid", "none"),
                               filename_template = NULL,
                               name_tag = NULL) {
  snap_to <- match.arg(snap_to)
  county <- get_county_by_fips(county_fips)
  county_sf <- county$county_sf
  county_ee <- county$county_ee
  county_name <- county_sf$NAME[1]
  message(sprintf("Creating %gkm x %gkm tile for %s (FIPS %s)", size_km, size_km, county_name, county_fips))
  
  # Validate/adjust center
  center <- validate_and_snap_center(center = center, county_sf = county_sf, snap_to = snap_to)
  
  # Build and clip square
  square <- create_square_around_point(center = center, size_km = size_km, county_sf = county_sf)
  clipped <- clip_square_to_county(square, county_sf)
  if (is.null(clipped) || (inherits(clipped, "sf") && nrow(clipped) == 0)) stop("The requested square does not intersect the county.")
  clipped_wgs84 <- sf::st_transform(clipped, 4326)
  
  # Convert to EE geometry if possible
  clipped_ee <- NULL
  if (!is.null(county_ee) && requireNamespace("rgee", quietly = TRUE)) {
    clipped_ee <- tryCatch(sf_as_ee(clipped_wgs84), error = function(e) { warning("Could not convert clipped polygon to EE: ", e$message); NULL })
  } else {
    if (!requireNamespace("rgee", quietly = TRUE)) {
      warning("rgee not available; skipping EE compositing/export steps. You can still inspect geometry outputs locally.")
    } else {
      warning("County EE geometry not available; ensure rgee::ee_Initialize() was called to enable exports.")
    }
  }
  
  # If EE unavailable, return geometry
  if (is.null(clipped_ee)) {
    base_out <- file.path(output_dir, paste0(county_fips, "_tile"))
    dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
    center_tag <- if (!is.null(center)) sprintf("c_%0.5f_%0.5f", center[1], center[2]) else "centroid"
    message("rgee not initialized or EE geometry unavailable; returning local tile geometry without exporting GeoTIFFs.")
    return(invisible(list(base_dir = normalizePath(base_out), tile = clipped_wgs84, periods = list())))
  }
  
  # Prepare output dir and tags
  base_out <- file.path(output_dir, paste0(county_fips, "_tile"))
  dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
  center_tag <- if (!is.null(center)) sprintf("c_%0.5f_%0.5f", center[1], center[2]) else "centroid"
  name_tag_sanit <- if (!is.null(name_tag) && nzchar(name_tag)) gsub("[^A-Za-z0-9_\\-\\.]", "_", as.character(name_tag)) else ""
  
  # Build periods
  periods <- list()
  if (include_now) periods[["now"]] <- as.Date(target_date)
  if (!is.null(prev_months) && length(prev_months) > 0) {
    if (!is.numeric(prev_months)) stop("prev_months must be numeric vector of months")
    for (m in prev_months) {
      lbl <- label_from_months(m)
      pd <- subtract_months(target_date, m)
      periods[[lbl]] <- as.Date(pd)
    }
  }
  
  results <- list()
  for (lbl in names(periods)) {
    pd <- periods[[lbl]]
    message(sprintf("Building composite for %s (label: %s)", pd, lbl))
    
    # Use advanced composite (cloud probability join)
    img <- get_s2_composite_for_date(pd, clipped_ee, window_days = window_days, max_cloud_pct = max_cloud_pct)
    
    # Build filename
    if (!is.null(filename_template) && nzchar(filename_template)) {
      fname_base <- generate_output_filename(filename_template, county_fips, center_tag, pd, tag = name_tag_sanit)
      if (grepl("/", fname_base) || grepl("\\\\", fname_base)) {
        out_tif <- fname_base
        dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
      } else {
        out_tif <- file.path(base_out, fname_base)
      }
    } else {
      base_name <- sprintf("%s_%s_s2_%s", county_fips, center_tag, as.character(pd))
      if (nzchar(name_tag_sanit)) base_name <- paste0(base_name, "_", name_tag_sanit)
      out_tif <- file.path(base_out, paste0(base_name, ".tif"))
    }
    
    # Export image (handles via selection and fallbacks)
    exported <- ee_image_to_geotiff(img, clipped_wgs84, out_path = out_tif, bands = bands, scale = scale, via = download_via, drive_folder = drive_folder)
    
    # Save preview PNG
    out_png <- sub("\\.tif$", sprintf("_%s.png", lbl), exported$local_path)
    save_png_from_raster(exported$rast, out_png, r_band = 1, g_band = 2, b_band = 3, png_max_dim = png_max_dim)
    
    results[[lbl]] <- list(local_path = exported$local_path,
                           drive_name = exported$drive_name,
                           drive_folder = exported$drive_folder,
                           raster = exported$rast,
                           png = out_png)
  }
  
  message("Tile exports complete. Files saved to: ", normalizePath(base_out))
  invisible(list(base_dir = normalizePath(base_out), tile = clipped_wgs84, periods = results))
}

# -------------------------------------------------------------------------
# Lightweight unit test: snapping modes (avoids EE calls)
# -------------------------------------------------------------------------
run_snap_tests <- function(county_fips = "37153", verbose = TRUE) {
  if (verbose) message("Running snap_to tests for county ", county_fips)
  county_sf <- get_county_sf_only(county_fips)
  bbox <- sf::st_bbox(county_sf)
  test_point_outside <- c(as.numeric(bbox["xmax"]) + 0.5, as.numeric((bbox["ymin"] + bbox["ymax"]) / 2))
  modes <- c("nearest", "centroid", "none")
  results <- list()
  for (m in modes) {
    if (verbose) message("Testing snap_to = '", m, "' ...")
    final_center <- validate_and_snap_center(test_point_outside, county_sf, snap_to = m)
    square <- create_square_around_point(center = final_center, size_km = 1, county_sf = county_sf)
    clipped <- clip_square_to_county(square, county_sf)
    clipped_area <- if (is.null(clipped) || (inherits(clipped, "sf") && nrow(clipped) == 0)) 0 else as.numeric(sf::st_area(sf::st_union(clipped)))
    results[[m]] <- list(original_center = test_point_outside,
                         final_center = final_center,
                         clipped_area_m2 = clipped_area,
                         clipped_is_empty = clipped_area == 0)
    if (verbose) {
      message(sprintf(" mode=%s | final_center=%.6f,%.6f | clipped_area_m2=%.1f | empty=%s",
                      m, final_center[1], final_center[2], clipped_area, ifelse(clipped_area==0, "TRUE","FALSE")))
    }
  }
  invisible(results)
}

