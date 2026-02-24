#library(rgee)
library(sf)
library(tigris)
library(lubridate)
library(dplyr)
library(terra)
library(png)



# --------------------------
# Utilities
# --------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper: convert an sf object to an Earth Engine object using rgee::sf_as_ee
# - This wrapper provides a clear error message if rgee is not installed or not initialized.
# - Returns an ee$FeatureCollection or ee$Geometry depending on the input.
sf_as_ee <- function(sf_obj) {
  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop("The 'rgee' package is required to convert sf objects to Earth Engine geometry. ",
         "Install it with install.packages('rgee') and follow rgee::ee_Authenticate() / rgee::ee_Initialize().")
  }
  # Attempt conversion and provide actionable error message on failure
  tryCatch(
    rgee::sf_as_ee(sf_obj),
    error = function(e) {
      stop("Failed to convert 'sf' object to Earth Engine object. Ensure rgee is installed and initialized:\n",
           "  1) rgee::ee_Authenticate() (if needed)\n",
           "  2) rgee::ee_Initialize(drive = TRUE)\n",
           "Original error: ", e$message)
    }
  )
}

# Get county sf by 5-digit FIPS (accepts numeric or string, e.g., 37085 or "37085").
# Returns list(county_sf, county_ee)
get_county_by_fips <- function(fips) {
  fips <- as.character(fips)
  if (nchar(fips) != 5) stop("Please provide a 5-digit FIPS code (state(2)+county(3)), e.g. '37085'.")
  state_fips <- substr(fips, 1, 2)
  options(tigris_use_cache = TRUE)
  counties_sf <- tigris::counties(state = state_fips, cb = TRUE, year = 2021) %>%
    sf::st_as_sf() %>%
    sf::st_transform(4326)
  county_sf <- counties_sf %>% filter(GEOID == fips)
  if (nrow(county_sf) != 1) stop(sprintf("County with FIPS %s not found.", fips))
  # convert to ee object (uses sf_as_ee wrapper above)
  county_ee <- sf_as_ee(county_sf)
  list(county_sf = county_sf, county_ee = county_ee)
}

mask_s2_sr <- function(image) {
  bnames <- image$bandNames()
  has_QA60 <- bnames$contains("QA60")
  has_SCL  <- bnames$contains("SCL")
  
  qa <- image$select("QA60")
  cloudBitMask <- bitwShiftL(1L, 10L)
  cirrusBitMask <- bitwShiftL(1L, 11L)
  mask_qa <- qa$bitwiseAnd(cloudBitMask)$eq(0)$And( qa$bitwiseAnd(cirrusBitMask)$eq(0) )
  img_qa_masked <- image$updateMask(mask_qa)
  
  scl <- image$select("SCL")
  mask_scl <- scl$neq(3)$And(scl$neq(8))$And(scl$neq(9))$And(scl$neq(10))
  img_scl_masked <- image$updateMask(mask_scl)
  
  masked_image <- ee$Image(
    ee$Algorithms$If(
      has_QA60,
      img_qa_masked,
      ee$Image(ee$Algorithms$If(has_SCL, img_scl_masked, image))
    )
  )
  
  # do NOT call select(...) here; leave band selection to the exporter
  masked_image$divide(10000)
}


# Build median Sentinel-2 SR composite for a date window
get_s2_composite_for_date <- function(target_date, aoi_ee, window_days = 90, max_cloud_pct = 60) {
  target_date <- as.Date(target_date)
  start_date <- as.character(target_date - window_days + 1)
  end_date <- as.character(target_date)
  message(sprintf("Building composite for %s (window %s to %s)", target_date, start_date, end_date))
  # Use the harmonized Sentinel-2 SR collection to avoid deprecated asset warnings
  col <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(start_date, end_date)$
    filterBounds(aoi_ee)$
    filter(ee$Filter$lte("CLOUDY_PIXEL_PERCENTAGE", max_cloud_pct))$
    map(mask_s2_sr)
  median_img <- col$median()$clip(aoi_ee)
  median_img
}

# Usage: generate_output_filename(template, county_fips, center_tag, date_str, tag = "warehouse")
generate_output_filename <- function(template, county_fips, center_tag, date_str, tag = NULL) {
  if (is.null(template) || template == "") return(NULL)
  fname <- template
  
  # sanitize tag for filenames
  if (!is.null(tag) && nzchar(tag)) {
    tag_sanit <- gsub("[^A-Za-z0-9_\\-\\.]", "_", as.character(tag))
  } else {
    tag_sanit <- ""
  }
  
  fname <- gsub("\\{fips\\}", county_fips, fname, fixed = FALSE)
  fname <- gsub("\\{center\\}", center_tag, fname, fixed = FALSE)
  fname <- gsub("\\{date\\}", as.character(date_str), fname, fixed = FALSE)
  # Replace {tag} placeholder (if present) or add nothing
  fname <- gsub("\\{tag\\}", tag_sanit, fname, fixed = FALSE)
  
  # If template doesn't include {tag} but a tag was provided, we don't inject it automatically here;
  # injection happens in export_county_tile for non-template default naming.
  # Ensure .tif extension
  if (!grepl("\\.tif$", fname, ignore.case = TRUE)) fname <- paste0(fname, ".tif")
  fname
}

ee_image_to_geotiff <- function(ee_image, region, out_path, bands = c("B4","B3","B2"),
                                scale = 10, via = "drive") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  
  # 1) Convert region to an ee Geometry (Rectangle) if it's an sf/sfc bbox;
  #    otherwise try to coerce existing ee objects to a geometry.
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
  
  # 2) Determine which of the requested bands actually exist on the image.
  #    Use a client-side getInfo() on bandNames (small metadata call).
  band_names_info <- tryCatch({
    ee_image$bandNames()$getInfo()
  }, error = function(e) {
    warning("Couldn't retrieve band names via getInfo(); proceeding assuming requested bands are present.")
    NULL
  })
  
  if (!is.null(band_names_info)) {
    # band_names_info is an R character vector of available band names.
    present_bands <- intersect(bands, band_names_info)
  } else {
    # If getInfo failed, fall back to attempting to use the requested bands directly.
    present_bands <- bands
  }
  
  if (length(present_bands) == 0) {
    stop("None of the requested bands are present on the image. Available bands: ",
         paste(band_names_info, collapse = ", "))
  }
  
  # 3) Select only the present bands (this avoids passing NULLs/lists to the Python EE API)
  ee_sel <- ee_image$select(present_bands)
  
  message("Exporting to GeoTIFF (via = ", via, ") -> ", out_path, "  (bands: ", paste(present_bands, collapse = ","), ")")
  
  # 4) Export via ee_as_rast (rgee's newer function). ee_as_rast typically returns a SpatRaster
  #    or (in some flows) a path to the written file. Handle both cases robustly.
  rast_obj <- rgee::ee_as_rast(
    image = ee_sel,
    region = region_ee,
    scale = scale,
    via = via
  )
  
  # 5) Normalize result to a terra::SpatRaster and ensure the GeoTIFF exists at out_path.
  #    If ee_as_rast returned a path, use that; otherwise write the SpatRaster to disk.
  if (is.character(rast_obj) && file.exists(rast_obj)) {
    # If a path was returned, prefer using it (could be a temporary file)
    r <- terra::rast(rast_obj)
    # If the returned path differs from desired out_path, write to out_path for consistent naming.
    if (normalizePath(rast_obj, winslash = "/", mustWork = FALSE) != normalizePath(out_path, winslash = "/", mustWork = FALSE)) {
      terra::writeRaster(r, filename = out_path, filetype = "GTiff", overwrite = TRUE)
      r <- terra::rast(out_path)
    }
  } else if (inherits(rast_obj, "SpatRaster")) {
    r <- rast_obj
    # write out to out_path so you have an on-disk GeoTIFF
    terra::writeRaster(r, filename = out_path, filetype = "GTiff", overwrite = TRUE)
    r <- terra::rast(out_path)
  } else {
    # try to coerce other raster-like objects
    r <- tryCatch({
      terra::rast(rast_obj)
    }, error = function(e) {
      stop("ee_as_rast returned an unrecognized object. Please inspect return value.")
    })
    terra::writeRaster(r, filename = out_path, filetype = "GTiff", overwrite = TRUE)
    r <- terra::rast(out_path)
  }
  
  r
}

# Save an RGB PNG from a terra::rast (3-band). Downsamples raster if too large for reasonable PNG size.
save_png_from_raster <- function(r, out_png, r_band = 1, g_band = 2, b_band = 3,
                                 png_max_dim = 1600, stretch = "lin", gamma = 1.0, overwrite = TRUE) {
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  if (terra::nlyr(r) < max(c(r_band, g_band, b_band))) stop("Raster does not contain requested bands.")
  ncol_r <- ncol(r)
  nrow_r <- nrow(r)
  max_dim <- max(ncol_r, nrow_r)
  agg_factor <- max(1, ceiling(max_dim / png_max_dim))
  if (agg_factor > 1) {
    message(sprintf("Downsampling raster by factor %d for PNG output (original %dx%d -> approx %dx%d).", agg_factor, ncol_r, nrow_r, ceiling(ncol_r/agg_factor), ceiling(nrow_r/agg_factor)))
    r_small <- terra::aggregate(r, fact = agg_factor, fun = mean)
  } else {
    r_small <- r
  }
  png(filename = out_png, width = ncol(r_small), height = nrow(r_small), units = "px", bg = "white", res = 72)
  par(mar = c(0,0,0,0))
  terra::plotRGB(r_small, r = r_band, g = g_band, b = b_band, scale = 1, stretch = stretch, gamma = gamma, interpolate = FALSE)
  dev.off()
  message("Saved PNG: ", out_png)
  invisible(out_png)
}

# --------------------------
# create square polygon (size_km x size_km) around center (lon,lat)
# Uses local UTM projected CRS determined from center lon to ensure distance units in meters.
# --------------------------
create_square_around_point <- function(center = NULL, size_km = 5, county_sf = NULL) {
  # center: numeric vector c(lon, lat) in WGS84; if NULL, will use centroid of county_sf (must be provided)
  if (is.null(center) && is.null(county_sf)) stop("Either center or county_sf must be provided to create square.")
  if (is.null(center)) {
    # compute centroid of county (geometry-only)
    cpt <- sf::st_centroid(county_sf$geometry)
    center <- sf::st_coordinates(cpt)[1, ]
  }
  lon <- as.numeric(center[1])
  lat <- as.numeric(center[2])
  # Determine appropriate UTM zone EPSG for northern hemisphere (US)
  zone <- floor((lon + 180) / 6) + 1
  utm_epsg <- 32600 + zone  # EPSG for WGS84 / UTM zone N
  # create point in lon/lat, transform to UTM
  pt <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  pt_utm <- sf::st_transform(pt, crs = utm_epsg)
  half <- (size_km * 1000) / 2
  # coordinates in UTM
  coords <- sf::st_coordinates(pt_utm)[1, ]
  x <- coords[1]; y <- coords[2]
  # square polygon coords (counter-clockwise)
  sq_coords <- matrix(c(
    x - half, y - half,
    x + half, y - half,
    x + half, y + half,
    x - half, y + half,
    x - half, y - half
  ), ncol = 2, byrow = TRUE)
  sq_poly_utm <- sf::st_sfc(sf::st_polygon(list(sq_coords)), crs = utm_epsg)
  sq_poly_wgs84 <- sf::st_transform(sq_poly_utm, crs = 4326)
  # Ensure geometry is valid and return as an sf data.frame for consistent downstream behavior
  sq_poly_wgs84 <- sf::st_make_valid(sq_poly_wgs84)
  sq_sf <- sf::st_sf(geometry = sq_poly_wgs84)
  sf::st_crs(sq_sf) <- 4326
  sq_sf
}

# Clip square to the county polygon (intersection); returns clipped polygon (possibly smaller if tile extends beyond county)
clip_square_to_county <- function(square_sf, county_sf) {
  # both simple features with crs 4326 expected
  sq <- sf::st_make_valid(square_sf)
  ct <- sf::st_make_valid(county_sf)
  inter <- sf::st_intersection(sq, ct)
  if (is.null(inter) || (inherits(inter, "sf") && nrow(inter) == 0)) {
    warning("Square does not intersect the county. Returned object may be empty.")
  }
  inter
}

# --------------------------
# Helper: validate and optionally snap a center to be inside a county
# --------------------------
validate_and_snap_center <- function(center, county_sf, snap_to = c("nearest", "centroid", "none")) {
  # center: numeric c(lon,lat) or NULL
  # county_sf: sf polygon of county (WGS84)
  # snap_to: "nearest" (default) | "centroid" | "none"
  snap_to <- match.arg(snap_to)
  if (is.null(center)) return(NULL)
  if (!is.numeric(center) || length(center) != 2) stop("center must be numeric c(lon, lat).")
  center_pt <- sf::st_sfc(sf::st_point(center), crs = 4326)
  inside <- tryCatch({
    sf::st_within(center_pt, county_sf, sparse = FALSE)[1, 1]
  }, error = function(e) FALSE)
  if (isTRUE(inside)) return(center)  # already inside
  # outside: apply policy
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
    # snap_to == "none": leave as-is and warn
    warning("Provided center lies outside the county and snap_to = 'none'. Proceeding with the provided center; clipped tile may be empty or small.")
    return(center)
  }
}

# --------------------------
# export a single tile (size_km x size_km) inside a county selected by FIPS
# Produces GeoTIFF and PNG for "now" and "one year ago" for that tile area
# Adds snap_to option to control how out-of-county centers are handled.
# snap_to: "nearest" (default) | "centroid" | "none"
# --------------------------
export_county_tile <- function(county_fips,
                               center = NULL,         # c(lon, lat) in WGS84; if NULL -> centroid
                               size_km = 5,
                               target_date = Sys.Date(),
                               window_days = 90,
                               bands = c("B4","B3","B2"),
                               scale = 10,            # meters
                               max_cloud_pct = 60,
                               output_dir = "county_tiles",
                               download_via = "drive",
                               png_max_dim = 1600,
                               snap_to = c("nearest", "centroid", "none")) {
  snap_to <- match.arg(snap_to)
  # get county and prepare county objects early
  county <- get_county_by_fips(county_fips)
  county_sf <- county$county_sf         # ensure county_sf is available for validation and later use
  county_ee <- county$county_ee
  county_name <- county_sf$NAME[1]
  message(sprintf("Creating %dkm x %dkm tile for %s (FIPS %s)", size_km, size_km, county_name, county_fips))
  
  # validate/adjust center (may modify center if outside county)
  center <- validate_and_snap_center(center = center, county_sf = county_sf, snap_to = snap_to)
  
  # build square around final (possibly snapped) center and clip to county
  square <- create_square_around_point(center = center, size_km = size_km, county_sf = county_sf)
  clipped <- clip_square_to_county(square, county_sf)
  if (is.null(clipped) || (inherits(clipped, "sf") && nrow(clipped) == 0)) stop("The requested square does not intersect the county.")
  clipped_wgs84 <- sf::st_transform(clipped, 4326)
  clipped_ee <- sf_as_ee(clipped_wgs84)
  # Build composites clipped to the tile geometry
  img_now <- get_s2_composite_for_date(target_date, clipped_ee, window_days, max_cloud_pct)
  img_prev <- get_s2_composite_for_date(as.Date(target_date) - years(1), clipped_ee, window_days, max_cloud_pct)
  # Output dir and filenames
  base_out <- file.path(output_dir, paste0(county_fips, "_tile"))
  dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
  center_tag <- if (!is.null(center)) sprintf("c_%0.5f_%0.5f", center[1], center[2]) else "centroid"
  now_tif <- file.path(base_out, sprintf("%s_%s_s2_%s.tif", county_fips, center_tag, as.character(target_date)))
  prev_tif <- file.path(base_out, sprintf("%s_%s_s2_%s.tif", county_fips, center_tag, as.character(as.Date(target_date) - years(1))))
  # Export GeoTIFFs for clipped tile region
  r_now <- ee_image_to_geotiff(img_now, clipped_wgs84, out_path = now_tif, bands = bands, scale = scale, via = download_via)
  r_prev <- ee_image_to_geotiff(img_prev, clipped_wgs84, out_path = prev_tif, bands = bands, scale = scale, via = download_via)
  # Save PNGs (downsampled when necessary)
  now_png <- sub("\\.tif$", ".png", now_tif)
  prev_png <- sub("\\.tif$", ".png", prev_tif)
  save_png_from_raster(r_now, now_png, r_band = 1, g_band = 2, b_band = 3, png_max_dim = png_max_dim)
  save_png_from_raster(r_prev, prev_png, r_band = 1, g_band = 2, b_band = 3, png_max_dim = png_max_dim)
  message("Tile exports complete. Files saved to: ", normalizePath(base_out))
  invisible(list(geo_now = now_tif, geo_prev = prev_tif, png_now = now_png, png_prev = prev_png, tile = clipped_wgs84))
}

# --------------------------
# Backward-compatible: export full-county GeoTIFF & PNG (unchanged interface)
# --------------------------
export_county_fullres <- function(county_fips,
                                  target_date = Sys.Date(),
                                  window_days = 90,
                                  bands = c("B4","B3","B2"),
                                  scale = 10,                     # 10 m full Sentinel-2
                                  max_cloud_pct = 60,
                                  output_dir = "county_outputs",
                                  download_via = "drive",
                                  png_max_dim = 1600,
                                  overwrite = FALSE) {
  county <- get_county_by_fips(county_fips)
  county_sf <- county$county_sf
  county_ee <- county$county_ee
  county_name <- county_sf$NAME[1]
  message(sprintf("Preparing exports for %s County (FIPS %s)", county_name, county_fips))
  img_now <- get_s2_composite_for_date(target_date, county_ee, window_days, max_cloud_pct)
  img_prev <- get_s2_composite_for_date(as.Date(target_date) - years(1), county_ee, window_days, max_cloud_pct)
  base_out <- file.path(output_dir, county_fips)
  dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
  now_tif <- file.path(base_out, sprintf("%s_s2_%s.tif", county_fips, as.character(target_date)))
  prev_tif <- file.path(base_out, sprintf("%s_s2_%s.tif", county_fips, as.character(as.Date(target_date) - years(1))))
  r_now <- ee_image_to_geotiff(img_now, county_sf, out_path = now_tif, bands = bands, scale = scale, via = download_via)
  r_prev <- ee_image_to_geotiff(img_prev, county_sf, out_path = prev_tif, bands = bands, scale = scale, via = download_via)
  now_png <- file.path(base_out, sprintf("%s_s2_%s.png", county_fips, as.character(target_date)))
  prev_png <- file.path(base_out, sprintf("%s_s2_%s.png", county_fips, as.character(as.Date(target_date) - years(1))))
  save_png_from_raster(r_now, now_png, r_band = 1, g_band = 2, b_band = 3, png_max_dim = png_max_dim)
  save_png_from_raster(r_prev, prev_png, r_band = 1, g_band = 2, b_band = 3, png_max_dim = png_max_dim)
  message("Exports complete. GeoTIFFs and PNGs saved to: ", normalizePath(base_out))
  invisible(list(geo_now = now_tif, geo_prev = prev_tif, png_now = now_png, png_prev = prev_png))
}