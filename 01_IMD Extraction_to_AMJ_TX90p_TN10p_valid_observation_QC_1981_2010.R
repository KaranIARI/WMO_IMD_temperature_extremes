## =========================================================================
## IMD 1x1 deg daily Tmax/Tmin .grd extraction (Mar-Jul, 1976-2025)
## Native binary read + STATE_BOUNDARY.shp masking + grid-wise QC
## =========================================================================

library(sf)
library(abind)

## ---- CONFIG ---------------------------------------------------------------
NLAT <- 31; NLON <- 31                      # IMD 1-deg India domain
LAT_SEQ <- seq(7.5, 37.5, by = 1)           # south -> north
LON_SEQ <- seq(67.5, 97.5, by = 1)          # west  -> east
MISSING_VAL <- 99.9
MISSING_TOL <- 0.01                          # float-compare tolerance
YEARS <- 1976:2025
SHAPEFILE_PATH <- "F:/WMO_IMD_R/WMO_IMD/shapefile/STATE_BOUNDARY.shp"


## File-name templates -- ADJUST to your actual .grd naming convention

tmax_path <- function(yr) sprintf("F:/WMO_IMD_R/WMO_IMD/data/raw/imdR_data/tmax/%d.grd", yr)
tmin_path <- function(yr) sprintf("F:/WMO_IMD_R/WMO_IMD/data/raw/imdR_data/tmin/%d.grd", yr)

is_leap <- function(yr) (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
ndays_year <- function(yr) if (is_leap(yr)) 366 else 365

## ---- 1. Spatial mask (grid-centre-within-boundary) ------------------------

build_land_mask <- function(shp_path, lon_seq, lat_seq) {
  poly <- st_read(shp_path, quiet = TRUE)
  poly <- st_make_valid(poly)
  
  # Convert shapefile from LCC_WGS84 to geographic WGS84
  poly <- st_transform(poly, 4326)
  
  grid_pts <- expand.grid(
    lon = lon_seq,
    lat = lat_seq
  )
  
  pts_sf <- st_as_sf(
    grid_pts,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  inside <- lengths(st_intersects(pts_sf, poly)) > 0
  
  mask_mat <- matrix(
    inside,
    nrow = length(lat_seq),
    ncol = length(lon_seq),
    byrow = TRUE
  )
  
  mask_mat  # [lat, lon], TRUE = retain
}

## ---- 2. Read one native .grd binary year file ------------------------------
## IMD convention: sequential daily records, each NLAT x NLON, 4-byte float,
## little-endian, row order = latitude ascending (south->north) within record.
read_grd_year <- function(path, yr, nlat = NLAT, nlon = NLON) {
  nday <- ndays_year(yr)
  n_expected <- nlat * nlon * nday
  con <- file(path, "rb")
  on.exit(close(con))
  raw <- readBin(con, what = "numeric", size = 4, n = n_expected, endian = "little")
  
  if (length(raw) != n_expected) {
    stop(sprintf("Year %d: read %d values, expected %d (nlat*nlon*nday). Check file/endianness.",
                 yr, length(raw), n_expected))
  }
  
  arr <- array(raw, dim = c(nlon, nlat, nday))   # fastest-varying = lon, per IMD record order
  arr <- aperm(arr, c(2, 1, 3))                  # -> [lat, lon, day]
  arr
}

## ---- 3. Missing-value + spatial masking (grid-wise, structure preserved) --
apply_qc_mask <- function(arr, land_mask, var_name = "T") {
  # 99.9 (within tolerance) -> NA
  arr[abs(arr - MISSING_VAL) < MISSING_TOL] <- NA_real_
  
  # physically implausible values -> NA (flagged, not silently dropped)
  implausible <- arr < -10 | arr > 55
  n_implausible <- sum(implausible, na.rm = TRUE)
  if (n_implausible > 0) {
    warning(sprintf("%s: %d implausible values (< -10 or > 55 C) set to NA", var_name, n_implausible))
    arr[implausible] <- NA_real_
  }
  
  # apply land mask: non-land cells -> NA, grid structure (31x31xday) preserved
  land_idx <- which(!land_mask, arr.ind = TRUE)   # [lat, lon] pairs to blank
  for (k in seq_len(nrow(land_idx))) {
    arr[land_idx[k, 1], land_idx[k, 2], ] <- NA_real_
  }
  arr
}

## ---- 4. Subset to March-July (day-of-year range, leap-year aware) ---------
subset_mar_jul <- function(arr, yr) {
  d1 <- as.Date(sprintf("%d-03-01", yr))
  d2 <- as.Date(sprintf("%d-07-31", yr))
  doy1 <- as.integer(format(d1, "%j"))
  doy2 <- as.integer(format(d2, "%j"))
  arr[, , doy1:doy2, drop = FALSE]
}

## ---- 5. QC summary ----------------------------------------------------------
## ---- 5. QC summary ----------------------------------------------------------

qc_summary <- function(arr, yr, var_name, land_mask) {
  
  n_land_cells <- sum(land_mask)
  
  # Only retained land-grid values
  land_vals <- arr[
    rep(as.vector(land_mask), dim(arr)[3])
  ]
  
  frac_na <- mean(is.na(land_vals))
  
  valid <- land_vals[!is.na(land_vals)]
  
  cat(sprintf(
    "[%s %d] dims=%s | land cells=%d | NA frac (land)=%.3f | range=[%.2f, %.2f]\n",
    var_name,
    yr,
    paste(dim(arr), collapse = "x"),
    n_land_cells,
    frac_na,
    min(valid),
    max(valid)
  ))
}

## ---- 6. Master extraction function (one variable, all years) --------------
extract_variable <- function(path_fn, var_name, land_mask, years = YEARS) {
  out <- vector("list", length(years))
  names(out) <- as.character(years)
  
  for (yr in years) {
    p <- path_fn(yr)
    if (!file.exists(p)) stop(sprintf("Missing file for %s %d: %s", var_name, yr, p))
    
    arr <- read_grd_year(p, yr)
    stopifnot(dim(arr)[1] == NLAT, dim(arr)[2] == NLON, dim(arr)[3] == ndays_year(yr))
    
    arr <- apply_qc_mask(arr, land_mask, var_name)
    arr_sub <- subset_mar_jul(arr, yr)
    
    dimnames(arr_sub) <- list(lat = LAT_SEQ, lon = LON_SEQ,
                              date = as.character(seq(as.Date(sprintf("%d-03-01", yr)),
                                                      as.Date(sprintf("%d-07-31", yr)), by = "day")))
    qc_summary(arr_sub, yr, var_name, land_mask)
    out[[as.character(yr)]] <- arr_sub
  }
  out
}

## ---- 7. RUN -----------------------------------------------------------------
land_mask <- build_land_mask(SHAPEFILE_PATH, LON_SEQ, LAT_SEQ)
n_land <- sum(land_mask)
cat(sprintf("Land grid cells retained: %d\n", n_land))
stopifnot(n_land == 283)
library(data.table)
library(arrow)

OUT <- "F:/WMO_IMD_R/WMO_IMD/data/IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet"

all_data <- rbindlist(lapply(YEARS, function(yr) {
  
  cat("Processing:", yr, "\n")
  
  tx <- apply_qc_mask(
    read_grd_year(tmax_path(yr), yr),
    land_mask,
    "Tmax"
  )
  
  tn <- apply_qc_mask(
    read_grd_year(tmin_path(yr), yr),
    land_mask,
    "Tmin"
  )
  
  tx <- subset_mar_jul(tx, yr)
  tn <- subset_mar_jul(tn, yr)
  
  dates <- seq(
    as.Date(sprintf("%d-03-01", yr)),
    as.Date(sprintf("%d-07-31", yr)),
    by = "day"
  )
  
  idx <- which(land_mask, arr.ind = TRUE)
  
  rbindlist(lapply(seq_len(nrow(idx)), function(i) {
    
    r <- idx[i, 1]
    c <- idx[i, 2]
    
    data.table(
      date = dates,
      grid_id = seq_len(nrow(idx))[i],
      lat = LAT_SEQ[r],
      lon = LON_SEQ[c],
      tmax = as.numeric(tx[r, c, ]),
      tmin = as.numeric(tn[r, c, ])
    )
  }))
}), use.names = TRUE)

setorder(all_data, date, grid_id)

write_parquet(
  all_data,
  OUT,
  compression = "zstd"
)

cat("\n================================\n")
cat("EXTRACTION COMPLETE\n")
cat("Rows:", nrow(all_data), "\n")
cat("Grids:", uniqueN(all_data$grid_id), "\n")
cat("Years:", uniqueN(year(all_data$date)), "\n")
cat("Output:", OUT, "\n")
cat("================================\n")

# ============================================================
# FINAL GRID-WISE QC
# ============================================================

cat("\n========================================\n")
cat("FINAL GRID-WISE QC\n")
cat("========================================\n")

# 1. Expected number of rows
expected_rows <- 50 * 153 * 283

cat("Expected rows :", expected_rows, "\n")
cat("Actual rows   :", nrow(all_data), "\n")

# ------------------------------------------------------------
# 2. Date × grid_id duplicates
# ------------------------------------------------------------

duplicate_n <- all_data[
  ,
  .N,
  by = .(date, grid_id)
][N > 1, .N]

cat("Duplicate date-grid:", duplicate_n, "\n")

# ------------------------------------------------------------
# 3. Check every date has exactly 283 grids
# ------------------------------------------------------------

date_qc <- all_data[
  ,
  .(
    n_grids = uniqueN(grid_id),
    rows = .N
  ),
  by = date
]

cat(
  "Dates with != 283 grids:",
  sum(date_qc$n_grids != 283),
  "\n"
)

# ------------------------------------------------------------
# 4. Check every grid has exactly 153 days per year
# ------------------------------------------------------------

all_data[
  ,
  year := year(date)
]

grid_year_qc <- all_data[
  ,
  .(
    n_days = uniqueN(date)
  ),
  by = .(year, grid_id)
]

cat(
  "Grid-years with != 153 days:",
  sum(grid_year_qc$n_days != 153),
  "\n"
)

# ------------------------------------------------------------
# 5. Missing Tmax / Tmin
# ------------------------------------------------------------

cat(
  "Tmax NA:",
  sum(is.na(all_data$tmax)),
  "\n"
)

cat(
  "Tmin NA:",
  sum(is.na(all_data$tmin)),
  "\n"
)

cat(
  "Both Tmax & Tmin NA:",
  sum(is.na(all_data$tmax) & is.na(all_data$tmin)),
  "\n"
)

# ------------------------------------------------------------
# 6. Grid-wise completeness
# ------------------------------------------------------------

grid_qc <- all_data[
  ,
  .(
    records = .N,
    tmax_valid = sum(!is.na(tmax)),
    tmin_valid = sum(!is.na(tmin)),
    tmax_NA = sum(is.na(tmax)),
    tmin_NA = sum(is.na(tmin))
  ),
  by = grid_id
]

cat(
  "Grids with records != 7650:",
  sum(grid_qc$records != 50 * 153),
  "\n"
)

# ------------------------------------------------------------
# 7. Tmax >= Tmin check
# ------------------------------------------------------------

bad_tmax_tmin <- all_data[
  !is.na(tmax) &
    !is.na(tmin) &
    tmax < tmin
]

cat(
  "Tmax < Tmin:",
  nrow(bad_tmax_tmin),
  "\n"
)

# ------------------------------------------------------------
# 8. Date range
# ------------------------------------------------------------

cat(
  "Date range:",
  as.character(min(all_data$date)),
  "to",
  as.character(max(all_data$date)),
  "\n"
)

# ------------------------------------------------------------
# 9. Grid count
# ------------------------------------------------------------

cat(
  "Unique grids:",
  uniqueN(all_data$grid_id),
  "\n"
)

# ------------------------------------------------------------
# 10. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("QC SUMMARY\n")
cat("========================================\n")

if (
  nrow(all_data) == expected_rows &&
  duplicate_n == 0 &&
  sum(date_qc$n_grids != 283) == 0 &&
  sum(grid_year_qc$n_days != 153) == 0 &&
  uniqueN(all_data$grid_id) == 283
) {
  cat("STRUCTURAL QC: PASS\n")
} else {
  cat("STRUCTURAL QC: FAIL\n")
}

cat("========================================\n")









# ============================================================
# COMPLETE FINAL QC — IMD GRD 283 GRIDS, 1976–2025
# ============================================================

library(data.table)

cat("\n========================================\n")
cat("COMPLETE FINAL DATA QC\n")
cat("========================================\n")

# ------------------------------------------------------------
# 0. Basic dimensions
# ------------------------------------------------------------

cat("\n[0] BASIC STRUCTURE\n")

cat("Rows        :", nrow(all_data), "\n")
cat("Columns     :", ncol(all_data), "\n")
cat("Grids       :", uniqueN(all_data$grid_id), "\n")
cat("Years       :", uniqueN(year(all_data$date)), "\n")
cat("Date range  :", as.character(min(all_data$date)),
    "to", as.character(max(all_data$date)), "\n")

# ------------------------------------------------------------
# 1. Expected row count
# ------------------------------------------------------------

expected_rows <- 50 * 153 * 283

cat("\n[1] ROW COUNT\n")
cat("Expected:", expected_rows, "\n")
cat("Actual  :", nrow(all_data), "\n")

# ------------------------------------------------------------
# 2. Duplicate date-grid
# ------------------------------------------------------------

dup <- all_data[
  ,
  .N,
  by = .(date, grid_id)
][N > 1]

cat("\n[2] DUPLICATES\n")
cat("Duplicate date-grid combinations:", nrow(dup), "\n")

# ------------------------------------------------------------
# 3. Every date must have 283 grids
# ------------------------------------------------------------

date_qc <- all_data[
  ,
  .(
    n_grids = uniqueN(grid_id),
    n_rows = .N
  ),
  by = date
]

bad_dates <- date_qc[n_grids != 283 | n_rows != 283]

cat("\n[3] DATE × GRID COVERAGE\n")
cat("Total dates:", nrow(date_qc), "\n")
cat("Dates with != 283 grids:", nrow(bad_dates), "\n")

if (nrow(bad_dates) > 0) {
  print(bad_dates)
}

# ------------------------------------------------------------
# 4. Every grid-year must have exactly 153 days
# ------------------------------------------------------------

grid_year_qc <- all_data[
  ,
  .(
    n_days = uniqueN(date)
  ),
  by = .(year, grid_id)
]

bad_grid_year <- grid_year_qc[n_days != 153]

cat("\n[4] GRID-YEAR COVERAGE\n")
cat("Expected grid-years:", 50 * 283, "\n")
cat("Actual grid-years  :", nrow(grid_year_qc), "\n")
cat("Grid-years != 153 days:", nrow(bad_grid_year), "\n")

if (nrow(bad_grid_year) > 0) {
  print(bad_grid_year)
}

# ------------------------------------------------------------
# 5. Check actual dates within every year
# ------------------------------------------------------------

date_seq_qc <- all_data[
  ,
  .(
    first_date = min(date),
    last_date = max(date),
    n_dates = uniqueN(date)
  ),
  by = year
]

cat("\n[5] YEARLY DATE RANGE\n")
print(date_seq_qc)

# ------------------------------------------------------------
# 6. Missing Tmax/Tmin — total
# ------------------------------------------------------------

cat("\n[6] TOTAL MISSING VALUES\n")

tmax_na <- sum(is.na(all_data$tmax))
tmin_na <- sum(is.na(all_data$tmin))

cat("Tmax NA:", tmax_na,
    sprintf("(%.4f%%)", 100 * tmax_na / nrow(all_data)), "\n")

cat("Tmin NA:", tmin_na,
    sprintf("(%.4f%%)", 100 * tmin_na / nrow(all_data)), "\n")

cat("Both NA:",
    sum(is.na(all_data$tmax) & is.na(all_data$tmin)),
    "\n")

# ------------------------------------------------------------
# 7. Missing values by YEAR
# ------------------------------------------------------------

year_missing <- all_data[
  ,
  .(
    total = .N,
    tmax_NA = sum(is.na(tmax)),
    tmin_NA = sum(is.na(tmin)),
    both_NA = sum(is.na(tmax) & is.na(tmin))
  ),
  by = year
]

cat("\n[7] MISSING BY YEAR\n")
print(year_missing)

# ------------------------------------------------------------
# 8. Missing values by GRID
# ------------------------------------------------------------

grid_missing <- all_data[
  ,
  .(
    total = .N,
    tmax_NA = sum(is.na(tmax)),
    tmin_NA = sum(is.na(tmin)),
    both_NA = sum(is.na(tmax) & is.na(tmin))
  ),
  by = grid_id
]

cat("\n[8] MISSING BY GRID\n")

cat(
  "Grids with Tmax missing:",
  sum(grid_missing$tmax_NA > 0), "\n"
)

cat(
  "Grids with Tmin missing:",
  sum(grid_missing$tmin_NA > 0), "\n"
)

cat(
  "Maximum Tmax NA in one grid:",
  max(grid_missing$tmax_NA), "\n"
)

cat(
  "Maximum Tmin NA in one grid:",
  max(grid_missing$tmin_NA), "\n"
)

# Worst grids
cat("\nTop 10 grids by Tmax missing:\n")
print(
  grid_missing[
    order(-tmax_NA)
  ][1:min(10, .N)]
)

cat("\nTop 10 grids by Tmin missing:\n")
print(
  grid_missing[
    order(-tmin_NA)
  ][1:min(10, .N)]
)

# ------------------------------------------------------------
# 9. Missing values by GRID × YEAR
# ------------------------------------------------------------

grid_year_missing <- all_data[
  ,
  .(
    tmax_NA = sum(is.na(tmax)),
    tmin_NA = sum(is.na(tmin))
  ),
  by = .(year, grid_id)
]

cat("\n[9] GRID-YEAR MISSINGNESS\n")

cat(
  "Grid-years with Tmax missing:",
  sum(grid_year_missing$tmax_NA > 0),
  "\n"
)

cat(
  "Grid-years with Tmin missing:",
  sum(grid_year_missing$tmin_NA > 0),
  "\n"
)

cat(
  "Maximum Tmax NA in one grid-year:",
  max(grid_year_missing$tmax_NA),
  "\n"
)

cat(
  "Maximum Tmin NA in one grid-year:",
  max(grid_year_missing$tmin_NA),
  "\n"
)

# ------------------------------------------------------------
# 10. Coordinates: one unique lat/lon per grid
# ------------------------------------------------------------

coord_qc <- all_data[
  ,
  .(
    n_lat = uniqueN(lat),
    n_lon = uniqueN(lon)
  ),
  by = grid_id
]

bad_coord <- coord_qc[
  n_lat != 1 | n_lon != 1
]

cat("\n[10] GRID COORDINATES\n")
cat("Grids with changing coordinates:", nrow(bad_coord), "\n")

if (nrow(bad_coord) > 0) {
  print(bad_coord)
}

# ------------------------------------------------------------
# 11. Grid IDs must be exactly 1:283
# ------------------------------------------------------------

expected_grid_ids <- 1:283
actual_grid_ids <- sort(unique(all_data$grid_id))

cat("\n[11] GRID IDs\n")

cat(
  "Missing grid IDs:",
  length(setdiff(expected_grid_ids, actual_grid_ids)),
  "\n"
)

cat(
  "Unexpected grid IDs:",
  length(setdiff(actual_grid_ids, expected_grid_ids)),
  "\n"
)

# ------------------------------------------------------------
# 12. Temperature ranges
# ------------------------------------------------------------

cat("\n[12] TEMPERATURE RANGES\n")

cat(
  "Tmax range:",
  range(all_data$tmax, na.rm = TRUE),
  "\n"
)

cat(
  "Tmin range:",
  range(all_data$tmin, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------
# 13. Tmax >= Tmin
# ------------------------------------------------------------

bad_temp_order <- all_data[
  !is.na(tmax) &
    !is.na(tmin) &
    tmax < tmin
]

cat("\n[13] Tmax < Tmin\n")
cat("Violations:", nrow(bad_temp_order), "\n")

# ------------------------------------------------------------
# 14. Values outside our QC physical range
# ------------------------------------------------------------

bad_tmax_range <- all_data[
  !is.na(tmax) &
    (tmax < -10 | tmax > 55)
]

bad_tmin_range <- all_data[
  !is.na(tmin) &
    (tmin < -10 | tmin > 55)
]

cat("\n[14] PHYSICAL RANGE CHECK\n")
cat("Tmax outside -10 to 55:", nrow(bad_tmax_range), "\n")
cat("Tmin outside -10 to 55:", nrow(bad_tmin_range), "\n")

# ------------------------------------------------------------
# 15. Final structural PASS/FAIL
# ------------------------------------------------------------

structural_pass <- (
  nrow(all_data) == expected_rows &&
    nrow(dup) == 0 &&
    nrow(bad_dates) == 0 &&
    nrow(bad_grid_year) == 0 &&
    uniqueN(all_data$grid_id) == 283 &&
    nrow(bad_coord) == 0 &&
    length(setdiff(expected_grid_ids, actual_grid_ids)) == 0 &&
    length(setdiff(actual_grid_ids, expected_grid_ids)) == 0 &&
    nrow(bad_temp_order) == 0 &&
    nrow(bad_tmax_range) == 0 &&
    nrow(bad_tmin_range) == 0
)

cat("\n========================================\n")
cat("FINAL STRUCTURAL QC:",
    ifelse(structural_pass, "PASS", "FAIL"), "\n")
cat("========================================\n")




# ============================================================
# AMJ 5-DAY WINDOW BUFFER QC — FINAL
#
# Purpose:
# Check data completeness required for TX90p / TN10p
# 5-day moving-window percentile calculation.
#
# Baseline : 1981–2010
# Target   : 1 April – 30 June (AMJ)
# Buffer   : 30 March – 2 July
# Grids    : 283
# ============================================================

library(data.table)

# ============================================================
# 1. CORRECT BUFFER SELECTION
# ============================================================

# IMPORTANT:
# March 30–July 2 is selected separately for EACH year.

buffer <- all_data[
  year(date) %in% 1981:2010 &
    format(date, "%m-%d") >= "03-30" &
    format(date, "%m-%d") <= "07-02"
]

# ============================================================
# 2. BUFFER STRUCTURE
# ============================================================

cat("\n========================================\n")
cat("AMJ BUFFER STRUCTURE QC\n")
cat("========================================\n")

cat("Rows:", nrow(buffer), "\n")
cat("Grids:", uniqueN(buffer$grid_id), "\n")
cat("Years:", uniqueN(year(buffer$date)), "\n")

expected_rows <- 30 * 95 * 283

cat("Expected rows:", expected_rows, "\n")
cat("Actual rows:", nrow(buffer), "\n")

# Every year should contain 95 days × 283 grids
buffer_year <- buffer[
  ,
  .(
    n_days = uniqueN(date),
    n_grids = uniqueN(grid_id),
    rows = .N
  ),
  by = year(date)
]

cat("\nYear-wise structure:\n")
print(buffer_year)

cat(
  "\nYears with != 95 days:",
  sum(buffer_year$n_days != 95),
  "\n"
)

cat(
  "Years with != 283 grids:",
  sum(buffer_year$n_grids != 283),
  "\n"
)

cat(
  "Years with != 26885 rows:",
  sum(buffer_year$rows != 95 * 283),
  "\n"
)

# ============================================================
# 3. DUPLICATE CHECK
# ============================================================

duplicate_buffer <- buffer[
  ,
  .N,
  by = .(date, grid_id)
][N > 1]

cat("\n========================================\n")
cat("DUPLICATE CHECK\n")
cat("========================================\n")

cat(
  "Duplicate date-grid records:",
  nrow(duplicate_buffer),
  "\n"
)

# ============================================================
# 4. OVERALL BUFFER MISSINGNESS
# ============================================================

total <- nrow(buffer)

Tmax_NA <- sum(is.na(buffer$tmax))
Tmin_NA <- sum(is.na(buffer$tmin))

Tmax_pct <- 100 * Tmax_NA / total
Tmin_pct <- 100 * Tmin_NA / total

cat("\n========================================\n")
cat("OVERALL AMJ BUFFER MISSINGNESS\n")
cat("========================================\n")

cat(
  "Tmax missing:",
  Tmax_NA,
  "/",
  total,
  sprintf(" = %.4f%%\n", Tmax_pct)
)

cat(
  "Tmin missing:",
  Tmin_NA,
  "/",
  total,
  sprintf(" = %.4f%%\n", Tmin_pct)
)

cat(
  "Both Tmax & Tmin missing:",
  sum(is.na(buffer$tmax) & is.na(buffer$tmin)),
  "\n"
)

# ============================================================
# 5. YEAR-WISE MISSINGNESS
# ============================================================

year_missing <- buffer[
  ,
  .(
    total = .N,
    Tmax_NA = sum(is.na(tmax)),
    Tmin_NA = sum(is.na(tmin))
  ),
  by = year(date)
]

year_missing[
  ,
  `:=`(
    Tmax_NA_pct = 100 * Tmax_NA / total,
    Tmin_NA_pct = 100 * Tmin_NA / total
  )
]

cat("\n========================================\n")
cat("YEAR-WISE BUFFER MISSINGNESS\n")
cat("========================================\n")

print(year_missing)

# ============================================================
# 6. 5-DAY MOVING WINDOWS
# ============================================================

window_qc <- rbindlist(
  lapply(1981:2010, function(y) {
    
    target_dates <- seq(
      as.IDate(paste0(y, "-04-01")),
      as.IDate(paste0(y, "-06-30")),
      by = "day"
    )
    
    rbindlist(
      lapply(target_dates, function(target_date) {
        
        window_dates <- seq(
          target_date - 2,
          target_date + 2,
          by = "day"
        )
        
        z <- buffer[
          date %in% window_dates
        ]
        
        data.table(
          year = y,
          target_date = target_date,
          
          expected = 5 * 283,
          
          actual = nrow(z),
          
          Tmax_NA = sum(is.na(z$tmax)),
          Tmin_NA = sum(is.na(z$tmin))
        )
      })
    )
  })
)

# ============================================================
# 7. 5-DAY WINDOW STRUCTURAL QC
# ============================================================

cat("\n========================================\n")
cat("5-DAY WINDOW STRUCTURE\n")
cat("========================================\n")

cat(
  "Total windows:",
  nrow(window_qc),
  "\n"
)

cat(
  "Expected windows:",
  30 * 91,
  "\n"
)

cat(
  "Windows with != 1415 records:",
  sum(window_qc$actual != 1415),
  "\n"
)

# ============================================================
# 8. WINDOW MISSINGNESS %
# ============================================================

window_qc[
  ,
  `:=`(
    Tmax_valid = expected - Tmax_NA,
    Tmin_valid = expected - Tmin_NA,
    
    Tmax_NA_pct = 100 * Tmax_NA / expected,
    Tmin_NA_pct = 100 * Tmin_NA / expected
  )
]

# ============================================================
# 9. WORST 5-DAY WINDOWS
# ============================================================

cat("\n========================================\n")
cat("WORST 20 — Tmax\n")
cat("========================================\n")

print(
  window_qc[
    order(-Tmax_NA_pct)
  ][1:min(20, .N)]
)

cat("\n========================================\n")
cat("WORST 20 — Tmin\n")
cat("========================================\n")

print(
  window_qc[
    order(-Tmin_NA_pct)
  ][1:min(20, .N)]
)

# ============================================================
# 10. FINAL 5-DAY WINDOW SUMMARY
# ============================================================

cat("\n========================================\n")
cat("FINAL 5-DAY WINDOW SUMMARY\n")
cat("========================================\n")

cat(
  "Minimum Tmax valid:",
  min(window_qc$Tmax_valid),
  "/ 1415\n"
)

cat(
  "Minimum Tmin valid:",
  min(window_qc$Tmin_valid),
  "/ 1415\n"
)

cat(
  "Maximum Tmax missing:",
  max(window_qc$Tmax_NA),
  "/ 1415\n"
)

cat(
  "Maximum Tmin missing:",
  max(window_qc$Tmin_NA),
  "/ 1415\n"
)

cat(
  "Maximum Tmax missing %:",
  sprintf("%.4f%%", max(window_qc$Tmax_NA_pct)),
  "\n"
)

cat(
  "Maximum Tmin missing %:",
  sprintf("%.4f%%", max(window_qc$Tmin_NA_pct)),
  "\n"
)

# ============================================================
# 11. FINAL STRUCTURAL PASS/FAIL
# ============================================================

structure_pass <-
  nrow(buffer) == 806550 &&
  all(buffer_year$n_days == 95) &&
  all(buffer_year$n_grids == 283) &&
  all(buffer_year$rows == 95 * 283) &&
  nrow(duplicate_buffer) == 0 &&
  nrow(window_qc) == 2730 &&
  all(window_qc$actual == 1415)

cat("\n========================================\n")
cat("FINAL QC STATUS\n")
cat("========================================\n")

if (structure_pass) {
  cat("BUFFER STRUCTURE: PASS\n")
  cat("5-DAY WINDOW STRUCTURE: PASS\n")
} else {
  cat("QC: FAIL — DO NOT PROCEED\n")
}

# ============================================================
# 12. SAVE FINAL QC TABLES
# ============================================================

write_parquet(
  buffer,
  "F:/WMO_IMD_R/WMO_IMD/data/AMJ_buffer_1981_2010_30Mar_02Jul.parquet"
)

write_parquet(
  year_missing,
  "F:/WMO_IMD_R/WMO_IMD/data/AMJ_buffer_year_QC_1981_2010.parquet"
)

write_parquet(
  window_qc,
  "F:/WMO_IMD_R/WMO_IMD/data/AMJ_5day_window_QC_1981_2010.parquet"
)

cat("\nQC TABLES SAVED.\n")

# ============================================================
# END
# ============================================================

# ============================================================
# VALID OBSERVATION QC FOR TX90p / TN10p
# Baseline: 1981–2010
# 5-day window
# 283 grids × 91 AMJ target days
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. Create target calendar days
# ------------------------------------------------------------

target_days <- seq(
  as.IDate("2001-04-01"),
  as.IDate("2001-06-30"),
  by = "day"
)

cat("Target AMJ days:", length(target_days), "\n")

# ------------------------------------------------------------
# 2. Calculate valid observations for every
#    grid × target day
# ------------------------------------------------------------

valid_qc <- rbindlist(
  lapply(target_days, function(td_ref) {
    
    month_day <- format(td_ref, "%m-%d")
    
    rbindlist(
      lapply(1:283, function(g) {
        
        # Collect the same calendar-day ±2 days
        # from every baseline year
        
        vals <- rbindlist(
          lapply(1981:2010, function(y) {
            
            target <- as.IDate(
              paste0(y, "-", month_day)
            )
            
            window_dates <- seq(
              target - 2,
              target + 2,
              by = "day"
            )
            
            buffer[
              grid_id == g &
                date %in% window_dates,
              .(
                tmax = tmax,
                tmin = tmin
              )
            ]
          })
        )
        
        data.table(
          target_day = month_day,
          grid_id = g,
          
          expected = 150,
          
          Tmax_valid = sum(!is.na(vals$tmax)),
          Tmin_valid = sum(!is.na(vals$tmin)),
          
          Tmax_NA = sum(is.na(vals$tmax)),
          Tmin_NA = sum(is.na(vals$tmin))
        )
      })
    )
  })
)

# ------------------------------------------------------------
# 3. Valid percentage
# ------------------------------------------------------------

valid_qc[
  ,
  `:=`(
    Tmax_valid_pct = 100 * Tmax_valid / expected,
    Tmin_valid_pct = 100 * Tmin_valid / expected
  )
]

# ------------------------------------------------------------
# 4. Basic structural check
# ------------------------------------------------------------

cat("\n========================================\n")
cat("VALID OBSERVATION QC\n")
cat("========================================\n")

cat(
  "Total grid-days:",
  nrow(valid_qc),
  "\n"
)

cat(
  "Expected grid-days:",
  91 * 283,
  "\n"
)

cat(
  "Rows with expected != 150:",
  sum(valid_qc$expected != 150),
  "\n"
)

# ------------------------------------------------------------
# 5. Overall distribution
# ------------------------------------------------------------

cat("\n========================================\n")
cat("VALID OBSERVATION DISTRIBUTION\n")
cat("========================================\n")

cat("\nTmax valid observations:\n")
print(
  summary(valid_qc$Tmax_valid)
)

cat("\nTmin valid observations:\n")
print(
  summary(valid_qc$Tmin_valid)
)

# ------------------------------------------------------------
# 6. Percentile distribution of valid counts
# ------------------------------------------------------------

cat("\n========================================\n")
cat("VALID COUNT PERCENTILES\n")
cat("========================================\n")

cat("\nTmax:\n")
print(
  quantile(
    valid_qc$Tmax_valid,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
)

cat("\nTmin:\n")
print(
  quantile(
    valid_qc$Tmin_valid,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 7. Valid percentage distribution
# ------------------------------------------------------------

cat("\n========================================\n")
cat("VALID PERCENTAGE DISTRIBUTION\n")
cat("========================================\n")

cat("\nTmax valid %:\n")
print(
  quantile(
    valid_qc$Tmax_valid_pct,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
)

cat("\nTmin valid %:\n")
print(
  quantile(
    valid_qc$Tmin_valid_pct,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 8. Count grid-days below different completeness levels
# ------------------------------------------------------------

cat("\n========================================\n")
cat("LOW-COMPLETENESS GRID-DAYS\n")
cat("========================================\n")

levels <- c(90, 80, 75, 70, 60, 50)

for (p in levels) {
  
  cat(
    "\nBelow", p, "% valid:\n"
  )
  
  cat(
    "Tmax:",
    sum(valid_qc$Tmax_valid_pct < p),
    "\n"
  )
  
  cat(
    "Tmin:",
    sum(valid_qc$Tmin_valid_pct < p),
    "\n"
  )
}

# ------------------------------------------------------------
# 9. Worst grid-days
# ------------------------------------------------------------

cat("\n========================================\n")
cat("WORST 20 Tmax GRID-DAYS\n")
cat("========================================\n")

print(
  valid_qc[
    order(Tmax_valid_pct)
  ][1:20]
)

cat("\n========================================\n")
cat("WORST 20 Tmin GRID-DAYS\n")
cat("========================================\n")

print(
  valid_qc[
    order(Tmin_valid_pct)
  ][1:20]
)

# ------------------------------------------------------------
# 10. Save QC
# ------------------------------------------------------------

write_parquet(
  valid_qc,
  "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_valid_observation_QC_1981_2010.parquet"
)

cat("\n========================================\n")
cat("VALID OBSERVATION QC COMPLETE\n")
cat("========================================\n")

































































































# ============================================================
# MISSING GRID-YEAR DIAGNOSIS
# ============================================================

# 1. Grid-year missingness table
grid_year_missing <- all_data[
  ,
  .(
    total_days = .N,
    Tmax_NA = sum(is.na(tmax)),
    Tmin_NA = sum(is.na(tmin)),
    Tmax_valid = sum(!is.na(tmax)),
    Tmin_valid = sum(!is.na(tmin))
  ),
  by = .(year, grid_id)
]

# ------------------------------------------------------------
# 2. Add coordinates
# ------------------------------------------------------------

grid_coords <- unique(
  all_data[, .(grid_id, lat, lon)]
)

grid_year_missing <- merge(
  grid_year_missing,
  grid_coords,
  by = "grid_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(grid_year_missing, year, grid_id)

# ------------------------------------------------------------
# 3. COMPLETE SEASON MISSING
# ------------------------------------------------------------

full_missing <- grid_year_missing[
  Tmax_NA == 153 | Tmin_NA == 153
]

cat("\n========================================\n")
cat("FULL SEASON MISSING GRID-YEARS\n")
cat("========================================\n")

cat(
  "Grid-years with complete Tmax/Tmin season missing:",
  nrow(full_missing),
  "\n"
)

print(full_missing)

# ------------------------------------------------------------
# 4. More than 10% missing
# ------------------------------------------------------------

partial_missing <- grid_year_missing[
  Tmax_NA > 15 | Tmin_NA > 15
]

cat("\n========================================\n")
cat(">10% MISSING GRID-YEARS\n")
cat("========================================\n")

cat(
  "Grid-years with >15 missing days:",
  nrow(partial_missing),
  "\n"
)

print(partial_missing)

# ------------------------------------------------------------
# 5. Missingness by grid — percentage
# ------------------------------------------------------------

grid_summary <- all_data[
  ,
  .(
    total = .N,
    Tmax_NA = sum(is.na(tmax)),
    Tmin_NA = sum(is.na(tmin))
  ),
  by = .(grid_id, lat, lon)
]

grid_summary[
  ,
  `:=`(
    Tmax_NA_pct = 100 * Tmax_NA / total,
    Tmin_NA_pct = 100 * Tmin_NA / total
  )
]

setorder(grid_summary, -Tmax_NA_pct)

cat("\n========================================\n")
cat("GRID-WISE MISSINGNESS\n")
cat("========================================\n")

print(grid_summary)

# ------------------------------------------------------------
# 6. Worst grids (>5% missing)
# ------------------------------------------------------------

worst_grids <- grid_summary[
  Tmax_NA_pct > 5 | Tmin_NA_pct > 5
]

cat("\n========================================\n")
cat("GRIDS WITH >5% MISSING\n")
cat("========================================\n")

cat(
  "Number of grids:",
  nrow(worst_grids),
  "\n"
)

print(worst_grids)

# ------------------------------------------------------------
# 7. Missingness by year
# ------------------------------------------------------------

year_summary <- all_data[
  ,
  .(
    total = .N,
    Tmax_NA = sum(is.na(tmax)),
    Tmin_NA = sum(is.na(tmin))
  ),
  by = year
]

year_summary[
  ,
  `:=`(
    Tmax_NA_pct = 100 * Tmax_NA / total,
    Tmin_NA_pct = 100 * Tmin_NA / total
  )
]

cat("\n========================================\n")
cat("YEAR-WISE MISSINGNESS\n")
cat("========================================\n")

print(year_summary)

# ------------------------------------------------------------
# 8. Exact locations of missing values
# ------------------------------------------------------------

missing_records <- all_data[
  is.na(tmax) | is.na(tmin),
  .(
    date,
    year,
    grid_id,
    lat,
    lon,
    tmax,
    tmin
  )
]

cat("\n========================================\n")
cat("TOTAL MISSING RECORDS\n")
cat("========================================\n")

cat(
  "Records with Tmax or Tmin missing:",
  nrow(missing_records),
  "\n"
)

# ------------------------------------------------------------
# 9. Save QC tables
# ------------------------------------------------------------

write_parquet(
  grid_year_missing,
  "F:/WMO_IMD_R/WMO_IMD/data/grid_year_missing_QC.parquet"
)

write_parquet(
  grid_summary,
  "F:/WMO_IMD_R/WMO_IMD/data/grid_missing_QC.parquet"
)

write_parquet(
  missing_records,
  "F:/WMO_IMD_R/WMO_IMD/data/missing_records_QC.parquet"
)

cat("\nQC TABLES SAVED.\n")































tmax_list <- extract_variable(tmax_path, "Tmax", land_mask, YEARS)
tmin_list <- extract_variable(tmin_path, "Tmin", land_mask, YEARS)

## Final cross-check: Tmax >= Tmin wherever both present
check_consistency <- function(tmax_list, tmin_list) {
  bad_total <- 0
  for (yr in names(tmax_list)) {
    diff <- tmax_list[[yr]] - tmin_list[[yr]]
    bad <- sum(diff < 0, na.rm = TRUE)
    bad_total <- bad_total + bad
    if (bad > 0) warning(sprintf("Year %s: %d cells with Tmax < Tmin", yr, bad))
  }
  cat(sprintf("Total Tmax < Tmin violations across all years: %d\n", bad_total))
}
check_consistency(tmax_list, tmin_list)

saveRDS(list(tmax = tmax_list, tmin = tmin_list, land_mask = land_mask,
             lat = LAT_SEQ, lon = LON_SEQ),
        "F:/WMO_IMD_R/WMO_IMD/data/IMD_grd_extracted_MarJul_1976_2025.rds")

cat("Done. Grid-wise Tmax/Tmin (lat x lon x day, per year) saved.\n")




#-----------------------------------------------
# ---- TEST: 1999 only -------------------------------------------------------

test_tmax <- extract_variable(
  tmax_path, "Tmax", land_mask, years = 1999
)

test_tmin <- extract_variable(
  tmin_path, "Tmin", land_mask, years = 1999
)

cat("\n1999 GRD TEST SUCCESSFUL\n")

## ---- Grid 105 validation: direct GRD vs known CSV --------------------------

gid <- 105
lat0 <- 25.5
lon0 <- 81.5

dates_check <- as.Date(c(
  "1999-03-30",
  "1999-03-31",
  "1999-04-01",
  "1999-04-02",
  "1999-04-03"
))

lat_i <- which(LAT_SEQ == lat0)
lon_i <- which(LON_SEQ == lon0)

date_i <- match(
  as.character(dates_check),
  dimnames(test_tmax[[1]])$date
)

check <- data.table(
  date = dates_check,
  lat = lat0,
  lon = lon0,
  GRD_Tmax = test_tmax[[1]][lat_i, lon_i, date_i],
  GRD_Tmin = test_tmin[[1]][lat_i, lon_i, date_i]
)

print(check)

# 1999 Tmax GRD: find known CSV value 33.09

p <- tmax_path(1999)

con <- file(p, "rb")
x <- readBin(
  con,
  what = "numeric",
  n = 31 * 31 * 366,
  size = 4,
  endian = "little"
)
close(con)

# 31x31 matrix for 31 March 1999
m31 <- matrix(
  x[((90 - 1) * 31 * 31 + 1):(90 * 31 * 31)],
  nrow = 31,
  ncol = 31,
  byrow = TRUE
)

# Find value corresponding to CSV = 33.09
pos <- which(
  abs(m31 - 33.09) < 0.01,
  arr.ind = TRUE
)

cat("Positions matching 33.09:\n")
print(pos)

if (nrow(pos) > 0) {
  print(
    data.table(
      row = pos[,1],
      col = pos[,2],
      value = m31[pos]
    )
  )
}



# 1999 Tmax - 31 March
p <- tmax_path(1999)

con <- file(p, "rb")
seek(con, where = (90 - 1) * 31 * 31 * 4, origin = "start")

v <- readBin(
  con,
  what = "numeric",
  n = 31 * 31,
  size = 4,
  endian = "little"
)
close(con)






## ---- Check GRD spatial orientation using land mask ------------------------

valid <- abs(m - 99.9) >= 0.01

scores <- c(
  original = sum(valid == land_mask),
  flip_lat = sum(valid == land_mask[31:1, ]),
  flip_lon = sum(valid == land_mask[, 31:1]),
  flip_both = sum(valid == land_mask[31:1, 31:1])
)

print(scores)

cat(
  "\nBest orientation:",
  names(which.max(scores)),
  "\n"
)

m <- matrix(v, nrow = 31, ncol = 31, byrow = TRUE)

cat("Original:\n")
print(m[1:3, 1:3])

cat("\nFlip latitude:\n")
print(m[31:29, 1:3])

cat("\nFlip longitude:\n")
print(m[1:3, 31:29])

cat("\nFlip both:\n")
print(m[31:29, 31:29])



cat("Latitude:", min(LAT_SEQ), "to", max(LAT_SEQ), "\n")
cat("Longitude:", min(LON_SEQ), "to", max(LON_SEQ), "\n")
cat("Grid spacing:", diff(LAT_SEQ)[1], "x", diff(LON_SEQ)[1], "degree\n")
cat("Total cells:", length(LAT_SEQ) * length(LON_SEQ), "\n")

# ============================================================
# VERIFY ORIGINAL GRD MISSINGNESS
# 1990 - Grid 36 - Tmax
# ============================================================

yr <- 1990
grid_id_test <- 36

# Grid 36 coordinates
coord <- all_data[
  grid_id == grid_id_test,
  .(grid_id, lat, lon)
][1]

print(coord)

# Grid position in our 31 x 31 IMD grid
row_id <- which(LAT_SEQ == coord$lat)
col_id <- which(LON_SEQ == coord$lon)

cat("Row:", row_id, "\n")
cat("Column:", col_id, "\n")

# Read ORIGINAL 1990 Tmax GRD
p <- tmax_path(yr)

con <- file(p, "rb")

raw <- readBin(
  con,
  what = "numeric",
  n = 31 * 31 * ndays_year(yr),
  size = 4,
  endian = "little"
)

close(con)

# Convert to [lat, lon, day]
arr <- array(
  raw,
  dim = c(NLON, NLAT, ndays_year(yr))
)

arr <- aperm(
  arr,
  c(2, 1, 3)
)

# March 1 - July 31
d1 <- as.Date("1990-03-01")
d2 <- as.Date("1990-07-31")

doy1 <- as.integer(format(d1, "%j"))
doy2 <- as.integer(format(d2, "%j"))

x <- arr[
  row_id,
  col_id,
  doy1:doy2
]

# Raw values BEFORE applying our QC
cat("\n========================================\n")
cat("ORIGINAL GRD RAW CHECK\n")
cat("========================================\n")

cat("Year:", yr, "\n")
cat("Grid:", grid_id_test, "\n")
cat("Lat:", coord$lat, "\n")
cat("Lon:", coord$lon, "\n")
cat("Days:", length(x), "\n")

cat("99.9 values:", sum(abs(x - 99.9) < 0.01), "\n")
cat("NA values:", sum(is.na(x)), "\n")

cat(
  "Non-99.9 values:",
  sum(abs(x - 99.9) >= 0.01, na.rm = TRUE),
  "\n"
)

cat(
  "Raw range excluding 99.9:",
  range(x[abs(x - 99.9) >= 0.01], na.rm = TRUE),
  "\n"
)

cat("\nFirst 20 RAW GRD values:\n")
print(x[1:20])
# ============================================================
# VERIFY ORIGINAL GRD MISSINGNESS
# 1990 - Grid 36 - Tmin
# ============================================================

yr <- 1990
grid_id_test <- 36

# Grid 36 coordinates
coord <- all_data[
  grid_id == grid_id_test,
  .(grid_id, lat, lon)
][1]

print(coord)

# Grid position
row_id <- which(LAT_SEQ == coord$lat)
col_id <- which(LON_SEQ == coord$lon)

cat("Row:", row_id, "\n")
cat("Column:", col_id, "\n")

# Read ORIGINAL 1990 Tmin GRD
p <- tmin_path(yr)

con <- file(p, "rb")

raw <- readBin(
  con,
  what = "numeric",
  n = 31 * 31 * ndays_year(yr),
  size = 4,
  endian = "little"
)

close(con)

# Convert to [lat, lon, day]
arr <- array(
  raw,
  dim = c(NLON, NLAT, ndays_year(yr))
)

arr <- aperm(
  arr,
  c(2, 1, 3)
)

# March 1 - July 31
d1 <- as.Date("1990-03-01")
d2 <- as.Date("1990-07-31")

doy1 <- as.integer(format(d1, "%j"))
doy2 <- as.integer(format(d2, "%j"))

x <- arr[
  row_id,
  col_id,
  doy1:doy2
]

# ------------------------------------------------------------
# RAW GRD CHECK — BEFORE OUR QC
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ORIGINAL GRD RAW CHECK — Tmin\n")
cat("========================================\n")

cat("Year:", yr, "\n")
cat("Grid:", grid_id_test, "\n")
cat("Lat:", coord$lat, "\n")
cat("Lon:", coord$lon, "\n")
cat("Days:", length(x), "\n")

cat(
  "99.9 values:",
  sum(abs(x - 99.9) < 0.01),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(x)),
  "\n"
)

cat(
  "Non-99.9 values:",
  sum(abs(x - 99.9) >= 0.01, na.rm = TRUE),
  "\n"
)

cat("\nFirst 20 RAW GRD values:\n")
print(x[1:20])








