###############################################################################
# FIGURE 1 (NEW, 2nd revision — Reviewer Comment 2)
# Study area and hot weather season (AMJ) maximum-temperature climatology
#
#   (a) Study domain: Indian state boundaries, all 961 nodes of the original
#       31 x 31 IMD 1-degree grid, and the 283 retained 1 x 1 degree land cells
#   (b) Hot weather season climatology: 1996-2025 cell-wise mean of the locked STEP 06
#       TXx, with an optional outline of cells whose mean TXx >= 40 degC
#       (indicative IMD plains threshold only)
#
# PROVENANCE (read-only inputs; nothing upstream is modified or recomputed)
#   * Grid definition and selection rule: copied verbatim from Script 01
#       LAT_SEQ = 7.5..37.5, LON_SEQ = 67.5..97.5 (31 x 31 = 961 nodes)
#       st_read -> st_make_valid -> st_transform(4326)
#       retained = lengths(st_intersects(grid_point, boundary)) > 0
#       grid_id  = order of which(land_mask, arr.ind = TRUE)
#   * Official 283-cell roster: STEP 06 and STEP 05 parquet outputs
#   * TXx values: STEP 06 parquet (TXx = max valid AMJ daily Tmax, degC)
#   * Raw IMD .grd folders: used ONLY to confirm the 31 x 31 record size
#     (file size = 31 x 31 x n_days x 4 bytes); no temperature is read.
#   * amj_indices.parquet is NOT used (it contains 284 cells).
#
# WHAT IS COMPUTED HERE
#   Only a descriptive 30-year mean of the existing STEP 06 TXx values per
#   cell (1996-2025). No index, threshold, trend or significance is
#   recalculated. No cell is removed: cells with missing TXx years are
#   reported, flagged in the CSV and outlined on the map.
#
# OUTPUTS
#   600-dpi PNG + TIFF, vector PDF
#   CSV: all 961 grid nodes with retained status
#   CSV: cell-wise mean TXx and valid-year counts (283 rows)
#   TXT: QC report (every check with PASS / WARN / FATAL)
###############################################################################

rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 0. PACKAGES (same stack as the existing pipeline / Figure 3)
# =============================================================================

required_packages <- c("data.table", "arrow", "sf", "ggplot2", "patchwork",
                       "grid", "scales", "tools")
for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package not installed: ", p,
         "\nInstall it with install.packages('", p, "') and re-run.", call. = FALSE)
}
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(sf)
  library(ggplot2);    library(patchwork); library(grid)
})

# =============================================================================
# 1. CONFIGURATION (edit paths here only)
# =============================================================================

# Use forward slashes on Windows. Spaces in folder names are fine.
BASE_DIR <- "H:/WMO_IMD_R/WMO_IMD"

# Survey of India state boundary. The script STOPS unless this file reproduces
# exactly the 283 cells (same grid_id and coordinates) used in STEP 05 / 06.
SHAPEFILE_PATH <- file.path(BASE_DIR, "shapefile/SOI/State_District_Subdistrict_PAN INDIA",
                            "State_District_Subdistrict_PAN INDIA/State Boundary/State Boundary.shp")

# Raw IMD yearly binaries (Script 01 naming: <var>/<year>.grd)
RAW_DIR   <- file.path(BASE_DIR, "data/raw/imdR_data")
tmax_path <- function(yr) file.path(RAW_DIR, "tmax", sprintf("%d.grd", yr))
tmin_path <- function(yr) file.path(RAW_DIR, "tmin", sprintf("%d.grd", yr))

# Official locked outputs (Scripts 05 and 06 write to <BASE_DIR>/data)
DATA_DIR    <- file.path(BASE_DIR, "data")
STEP06_FILE <- file.path(DATA_DIR, "AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet")
STEP05_FILE <- file.path(DATA_DIR, "AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet")

OUT_DIR <- file.path(BASE_DIR, "STEP08/12_figures/FIG01_STUDY_DOMAIN")

# ---- grid definition: verbatim from Script 01 --------------------------------
NLAT <- 31L; NLON <- 31L
LAT_SEQ <- seq(7.5, 37.5, by = 1)
LON_SEQ <- seq(67.5, 97.5, by = 1)
EXPECTED_NODES    <- 961L
EXPECTED_RETAINED <- 283L

# ---- STEP 06 structure (from Script 06) -------------------------------------
STEP06_YR_START <- 1976L; STEP06_YR_END <- 2025L
EXPECTED_STEP06_ROWS <- 14150L                    # 283 x 50

# ---- Panel (b) settings ------------------------------------------------------
TXX_START <- 1996L; TXX_END <- 2025L              # primary trend period
N_TXX_YEARS <- TXX_END - TXX_START + 1L           # 30
MIN_VALID_FRAC <- 0.80                            # same 80% rule as STEP 07; FLAG ONLY
MIN_VALID_YEARS <- as.integer(ceiling(MIN_VALID_FRAC * N_TXX_YEARS))   # 24
SHOW_HEAT_OUTLINE <- TRUE
HEAT_THRESHOLD_C  <- 40                           # indicative IMD plains threshold
TXX_PLAUSIBLE <- c(0, 60)                         # degC envelope used in Script 06

COORD_TOL <- 1e-6

# ---- Figure 3 style constants -----------------------------------------------
MAP_XLIM <- c(66.5, 98.5)
MAP_YLIM <- c(6.0, 39.6)
NORTH_X  <- 96.30; NORTH_Y <- 35.40
COL_BOUNDARY <- "grey10"; LW_BOUNDARY <- 0.42
COL_CELL_EDGE <- "grey55"; LW_CELL_EDGE <- 0.08

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

F_PNG   <- file.path(OUT_DIR, "FIG01_STUDY_AREA_AMJ_TXx_600dpi.png")
F_TIFF  <- file.path(OUT_DIR, "FIG01_STUDY_AREA_AMJ_TXx_600dpi.tiff")
F_PDF   <- file.path(OUT_DIR, "FIG01_STUDY_AREA_AMJ_TXx.pdf")
F_GRID  <- file.path(OUT_DIR, "FIG01_GRID_961_NODES_RETAINED_STATUS.csv")
F_TXX   <- file.path(OUT_DIR, "FIG01_CELLWISE_MEAN_TXx_1996_2025.csv")
F_QC    <- file.path(OUT_DIR, "FIG01_QC_REPORT.txt")

# =============================================================================
# 2. QC HELPERS (same PASS / WARN / FATAL convention as the pipeline)
# =============================================================================

.REPORT <- character(0)
say <- function(...) {
  line <- paste0(..., collapse = "")
  cat(line, "\n", sep = ""); .REPORT <<- c(.REPORT, line); invisible(line)
}
head_say <- function(txt) {
  bar <- strrep("-", 78); say(""); say(bar); say(txt); say(bar)
}
CHECKS <- list()
write_report <- function() {
  qc <- if (length(CHECKS)) rbindlist(CHECKS) else data.table()
  out <- c(.REPORT, "", "SUMMARY",
           sprintf("  PASS  : %d", if (nrow(qc)) qc[status == "PASS",  .N] else 0L),
           sprintf("  WARN  : %d", if (nrow(qc)) qc[status == "WARN",  .N] else 0L),
           sprintf("  FATAL : %d", if (nrow(qc)) qc[status == "FATAL", .N] else 0L),
           "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
  writeLines(out, F_QC)
}
check <- function(name, passed, detail, severity = c("FATAL", "WARN")) {
  severity <- match.arg(severity)
  ok <- length(passed) == 1L && isTRUE(passed)
  status <- if (ok) "PASS" else severity
  CHECKS[[length(CHECKS) + 1L]] <<- data.table(check = name, status = status, detail = detail)
  say(sprintf("  %-5s | %-44s | %s", status, name, detail))
  if (!ok && severity == "FATAL") {
    say(""); say("RUN STOPPED: FATAL QC failure in '", name, "'. No figure written.")
    write_report()
    stop("FATAL QC failure: ", name, " -- ", detail,
         "\nSee QC report: ", F_QC, call. = FALSE)
  }
  invisible(ok)
}
need_file <- function(path, label) {
  if (!file.exists(path)) {
    say("MISSING INPUT: ", label, " -> ", path)
    write_report()
    stop("Required input not found (", label, "): ", path,
         "\nEdit the path in Section 1. Nothing was guessed or substituted.", call. = FALSE)
  }
}

head_say("FIGURE 1 | Study area and hot weather season (AMJ) TXx climatology")
say("Run time         : ", format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"))
say("R version        : ", R.version.string)
say("sf / GEOS / GDAL : ", as.character(packageVersion("sf")), " / ",
    sf::sf_extSoftVersion()[["GEOS"]], " / ", sf::sf_extSoftVersion()[["GDAL"]])
say("sf_use_s2()      : ", sf::sf_use_s2(),
    "  (left at the session default, as in Script 01)")
say("Shapefile        : ", SHAPEFILE_PATH)
say("Raw IMD .grd     : ", RAW_DIR)
say("STEP 06          : ", STEP06_FILE)
say("STEP 05          : ", STEP05_FILE)
say("Output dir       : ", OUT_DIR)
say("NOT USED         : amj_indices.parquet (284 cells; not the official roster)")

need_file(SHAPEFILE_PATH, "state boundary shapefile")
need_file(STEP06_FILE,    "STEP 06 parquet")
need_file(STEP05_FILE,    "STEP 05 parquet")

INPUTS <- c(SHAPEFILE_PATH, STEP06_FILE, STEP05_FILE)
.md5_before   <- tools::md5sum(INPUTS)
.mtime_before <- file.mtime(INPUTS)

# =============================================================================
# 3. REBUILD THE ORIGINAL GRID AND SELECTION RULE (verbatim from Script 01)
# =============================================================================
head_say("3. ORIGINAL 31 x 31 IMD GRID AND st_intersects SELECTION (Script 01)")

poly_raw <- st_read(SHAPEFILE_PATH, quiet = TRUE)
say("Boundary features: ", nrow(poly_raw))
check("boundary_has_crs", !is.na(st_crs(poly_raw)),
      sprintf("input CRS: %s", if (is.na(st_crs(poly_raw))) "NA" else st_crs(poly_raw)$Name))

# --- build_land_mask() exactly as in Script 01 ---------------------------------
build_land_mask <- function(shp_path, lon_seq, lat_seq) {
  poly <- st_read(shp_path, quiet = TRUE)
  poly <- st_make_valid(poly)
  # Convert shapefile from LCC_WGS84 to geographic WGS84
  poly <- st_transform(poly, 4326)
  grid_pts <- expand.grid(lon = lon_seq, lat = lat_seq)
  pts_sf <- st_as_sf(grid_pts, coords = c("lon", "lat"), crs = 4326)
  inside <- lengths(st_intersects(pts_sf, poly)) > 0
  mask_mat <- matrix(inside, nrow = length(lat_seq), ncol = length(lon_seq), byrow = TRUE)
  mask_mat  # [lat, lon], TRUE = retain
}
land_mask <- build_land_mask(SHAPEFILE_PATH, LON_SEQ, LAT_SEQ)

# Boundary for plotting: identical processing chain
bnd <- st_transform(st_make_valid(poly_raw), 4326)

# --- all 961 nodes, with the same orientation as the mask ---------------------
nodes <- data.table(expand.grid(lon = LON_SEQ, lat = LAT_SEQ))
nodes[, `:=`(row = match(lat, LAT_SEQ), col = match(lon, LON_SEQ))]
nodes[, inside_boundary := land_mask[cbind(row, col)]]

check("grid_31x31", length(LAT_SEQ) == NLAT && length(LON_SEQ) == NLON,
      sprintf("%d lat x %d lon", length(LAT_SEQ), length(LON_SEQ)))
check("grid_961_nodes", nrow(nodes) == EXPECTED_NODES,
      sprintf("observed %d vs expected %d", nrow(nodes), EXPECTED_NODES))
check("grid_nodes_unique", !anyDuplicated(nodes, by = c("lat", "lon")),
      sprintf("%d duplicated (lat, lon)", sum(duplicated(nodes, by = c("lat", "lon")))))
check("grid_spacing_1deg", all(abs(diff(LAT_SEQ) - 1) < COORD_TOL) &&
        all(abs(diff(LON_SEQ) - 1) < COORD_TOL), "uniform 1-degree spacing")

# --- raw IMD binaries: record size must be 31 x 31 x n_days x 4 bytes --------
is_leap <- function(yr) (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
if (dir.exists(RAW_DIR)) {
  yrs <- STEP06_YR_START:STEP06_YR_END
  raw <- rbindlist(lapply(c("tmax", "tmin"), function(v) {
    f <- if (v == "tmax") tmax_path(yrs) else tmin_path(yrs)
    data.table(var = v, year = yrs, file = f, exists = file.exists(f),
               size = ifelse(file.exists(f), file.size(f), NA_real_),
               expected = NLAT * NLON * ifelse(is_leap(yrs), 366, 365) * 4)
  }))
  n_missing_raw <- raw[exists == FALSE, .N]
  n_bad_size    <- raw[exists == TRUE & size != expected, .N]
  check("raw_grd_files_present", n_missing_raw == 0L,
        sprintf("%d of %d tmax/tmin yearly files missing", n_missing_raw, nrow(raw)),
        severity = "WARN")
  check("raw_grd_record_is_31x31", n_bad_size == 0L,
        sprintf("%d files whose size != 31 x 31 x n_days x 4 bytes (checked %d)",
                n_bad_size, raw[exists == TRUE, .N]))
  if (n_bad_size > 0L) print(raw[exists == TRUE & size != expected])
} else {
  check("raw_grd_dir_found", FALSE,
        sprintf("RAW_DIR not found (%s); 31 x 31 record-size check skipped", RAW_DIR),
        severity = "WARN")
}

n_mask <- sum(land_mask)
check("mask_reproduces_283", n_mask == EXPECTED_RETAINED,
      sprintf("st_intersects mask retains %d nodes vs expected %d", n_mask, EXPECTED_RETAINED))

# --- grid_id exactly as assigned in Script 01 ---------------------------------
idx <- which(land_mask, arr.ind = TRUE)
roster_mask <- data.table(grid_id = as.character(seq_len(nrow(idx))),
                          lat = LAT_SEQ[idx[, 1]], lon = LON_SEQ[idx[, 2]])
nodes <- merge(nodes, roster_mask, by = c("lat", "lon"), all.x = TRUE, sort = FALSE)

# =============================================================================
# 4. OFFICIAL ROSTER: STEP 06 AND STEP 05
# =============================================================================
head_say("4. OFFICIAL 283-CELL ROSTER (STEP 06, STEP 05)")

s6 <- as.data.table(arrow::read_parquet(STEP06_FILE))
say("STEP 06 rows x cols: ", nrow(s6), " x ", ncol(s6))
need6 <- c("grid_id", "year", "lat", "lon", "TXx", "usable_TXx", "n_obs_tmax")
miss6 <- setdiff(need6, names(s6))
if (length(miss6)) {
  say("MISSING STEP 06 COLUMN(S): ", paste(miss6, collapse = ", "))
  say("Columns present: ", paste(names(s6), collapse = ", "))
  write_report()
  stop("STEP 06 schema error. Missing: ", paste(miss6, collapse = ", "),
       "\nTXx cannot be verified, so panel (b) is not drawn. Nothing was substituted.",
       call. = FALSE)
}
check("step06_schema", TRUE, sprintf("required columns present: %s", paste(need6, collapse = ", ")))
s6 <- s6[, ..need6]
s6[, `:=`(grid_id = as.character(grid_id), year = as.integer(year),
          lat = as.numeric(lat), lon = as.numeric(lon), TXx = as.numeric(TXx))]

check("step06_row_count", nrow(s6) == EXPECTED_STEP06_ROWS,
      sprintf("observed %d vs expected %d (%d x %d)", nrow(s6), EXPECTED_STEP06_ROWS,
              EXPECTED_RETAINED, STEP06_YR_END - STEP06_YR_START + 1L))
check("step06_year_coverage", identical(sort(unique(s6$year)), STEP06_YR_START:STEP06_YR_END),
      sprintf("%d years [%d-%d]", uniqueN(s6$year), min(s6$year), max(s6$year)))
check("step06_no_duplicate_grid_year", !anyDuplicated(s6, by = c("grid_id", "year")),
      sprintf("%d duplicated (grid_id, year)", sum(duplicated(s6, by = c("grid_id", "year")))))

roster6 <- unique(s6[, .(grid_id, lat, lon)])
check("step06_283_grids", uniqueN(s6$grid_id) == EXPECTED_RETAINED,
      sprintf("observed %d vs expected %d", uniqueN(s6$grid_id), EXPECTED_RETAINED))
check("step06_one_coordinate_per_grid", nrow(roster6) == uniqueN(s6$grid_id),
      sprintf("%d grid/lat/lon triplets for %d grids", nrow(roster6), uniqueN(s6$grid_id)))
check("step06_no_duplicate_coordinates", !anyDuplicated(roster6, by = c("lat", "lon")),
      sprintf("%d grids share a coordinate", sum(duplicated(roster6, by = c("lat", "lon")))))
check("step06_coords_on_imd_grid",
      all(vapply(roster6$lat, function(v) any(abs(v - LAT_SEQ) < COORD_TOL), logical(1))) &&
        all(vapply(roster6$lon, function(v) any(abs(v - LON_SEQ) < COORD_TOL), logical(1))),
      "every STEP 06 lat/lon is a node of the 31 x 31 grid")

s5 <- as.data.table(arrow::read_parquet(STEP05_FILE))
need5 <- c("grid_id", "lat", "lon")
miss5 <- setdiff(need5, names(s5))
check("step05_schema", length(miss5) == 0L,
      if (length(miss5)) paste("missing:", paste(miss5, collapse = ", ")) else "grid_id, lat, lon present")
roster5 <- unique(s5[, .(grid_id = as.character(grid_id),
                         lat = as.numeric(lat), lon = as.numeric(lon))])
rm(s5)
check("step05_283_grids", uniqueN(roster5$grid_id) == EXPECTED_RETAINED &&
        nrow(roster5) == EXPECTED_RETAINED,
      sprintf("%d grids, %d triplets", uniqueN(roster5$grid_id), nrow(roster5)))

# --- exact coordinate + grid_id matching, three ways --------------------------
cmp_roster <- function(a, b, la, lb) {
  m <- merge(a, b, by = "grid_id", all = TRUE, suffixes = c(".a", ".b"))
  n_only <- m[is.na(lat.a) | is.na(lat.b), .N]
  n_diff <- m[!is.na(lat.a) & !is.na(lat.b) &
                (abs(lat.a - lat.b) > COORD_TOL | abs(lon.a - lon.b) > COORD_TOL), .N]
  check(sprintf("roster_match_%s_vs_%s", la, lb), n_only == 0L && n_diff == 0L,
        sprintf("%d grid_ids in only one source, %d coordinate mismatches", n_only, n_diff))
}
cmp_roster(roster6, roster5, "STEP06", "STEP05")
cmp_roster(roster6, roster_mask, "STEP06", "Script01mask")

nodes[, retained_official := !is.na(grid_id) & grid_id %in% roster6$grid_id]
check("nodes_retained_eq_283", nodes[retained_official == TRUE, .N] == EXPECTED_RETAINED,
      sprintf("%d of 961 nodes retained", nodes[retained_official == TRUE, .N]))
check("retained_iff_inside_boundary", nodes[retained_official != inside_boundary, .N] == 0L,
      sprintf("%d nodes where official status disagrees with st_intersects",
              nodes[retained_official != inside_boundary, .N]))

# =============================================================================
# 5. TXx (STEP 06) — VERIFY, THEN 1996-2025 CELL-WISE MEAN
# =============================================================================
head_say("5. STEP 06 TXx: VERIFICATION AND 1996-2025 MEAN")

say("Units: STEP 06 defines TXx = max of valid AMJ daily Tmax from the IMD")
say("       gridded Tmax (degC). The parquet carries no unit metadata, so the")
say("       unit is verified here by the Script 06 plausibility envelope only.")

check("TXx_usable_flag_matches_NA", s6[usable_TXx != !is.na(TXx), .N] == 0L,
      sprintf("%d rows where usable_TXx disagrees with is.na(TXx)", s6[usable_TXx != !is.na(TXx), .N]))
check("TXx_NA_iff_no_tmax_obs",
      s6[(is.na(TXx) & n_obs_tmax > 0L) | (!is.na(TXx) & n_obs_tmax == 0L), .N] == 0L,
      "TXx is NA exactly where the grid-year has zero valid Tmax days")
rng <- range(s6$TXx, na.rm = TRUE)
check("TXx_within_degC_envelope",
      s6[!is.na(TXx) & (TXx < TXX_PLAUSIBLE[1] | TXx > TXX_PLAUSIBLE[2]), .N] == 0L,
      sprintf("observed range [%.2f, %.2f]; envelope [%g, %g] degC",
              rng[1], rng[2], TXX_PLAUSIBLE[1], TXX_PLAUSIBLE[2]))

tx <- s6[year >= TXX_START & year <= TXX_END]
check("TXx_period_rows", nrow(tx) == EXPECTED_RETAINED * N_TXX_YEARS,
      sprintf("%d rows vs expected %d (%d x %d)", nrow(tx), EXPECTED_RETAINED * N_TXX_YEARS,
              EXPECTED_RETAINED, N_TXX_YEARS))
check("TXx_period_30_years_every_grid", nrow(tx[, .N, by = grid_id][N != N_TXX_YEARS]) == 0L,
      sprintf("%d grids without exactly 30 year-rows", nrow(tx[, .N, by = grid_id][N != N_TXX_YEARS])))

cell <- tx[, .(
  n_years_period = .N,
  n_valid_TXx    = sum(!is.na(TXx)),
  n_missing_TXx  = sum(is.na(TXx)),
  missing_years  = paste(sort(year[is.na(TXx)]), collapse = ";"),
  mean_TXx_C     = if (any(!is.na(TXx))) mean(TXx, na.rm = TRUE) else NA_real_,
  sd_TXx_C       = if (sum(!is.na(TXx)) > 1L) stats::sd(TXx, na.rm = TRUE) else NA_real_,
  min_TXx_C      = if (any(!is.na(TXx))) min(TXx, na.rm = TRUE) else NA_real_,
  max_TXx_C      = if (any(!is.na(TXx))) max(TXx, na.rm = TRUE) else NA_real_
), by = .(grid_id, lat, lon)]
cell[, below_min_valid_years := n_valid_TXx < MIN_VALID_YEARS]
cell[, mean_TXx_ge_threshold := !is.na(mean_TXx_C) & mean_TXx_C >= HEAT_THRESHOLD_C]
setorder(cell, lat, lon)

n_gap_cells <- cell[n_missing_TXx > 0L, .N]
n_low_cells <- cell[below_min_valid_years == TRUE, .N]
n_na_cells  <- cell[is.na(mean_TXx_C), .N]
check("TXx_cells_283", nrow(cell) == EXPECTED_RETAINED,
      sprintf("%d cells summarised", nrow(cell)))
check("TXx_no_missing_years", n_gap_cells == 0L,
      sprintf("%d cells with >= 1 missing TXx year (kept, listed in CSV)", n_gap_cells),
      severity = "WARN")
check("TXx_min_valid_years", n_low_cells == 0L,
      sprintf("%d cells with < %d of %d valid years (kept, outlined dotted)",
              n_low_cells, MIN_VALID_YEARS, N_TXX_YEARS), severity = "WARN")
check("TXx_mean_defined_all_cells", n_na_cells == 0L,
      sprintf("%d cells with no valid TXx year (kept, mapped as empty)", n_na_cells),
      severity = "WARN")
if (n_gap_cells > 0L) {
  say("Cells with missing TXx years (grid_id | lat | lon | n_valid | missing years):")
  for (i in which(cell$n_missing_TXx > 0L))
    say(sprintf("    %-5s | %5.1f | %5.1f | %2d | %s", cell$grid_id[i], cell$lat[i],
                cell$lon[i], cell$n_valid_TXx[i], cell$missing_years[i]))
}

say(sprintf("Mean TXx 1996-2025 across cells: min %.2f | median %.2f | max %.2f degC",
            min(cell$mean_TXx_C, na.rm = TRUE), stats::median(cell$mean_TXx_C, na.rm = TRUE),
            max(cell$mean_TXx_C, na.rm = TRUE)))
n_hot <- cell[mean_TXx_ge_threshold == TRUE, .N]
say(sprintf("Cells with mean TXx >= %g degC: %d of %d (indicative plains threshold)",
            HEAT_THRESHOLD_C, n_hot, nrow(cell)))

# =============================================================================
# 6. FIGURE
# =============================================================================
head_say("6. BUILD FIGURE")

theme_map_axes <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.65),
      panel.grid.major  = element_line(colour = "grey90", linewidth = 0.25),
      axis.text         = element_text(size = base_size - 2, colour = "grey20"),
      axis.text.y       = element_text(margin = margin(r = 3)),
      axis.text.x       = element_text(margin = margin(t = 3)),
      axis.ticks        = element_line(colour = "grey40", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      plot.title        = element_text(face = "bold", size = base_size + 1.5,
                                       hjust = 0, margin = margin(b = 4)),
      legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.25),
      legend.key        = element_rect(fill = "white", colour = NA),
      legend.title      = element_text(face = "bold", size = base_size - 1.5),
      legend.text       = element_text(size = base_size - 2),
      plot.margin       = margin(6, 8, 6, 6)
    )
}

# North arrow: same dart as Figure 3
north_arrow <- function(p, x = NORTH_X, y = NORTH_Y, h = 2.10, w = 0.62, lab_size = 4.4) {
  p +
    annotate("text", x = x, y = y + h + 0.72, label = "N",
             fontface = "bold", size = lab_size, colour = "black") +
    annotate("polygon", x = c(x, x - w, x, x + w), y = c(y + h, y, y + 0.30 * h, y),
             fill = "black", colour = "black", linewidth = 0.30)
}

# Scale bar: degrees of longitude are converted to km at a stated latitude,
# because the map is drawn in geographic coordinates (EPSG:4326).
SCALE_LAT <- 20; SCALE_KM <- 500
scale_bar <- function(p, x0 = 67.6, y0 = 7.0) {
  L <- SCALE_KM / (111.32 * cos(SCALE_LAT * pi / 180))
  p +
    annotate("rect", xmin = x0, xmax = x0 + L / 2, ymin = y0, ymax = y0 + 0.28,
             fill = "black", colour = "black", linewidth = 0.25) +
    annotate("rect", xmin = x0 + L / 2, xmax = x0 + L, ymin = y0, ymax = y0 + 0.28,
             fill = "white", colour = "black", linewidth = 0.25) +
    annotate("text", x = c(x0, x0 + L / 2, x0 + L), y = y0 + 0.75,
             label = c("0", SCALE_KM / 2, paste(SCALE_KM, "km")), size = 2.6) +
    annotate("text", x = x0 + L / 2, y = y0 - 0.45,
             label = sprintf("(scale at %d\u00B0N)", SCALE_LAT), size = 2.2, colour = "grey30")
}

lon_lab <- function(x) paste0(x, "\u00B0E")
lat_lab <- function(x) paste0(x, "\u00B0N")
map_coord <- coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE,
                      crs = 4326, default_crs = 4326, datum = sf::st_crs(4326))
map_scales <- list(
  scale_x_continuous(breaks = seq(70, 95, 5), labels = lon_lab),
  scale_y_continuous(breaks = seq(10, 35, 5), labels = lat_lab)
)

retained <- cell[, .(grid_id, lat, lon)]
LAB_CELL <- sprintf("Retained 1\u00B0 \u00D7 1\u00B0 land cell (n = %d)", EXPECTED_RETAINED)
LAB_NODE <- sprintf("IMD 1\u00B0 grid node (31 \u00D7 31 = %d)", EXPECTED_NODES)

# ---- Panel (a) -----------------------------------------------------------------
pA <- ggplot() +
  geom_tile(data = retained, aes(x = lon, y = lat, fill = LAB_CELL),
            width = 1, height = 1, colour = "#08519C", linewidth = 0.15) +
  geom_sf(data = bnd, inherit.aes = FALSE, fill = NA,
          colour = COL_BOUNDARY, linewidth = LW_BOUNDARY) +
  geom_point(data = nodes, aes(x = lon, y = lat, colour = LAB_NODE),
             size = 0.55, stroke = 0, shape = 16) +
  scale_fill_manual(values = setNames("#9ECAE1", LAB_CELL), name = NULL) +
  scale_colour_manual(values = setNames("grey30", LAB_NODE), name = NULL) +
  guides(fill   = guide_legend(order = 1, override.aes = list(colour = "#08519C")),
         colour = guide_legend(order = 2, override.aes = list(size = 1.6))) +
  map_scales + map_coord + theme_map_axes(10) +
  theme(legend.position = c(0.985, 0.015), legend.justification = c(1, 0),
        legend.margin = margin(1, 4, 1, 3), legend.spacing.y = unit(0, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.25),
        legend.box.margin = margin(2, 2, 2, 2)) +
  labs(title = "(a) Study domain and grid selection", x = NULL, y = NULL) +
  annotate("label", x = 89.6, y = 12.2, hjust = 0.5,
           label = sprintf("n = %d", EXPECTED_RETAINED),
           fontface = "bold", size = 4.2, label.size = 0.3, fill = "white")
pA <- scale_bar(north_arrow(pA))

# ---- Panel (b) -----------------------------------------------------------------
# Colour bar spans the map height; threshold key sits below it, level with
# the bottom of the map (layout requested in the 2nd-revision mark-up).
COLOURBAR_CM  <- 12.2   # adjust if the bar does not reach the map bottom
LEGEND_GAP_CM <- 0.15

br_lo <- floor(min(cell$mean_TXx_C, na.rm = TRUE) / 2) * 2
br_hi <- ceiling(max(cell$mean_TXx_C, na.rm = TRUE) / 2) * 2
TXX_BREAKS <- seq(br_lo, br_hi, by = 2)
TXX_COLS <- scales::brewer_pal(palette = "YlOrRd")(9)

pB <- ggplot() +
  # Cell outline drawn in each cell's own fill colour: no visible grid lines,
  # and no anti-aliasing seams between adjacent cells. Colours are unchanged.
  geom_tile(data = cell[!is.na(mean_TXx_C)],
            aes(x = lon, y = lat, fill = mean_TXx_C, colour = mean_TXx_C),
            width = 1, height = 1, linewidth = 0.25) +
  scale_colour_stepsn(colours = TXX_COLS, breaks = TXX_BREAKS,
                      limits = range(TXX_BREAKS), guide = "none") +
  scale_fill_stepsn(colours = TXX_COLS, breaks = TXX_BREAKS,
                    limits = range(TXX_BREAKS), na.value = "white",
                    name = sprintf("Mean AMJ TXx\n%d\u2013%d (\u00B0C)", TXX_START, TXX_END),
                    guide = guide_coloursteps(order = 1,
                                              title.position = "top", title.hjust = 0,
                                              barheight = unit(COLOURBAR_CM, "cm"),
                                              barwidth  = unit(0.42, "cm"),
                                              frame.colour = "grey20",
                                              show.limits = TRUE))

# cells kept but with too few valid years: white fill + dotted outline (Fig 3 convention)
if (n_low_cells > 0L || n_na_cells > 0L) {
  pB <- pB + geom_tile(data = cell[below_min_valid_years == TRUE | is.na(mean_TXx_C)],
                       aes(x = lon, y = lat), width = 1, height = 1, fill = NA,
                       colour = "grey25", linewidth = 0.38, linetype = "dotted")
}

heat_outline <- NULL
if (isTRUE(SHOW_HEAT_OUTLINE) && n_hot > 0L) {
  # Union of the 1-degree squares in planar lon/lat space (no s2 curvature),
  # then tagged EPSG:4326 for plotting. Pure geometry; no values change.
  hot <- cell[mean_TXx_ge_threshold == TRUE]
  sq <- lapply(seq_len(nrow(hot)), function(i) {
    x <- hot$lon[i]; y <- hot$lat[i]
    st_polygon(list(matrix(c(x - .5, y - .5, x + .5, y - .5, x + .5, y + .5,
                             x - .5, y + .5, x - .5, y - .5), ncol = 2, byrow = TRUE)))
  })
  heat_outline <- st_set_crs(st_union(st_sfc(sq, crs = NA_crs_)), 4326)
  LAB_HEAT <- sprintf("Mean TXx \u2265 %g \u00B0C\n(%d cells; indicative\nIMD plains threshold)",
                      HEAT_THRESHOLD_C, n_hot)
  pB <- pB +
    geom_sf(data = st_sf(lab = LAB_HEAT, geometry = heat_outline),
            aes(linetype = lab), inherit.aes = FALSE, fill = NA,
            colour = "black", linewidth = 0.75) +
    scale_linetype_manual(values = setNames("solid", LAB_HEAT), name = NULL) +
    guides(linetype = guide_legend(order = 2, override.aes = list(fill = NA),
                                   title = "\n", title.theme = element_text(size = 14),
                                   keywidth = unit(0.55, "cm"), keyheight = unit(0.55, "cm")))
}

pB <- pB +
  geom_sf(data = bnd, inherit.aes = FALSE, fill = NA,
          colour = COL_BOUNDARY, linewidth = LW_BOUNDARY) +
  map_scales + map_coord + theme_map_axes(10) +
  theme(legend.position = "right",
        legend.box = "vertical", legend.box.just = "left",
        legend.spacing.y = unit(LEGEND_GAP_CM, "cm"),
        legend.background = element_rect(fill = "white", colour = NA)) +
  labs(title = sprintf("(b) Mean hot weather season TXx (AMJ), %d\u2013%d", TXX_START, TXX_END),
       x = NULL, y = NULL)
pB <- scale_bar(north_arrow(pB))

FOOTNOTE <- paste0(
  "(a) The ", EXPECTED_RETAINED, " retained land cells are those of the 31 \u00D7 31 IMD 1\u00B0 grid ",
  "(", EXPECTED_NODES, " nodes, 7.5\u201337.5\u00B0N, 67.5\u201397.5\u00B0E) whose centres fall within ",
  "the Indian boundary.\n",
  "(b) TXx is the highest daily maximum temperature recorded during each April\u2013June (AMJ) hot weather ",
  "season, from IMD gridded daily maximum temperature; each cell shows its ", TXX_START, "\u2013", TXX_END,
  " average.\n",
  if (isTRUE(SHOW_HEAT_OUTLINE) && n_hot > 0L)
    paste0("The black outline encloses cells with mean TXx \u2265 ", HEAT_THRESHOLD_C, " \u00B0C, the IMD ",
           "heatwave criterion for the plains, shown as an indicative reference only; IMD applies different ",
           "criteria for hilly and coastal regions.\n") else "",
  "Native 1\u00B0 cells, no interpolation. Scale bar valid at ", SCALE_LAT, "\u00B0N.")

FIG1 <- (pA | pB) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title   = "Study area and hot weather season (AMJ) maximum-temperature climatology",
    caption = FOOTNOTE,
    theme = theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 16,
                                     margin = margin(b = 6)),
      plot.caption    = element_text(hjust = 0, size = 8.8, colour = "grey20",
                                     lineheight = 1.15, margin = margin(t = 8)),
      plot.caption.position = "plot",
      plot.margin     = margin(10, 12, 8, 12),
      plot.background = element_rect(fill = "white", colour = "grey35", linewidth = 0.5)))

# =============================================================================
# 7. WRITE OUTPUTS
# =============================================================================
head_say("7. WRITE OUTPUTS")

FIG_W <- 15; FIG_H <- 8.2
ggsave(F_PNG,  FIG1, width = FIG_W, height = FIG_H, units = "in", dpi = 600, bg = "white")
ggsave(F_TIFF, FIG1, width = FIG_W, height = FIG_H, units = "in", dpi = 600, bg = "white",
       compression = "lzw")
has_cairo <- isTRUE(capabilities("cairo"))
check("cairo_pdf_available", has_cairo,
      if (has_cairo) "vector PDF via cairo_pdf (as in Figures 1-6)" else
        "cairo not available: base pdf() used; en-dash / >= glyphs may be replaced by dots in the PDF only",
      severity = "WARN")
pdf_dev <- if (has_cairo) cairo_pdf else "pdf"
ggsave(F_PDF,  FIG1, width = FIG_W, height = FIG_H, units = "in", device = pdf_dev, bg = "white")

grid_out <- nodes[, .(row_lat_index = row, col_lon_index = col, lat, lon,
                      inside_boundary_st_intersects = inside_boundary,
                      retained_official = retained_official,
                      grid_id = fifelse(is.na(grid_id), NA_character_, grid_id))]
setorder(grid_out, row_lat_index, col_lon_index)
fwrite(grid_out, F_GRID)

txx_out <- copy(cell)
txx_out[, `:=`(period = sprintf("%d-%d", TXX_START, TXX_END), units = "degC",
               source = basename(STEP06_FILE),
               threshold_C = HEAT_THRESHOLD_C, min_valid_years_flag = MIN_VALID_YEARS)]
setcolorder(txx_out, c("grid_id", "lat", "lon", "period", "n_years_period",
                       "n_valid_TXx", "n_missing_TXx", "missing_years",
                       "mean_TXx_C", "sd_TXx_C", "min_TXx_C", "max_TXx_C", "units",
                       "below_min_valid_years", "mean_TXx_ge_threshold",
                       "threshold_C", "min_valid_years_flag", "source"))
fwrite(txx_out, F_TXX)

for (f in c(F_PNG, F_TIFF, F_PDF, F_GRID, F_TXX))
  check(paste0("written_", tools::file_ext(f), "_", substr(basename(f), 7, 18)),
        file.exists(f) && file.size(f) > 0,
        sprintf("%s (%.2f MB)", basename(f), file.size(f) / 1024^2))

# read-back
rb_g <- fread(F_GRID); rb_t <- fread(F_TXX)
check("readback_grid_csv", nrow(rb_g) == EXPECTED_NODES &&
        sum(rb_g$retained_official) == EXPECTED_RETAINED,
      sprintf("%d rows, %d retained", nrow(rb_g), sum(rb_g$retained_official)))
check("readback_txx_csv", nrow(rb_t) == EXPECTED_RETAINED &&
        isTRUE(all.equal(rb_t$mean_TXx_C[order(rb_t$lat, rb_t$lon)],
                         cell$mean_TXx_C, tolerance = 1e-9)),
      sprintf("%d rows, mean TXx identical after round-trip", nrow(rb_t)))

# inputs untouched
check("inputs_unmodified",
      identical(unname(tools::md5sum(INPUTS)), unname(.md5_before)) &&
        identical(file.mtime(INPUTS), .mtime_before),
      "shapefile, STEP 06 and STEP 05: md5 and modification time unchanged")

head_say("RUN COMPLETE")
say("Retained cells plotted : ", nrow(retained), " of ", EXPECTED_NODES, " grid nodes")
say("Statistics recomputed  : NONE (descriptive 30-year mean of locked STEP 06 TXx only)")
say("Cells removed          : NONE")
say("Heat outline           : ", if (isTRUE(SHOW_HEAT_OUTLINE))
  sprintf("mean TXx >= %g degC, %d cells (indicative plains threshold)", HEAT_THRESHOLD_C, n_hot)
  else "disabled")
write_report()
cat("\nQC report written to:", F_QC, "\n")
if (interactive()) print(FIG1)   # on-screen preview in RStudio only