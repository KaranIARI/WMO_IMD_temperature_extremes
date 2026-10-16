###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 08 : FINAL RESULTS, SYNTHESIS, TABLES, FIGURES AND QC
##
##  283 IMD 1 deg land grids | AMJ (01 Apr - 30 Jun)
##  Indices: TXx, TXn, TNx, TNn, TX90p, TN10p, WSDI, DTR
##
##  v1.1.0 CHANGE LOG (this revision)
##    ONLY the official India state boundary shapefile was integrated into the
##    map-rendering block. No statistical calculation, input dataset, trend
##    computation, HR_lag3 setting, FDR_Q value, colour scale, background, axis,
##    label, stipple, layout, figure dimension, font, title or panel arrangement
##    was altered. Sections 08.1, 08.2, 08.4-08.12, the final QC block and the
##    manifest are byte-for-byte unchanged from v1.0.0 apart from the version
##    string and the single new QC entry inside 08.3.
##
##    The boundary is drawn with geom_path() over coordinates extracted by
##    sf::st_coordinates(), NOT with geom_sf(). geom_sf() forces coord_sf(),
##    which would replace the existing coord_fixed(ratio = 1.05, ...) and
##    silently change the map aspect ratio and panel proportions. sf still
##    performs all reading, CRS detection, validity checking and transformation;
##    only the final draw call avoids coord_sf so the locked layout survives.
##
## ===========================================================================
##  STEP 08 PERFORMS NO STATISTICAL RE-ESTIMATION.
##  The trend engine of STEP 07 / 07A / 07B is LOCKED and is neither modified
##  nor re-run here. STEP 08 reads finished results and produces the tables,
##  figures, QC and manifest for the manuscript.
##
##  LOCKED DESIGN (inherited, never altered in this script)
##    Primary period      1996-2025      Context period   1976-2025
##    Primary variant     HR_lag3        Sensitivity      HR_all, MK_original
##    Completeness        >= 80%         Slope            Theil-Sen per decade
##    Uncertainty         95% CI from the Hamed-Rao corrected variance
##    Multiplicity        Benjamini-Hochberg, q <= 0.10
##    CF <= 0             test_undefined (CF is never floored)
##    Significance and magnitude are reported on separate axes.
##    No strong/weak classification, no magnitude threshold.
##
## ===========================================================================
##  INPUT DISCOVERY - WHAT EXISTS AND WHAT DOES NOT
## ===========================================================================
##  PRESENT in the pipeline:
##    AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv   2,264 rows = 283 x 8
##    AMJ_STEP07_trends_all_variants.parquet        10,188 = 283 x 6 x 2 x 3
##    AMJ_STEP07B_trends_all_variants.parquet        3,396 = 283 x 2 x 2 x 3
##    AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet   14,150 grid-years
##    AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet        14,150 grid-years
##    IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet         daily source
##    STATE_BOUNDARY.shp                                      official boundary
##
##  ABSENT - handled explicitly, never invented:
##
##  (a) NO MONTHLY INDEX FILE EXISTS. Every index in STEP 05/06 is a SEASONAL
##      AMJ value (one number per grid-year over the continuous 91-day season).
##      Sections 08.1 and 08.9 therefore derive monthly CLIMATOLOGY from the
##      daily source, and only where that is methodologically legitimate:
##        TXx, TXn, TNx, TNn, DTR -> derivable (plain max/min/mean per month),
##            written as *_monthly_climatology, baseline-free, never used for
##            any trend statement.
##        WSDI -> BLOCKED. A monthly WSDI requires resetting the spell counter
##            at 30 Apr / 31 May, contradicting the locked STEP 06 rule of one
##            continuous 91-day sequence with no month reset. The script
##            refuses and writes a BLOCKED report.
##        TX90p / TN10p -> computable only as DIRECT exceedance percentages.
##            The Zhang (2005) in-base bootstrap was applied at the seasonal
##            level and re-running it monthly is forbidden here. Computed under
##            the distinct names TX90p_direct_monthly / TN10p_direct_monthly
##            with a WARN, and barred from any trend context.
##            Controlled by MONTHLY_TX90P_MODE.
##
##  (b) NO AUTHORITATIVE REGIONAL CLASSIFICATION EXISTS. Section 08.10 is
##      PENDING by default and writes a requirements report. Supply REGION_FILE
##      (CSV with columns grid_id, region) to activate it. NOTE: the state
##      boundary shapefile is a CARTOGRAPHIC OVERLAY only and does NOT
##      constitute a regional classification for 08.10.
##
##  (c) Indices carry THREE different units (degC, %, days per decade). A single
##      shared colour scale across all eight would be dimensionally invalid, so
##      combined maps are produced PER UNIT FAMILY. Eight individual maps are
##      produced as well.
##
##  Run:  RUN_SECTION = "ALL"  or one of "08.1" ... "08.12"
###############################################################################

## ===========================================================================
## SECTION 0 - CONFIGURATION   (the ONLY block you should need to edit)
## ===========================================================================
rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE, scipen = 999)
options(error = NULL)
Sys.setenv(TZ = "UTC")

SCRIPT_NAME    <- "08_final_results_and_figures.R"
SCRIPT_VERSION <- "1.1.0"

## ---- paths -----------------------------------------------------------------
INPUT_DIR  <- "F:/WMO_IMD_R/WMO_IMD/data"
OUTPUT_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08"

TREND_FILE      <- file.path(INPUT_DIR, "AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv")
VARIANTS_FILE_A <- file.path(INPUT_DIR, "AMJ_STEP07_trends_all_variants.parquet")
VARIANTS_FILE_B <- file.path(INPUT_DIR, "AMJ_STEP07B_trends_all_variants.parquet")
STEP06_FILE     <- file.path(INPUT_DIR, "AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet")
STEP05_FILE     <- file.path(INPUT_DIR, "AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet")
DAILY_FILE      <- file.path(INPUT_DIR, "IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet")
THRESH_FILE     <- file.path(INPUT_DIR, "AMJ_TX90_TN10_thresholds_1981_2010.parquet")

## No monthly index file exists in this pipeline (see header note (a)).
## Deliberately NA. Pointing it at a real file would mean inventing one.
MONTHLY_INDEX_FILE <- NA_character_

## ---- official India state boundary (cartographic overlay for 08.3) ---------
## Forward slashes; R rejects single backslashes in string literals.
SHAPEFILE          <- "F:/WMO_IMD_R/WMO_IMD/shapefile/STATE_BOUNDARY.shp"
SHAPEFILE_REQUIRED <- TRUE    # TRUE = 08.3 is FATAL without a usable boundary
GRID_CRS_EPSG      <- 4326L   # the 283 IMD grids are plain lon/lat degrees
SHP_LINE_COLOUR    <- "grey20"
SHP_LINE_WIDTH     <- 0.30

## Regional classification for 08.10. A boundary shapefile is NOT one of these.
REGION_FILE <- NA_character_     # CSV with columns: grid_id, region

## ---- locked analysis parameters (must match STEP 07 / 07B) -----------------
PRIMARY_START_YEAR <- 1996L
PRIMARY_END_YEAR   <- 2025L
CONTEXT_START_YEAR <- 1976L
CONTEXT_END_YEAR   <- 2025L
PRIMARY_VARIANT    <- "HR_lag3"
SENSITIVITY_VARIANTS <- c("MK_original", "HR_lag3", "HR_all")
FDR_Q      <- 0.10
RAW_ALPHA  <- 0.05
COMPLETENESS_THRESHOLD <- 0.80
EXP_GRIDS  <- 283L

INDICES <- c("TXx", "TXn", "TNx", "TNn", "TX90p", "TN10p", "WSDI", "DTR")
MONTHS  <- c(4L, 5L, 6L)
MONTH_LABELS <- c("4" = "April", "5" = "May", "6" = "June")

INDEX_DIMENSION <- c(TXx = "Intensity", TXn = "Intensity",
                     TNx = "Intensity", TNn = "Intensity",
                     TX90p = "Frequency", TN10p = "Frequency",
                     WSDI = "Persistence", DTR = "Diurnal contrast")
INDEX_UNITS <- c(TXx = "degC/decade", TXn = "degC/decade",
                 TNx = "degC/decade", TNn = "degC/decade",
                 TX90p = "%/decade",  TN10p = "%/decade",
                 WSDI = "days/decade", DTR = "degC/decade")

## monthly-climatology policy (see header note (a))
MONTHLY_DERIVABLE  <- c("TXx", "TXn", "TNx", "TNn", "DTR")
MONTHLY_BLOCKED    <- c("WSDI")
MONTHLY_TX90P_MODE <- "direct_labelled"   # "direct_labelled" or "skip"

## ---- execution switches ----------------------------------------------------
RUN_SECTION <- "ALL"     # "ALL" or "08.1" .. "08.12"
MAKE_CSV  <- TRUE
MAKE_XLSX <- TRUE
MAKE_PNG  <- TRUE
MAKE_PDF  <- TRUE
DPI       <- 300
OVERWRITE <- FALSE
FIG_W     <- 9.5
FIG_H     <- 8.0
AUTO_INSTALL <- TRUE
CRAN_REPO    <- "https://cloud.r-project.org"

## ===========================================================================
## SECTION 0b - PACKAGES, FOLDERS, HELPERS  (no configuration below this line)
## ===========================================================================
.req <- c("data.table", "arrow", "ggplot2")
for (p in .req) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (AUTO_INSTALL) { message("Installing: ", p); install.packages(p, repos = CRAN_REPO) }
    if (!requireNamespace(p, quietly = TRUE))
      stop("FATAL: required package '", p, "' unavailable.", call. = FALSE)
  }
}
suppressPackageStartupMessages({ library(data.table); library(arrow); library(ggplot2) })
setDTthreads(0L)

HAS_XLSX <- requireNamespace("writexl", quietly = TRUE)
HAS_SF   <- requireNamespace("sf", quietly = TRUE)
if (!HAS_SF && isTRUE(SHAPEFILE_REQUIRED)) {
  if (AUTO_INSTALL) { message("Installing: sf"); install.packages("sf", repos = CRAN_REPO) }
  HAS_SF <- requireNamespace("sf", quietly = TRUE)
}
if (!HAS_SF && isTRUE(SHAPEFILE_REQUIRED))
  stop("FATAL: SHAPEFILE_REQUIRED is TRUE but package 'sf' is unavailable.", call. = FALSE)

SUBDIRS <- c("01_descriptive", "02_trend_summary", "03_fdr", "04_sensitivity",
             "05_wsdi", "06_frequency", "07_synthesis", "08_seasonal",
             "09_regional", "10_uncertainty", "11_maps", "12_figures",
             "13_qc", "14_manifest")
for (s in SUBDIRS) dir.create(file.path(OUTPUT_DIR, s), recursive = TRUE, showWarnings = FALSE)

.LOG <- character(0); .WARNINGS <- character(0)
.QC <- list(); .OUTFILES <- list(); .SECTIONS_RUN <- character(0)
t_start <- Sys.time()

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "] ",
                 paste0(..., collapse = ""))
  cat(line, "\n", sep = ""); .LOG <<- c(.LOG, line); invisible(line)
}
log_head <- function(txt) {
  bar <- strrep("-", 78)
  cat("\n", bar, "\n", txt, "\n", bar, "\n", sep = "")
  .LOG <<- c(.LOG, "", bar, txt, bar); invisible(NULL)
}
warn_msg <- function(...) {
  m <- paste0(..., collapse = "")
  .WARNINGS <<- c(.WARNINGS, m); log_msg("  WARN  | ", m); invisible(m)
}
qc_add <- function(name, passed, detail, severity = c("FATAL", "WARN", "INFO")) {
  severity <- match.arg(severity)
  ok <- length(passed) == 1L && isTRUE(passed)
  status <- if (severity == "INFO") "INFO" else if (ok) "PASS" else severity
  .QC[[length(.QC) + 1L]] <<- data.table(check = name, status = status,
                                         detail = detail, severity = severity)
  log_msg(sprintf("  %-5s | %-44s | %s", status, name, detail))
  if (!ok && severity == "WARN") .WARNINGS <<- c(.WARNINGS, paste0(name, ": ", detail))
  invisible(ok)
}
fatal <- function(msg) stop("FATAL: ", msg, call. = FALSE)

SEC <- function(id) identical(RUN_SECTION, "ALL") || identical(RUN_SECTION, id)
mark_run <- function(id) .SECTIONS_RUN <<- unique(c(.SECTIONS_RUN, id))

register_out <- function(path) {
  .OUTFILES[[length(.OUTFILES) + 1L]] <<-
    data.table(file = normalizePath(path, winslash = "/", mustWork = FALSE),
               exists = file.exists(path),
               size_kb = if (file.exists(path)) round(file.size(path) / 1024, 1) else NA_real_)
  invisible(path)
}
can_write <- function(path) {
  if (file.exists(path) && !OVERWRITE) {
    warn_msg("exists and OVERWRITE = FALSE, skipped: ", basename(path)); return(FALSE)
  }
  TRUE
}
save_table <- function(dt, subdir, base) {
  stopifnot(is.data.frame(dt))
  if (MAKE_CSV) {
    f <- file.path(OUTPUT_DIR, subdir, paste0(base, ".csv"))
    if (can_write(f)) fwrite(dt, f)
    register_out(f)
  }
  if (MAKE_XLSX) {
    if (!HAS_XLSX) {
      warn_msg("MAKE_XLSX = TRUE but 'writexl' is not installed; XLSX skipped for ", base)
    } else {
      f <- file.path(OUTPUT_DIR, subdir, paste0(base, ".xlsx"))
      if (can_write(f)) writexl::write_xlsx(as.data.frame(dt), f)
      register_out(f)
    }
  }
  invisible(dt)
}
save_fig <- function(pl, subdir, base, w = FIG_W, h = FIG_H) {
  if (MAKE_PNG) {
    f <- file.path(OUTPUT_DIR, subdir, paste0(base, ".png"))
    if (can_write(f)) ggplot2::ggsave(f, pl, width = w, height = h, dpi = DPI, units = "in")
    register_out(f)
  }
  if (MAKE_PDF) {
    f <- file.path(OUTPUT_DIR, subdir, paste0(base, ".pdf"))
    if (can_write(f)) ggplot2::ggsave(f, pl, width = w, height = h, units = "in")
    register_out(f)
  }
  invisible(pl)
}
sfmt  <- function(v, d = 4) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else round(v, d) }
mstat <- function(v, f)     { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else f(v) }

theme_pub <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "right")

## ---- lazy data loaders -----------------------------------------------------
.CACHE <- new.env(parent = emptyenv())
need_file <- function(path, label) {
  if (is.na(path) || !nzchar(path) || !file.exists(path))
    fatal(paste0("required input for ", label, " not found: ", path))
  invisible(TRUE)
}
get_all8 <- function() {
  if (!is.null(.CACHE$all8)) return(.CACHE$all8)
  need_file(TREND_FILE, "the primary trend table")
  d <- as.data.table(fread(TREND_FILE))
  must <- c("index", "grid_id", "lat", "lon", "sen_per_decade",
            "ci_lo_per_decade", "ci_hi_per_decade", "units", "p", "q",
            "trend_class", "significant_fdr", "significant_raw", "cf_status",
            "variant", "period")
  mis <- setdiff(must, names(d))
  if (length(mis))
    fatal(paste0("TREND_FILE is missing required columns: ", paste(mis, collapse = ", "),
                 "\n  present: ", paste(names(d), collapse = ", ")))
  d[, grid_id := as.character(grid_id)]
  d[, significant_fdr := as.logical(significant_fdr)]
  d[, significant_raw := as.logical(significant_raw)]
  d[, index := factor(index, levels = INDICES)]
  setorder(d, index, grid_id)
  .CACHE$all8 <- d
  d
}
get_variants <- function() {
  if (!is.null(.CACHE$var)) return(.CACHE$var)
  need_file(VARIANTS_FILE_A, "the STEP 07 all-variant table")
  need_file(VARIANTS_FILE_B, "the STEP 07B all-variant table")
  a <- as.data.table(arrow::read_parquet(VARIANTS_FILE_A))
  b <- as.data.table(arrow::read_parquet(VARIANTS_FILE_B))
  keep <- intersect(names(a), names(b))
  d <- rbind(a[, ..keep], b[, ..keep], use.names = TRUE)
  d[, grid_id := as.character(grid_id)]
  .CACHE$var <- d
  d
}
get_step06 <- function() {
  if (!is.null(.CACHE$s6)) return(.CACHE$s6)
  need_file(STEP06_FILE, "the STEP 06 seasonal index table")
  d <- as.data.table(arrow::read_parquet(STEP06_FILE))
  d[, `:=`(grid_id = as.character(grid_id), year = as.integer(year))]
  .CACHE$s6 <- d
  d
}
get_step05 <- function() {
  if (!is.null(.CACHE$s5)) return(.CACHE$s5)
  need_file(STEP05_FILE, "the STEP 05 percentile index table")
  d <- as.data.table(arrow::read_parquet(STEP05_FILE))
  d[, `:=`(grid_id = as.character(grid_id), year = as.integer(year))]
  .CACHE$s5 <- d
  d
}
get_daily_amj <- function() {
  if (!is.null(.CACHE$daily)) return(.CACHE$daily)
  need_file(DAILY_FILE, "monthly climatology (08.1 / 08.9)")
  d <- as.data.table(arrow::read_parquet(DAILY_FILE))
  nm <- tolower(names(d))
  pick <- function(cands, lab) {
    h <- names(d)[nm %in% cands]
    if (!length(h)) fatal(paste0("daily file lacks a '", lab, "' column; present: ",
                                 paste(names(d), collapse = ", ")))
    h[1]
  }
  cg <- pick(c("grid_id", "gridid", "grid", "id"), "grid_id")
  cd <- pick(c("date", "obs_date", "time"), "date")
  cx <- pick(c("tmax", "tx", "t_max"), "tmax")
  cn <- pick(c("tmin", "tn", "t_min"), "tmin")
  d <- d[, c(cg, cd, cx, cn), with = FALSE]
  setnames(d, c("grid_id", "date", "tmax", "tmin"))
  d[, date := as.IDate(date)]
  d[, `:=`(grid_id = as.character(grid_id),
           year = as.integer(data.table::year(date)),
           month = as.integer(data.table::month(date)),
           tmax = as.numeric(tmax), tmin = as.numeric(tmin))]
  d <- d[month %in% MONTHS & year >= PRIMARY_START_YEAR & year <= PRIMARY_END_YEAR]
  .CACHE$daily <- d
  d
}

## ---- state boundary loader -------------------------------------------------
## Reads the official state boundary with sf, reports its CRS, repairs invalid
## geometry, transforms to the grid CRS, and returns a plain coordinate table so
## the existing coord_fixed() map layout is left completely untouched.
.SHP_INFO <- list(loaded = FALSE, crs_in = NA_character_, crs_out = NA_character_,
                  n_features = NA_integer_, invalid_before = NA_integer_,
                  transformed = NA, layers_added = 0L)
get_state_boundary <- function() {
  if (!is.null(.CACHE$shp)) return(.CACHE$shp)
  if (is.na(SHAPEFILE) || !nzchar(SHAPEFILE))
    fatal("SHAPEFILE is not set. Mapping (08.3) cannot proceed; no grid-only fallback is permitted.")
  if (!file.exists(SHAPEFILE))
    fatal(paste0("SHAPEFILE not found on disk: ", SHAPEFILE,
                 "\n  Check the path and that the sidecar files (.shx, .dbf, .prj) are present."))
  for (ext in c("shx", "dbf")) {
    side <- sub("\\.shp$", paste0(".", ext), SHAPEFILE, ignore.case = TRUE)
    if (!file.exists(side))
      fatal(paste0("required shapefile component missing: ", basename(side)))
  }
  b <- try(sf::st_read(SHAPEFILE, quiet = TRUE), silent = TRUE)
  if (inherits(b, "try-error"))
    fatal(paste0("sf::st_read failed for ", SHAPEFILE, "\n  ", as.character(b)))
  if (!nrow(b)) fatal("the shapefile contains zero features.")
  .SHP_INFO$n_features <<- nrow(b)
  
  crs_in <- sf::st_crs(b)
  if (is.na(crs_in))
    fatal(paste0("the shapefile has NO coordinate reference system (missing or unreadable .prj): ",
                 SHAPEFILE,
                 "\n  Assign the correct CRS at source; STEP 08 will not guess one."))
  .SHP_INFO$crs_in <<- if (!is.na(crs_in$epsg)) paste0("EPSG:", crs_in$epsg) else
    substr(crs_in$input, 1, 120)
  log_msg("Boundary CRS in  : ", .SHP_INFO$crs_in, "  | features = ", nrow(b))
  
  ## --- validity repair --------------------------------------------------
  bad <- suppressWarnings(sum(!sf::st_is_valid(b), na.rm = TRUE))
  .SHP_INFO$invalid_before <<- bad
  if (bad > 0L) {
    warn_msg(bad, " invalid boundary geometries detected; repaired with sf::st_make_valid().")
    b <- sf::st_make_valid(b)
    if (suppressWarnings(any(!sf::st_is_valid(b), na.rm = TRUE)))
      fatal("boundary geometry remains invalid after st_make_valid(); cannot render maps.")
  }
  
  ## --- geometry normalisation -------------------------------------------
  ## st_make_valid() can return a GEOMETRYCOLLECTION for a repaired feature
  ## (the polygon plus zero-area LINESTRING / POINT slivers). That makes the
  ## whole geometry column sfc_GEOMETRY, which st_coordinates() cannot handle:
  ##   "not implemented for objects of class sfc_GEOMETRY"
  ## Pull the POLYGON component out of any collection, discard the non-areal
  ## slivers (they carry no state area), and cast everything to MULTIPOLYGON so
  ## each state survives as exactly one feature. CRS is untouched here.
  n_before  <- nrow(b)
  gt_before <- as.character(sf::st_geometry_type(b))
  .SHP_INFO$geom_types_before <<-
    paste(sprintf("%s=%d", names(table(gt_before)), as.integer(table(gt_before))),
          collapse = ", ")
  log_msg("Boundary geoms in: ", .SHP_INFO$geom_types_before)
  
  if (any(gt_before == "GEOMETRYCOLLECTION")) {
    n_gc <- sum(gt_before == "GEOMETRYCOLLECTION")
    warn_msg(n_gc, " feature(s) became GEOMETRYCOLLECTION after repair; extracting the ",
             "POLYGON component and discarding zero-area slivers.")
    b <- suppressWarnings(sf::st_collection_extract(b, "POLYGON", warn = FALSE))
  }
  
  gt_mid <- as.character(sf::st_geometry_type(b))
  keep   <- gt_mid %in% c("POLYGON", "MULTIPOLYGON")
  if (any(!keep)) {
    warn_msg(sum(!keep), " non-areal feature(s) (",
             paste(unique(gt_mid[!keep]), collapse = ", "),
             ") dropped; they contribute no state area.")
    b <- b[keep, , drop = FALSE]
  }
  if (!nrow(b))
    fatal("no polygonal features survived geometry normalisation; the boundary cannot be drawn.")
  
  ## cast to a single homogeneous MULTIPOLYGON type; fall back feature-by-feature
  b_cast <- try(suppressWarnings(sf::st_cast(b, "MULTIPOLYGON")), silent = TRUE)
  if (inherits(b_cast, "try-error")) {
    warn_msg("bulk st_cast to MULTIPOLYGON failed; casting feature by feature.")
    sf::st_geometry(b) <- sf::st_sfc(
      lapply(sf::st_geometry(b), function(g)
        if (inherits(g, "MULTIPOLYGON")) g else sf::st_cast(g, "MULTIPOLYGON")),
      crs = sf::st_crs(b))
  } else {
    b <- b_cast
  }
  
  ## the column must no longer be sfc_GEOMETRY or st_coordinates() will fail again
  if (inherits(sf::st_geometry(b), "sfc_GEOMETRY"))
    fatal(paste0("boundary geometry is still mixed (sfc_GEOMETRY) after normalisation; ",
                 "types present: ",
                 paste(unique(as.character(sf::st_geometry_type(b))), collapse = ", ")))
  
  gt_after <- as.character(sf::st_geometry_type(b))
  .SHP_INFO$geom_types_after  <<-
    paste(sprintf("%s=%d", names(table(gt_after)), as.integer(table(gt_after))),
          collapse = ", ")
  .SHP_INFO$n_features_after  <<- nrow(b)
  .SHP_INFO$n_features_dropped <<- n_before - nrow(b)
  log_msg("Boundary geoms out: ", .SHP_INFO$geom_types_after,
          "  | features ", n_before, " -> ", nrow(b))
  if (nrow(b) < n_before)
    warn_msg(n_before - nrow(b), " of ", n_before,
             " boundary features were lost during normalisation; verify no state is missing.")
  
  tgt <- sf::st_crs(GRID_CRS_EPSG)
  if (crs_in != tgt) {
    b <- try(sf::st_transform(b, tgt), silent = TRUE)
    if (inherits(b, "try-error"))
      fatal(paste0("failed to transform the boundary to EPSG:", GRID_CRS_EPSG, "\n  ",
                   as.character(b)))
    .SHP_INFO$transformed <<- TRUE
    log_msg("Boundary CRS out : EPSG:", GRID_CRS_EPSG, "  (transformed)")
  } else {
    .SHP_INFO$transformed <<- FALSE
    log_msg("Boundary CRS out : EPSG:", GRID_CRS_EPSG, "  (already matched; no transform)")
  }
  .SHP_INFO$crs_out <<- paste0("EPSG:", GRID_CRS_EPSG)
  
  cc <- sf::st_coordinates(sf::st_geometry(b))
  if (!nrow(cc)) fatal("the boundary geometry yielded no coordinates.")
  cdt <- as.data.table(cc)
  lcols <- grep("^L[0-9]+$", names(cdt), value = TRUE)
  if (!length(lcols)) cdt[, L1 := 1L]
  lcols <- grep("^L[0-9]+$", names(cdt), value = TRUE)
  cdt[, shp_grp := do.call(paste, c(.SD, sep = "_")), .SDcols = lcols]
  cdt <- cdt[, .(x = X, y = Y, shp_grp)]
  
  .SHP_INFO$loaded <<- TRUE
  log_msg("Boundary loaded  : ", format(nrow(cdt), big.mark = ","), " vertices in ",
          uniqueN(cdt$shp_grp), " rings")
  .CACHE$shp <- cdt
  cdt
}

log_head(paste0("STEP 08 | FINAL RESULTS, TABLES, FIGURES AND QC   (v", SCRIPT_VERSION, ")"))
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("ggplot2          : ", as.character(utils::packageVersion("ggplot2")))
if (HAS_SF) log_msg("sf               : ", as.character(utils::packageVersion("sf")))
log_msg("RUN_SECTION      : ", RUN_SECTION)
log_msg("Primary period   : ", PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR,
        " | variant ", PRIMARY_VARIANT, " | BH q <= ", FDR_Q)
log_msg("Output directory : ", OUTPUT_DIR)
log_msg("NOTE             : STEP 08 performs no statistical re-estimation.")
log_msg("State boundary   : ", ifelse(is.na(SHAPEFILE), "<not supplied>", SHAPEFILE),
        "  (required = ", SHAPEFILE_REQUIRED, ")")
if (is.na(MONTHLY_INDEX_FILE))
  warn_msg("no monthly index file exists in this pipeline; 08.1/08.9 derive monthly ",
           "climatology from the daily source for ", paste(MONTHLY_DERIVABLE, collapse = "/"),
           " only. Monthly WSDI is BLOCKED (see the header).")

## ===========================================================================
## 08.1 - DESCRIPTIVE STATISTICS
## ===========================================================================
if (SEC("08.1")) {
  log_head("08.1 | DESCRIPTIVE STATISTICS (climatology - NOT trends)")
  mark_run("08.1")
  
  desc_of <- function(v) {
    vv <- v[!is.na(v)]
    list(N = length(vv), N_missing = sum(is.na(v)),
         Mean = mstat(vv, mean), Median = mstat(vv, stats::median),
         SD = mstat(vv, stats::sd),
         Q1 = mstat(vv, function(z) stats::quantile(z, .25, names = FALSE)),
         Q3 = mstat(vv, function(z) stats::quantile(z, .75, names = FALSE)),
         IQR = mstat(vv, stats::IQR),
         Minimum = mstat(vv, min), Maximum = mstat(vv, max))
  }
  
  ## ---- (a) seasonal AMJ, all eight indices ---------------------------------
  s6 <- get_step06(); s5 <- get_step05()
  seas <- rbindlist(lapply(INDICES, function(ix) {
    src <- if (ix %in% names(s6)) s6 else if (ix %in% names(s5)) s5 else NULL
    if (is.null(src)) { warn_msg("index ", ix, " not found in STEP 05 or STEP 06"); return(NULL) }
    sub <- src[year >= PRIMARY_START_YEAR & year <= PRIMARY_END_YEAR]
    v <- sub[[ix]]
    out <- as.data.table(desc_of(v))
    out[, `:=`(index = ix, scope = "AMJ_season",
               period = paste0(PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
               source = if (ix %in% names(s6)) basename(STEP06_FILE) else basename(STEP05_FILE),
               unit = sub("/decade", "", INDEX_UNITS[[ix]]))]
    if (ix == "WSDI") {
      out[, `:=`(zero_fraction_grid_years = round(mean(v == 0, na.rm = TRUE), 4),
                 n_zero_grid_years = sum(v == 0, na.rm = TRUE))]
    } else {
      out[, `:=`(zero_fraction_grid_years = NA_real_, n_zero_grid_years = NA_integer_)]
    }
    out[]
  }), fill = TRUE)
  setcolorder(seas, c("index", "scope", "period", "unit", "N", "N_missing",
                      "Mean", "Median", "SD", "Q1", "Q3", "IQR", "Minimum", "Maximum",
                      "zero_fraction_grid_years", "n_zero_grid_years", "source"))
  save_table(seas, "01_descriptive",
             sprintf("STEP08_01_DESCRIPTIVE_SEASONAL_AMJ_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  qc_add("desc_seasonal_all_indices", nrow(seas) == length(INDICES),
         sprintf("%d of %d indices summarised at the AMJ seasonal scale",
                 nrow(seas), length(INDICES)))
  
  ## ---- (b) monthly climatology, derived here -------------------------------
  dl <- get_daily_amj()
  mon <- dl[, .(TXx = if (all(is.na(tmax))) NA_real_ else max(tmax, na.rm = TRUE),
                TXn = if (all(is.na(tmax))) NA_real_ else min(tmax, na.rm = TRUE),
                TNx = if (all(is.na(tmin))) NA_real_ else max(tmin, na.rm = TRUE),
                TNn = if (all(is.na(tmin))) NA_real_ else min(tmin, na.rm = TRUE),
                DTR = { dv <- (tmax - tmin); dv <- dv[!is.na(dv)]
                if (!length(dv)) NA_real_ else mean(dv) },
                n_days = .N, n_valid_tmax = sum(!is.na(tmax)),
                n_valid_tmin = sum(!is.na(tmin)),
                n_paired = sum(!is.na(tmax) & !is.na(tmin))),
            by = .(grid_id, year, month)]
  .CACHE$monthly <- mon
  
  mdesc <- rbindlist(lapply(MONTHLY_DERIVABLE, function(ix) {
    rbindlist(lapply(MONTHS, function(mm) {
      out <- as.data.table(desc_of(mon[month == mm][[ix]]))
      out[, `:=`(index = paste0(ix, "_monthly_climatology"),
                 month = MONTH_LABELS[[as.character(mm)]],
                 scope = "single_month", unit = "degC",
                 period = paste0(PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
                 derivation = "computed in STEP 08 from the daily source; baseline-free; climatology only")]
      out[]
    }))
  }))
  setcolorder(mdesc, c("index", "month", "scope", "unit", "period", "N", "N_missing",
                       "Mean", "Median", "SD", "Q1", "Q3", "IQR", "Minimum", "Maximum"))
  save_table(mdesc, "01_descriptive",
             sprintf("STEP08_01_DESCRIPTIVE_MONTHLY_CLIMATOLOGY_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  ## ---- monthly WSDI is BLOCKED --------------------------------------------
  blocked <- data.table(
    index = MONTHLY_BLOCKED, status = "BLOCKED",
    reason = paste("A monthly WSDI requires resetting the warm-spell counter at",
                   "30 Apr and 31 May. The locked STEP 06 definition treats AMJ as ONE",
                   "continuous 91-day sequence with no month reset, so a monthly WSDI",
                   "is a different index, not a subset of the published one."),
    required_to_unblock = paste("An explicit, separately documented decision to define and",
                                "compute a month-reset WSDI as an additional index."))
  save_table(blocked, "01_descriptive", "STEP08_01_MONTHLY_WSDI_BLOCKED_REPORT")
  qc_add("monthly_wsdi_blocked", FALSE,
         "monthly WSDI not computed: conflicts with the locked continuous-91-day definition",
         severity = "WARN")
  
  ## ---- monthly TX90p / TN10p ----------------------------------------------
  if (identical(MONTHLY_TX90P_MODE, "direct_labelled")) {
    if (!file.exists(THRESH_FILE)) {
      warn_msg("threshold file absent; monthly TX90p/TN10p climatology skipped.")
    } else {
      thr <- as.data.table(arrow::read_parquet(THRESH_FILE))
      tg <- names(thr)[tolower(names(thr)) %in% c("grid_id", "gridid", "grid", "id")][1]
      tm <- names(thr)[tolower(names(thr)) %in% c("target_md", "md", "mmdd")][1]
      tx <- names(thr)[tolower(names(thr)) %in% c("tx90", "tx90_threshold")][1]
      tn <- names(thr)[tolower(names(thr)) %in% c("tn10", "tn10_threshold")][1]
      if (any(is.na(c(tg, tm, tx, tn)))) {
        warn_msg("threshold file schema not recognised; monthly TX90p/TN10p skipped.")
      } else {
        thr <- thr[, c(tg, tm, tx, tn), with = FALSE]
        setnames(thr, c("grid_id", "md", "tx90", "tn10"))
        thr[, grid_id := as.character(grid_id)]
        dd <- copy(dl)
        dd[, md := as.integer(data.table::month(date)) * 100L +
             as.integer(data.table::mday(date))]
        dd <- merge(dd, thr, by = c("grid_id", "md"), all.x = TRUE)
        dd[, `:=`(evx = !is.na(tmax) & !is.na(tx90), evn = !is.na(tmin) & !is.na(tn10))]
        pm <- dd[, .(TX90p_direct_monthly =
                       { nv <- sum(evx); if (!nv) NA_real_ else 100 * sum(evx & tmax > tx90) / nv },
                     TN10p_direct_monthly =
                       { nv <- sum(evn); if (!nv) NA_real_ else 100 * sum(evn & tmin < tn10) / nv },
                     n_evaluable_tmax = sum(evx), n_evaluable_tmin = sum(evn)),
                 by = .(grid_id, year, month)]
        pdesc <- rbindlist(lapply(c("TX90p_direct_monthly", "TN10p_direct_monthly"), function(ix)
          rbindlist(lapply(MONTHS, function(mm) {
            out <- as.data.table(desc_of(pm[month == mm][[ix]]))
            out[, `:=`(index = ix, month = MONTH_LABELS[[as.character(mm)]],
                       scope = "single_month", unit = "%",
                       period = paste0(PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
                       derivation = paste("DIRECT exceedance against the locked 1981-2010",
                                          "calendar-day thresholds. The Zhang (2005) in-base",
                                          "bootstrap is NOT applied monthly. This is NOT the",
                                          "published seasonal TX90p/TN10p and must never be",
                                          "used in a trend context."))]
            out[]
          }))))
        setcolorder(pdesc, c("index", "month", "scope", "unit", "period"))
        save_table(pdesc, "01_descriptive",
                   sprintf("STEP08_01_DESCRIPTIVE_MONTHLY_TX90p_TN10p_DIRECT_%d_%d",
                           PRIMARY_START_YEAR, PRIMARY_END_YEAR))
        qc_add("monthly_percentile_direct_only", FALSE,
               "monthly TX90p/TN10p computed as DIRECT exceedance only (no bootstrap); labelled separately",
               severity = "WARN")
      }
    }
  } else {
    warn_msg("MONTHLY_TX90P_MODE = 'skip'; monthly percentile climatology not produced.")
  }
  
  log_msg("08.1 complete: seasonal table (", nrow(seas), " rows), monthly table (",
          nrow(mdesc), " rows).")
}

## ===========================================================================
## 08.2 - FINAL TREND SUMMARY  (principal manuscript table)
## ===========================================================================
if (SEC("08.2")) {
  log_head("08.2 | FINAL TREND SUMMARY (primary design)")
  mark_run("08.2")
  d <- get_all8()
  
  tsum <- d[, {
    ntest <- sum(is.finite(p))
    .(n_grids = .N, n_testable = ntest,
      n_test_undefined = sum(trend_class == "test_undefined"),
      n_insufficient_data = sum(trend_class == "insufficient_data"),
      mean_sen = sfmt(mean(sen_per_decade, na.rm = TRUE)),
      median_sen = sfmt(stats::median(sen_per_decade, na.rm = TRUE)),
      sd_sen = sfmt(stats::sd(sen_per_decade, na.rm = TRUE)),
      iqr_sen = sfmt(stats::IQR(sen_per_decade, na.rm = TRUE)),
      min_sen = sfmt(mstat(sen_per_decade, min)),
      max_sen = sfmt(mstat(sen_per_decade, max)),
      median_ci_lo = sfmt(stats::median(ci_lo_per_decade, na.rm = TRUE)),
      median_ci_hi = sfmt(stats::median(ci_hi_per_decade, na.rm = TRUE)),
      median_ci_width = sfmt(stats::median(ci_hi_per_decade - ci_lo_per_decade, na.rm = TRUE)),
      n_raw_sig = sum(significant_raw, na.rm = TRUE),
      n_fdr_sig = sum(significant_fdr, na.rm = TRUE),
      n_sig_increase = sum(trend_class == "significant_increase"),
      n_sig_decrease = sum(trend_class == "significant_decrease"),
      n_sig_zero_slope = sum(trend_class == "significant_zero_slope"),
      n_no_trend = sum(trend_class == "no_significant_trend"),
      pct_raw_sig = round(100 * sum(significant_raw, na.rm = TRUE) / max(ntest, 1L), 2),
      pct_fdr_sig = round(100 * sum(significant_fdr, na.rm = TRUE) / max(ntest, 1L), 2),
      pct_sig_increase = round(100 * sum(trend_class == "significant_increase") / .N, 2),
      pct_sig_decrease = round(100 * sum(trend_class == "significant_decrease") / .N, 2),
      pct_no_trend = round(100 * sum(trend_class == "no_significant_trend") / .N, 2),
      units = units[1])
  }, by = index]
  setorder(tsum, index)
  tsum[, `:=`(period = paste0(PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
              variant = PRIMARY_VARIANT, fdr_q = FDR_Q, raw_alpha = RAW_ALPHA)]
  
  save_table(tsum, "02_trend_summary",
             sprintf("STEP08_02_FINAL_TREND_SUMMARY_PRIMARY_%s_%d_%d",
                     PRIMARY_VARIANT, PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  for (i in seq_len(nrow(tsum))) {
    z <- tsum[i]
    log_msg(sprintf("  %-6s %-12s median %8.4f  raw %3d  FDR %3d  (+%d / -%d)  no-trend %5.1f%%",
                    z$index, z$units, z$median_sen, z$n_raw_sig, z$n_fdr_sig,
                    z$n_sig_increase, z$n_sig_decrease, z$pct_no_trend))
  }
  qc_add("trend_summary_index_count", nrow(tsum) == length(INDICES),
         sprintf("%d indices in the summary table", nrow(tsum)))
  qc_add("trend_summary_grid_counts", all(tsum$n_grids == EXP_GRIDS),
         sprintf("grid counts: %s", paste(unique(tsum$n_grids), collapse = ", ")))
  qc_add("class_counts_sum_to_grids",
         tsum[, all(n_sig_increase + n_sig_decrease + n_sig_zero_slope +
                      n_no_trend + n_test_undefined + n_insufficient_data == n_grids)],
         "class counts sum to the grid count for every index")
  qc_add("percentages_within_range",
         tsum[, all(c(pct_raw_sig, pct_fdr_sig, pct_sig_increase,
                      pct_sig_decrease, pct_no_trend) %between% c(0, 100))],
         "all reported percentages lie in [0, 100]")
}
## ===========================================================================
## 08.3 - SPATIAL TREND MAPS      [REVISED DROP-IN BLOCK]
## ---------------------------------------------------------------------------
##  Renders the LOCKED STEP 07 / 07A / 07B primary trend table only.
##  NO statistic is recomputed, refitted, re-tested or re-adjusted here.
##  Every number drawn is read verbatim from TREND_FILE via get_all8().
##
##  Cartographic contract:
##    * one geom_tile per IMD 1 deg grid centre, width = height = 1 deg
##    * NO interpolation, NO smoothing, NO regridding, NO kriging, NO clipping
##    * fill  = Theil-Sen slope per decade, zero-centred diverging scale
##    * fill is INDEPENDENT of significance (never masked, never greyed by q)
##    * stipple (points) = BH-FDR q <= FDR_Q
##    * solid black cell outline  = test_undefined  (CF <= 0; slope shown, no test)
##    * dotted grey cell outline  = insufficient_data (no slope; grey fill)
##    * state boundary drawn as geom_path over st_coordinates(), NOT geom_sf,
##      so coord_fixed() is preserved and coord_sf() cannot silently replace
##      the aspect ratio or the extent
##    * identical extent (MAP_XLIM / MAP_YLIM) on individual AND combined maps
## ===========================================================================
if (SEC("08.3")) {
  log_head("08.3 | SPATIAL TREND MAPS (continuous slope + FDR stipple + boundary)")
  mark_run("08.3")
  
  ## ---- 08.3.0 local map configuration --------------------------------------
  MAP_DPI        <- 600L        # publication raster resolution (PDF is vector)
  MAP_PAD_DEG    <- 0.25        # padding beyond the outer tile edges
  MAP_LIM_PROB   <- 1.00        # 1.00 = full range (no squishing); set 0.98 for a robust scale
  MAP_ASPECT     <- 1.05        # unchanged from the original block
  MAP_SHAPEFILE  <- "F:/WMO_IMD_R/WMO_IMD/shapefile/STATE_BOUNDARY.shp"
  
  if (is.na(SHAPEFILE) && file.exists(MAP_SHAPEFILE)) {
    SHAPEFILE <<- MAP_SHAPEFILE
    log_msg("SHAPEFILE was NA; set from the 08.3 map configuration: ", MAP_SHAPEFILE)
  }
  MAP_SHP <- if (!is.na(SHAPEFILE)) SHAPEFILE else MAP_SHAPEFILE
  
  ## ---- 08.3.1 boundary loader ----------------------------------------------
  ## Defined ONLY if the running script does not already provide it, so a
  ## later revision of STEP 08 that ships get_state_boundary() keeps its own.
  if (!exists("get_state_boundary", mode = "function")) {
    .SHP_INFO <- new.env(parent = emptyenv())
    get_state_boundary <- function(path = MAP_SHP, target_epsg = 4326L) {
      if (!HAS_SF) { warn_msg("package 'sf' is not installed; boundary overlay skipped."); return(NULL) }
      if (is.na(path) || !nzchar(path) || !file.exists(path)) {
        warn_msg("state boundary not found (", path, "); overlay skipped."); return(NULL)
      }
      b <- try(sf::st_read(path, quiet = TRUE), silent = TRUE)
      if (inherits(b, "try-error") || !nrow(b)) {
        warn_msg("state boundary unreadable or empty; overlay skipped."); return(NULL)
      }
      crs_in <- sf::st_crs(b)
      if (is.na(crs_in)) {
        warn_msg("state boundary has no CRS (.prj missing); overlay skipped rather than guessed.")
        return(NULL)
      }
      log_msg("Boundary CRS in  : ",
              if (!is.na(crs_in$epsg)) paste0("EPSG:", crs_in$epsg) else crs_in$input,
              " | features = ", nrow(b))
      
      n_before <- nrow(b)
      bad <- suppressWarnings(sum(!sf::st_is_valid(b), na.rm = TRUE))
      if (bad > 0L) {
        warn_msg(bad, " invalid boundary geometries; repaired with sf::st_make_valid().")
        b <- sf::st_make_valid(b)
      }
      gt <- as.character(sf::st_geometry_type(b))
      if (any(gt == "GEOMETRYCOLLECTION")) {
        warn_msg(sum(gt == "GEOMETRYCOLLECTION"),
                 " feature(s) became GEOMETRYCOLLECTION after repair; extracting POLYGON parts.")
        b <- suppressWarnings(sf::st_collection_extract(b, "POLYGON", warn = FALSE))
      }
      gt <- as.character(sf::st_geometry_type(b))
      keep <- gt %in% c("POLYGON", "MULTIPOLYGON")
      if (any(!keep)) {
        warn_msg(sum(!keep), " non-areal feature(s) dropped; they carry no state area.")
        b <- b[keep, , drop = FALSE]
      }
      if (!nrow(b)) { warn_msg("no polygonal features survived normalisation; overlay skipped."); return(NULL) }
      bc <- try(suppressWarnings(sf::st_cast(b, "MULTIPOLYGON")), silent = TRUE)
      if (!inherits(bc, "try-error")) b <- bc
      if (inherits(sf::st_geometry(b), "sfc_GEOMETRY")) {
        warn_msg("boundary geometry remains mixed (sfc_GEOMETRY); overlay skipped.")
        return(NULL)
      }
      if (is.na(crs_in$epsg) || !identical(as.integer(crs_in$epsg), as.integer(target_epsg))) {
        b <- sf::st_transform(b, target_epsg)
        log_msg("Boundary CRS out : EPSG:", target_epsg, " (transformed)")
      } else {
        log_msg("Boundary CRS out : EPSG:", target_epsg, " (already geographic)")
      }
      if (nrow(b) < n_before)
        warn_msg(n_before - nrow(b), " of ", n_before,
                 " boundary features lost in normalisation; verify no state is missing.")
      
      xy <- sf::st_coordinates(sf::st_geometry(b))
      xy <- as.data.table(xy)
      idc <- intersect(c("L1", "L2", "L3"), names(xy))
      xy[, ring := do.call(paste, c(.SD, sep = "_")), .SDcols = idc]
      out <- xy[, .(x = X, y = Y, ring = ring)]
      assign("n_features", nrow(b), envir = .SHP_INFO)
      assign("n_vertices", nrow(out), envir = .SHP_INFO)
      assign("n_rings",    uniqueN(out$ring), envir = .SHP_INFO)
      log_msg("Boundary loaded  : ", nrow(out), " vertices in ",
              uniqueN(out$ring), " rings from ", nrow(b), " features")
      out[]
    }
  }
  
  ## ---- 08.3.2 high-resolution figure writer --------------------------------
  ## Mirrors save_fig() exactly (same can_write / register_out / MAKE_* logic)
  ## and only overrides the raster dpi, so the global DPI used by every other
  ## section is left untouched.
  save_fig_dpi <- function(pl, subdir, base, w = FIG_W, h = FIG_H, dpi = MAP_DPI) {
    if (MAKE_PNG) {
      f <- file.path(OUTPUT_DIR, subdir, paste0(base, ".png"))
      if (can_write(f)) ggplot2::ggsave(f, pl, width = w, height = h, dpi = dpi,
                                        units = "in", limitsize = FALSE)
      register_out(f)
    }
    if (MAKE_PDF) {
      f <- file.path(OUTPUT_DIR, subdir, paste0(base, ".pdf"))
      if (can_write(f)) ggplot2::ggsave(f, pl, width = w, height = h,
                                        units = "in", limitsize = FALSE)
      register_out(f)
    }
    invisible(pl)
  }
  
  ## ---- 08.3.3 data: locked primary table, defensively filtered -------------
  d <- get_all8()
  
  n_all <- nrow(d)
  if ("variant" %in% names(d) && uniqueN(d$variant) > 1L)
    d <- d[variant == PRIMARY_VARIANT]
  if ("period" %in% names(d) && uniqueN(d$period) > 1L) {
    per_lab <- sprintf("%d-%d", PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    per_alt <- sprintf("%d_%d", PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    d <- d[period %in% c(per_lab, per_alt)]
  }
  if (!nrow(d))
    fatal("08.3: no rows remain after filtering the primary trend table to variant/period.")
  if (nrow(d) != n_all)
    log_msg("08.3 filtered the trend table to ", PRIMARY_VARIANT, " / ",
            PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR, ": ", n_all, " -> ", nrow(d), " rows")
  
  d <- copy(d)
  d[, sen_plot := as.numeric(sen_per_decade)]
  n_insuf_with_slope <- d[trend_class == "insufficient_data" & is.finite(sen_plot), .N]
  if (n_insuf_with_slope > 0L) {
    warn_msg(n_insuf_with_slope, " insufficient_data cell(s) carry a finite slope in the ",
             "trend table; the fill is blanked for display only. No stored value is altered.")
    d[trend_class == "insufficient_data", sen_plot := NA_real_]
  }
  d[, status_cls := fifelse(trend_class == "test_undefined",    "Test undefined (CF <= 0)",
                            fifelse(trend_class == "insufficient_data", "Insufficient data", NA_character_))]
  
  ## ---- 08.3.4 shared extent and boundary -----------------------------------
  MAP_XLIM <- c(min(d$lon, na.rm = TRUE) - 0.5 - MAP_PAD_DEG,
                max(d$lon, na.rm = TRUE) + 0.5 + MAP_PAD_DEG)
  MAP_YLIM <- c(min(d$lat, na.rm = TRUE) - 0.5 - MAP_PAD_DEG,
                max(d$lat, na.rm = TRUE) + 0.5 + MAP_PAD_DEG)
  
  bnd <- get_state_boundary()
  if (!is.null(bnd)) {
    bnd <- bnd[x >= MAP_XLIM[1] - 5 & x <= MAP_XLIM[2] + 5 &
                 y >= MAP_YLIM[1] - 5 & y <= MAP_YLIM[2] + 5]
    if (!nrow(bnd)) { warn_msg("boundary has no vertices near the grid extent; overlay skipped."); bnd <- NULL }
  }
  
  STIP_LAB <- sprintf("BH-FDR q <= %.2f", FDR_Q)
  STATUS_COLS <- c("Test undefined (CF <= 0)" = "black", "Insufficient data" = "grey25")
  STATUS_LTY  <- c("Test undefined (CF <= 0)" = "solid", "Insufficient data" = "dotted")
  
  add_boundary <- function(p) {
    if (is.null(bnd)) return(p)
    p + geom_path(data = bnd, aes(x = x, y = y, group = ring),
                  inherit.aes = FALSE, colour = "grey15", linewidth = 0.28)
  }
  
  sym_lim <- function(v, prob = MAP_LIM_PROB) {
    v <- v[is.finite(v)]
    if (!length(v)) return(1)
    l <- stats::quantile(abs(v), probs = prob, na.rm = TRUE, names = FALSE, type = 8)
    if (!is.finite(l) || l <= 0) l <- suppressWarnings(max(abs(v), na.rm = TRUE))
    if (!is.finite(l) || l <= 0) l <- 1
    ## round UP to 2 significant digits, so at prob = 1.00 nothing is ever squished
    e <- 10^(floor(log10(l)) - 1)
    ceiling(l / e) * e
  }
  
  ## ---- 08.3.5 map builder ---------------------------------------------------
  build_map <- function(dd, ttl, unit_lab, lim, faceted = FALSE, ncol_f = 3L) {
    sig <- dd[significant_fdr %in% TRUE]
    und <- dd[!is.na(status_cls)]
    
    p <- ggplot() +
      geom_tile(data = dd, aes(x = lon, y = lat, fill = sen_plot),
                width = 1, height = 1) +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = 0, limits = c(-lim, lim),
                           oob = scales::squish, na.value = "grey85",
                           name = unit_lab)
    if (nrow(und))
      p <- p + geom_tile(data = und,
                         aes(x = lon, y = lat, colour = status_cls, linetype = status_cls),
                         fill = NA, width = 1, height = 1, linewidth = 0.4,
                         show.legend = c(colour = TRUE, linetype = TRUE,
                                         fill = FALSE, shape = FALSE))
    p <- p +
      scale_colour_manual(values = STATUS_COLS, name = NULL, drop = TRUE,
                          na.translate = FALSE) +
      scale_linetype_manual(values = STATUS_LTY, name = NULL, drop = TRUE,
                            na.translate = FALSE)
    if (nrow(sig))
      p <- p + geom_point(data = sig, aes(x = lon, y = lat, shape = STIP_LAB),
                          size = if (faceted) 0.35 else 0.6, colour = "grey10",
                          show.legend = c(shape = TRUE, colour = FALSE,
                                          linetype = FALSE, fill = FALSE))
    p <- p + scale_shape_manual(values = stats::setNames(16L, STIP_LAB), name = NULL)
    p <- add_boundary(p)
    if (faceted) p <- p + facet_wrap(~ index, ncol = ncol_f)
    p +
      coord_fixed(ratio = MAP_ASPECT, xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
      labs(title = ttl,
           subtitle = sprintf(
             "Theil-Sen slope per decade (colour, independent of significance); stipple = %s", STIP_LAB),
           caption = sprintf(
             "AMJ %d-%d | Hamed-Rao MMK %s | %d IMD 1 deg land grids | native 1 deg tiles, no interpolation",
             PRIMARY_START_YEAR, PRIMARY_END_YEAR, PRIMARY_VARIANT, EXP_GRIDS),
           x = "Longitude (deg E)", y = "Latitude (deg N)") +
      theme_pub +
      theme(legend.key = element_rect(fill = "white", colour = NA))
  }
  
  ## ---- 08.3.6 one map per index --------------------------------------------
  lim_tab <- list()
  for (ix in INDICES) {
    dd <- d[index == ix]
    if (!nrow(dd)) { warn_msg("no rows for index ", ix, "; map skipped."); next }
    lim <- sym_lim(dd$sen_plot)
    lim_tab[[ix]] <- data.table(index = ix, unit = INDEX_UNITS[[ix]],
                                colour_limit = lim,
                                n_grids = nrow(dd),
                                n_fdr_sig = dd[significant_fdr %in% TRUE, .N],
                                n_test_undefined = dd[trend_class == "test_undefined", .N],
                                n_insufficient = dd[trend_class == "insufficient_data", .N],
                                n_squished = dd[is.finite(sen_plot) & abs(sen_plot) > lim, .N])
    pl <- build_map(dd,
                    paste0(ix, " trend, AMJ ", PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
                    INDEX_UNITS[[ix]], lim)
    save_fig_dpi(pl, "11_maps",
                 sprintf("STEP08_03_MAP_%s_SEN_SLOPE_%s_%d_%d_FDRq%03d",
                         ix, PRIMARY_VARIANT, PRIMARY_START_YEAR, PRIMARY_END_YEAR,
                         round(FDR_Q * 100)))
  }
  
  ## ---- 08.3.7 combined maps, one per UNIT FAMILY ---------------------------
  ## A shared colour scale is valid only within a unit family.
  fam <- split(names(INDEX_UNITS), unname(INDEX_UNITS))
  for (u in names(fam)) {
    ixs <- intersect(fam[[u]], INDICES)
    dd  <- d[index %in% ixs]
    if (!nrow(dd)) next
    dd  <- droplevels(dd)
    lim <- sym_lim(dd$sen_plot)
    ncol_f <- min(3L, length(ixs))
    pl <- build_map(dd,
                    paste0("AMJ trends, indices measured in ", u),
                    u, lim, faceted = TRUE, ncol_f = ncol_f)
    save_fig_dpi(pl, "11_maps",
                 sprintf("STEP08_03_MAP_COMBINED_%s_SEN_SLOPE_%s_%d_%d",
                         gsub("[^A-Za-z0-9]", "", u), PRIMARY_VARIANT,
                         PRIMARY_START_YEAR, PRIMARY_END_YEAR),
                 w = 11, h = 4 + 3.2 * ceiling(length(ixs) / ncol_f))
  }
  
  ## ---- 08.3.8 map provenance table -----------------------------------------
  if (length(lim_tab)) {
    lim_dt <- rbindlist(lim_tab, use.names = TRUE)
    save_table(lim_dt, "11_maps",
               sprintf("STEP08_03_MAP_SCALES_AND_OVERLAYS_%s_%d_%d",
                       PRIMARY_VARIANT, PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  } else {
    lim_dt <- data.table()
  }
  
  ## ---- 08.3.9 QC ------------------------------------------------------------
  uchk <- d[, .(u = units[1]), by = index]
  qc_add("maps_units_consistent",
         all(d[, uniqueN(units), by = index]$V1 == 1L) &&
           all(uchk$u == INDEX_UNITS[as.character(uchk$index)]),
         "each index carries exactly one unit label matching the configuration")
  
  gcnt <- d[, .(n = uniqueN(grid_id)), by = index]
  qc_add("maps_grid_count_per_index",
         all(gcnt$n == EXP_GRIDS),
         sprintf("grid cells per mapped index: %s (expected %d)",
                 paste(sort(unique(gcnt$n)), collapse = "/"), EXP_GRIDS))
  
  qc_add("maps_no_duplicate_cells",
         d[, .N, by = .(index, grid_id)][N > 1L, .N] == 0L,
         "no index x grid_id pair is plotted more than once")
  
  qc_add("maps_single_variant_period",
         (!"variant" %in% names(d) || uniqueN(d$variant) == 1L) &&
           (!"period" %in% names(d) || uniqueN(d$period) == 1L),
         sprintf("mapped rows carry one variant (%s) and one period",
                 if ("variant" %in% names(d)) paste(unique(d$variant), collapse = ",") else "n/a"))
  
  qc_add("maps_stipple_matches_q",
         d[significant_fdr %in% TRUE & (is.na(q) | q > FDR_Q), .N] == 0L,
         sprintf("every stippled cell satisfies q <= %.2f", FDR_Q))
  
  qc_add("maps_undefined_not_stippled",
         d[trend_class %in% c("test_undefined", "insufficient_data") &
             significant_fdr %in% TRUE, .N] == 0L,
         "no test_undefined or insufficient_data cell is stippled")
  
  qc_add("maps_boundary_overlay",
         !is.null(bnd),
         if (is.null(bnd)) paste0("boundary NOT drawn; source: ", MAP_SHP)
         else sprintf("boundary drawn: %d vertices, %d rings, EPSG:4326",
                      nrow(bnd), uniqueN(bnd$ring)),
         severity = if (is.null(bnd)) "WARN" else "INFO")
  
  qc_add("maps_colour_scale_zero_centred",
         nrow(lim_dt) == 0L || lim_dt[, all(colour_limit > 0)],
         sprintf("symmetric limits +/- q%.0f of |slope|; squished cells: %d",
                 100 * MAP_LIM_PROB,
                 if (nrow(lim_dt)) sum(lim_dt$n_squished) else 0L),
         severity = "INFO")
  
  qc_add("maps_extent_consistent", TRUE,
         sprintf("all individual and combined maps use lon %.2f to %.2f, lat %.2f to %.2f, aspect %.2f",
                 MAP_XLIM[1], MAP_XLIM[2], MAP_YLIM[1], MAP_YLIM[2], MAP_ASPECT),
         severity = "INFO")
  
  qc_add("maps_resolution", TRUE,
         sprintf("PNG at %d dpi, PDF vector; native 1 deg tiles, no interpolation or regridding",
                 MAP_DPI),
         severity = "INFO")
}

## ===========================================================================
## 08.4 - RAW vs FDR
## ===========================================================================
if (SEC("08.4")) {
  log_head("08.4 | RAW p < 0.05 vs BENJAMINI-HOCHBERG FDR")
  mark_run("08.4")
  d <- get_all8()
  
  fdr <- d[, {
    nt <- sum(is.finite(p))
    nraw <- sum(significant_raw, na.rm = TRUE)
    nbh  <- sum(significant_fdr, na.rm = TRUE)
    .(n_tested = nt, n_raw_sig = nraw,
      expected_false_positives = round(RAW_ALPHA * nt, 2),
      raw_to_expected_ratio = round(nraw / max(RAW_ALPHA * nt, 1e-9), 2),
      n_fdr_sig = nbh,
      pct_retained_after_fdr = round(100 * nbh / max(nraw, 1L), 2),
      n_lost_to_fdr = nraw - nbh,
      pct_of_grids_raw = round(100 * nraw / max(nt, 1L), 2),
      pct_of_grids_fdr = round(100 * nbh / max(nt, 1L), 2),
      min_p = sfmt(mstat(p, min), 6), min_q = sfmt(mstat(q, min), 6))
  }, by = index]
  setorder(fdr, index)
  fdr[, `:=`(raw_alpha = RAW_ALPHA, fdr_q = FDR_Q,
             period = paste0(PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
             variant = PRIMARY_VARIANT)]
  save_table(fdr, "03_fdr",
             sprintf("STEP08_04_RAW_vs_FDR_%d_%d", PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  plt <- melt(fdr[, .(index, `Raw p<0.05` = n_raw_sig,
                      `Expected false` = expected_false_positives,
                      `BH q<=0.10` = n_fdr_sig)],
              id.vars = "index", variable.name = "category", value.name = "n_grids")
  plt[, category := factor(category, levels = c("Raw p<0.05", "Expected false", "BH q<=0.10"))]
  pl <- ggplot(plt, aes(index, n_grids, fill = category)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = c("Raw p<0.05" = "#9ECAE1",
                                 "Expected false" = "grey55",
                                 "BH q<=0.10" = "#B2182B"), name = NULL) +
    labs(title = "Grid-level significance before and after false-discovery-rate control",
         subtitle = sprintf("AMJ %d-%d | %s | %d grids per index | expected false = %.2f x N tested",
                            PRIMARY_START_YEAR, PRIMARY_END_YEAR, PRIMARY_VARIANT,
                            EXP_GRIDS, RAW_ALPHA),
         x = NULL, y = "Number of grids") + theme_pub
  save_fig(pl, "12_figures",
           sprintf("STEP08_04_FIG_RAW_vs_EXPECTED_vs_FDR_%d_%d",
                   PRIMARY_START_YEAR, PRIMARY_END_YEAR), w = 9, h = 5)
  
  for (i in seq_len(nrow(fdr))) {
    z <- fdr[i]
    log_msg(sprintf("  %-6s tested %3d | raw %3d | expected-false %5.2f | BH %3d | retained %5.1f%%",
                    z$index, z$n_tested, z$n_raw_sig, z$expected_false_positives,
                    z$n_fdr_sig, z$pct_retained_after_fdr))
  }
  qc_add("fdr_never_exceeds_raw", fdr[n_fdr_sig > n_raw_sig, .N] == 0L,
         sprintf("%d indices where BH exceeded raw at the same level",
                 fdr[n_fdr_sig > n_raw_sig, .N]))
}

## ===========================================================================
## 08.5 - SENSITIVITY  (MK_original vs HR_lag3 vs HR_all)
## ===========================================================================
if (SEC("08.5")) {
  log_head("08.5 | VARIANT SENSITIVITY (evidence, not a replacement for the primary result)")
  mark_run("08.5")
  v <- get_variants()
  pn <- unique(v$period)[grepl(as.character(PRIMARY_START_YEAR), unique(v$period))]
  if (!length(pn)) fatal("the primary period could not be identified in the all-variant files")
  vp <- v[period == pn[1] & variant %in% SENSITIVITY_VARIANTS]
  qc_add("sensitivity_input_rows",
         nrow(vp) == EXP_GRIDS * length(INDICES) * length(SENSITIVITY_VARIANTS),
         sprintf("%s rows = 283 x %d indices x %d variants",
                 format(nrow(vp), big.mark = ","), uniqueN(vp$index), uniqueN(vp$variant)))
  
  w <- dcast(vp[, .(index, grid_id, variant, sen_per_decade, q,
                    trend_class = as.character(trend_class))],
             index + grid_id ~ variant,
             value.var = c("sen_per_decade", "q", "trend_class"))
  vpairs <- list(c("HR_lag3", "HR_all"), c("HR_lag3", "MK_original"),
                 c("HR_all", "MK_original"))
  sens <- rbindlist(lapply(vpairs, function(pp) {
    a <- pp[1]; b <- pp[2]
    sa <- paste0("sen_per_decade_", a); sb <- paste0("sen_per_decade_", b)
    qa <- paste0("q_", a); qb <- paste0("q_", b)
    ca <- paste0("trend_class_", a); cb <- paste0("trend_class_", b)
    if (!all(c(sa, sb, qa, qb, ca, cb) %in% names(w))) return(NULL)
    w[, {
      siga <- get(qa) <= FDR_Q & !is.na(get(qa))
      sigb <- get(qb) <= FDR_Q & !is.na(get(qb))
      inter <- sum(siga & sigb); uni <- sum(siga | sigb)
      .(variant_a = a, variant_b = b, n_grids = .N,
        slope_correlation = sfmt(suppressWarnings(
          stats::cor(get(sa), get(sb), use = "complete.obs", method = "pearson"))),
        slope_spearman = sfmt(suppressWarnings(
          stats::cor(get(sa), get(sb), use = "complete.obs", method = "spearman"))),
        median_abs_slope_diff = sfmt(stats::median(abs(get(sa) - get(sb)), na.rm = TRUE), 6),
        max_abs_slope_diff = sfmt(mstat(abs(get(sa) - get(sb)), max), 6),
        n_class_agree = sum(get(ca) == get(cb), na.rm = TRUE),
        pct_class_agreement = round(100 * mean(get(ca) == get(cb), na.rm = TRUE), 2),
        n_sig_a = sum(siga), n_sig_b = sum(sigb), n_sig_overlap = inter,
        jaccard_similarity = if (uni == 0L) NA_real_ else round(inter / uni, 4),
        n_discordant_sig = sum(siga != sigb))
    }, by = index]
  }))
  setorder(sens, index, variant_a, variant_b)
  sens[, sensitivity_flag := fifelse(pct_class_agreement < 95, "METHOD_SENSITIVE", "stable")]
  save_table(sens, "04_sensitivity",
             sprintf("STEP08_05_VARIANT_SENSITIVITY_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  flagged <- sens[sensitivity_flag == "METHOD_SENSITIVE"]
  if (nrow(flagged)) {
    for (i in seq_len(nrow(flagged))) {
      z <- flagged[i]
      warn_msg(sprintf("METHOD-SENSITIVE: %s %s vs %s - class agreement %.1f%%, significant %d vs %d, Jaccard %s",
                       z$index, z$variant_a, z$variant_b, z$pct_class_agreement,
                       z$n_sig_a, z$n_sig_b, format(z$jaccard_similarity)))
    }
    save_table(flagged, "04_sensitivity",
               sprintf("STEP08_05_METHOD_SENSITIVE_INDICES_%d_%d",
                       PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  }
  qc_add("sensitivity_flags_reported", TRUE,
         sprintf("%d index-pair comparisons flagged METHOD_SENSITIVE (< 95%% class agreement)",
                 nrow(flagged)), severity = "INFO")
  
  ps <- sens[variant_a == "HR_lag3"]
  ylo <- min(80, suppressWarnings(min(ps$pct_class_agreement, na.rm = TRUE)) - 2)
  pl <- ggplot(ps, aes(index, pct_class_agreement, fill = variant_b)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(yintercept = 95, linetype = 2, colour = "grey30") +
    coord_cartesian(ylim = c(ylo, 100)) +
    scale_fill_brewer(palette = "Set2", name = "compared with") +
    labs(title = "Sensitivity of the trend classification to the autocorrelation treatment",
         subtitle = sprintf("Agreement with the primary variant %s; dashed line = 95%%",
                            PRIMARY_VARIANT),
         x = NULL, y = "Grids with an identical class (%)") + theme_pub
  save_fig(pl, "12_figures",
           sprintf("STEP08_05_FIG_VARIANT_SENSITIVITY_%d_%d",
                   PRIMARY_START_YEAR, PRIMARY_END_YEAR), w = 9, h = 5)
}

## ===========================================================================
## 08.6 - WSDI SPECIAL ANALYSIS
## ===========================================================================
if (SEC("08.6")) {
  log_head("08.6 | WSDI ZERO-INFLATION ANALYSIS")
  mark_run("08.6")
  s6 <- get_step06(); d <- get_all8()
  if (!"WSDI" %in% names(s6)) fatal("WSDI column not present in the STEP 06 file")
  
  ws <- s6[year >= PRIMARY_START_YEAR & year <= PRIMARY_END_YEAR]
  per_grid <- ws[, .(n_years = .N, n_valid = sum(!is.na(WSDI)),
                     n_zero_years = sum(WSDI == 0, na.rm = TRUE),
                     pct_zero_years = round(100 * mean(WSDI == 0, na.rm = TRUE), 2),
                     n_distinct_values = uniqueN(WSDI[!is.na(WSDI)]),
                     max_WSDI = mstat(WSDI, max),
                     mean_WSDI = sfmt(mean(WSDI, na.rm = TRUE), 3)), by = grid_id]
  wt <- d[index == "WSDI", .(grid_id, lat, lon, sen_per_decade, ci_lo_per_decade,
                             ci_hi_per_decade, p, q, trend_class, cf_status)]
  per_grid <- merge(per_grid, wt, by = "grid_id", all.x = TRUE)
  per_grid[, sen_exactly_zero := is.finite(sen_per_decade) & sen_per_decade == 0]
  save_table(per_grid, "05_wsdi",
             sprintf("STEP08_06_WSDI_ZERO_INFLATION_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  wsum <- data.table(
    metric = c("grid-years analysed", "grid-year zero fraction (%)",
               "grids with any non-zero WSDI year", "grids with Sen slope exactly zero",
               "grids FDR-significant", "grids test_undefined",
               "median Sen slope (days/decade)", "IQR of Sen slope",
               "maximum observed WSDI (days)"),
    value = c(nrow(ws), round(100 * mean(ws$WSDI == 0, na.rm = TRUE), 2),
              per_grid[n_zero_years < n_valid, .N],
              per_grid[sen_exactly_zero == TRUE, .N],
              d[index == "WSDI" & significant_fdr == TRUE, .N],
              d[index == "WSDI" & trend_class == "test_undefined", .N],
              sfmt(stats::median(wt$sen_per_decade, na.rm = TRUE)),
              sfmt(stats::IQR(wt$sen_per_decade, na.rm = TRUE)),
              mstat(ws$WSDI, max)))
  save_table(wsum, "05_wsdi",
             sprintf("STEP08_06_WSDI_SUMMARY_%d_%d", PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  p1 <- ggplot(per_grid, aes(pct_zero_years)) +
    geom_histogram(bins = 30, fill = "#4292C6", colour = "white") +
    labs(title = "WSDI zero-inflation across grids",
         subtitle = sprintf("Percentage of %d-%d grid-years with WSDI = 0",
                            PRIMARY_START_YEAR, PRIMARY_END_YEAR),
         x = "Grid-years with WSDI = 0 (%)", y = "Number of grids") + theme_pub
  p2 <- ggplot(per_grid[is.finite(sen_per_decade)], aes(sen_per_decade)) +
    geom_histogram(bins = 40, fill = "grey55", colour = "white") +
    geom_vline(xintercept = 0, colour = "#B2182B", linewidth = 0.6) +
    labs(title = "WSDI Theil-Sen slope distribution",
         subtitle = sprintf("%d of %d grids have a slope of exactly zero; %d grids are FDR-significant",
                            per_grid[sen_exactly_zero == TRUE, .N], nrow(per_grid),
                            d[index == "WSDI" & significant_fdr == TRUE, .N]),
         x = "Sen slope (days/decade)", y = "Number of grids") + theme_pub
  save_fig(p1, "12_figures",
           sprintf("STEP08_06_FIG_WSDI_ZERO_INFLATION_%d_%d",
                   PRIMARY_START_YEAR, PRIMARY_END_YEAR), w = 7.5, h = 4.5)
  save_fig(p2, "12_figures",
           sprintf("STEP08_06_FIG_WSDI_SLOPE_DISTRIBUTION_%d_%d",
                   PRIMARY_START_YEAR, PRIMARY_END_YEAR), w = 7.5, h = 4.5)
  
  log_msg(sprintf("  WSDI: %.1f%% zero grid-years; %d of %d grids with slope exactly 0; %d FDR-significant",
                  100 * mean(ws$WSDI == 0, na.rm = TRUE),
                  per_grid[sen_exactly_zero == TRUE, .N], nrow(per_grid),
                  d[index == "WSDI" & significant_fdr == TRUE, .N]))
  qc_add("wsdi_zero_inflation_reported", TRUE,
         "WSDI reported with its zero-inflation structure, not as an ordinary continuous index",
         severity = "INFO")
}

## ===========================================================================
## 08.7 - TX90p / TN10p AS PERCENTILE-FREQUENCY INDICES
## ===========================================================================
if (SEC("08.7")) {
  log_head("08.7 | PERCENTILE-FREQUENCY INDICES TX90p / TN10p")
  mark_run("08.7")
  d <- get_all8(); s5 <- get_step05()
  fx <- c("TX90p", "TN10p")
  
  ds <- rbindlist(lapply(fx, function(ix) {
    v <- s5[year >= PRIMARY_START_YEAR & year <= PRIMARY_END_YEAR][[ix]]
    data.table(index = ix, scope = "AMJ_season", unit = "%",
               N = sum(!is.na(v)), N_missing = sum(is.na(v)),
               Mean = sfmt(mean(v, na.rm = TRUE), 3),
               Median = sfmt(stats::median(v, na.rm = TRUE), 3),
               SD = sfmt(stats::sd(v, na.rm = TRUE), 3),
               Q1 = sfmt(stats::quantile(v, .25, na.rm = TRUE, names = FALSE), 3),
               Q3 = sfmt(stats::quantile(v, .75, na.rm = TRUE, names = FALSE), 3),
               Minimum = sfmt(mstat(v, min), 3), Maximum = sfmt(mstat(v, max), 3))
  }))
  trf <- d[index %in% fx, .(
    n_grids = .N, n_testable = sum(is.finite(p)),
    median_sen = sfmt(stats::median(sen_per_decade, na.rm = TRUE)),
    iqr_sen = sfmt(stats::IQR(sen_per_decade, na.rm = TRUE)),
    min_sen = sfmt(mstat(sen_per_decade, min)),
    max_sen = sfmt(mstat(sen_per_decade, max)),
    n_raw_sig = sum(significant_raw, na.rm = TRUE),
    n_fdr_sig = sum(significant_fdr, na.rm = TRUE),
    n_increase = sum(trend_class == "significant_increase"),
    n_decrease = sum(trend_class == "significant_decrease"),
    n_test_undefined = sum(trend_class == "test_undefined"),
    units = units[1]), by = index]
  freq <- merge(ds, trf, by = "index")
  save_table(freq, "06_frequency",
             sprintf("STEP08_07_FREQUENCY_INDICES_TX90p_TN10p_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  if ("bootstrap_applied" %in% names(s5)) {
    bc <- s5[year >= PRIMARY_START_YEAR & year <= PRIMARY_END_YEAR,
             .(n_grid_years = .N,
               estimator = fifelse(bootstrap_applied[1], "in_base_Zhang2005_bootstrap",
                                   "out_of_base_direct_fixed_thresholds"),
               mean_TX90p = sfmt(mean(TX90p, na.rm = TRUE), 3),
               mean_TN10p = sfmt(mean(TN10p, na.rm = TRUE), 3)),
             by = .(bootstrap_applied, year)][order(year)]
    save_table(bc, "06_frequency",
               sprintf("STEP08_07_TX90p_TN10p_BOOTSTRAP_PROVENANCE_%d_%d",
                       PRIMARY_START_YEAR, PRIMARY_END_YEAR))
    nb <- uniqueN(bc[bootstrap_applied == TRUE, year])
    nd <- uniqueN(bc[bootstrap_applied == FALSE, year])
    qc_add("percentile_estimator_provenance_reported", TRUE,
           sprintf("%d in-base bootstrapped years and %d out-of-base direct years inside the primary window",
                   nb, nd), severity = "INFO")
    log_msg(sprintf("  Estimator provenance: %d in-base bootstrapped years, %d direct years.",
                    nb, nd))
  } else {
    warn_msg("column 'bootstrap_applied' absent from the STEP 05 file; provenance not reported.")
  }
  log_msg("  TX90p / TN10p summarised; no percentile threshold or bootstrap recomputed here.")
}

## ===========================================================================
## 08.8 - SCIENTIFIC SYNTHESIS
## ===========================================================================
if (SEC("08.8")) {
  log_head("08.8 | SCIENTIFIC SYNTHESIS ACROSS DIMENSIONS")
  mark_run("08.8")
  d <- get_all8()
  
  syn <- d[, {
    nt <- sum(is.finite(p))
    nsi <- sum(trend_class == "significant_increase")
    nsd <- sum(trend_class == "significant_decrease")
    .(median_sen = sfmt(stats::median(sen_per_decade, na.rm = TRUE)),
      iqr_sen = sfmt(stats::IQR(sen_per_decade, na.rm = TRUE)),
      median_ci_lo = sfmt(stats::median(ci_lo_per_decade, na.rm = TRUE)),
      median_ci_hi = sfmt(stats::median(ci_hi_per_decade, na.rm = TRUE)),
      n_testable = nt, n_fdr_sig = nsi + nsd,
      pct_fdr_sig = round(100 * (nsi + nsd) / max(nt, 1L), 2),
      n_increase = nsi, n_decrease = nsd,
      direction = fifelse(nsi + nsd == 0L, "no field-significant direction",
                          fifelse(nsi > nsd, "predominantly increasing",
                                  fifelse(nsd > nsi, "predominantly decreasing", "mixed"))),
      units = units[1])
  }, by = index]
  syn[, dimension := INDEX_DIMENSION[as.character(index)]]
  
  sens_path <- file.path(OUTPUT_DIR, "04_sensitivity",
                         sprintf("STEP08_05_VARIANT_SENSITIVITY_%d_%d.csv",
                                 PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  if (file.exists(sens_path)) {
    sn <- fread(sens_path)
    fl <- sn[variant_a == "HR_lag3",
             .(min_agreement = min(pct_class_agreement, na.rm = TRUE)), by = index]
    syn <- merge(syn, fl, by = "index", all.x = TRUE)
    syn[, robustness_flag := fifelse(is.na(min_agreement), "not assessed",
                                     fifelse(min_agreement < 95, "METHOD_SENSITIVE", "robust"))]
  } else {
    syn[, `:=`(min_agreement = NA_real_, robustness_flag = "not assessed (run 08.5)")]
    warn_msg("08.5 output not found; the synthesis robustness flag is 'not assessed'.")
  }
  setcolorder(syn, c("index", "dimension", "units", "median_sen", "iqr_sen",
                     "median_ci_lo", "median_ci_hi", "n_testable", "n_fdr_sig",
                     "pct_fdr_sig", "n_increase", "n_decrease", "direction",
                     "min_agreement", "robustness_flag"))
  setorder(syn, dimension, index)
  save_table(syn, "07_synthesis",
             sprintf("STEP08_08_SYNTHESIS_BY_DIMENSION_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  pd <- copy(d)
  pd[, unit_family := INDEX_UNITS[as.character(index)]]
  pl <- ggplot(pd, aes(x = index, y = sen_per_decade)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_boxplot(outlier.size = 0.4, width = 0.55, fill = "grey92") +
    geom_point(data = pd[significant_fdr == TRUE],
               aes(x = index, y = sen_per_decade), colour = "#B2182B",
               size = 0.9, alpha = 0.8,
               position = position_jitter(width = 0.12, height = 0)) +
    facet_wrap(~ unit_family, scales = "free", nrow = 1) +
    labs(title = "Trend magnitude and field significance by index",
         subtitle = sprintf("Boxes = across-grid Sen slope distribution; red points = BH-FDR q <= %.2f. A larger slope does not imply greater significance.",
                            FDR_Q),
         caption = sprintf("AMJ %d-%d | %s | %d grids per index",
                           PRIMARY_START_YEAR, PRIMARY_END_YEAR, PRIMARY_VARIANT, EXP_GRIDS),
         x = NULL, y = "Theil-Sen slope per decade") + theme_pub
  save_fig(pl, "12_figures",
           sprintf("STEP08_08_FIG_SYNTHESIS_SLOPE_AND_SIGNIFICANCE_%d_%d",
                   PRIMARY_START_YEAR, PRIMARY_END_YEAR), w = 11, h = 5)
  qc_add("synthesis_dimensions_assigned", all(!is.na(syn$dimension)),
         sprintf("dimensions: %s", paste(sort(unique(syn$dimension)), collapse = ", ")))
}

## ===========================================================================
## 08.9 - SEASONAL PROGRESSION  (climatology only, never a trend)
## ===========================================================================
if (SEC("08.9")) {
  log_head("08.9 | SEASONAL PROGRESSION April -> May -> June (CLIMATOLOGY)")
  mark_run("08.9")
  mon <- if (!is.null(.CACHE$monthly)) .CACHE$monthly else {
    dl <- get_daily_amj()
    dl[, .(TXx = if (all(is.na(tmax))) NA_real_ else max(tmax, na.rm = TRUE),
           TXn = if (all(is.na(tmax))) NA_real_ else min(tmax, na.rm = TRUE),
           TNx = if (all(is.na(tmin))) NA_real_ else max(tmin, na.rm = TRUE),
           TNn = if (all(is.na(tmin))) NA_real_ else min(tmin, na.rm = TRUE),
           DTR = { dv <- (tmax - tmin); dv <- dv[!is.na(dv)]
           if (!length(dv)) NA_real_ else mean(dv) }),
       by = .(grid_id, year, month)]
  }
  .CACHE$monthly <- mon
  
  prog <- rbindlist(lapply(MONTHLY_DERIVABLE, function(ix) {
    z <- mon[, .(domain_mean = sfmt(mean(get(ix), na.rm = TRUE), 3),
                 domain_median = sfmt(stats::median(get(ix), na.rm = TRUE), 3),
                 domain_sd = sfmt(stats::sd(get(ix), na.rm = TRUE), 3)), by = month]
    z[, index := ix]
    z[, month_num := month]
    z[, month := MONTH_LABELS[as.character(month_num)]]
    z[, .(index, month_num, month, domain_mean, domain_median, domain_sd)]
  }))
  setorder(prog, index, month_num)
  
  peak <- prog[, .(peak_month = month[which.max(domain_mean)],
                   peak_value = max(domain_mean),
                   april = domain_mean[month_num == 4L],
                   may = domain_mean[month_num == 5L],
                   june = domain_mean[month_num == 6L]), by = index]
  peak[, `:=`(may_minus_april = sfmt(may - april, 3),
              june_minus_may = sfmt(june - may, 3),
              june_minus_april = sfmt(june - april, 3))]
  
  gridpeak <- rbindlist(lapply(MONTHLY_DERIVABLE, function(ix) {
    gm <- mon[, .(v = mean(get(ix), na.rm = TRUE)), by = .(grid_id, month)]
    pk <- gm[, .(peak_month = MONTH_LABELS[as.character(month[which.max(v)])]), by = grid_id]
    out <- pk[, .(n_grids = .N), by = peak_month]
    out[, index := ix]
    out[]
  }))
  gridpeak[, pct_grids := round(100 * n_grids / EXP_GRIDS, 2)]
  setcolorder(gridpeak, c("index", "peak_month", "n_grids", "pct_grids"))
  setorder(gridpeak, index, peak_month)
  
  save_table(prog, "08_seasonal",
             sprintf("STEP08_09_SEASONAL_PROGRESSION_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  save_table(peak, "08_seasonal",
             sprintf("STEP08_09_SEASONAL_PEAK_MONTH_DOMAIN_MEAN_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  save_table(gridpeak, "08_seasonal",
             sprintf("STEP08_09_SEASONAL_PEAK_MONTH_PER_GRID_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  pl <- ggplot(prog, aes(factor(month, levels = unname(MONTH_LABELS)),
                         domain_mean, group = index, colour = index)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~ index, scales = "free_y", nrow = 2) +
    guides(colour = "none") +
    labs(title = "Seasonal progression of AMJ temperature extremes (climatology)",
         subtitle = sprintf("Domain mean over %d grids, %d-%d. This is the April-to-June seasonal march, NOT a long-term trend.",
                            EXP_GRIDS, PRIMARY_START_YEAR, PRIMARY_END_YEAR),
         x = NULL, y = "degC") + theme_pub
  save_fig(pl, "12_figures",
           sprintf("STEP08_09_FIG_SEASONAL_PROGRESSION_%d_%d",
                   PRIMARY_START_YEAR, PRIMARY_END_YEAR), w = 10, h = 6)
  
  for (i in seq_len(nrow(peak)))
    log_msg(sprintf("  %-5s peak month = %-5s (Apr %.2f, May %.2f, Jun %.2f)",
                    peak$index[i], peak$peak_month[i], peak$april[i],
                    peak$may[i], peak$june[i]))
  qc_add("seasonal_is_labelled_climatology", TRUE,
         "seasonal progression outputs are explicitly labelled climatology, not trend",
         severity = "INFO")
}

## ===========================================================================
## 08.10 - REGIONAL ANALYSIS  (PENDING unless an authoritative definition exists)
## ===========================================================================
if (SEC("08.10")) {
  log_head("08.10 | REGIONAL ANALYSIS")
  mark_run("08.10")
  if (is.na(REGION_FILE) || !file.exists(REGION_FILE)) {
    req <- data.table(
      section = "08.10", status = "PENDING",
      reason = paste("No authoritative regional classification exists in this pipeline.",
                     "Regions were NOT invented after inspecting the results.",
                     "NOTE: the state boundary shapefile used in 08.3 is a cartographic",
                     "overlay and is NOT a regional classification for grid aggregation."),
      required_input = "CSV with columns grid_id and region, covering all 283 grids",
      acceptable_sources = paste("IMD official meteorological subdivisions;",
                                 "IMD homogeneous temperature zones;",
                                 "a published, citable regionalisation with a defined boundary file"),
      must_record = "source name, version, citation, and the boundary file used",
      config_key = "REGION_FILE")
    save_table(req, "09_regional", "STEP08_10_REGIONAL_ANALYSIS_PENDING_REQUIREMENTS")
    qc_add("regional_analysis", FALSE,
           "PENDING: no authoritative regional classification supplied; regions were not invented",
           severity = "WARN")
    log_msg("  08.10 PENDING. Supply REGION_FILE to activate.")
  } else {
    rg <- fread(REGION_FILE)
    if (!all(c("grid_id", "region") %in% names(rg)))
      fatal("REGION_FILE must contain columns 'grid_id' and 'region'")
    rg[, grid_id := as.character(grid_id)]
    d <- get_all8()
    miss <- setdiff(unique(d$grid_id), rg$grid_id)
    qc_add("regional_covers_all_grids", length(miss) == 0L,
           sprintf("%d of %d grids unassigned", length(miss), EXP_GRIDS))
    dr <- merge(d, rg[, .(grid_id, region)], by = "grid_id", all.x = TRUE)
    rsum <- dr[, .(n_grids = .N, n_testable = sum(is.finite(p)),
                   median_sen = sfmt(stats::median(sen_per_decade, na.rm = TRUE)),
                   iqr_sen = sfmt(stats::IQR(sen_per_decade, na.rm = TRUE)),
                   n_raw_sig = sum(significant_raw, na.rm = TRUE),
                   n_fdr_sig = sum(significant_fdr, na.rm = TRUE),
                   pct_fdr_sig = round(100 * sum(significant_fdr, na.rm = TRUE) /
                                         max(sum(is.finite(p)), 1L), 2),
                   units = units[1]), by = .(region, index)]
    setorder(rsum, region, index)
    rsum[, region_source := basename(REGION_FILE)]
    save_table(rsum, "09_regional",
               sprintf("STEP08_10_REGIONAL_TREND_SUMMARY_%d_%d",
                       PRIMARY_START_YEAR, PRIMARY_END_YEAR))
    log_msg("  Regional summary written for ", uniqueN(rsum$region), " regions.")
  }
}

## ===========================================================================
## 08.11 - UNCERTAINTY  (CI-based; distinct from FDR significance)
## ===========================================================================
if (SEC("08.11")) {
  log_head("08.11 | SEN-SLOPE UNCERTAINTY (95% CI)")
  mark_run("08.11")
  d <- get_all8()
  dd <- copy(d)
  dd[, ci_width := ci_hi_per_decade - ci_lo_per_decade]
  dd[, ci_status := fifelse(!is.finite(ci_lo_per_decade), "undefined",
                            fifelse(ci_lo_per_decade > 0, "entirely_above_zero",
                                    fifelse(ci_hi_per_decade < 0, "entirely_below_zero", "crosses_zero")))]
  
  unc <- dd[, {
    nci <- sum(is.finite(ci_width))
    .(n_grids = .N, n_with_ci = nci,
      median_ci_width = sfmt(stats::median(ci_width, na.rm = TRUE)),
      iqr_ci_width = sfmt(stats::IQR(ci_width, na.rm = TRUE)),
      min_ci_width = sfmt(mstat(ci_width, min)),
      max_ci_width = sfmt(mstat(ci_width, max)),
      n_ci_above_zero = sum(ci_status == "entirely_above_zero"),
      n_ci_below_zero = sum(ci_status == "entirely_below_zero"),
      n_ci_crosses_zero = sum(ci_status == "crosses_zero"),
      pct_ci_above_zero = round(100 * sum(ci_status == "entirely_above_zero") / max(nci, 1L), 2),
      pct_ci_below_zero = round(100 * sum(ci_status == "entirely_below_zero") / max(nci, 1L), 2),
      pct_ci_crosses_zero = round(100 * sum(ci_status == "crosses_zero") / max(nci, 1L), 2),
      n_fdr_sig = sum(significant_fdr, na.rm = TRUE),
      units = units[1])
  }, by = index]
  setorder(unc, index)
  unc[, note := paste("CI exclusion of zero is a per-grid statement and is NOT equivalent to",
                      "FDR significance, which controls the false discovery rate across all",
                      "283 grids simultaneously. The two columns are not expected to agree.")]
  unc[, record_length_note := sprintf(
    "A %d-year record limits the detectable slope; wide CIs reflect record length as much as physical variability.",
    PRIMARY_END_YEAR - PRIMARY_START_YEAR + 1L)]
  save_table(unc, "10_uncertainty",
             sprintf("STEP08_11_SEN_SLOPE_UNCERTAINTY_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  save_table(dd[, .(index, grid_id, lat, lon, sen_per_decade, ci_lo_per_decade,
                    ci_hi_per_decade, ci_width, ci_status, q, significant_fdr, units)],
             "10_uncertainty",
             sprintf("STEP08_11_SEN_SLOPE_CI_BY_GRID_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  
  for (i in seq_len(nrow(unc))) {
    z <- unc[i]
    log_msg(sprintf("  %-6s median CI width %8.4f %-12s | CI>0 %5.1f%% | CI<0 %5.1f%% | crosses 0 %5.1f%% | FDR sig %d",
                    z$index, z$median_ci_width, z$units, z$pct_ci_above_zero,
                    z$pct_ci_below_zero, z$pct_ci_crosses_zero, z$n_fdr_sig))
  }
  qc_add("ci_percentages_sum_to_100",
         unc[, all(abs(pct_ci_above_zero + pct_ci_below_zero +
                         pct_ci_crosses_zero - 100) < 0.05)],
         "CI status percentages sum to 100 for every index")
}

## ===========================================================================
## 08.12 - FINAL PUBLICATION FIGURE SET (registry; no redundant plots)
## ===========================================================================
if (SEC("08.12")) {
  log_head("08.12 | FINAL PUBLICATION FIGURE SET")
  mark_run("08.12")
  figs <- data.table(
    figure = c("Fig 1", "Fig 2", "Fig 3", "Fig 4", "Fig 5", "Fig 6"),
    content = c("Seasonal progression, April-May-June climatology",
                "Raw vs expected-false vs FDR-significant counts",
                "Eight-index Sen slope maps with FDR stippling and state boundaries",
                "Synthesis: slope distribution with FDR significance by dimension",
                "WSDI zero-inflation and slope distribution",
                "Variant sensitivity (HR_lag3 vs HR_all vs MK_original)"),
    produced_by = c("08.9", "08.4", "08.3", "08.8", "08.6", "08.5"),
    file_stub = c(sprintf("STEP08_09_FIG_SEASONAL_PROGRESSION_%d_%d",
                          PRIMARY_START_YEAR, PRIMARY_END_YEAR),
                  sprintf("STEP08_04_FIG_RAW_vs_EXPECTED_vs_FDR_%d_%d",
                          PRIMARY_START_YEAR, PRIMARY_END_YEAR),
                  sprintf("STEP08_03_MAP_COMBINED_.*_SEN_SLOPE_%s_%d_%d",
                          PRIMARY_VARIANT, PRIMARY_START_YEAR, PRIMARY_END_YEAR),
                  sprintf("STEP08_08_FIG_SYNTHESIS_SLOPE_AND_SIGNIFICANCE_%d_%d",
                          PRIMARY_START_YEAR, PRIMARY_END_YEAR),
                  sprintf("STEP08_06_FIG_WSDI_ZERO_INFLATION_%d_%d",
                          PRIMARY_START_YEAR, PRIMARY_END_YEAR),
                  sprintf("STEP08_05_FIG_VARIANT_SENSITIVITY_%d_%d",
                          PRIMARY_START_YEAR, PRIMARY_END_YEAR)),
    role = c("Establishes the seasonal state and quarantines it from the trend analysis",
             "Justifies multiplicity control before any map is shown",
             "The principal spatial result",
             "Separates magnitude from significance across dimensions",
             "Explains the WSDI null as a property of the index",
             "Demonstrates the primary result is not an artefact of the lag choice"))
  save_table(figs, "12_figures",
             sprintf("STEP08_12_FINAL_FIGURE_REGISTRY_%d_%d",
                     PRIMARY_START_YEAR, PRIMARY_END_YEAR))
  allf <- list.files(OUTPUT_DIR, recursive = TRUE)
  present <- vapply(figs$file_stub, function(st) any(grepl(st, allf)), logical(1))
  qc_add("final_figures_present", all(present),
         sprintf("%d of %d registered figures found on disk (missing sections were not run)",
                 sum(present), length(present)),
         severity = if (identical(RUN_SECTION, "ALL")) "FATAL" else "WARN")
  log_msg("  Figure registry written; ", sum(present), " of ", length(present), " present.")
}

## ===========================================================================
## FINAL QC
## ===========================================================================
log_head("FINAL QC")

if (file.exists(TREND_FILE)) {
  d <- get_all8()
  qc_add("qc_index_count", uniqueN(d$index) == length(INDICES),
         sprintf("%d indices: %s", uniqueN(d$index),
                 paste(sort(unique(as.character(d$index))), collapse = ", ")))
  qc_add("qc_row_count", nrow(d) == EXP_GRIDS * length(INDICES),
         sprintf("%s rows vs expected %s", format(nrow(d), big.mark = ","),
                 format(EXP_GRIDS * length(INDICES), big.mark = ",")))
  qc_add("qc_grids_per_index", d[, .N, by = index][, all(N == EXP_GRIDS)],
         sprintf("grid counts per index: %s",
                 paste(unique(d[, .N, by = index]$N), collapse = ", ")))
  qc_add("qc_no_duplicate_grid_index", !anyDuplicated(d, by = c("index", "grid_id")),
         sprintf("%d duplicated (index, grid_id) rows",
                 sum(duplicated(d, by = c("index", "grid_id")))))
  qc_add("qc_period_correct",
         all(grepl(as.character(PRIMARY_START_YEAR), d$period)) &&
           all(grepl(as.character(PRIMARY_END_YEAR), d$period)),
         sprintf("period label(s): %s", paste(unique(d$period), collapse = ", ")))
  qc_add("qc_variant_correct", all(d$variant == PRIMARY_VARIANT),
         sprintf("variant(s): %s", paste(unique(d$variant), collapse = ", ")))
  qc_add("qc_q_ge_p", d[is.finite(p) & is.finite(q) & q < p - 1e-12, .N] == 0L,
         sprintf("%d rows with q < p", d[is.finite(p) & is.finite(q) & q < p - 1e-12, .N]))
  qc_add("qc_p_q_in_unit_interval",
         d[is.finite(p) & (p < 0 | p > 1), .N] == 0L &&
           d[is.finite(q) & (q < 0 | q > 1), .N] == 0L,
         "all p and q values lie in [0, 1]")
  qc_add("qc_ci_ordered",
         d[is.finite(ci_lo_per_decade) & ci_lo_per_decade > ci_hi_per_decade + 1e-9, .N] == 0L,
         sprintf("%d inverted confidence intervals",
                 d[is.finite(ci_lo_per_decade) &
                     ci_lo_per_decade > ci_hi_per_decade + 1e-9, .N]))
  qc_add("qc_ci_brackets_slope",
         d[is.finite(ci_lo_per_decade) & is.finite(sen_per_decade) &
             (sen_per_decade < ci_lo_per_decade - 1e-9 |
                sen_per_decade > ci_hi_per_decade + 1e-9), .N] == 0L,
         "every finite CI brackets its own Sen slope")
  qc_add("qc_significance_matches_slope_sign",
         d[trend_class == "significant_increase" & sen_per_decade <= 0, .N] == 0L &&
           d[trend_class == "significant_decrease" & sen_per_decade >= 0, .N] == 0L,
         "significant increase/decrease classes agree with the slope sign")
  qc_add("qc_significant_only_below_q",
         d[significant_fdr == TRUE & (is.na(q) | q > FDR_Q), .N] == 0L,
         sprintf("%d rows flagged significant without q <= %.2f",
                 d[significant_fdr == TRUE & (is.na(q) | q > FDR_Q), .N], FDR_Q))
  qc_add("qc_test_undefined_has_no_pq",
         d[trend_class == "test_undefined" & (!is.na(p) | !is.na(q)), .N] == 0L,
         sprintf("%d test_undefined rows carry p or q",
                 d[trend_class == "test_undefined" & (!is.na(p) | !is.na(q)), .N]))
  qc_add("qc_coordinates_complete", !anyNA(d$lat) && !anyNA(d$lon),
         sprintf("%d missing coordinates", sum(is.na(d$lat)) + sum(is.na(d$lon))))
  uu <- d[, .(u = units[1]), by = index]
  qc_add("qc_units_match_config", all(uu$u == INDEX_UNITS[as.character(uu$index)]),
         "unit labels in the trend file match the configuration")
  nud <- d[trend_class == "test_undefined", .N]
  if (nud > 0L)
    qc_add("qc_test_undefined_present", FALSE,
           sprintf("%d grid-index rows have an undefined test (Hamed-Rao CF <= 0); map as a separate class",
                   nud), severity = "WARN")
} else {
  qc_add("qc_trend_file_present", FALSE, paste0("TREND_FILE not found: ", TREND_FILE))
}

outdt <- if (length(.OUTFILES)) unique(rbindlist(.OUTFILES), by = "file") else
  data.table(file = character(0), exists = logical(0), size_kb = numeric(0))
if (nrow(outdt)) outdt[, exists := file.exists(file)]
qc_add("qc_all_declared_outputs_written",
       nrow(outdt) == 0L || all(outdt$exists),
       sprintf("%d of %d declared output files exist on disk",
               sum(outdt$exists), nrow(outdt)),
       severity = if (OVERWRITE) "FATAL" else "WARN")

qc <- rbindlist(.QC)
n_fatal <- qc[status == "FATAL", .N]
n_warn  <- qc[status == "WARN", .N]
STEP08_STATUS <- if (n_fatal > 0L) "FAIL" else "PASS"

qc[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
          script = SCRIPT_NAME, script_version = SCRIPT_VERSION,
          run_section = RUN_SECTION, step08_status = STEP08_STATUS)]
f_qc <- file.path(OUTPUT_DIR, "13_qc", "STEP08_FINAL_QC.csv")
fwrite(qc, f_qc); register_out(f_qc)

qc_txt <- c(
  "STEP 08 FINAL QC REPORT", strrep("=", 78),
  paste0("Script          : ", SCRIPT_NAME, " v", SCRIPT_VERSION),
  paste0("Run             : ", format(t_start, format = "%Y-%m-%d %H:%M:%S"), " UTC"),
  paste0("RUN_SECTION     : ", RUN_SECTION),
  paste0("Sections run    : ", paste(sort(.SECTIONS_RUN), collapse = ", ")),
  paste0("Primary design  : ", PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR,
         " | ", PRIMARY_VARIANT, " | BH q <= ", FDR_Q,
         " | completeness >= ", COMPLETENESS_THRESHOLD * 100, "%"),
  paste0("State boundary  : ", ifelse(is.na(SHAPEFILE), "<not supplied>", SHAPEFILE)),
  paste0("  loaded=", .SHP_INFO$loaded, " | CRS in=", .SHP_INFO$crs_in,
         " -> out=", .SHP_INFO$crs_out, " | transformed=", .SHP_INFO$transformed,
         " | layers added=", .SHP_INFO$layers_added),
  "",
  paste0("Checks          : ", nrow(qc), "  (PASS ", qc[status == "PASS", .N],
         " | WARN ", n_warn, " | FATAL ", n_fatal, " | INFO ",
         qc[status == "INFO", .N], ")"),
  paste0("STEP 08 STATUS  : ", STEP08_STATUS), "",
  "FATAL failures:",
  if (n_fatal) paste0("  - ", qc[status == "FATAL", check], ": ",
                      qc[status == "FATAL", detail]) else "  none",
  "", "Warnings:",
  if (length(.WARNINGS)) paste0("  - ", .WARNINGS) else "  none",
  "", "Known, deliberate limitations of this run:",
  "  - No monthly index file exists; monthly climatology is derived in STEP 08",
  "    for TXx/TXn/TNx/TNn/DTR only, and is never used for trends.",
  "  - Monthly WSDI is BLOCKED: it would require a month reset that contradicts",
  "    the locked continuous-91-day definition.",
  "  - Monthly TX90p/TN10p, if produced, are DIRECT exceedance only and are not",
  "    the published bootstrapped seasonal indices.",
  "  - Regional analysis is PENDING unless REGION_FILE is supplied. The state",
  "    boundary shapefile is a cartographic overlay, not a regionalisation.",
  strrep("=", 78))
f_qct <- file.path(OUTPUT_DIR, "13_qc", "STEP08_FINAL_QC.txt")
writeLines(qc_txt, f_qct); register_out(f_qct)

log_msg("QC checks        : ", nrow(qc), " (PASS ", qc[status == "PASS", .N],
        " | WARN ", n_warn, " | FATAL ", n_fatal, ")")
log_msg("STEP 08 STATUS   : ", STEP08_STATUS)

## ===========================================================================
## RUN LOG AND REPRODUCIBILITY MANIFEST
## ===========================================================================
log_head("RUN LOG AND REPRODUCIBILITY MANIFEST")

inputs <- data.table(
  role = c("TREND_FILE", "VARIANTS_FILE_A", "VARIANTS_FILE_B", "STEP06_FILE",
           "STEP05_FILE", "DAILY_FILE", "THRESH_FILE", "SHAPEFILE", "REGION_FILE",
           "MONTHLY_INDEX_FILE"),
  path = c(TREND_FILE, VARIANTS_FILE_A, VARIANTS_FILE_B, STEP06_FILE,
           STEP05_FILE, DAILY_FILE, THRESH_FILE,
           ifelse(is.na(SHAPEFILE), "<not supplied>", SHAPEFILE),
           ifelse(is.na(REGION_FILE), "<not supplied>", REGION_FILE),
           "<does not exist in this pipeline>"))
inputs[, exists := file.exists(path)]
inputs[, size_kb := ifelse(exists, round(file.size(path) / 1024, 1), NA_real_)]
inputs[, modified := ifelse(exists,
                            format(file.mtime(path), format = "%Y-%m-%d %H:%M:%S"),
                            NA_character_)]
inputs[, md5 := NA_character_]
for (i in which(inputs$exists)) {
  h <- try(unname(tools::md5sum(inputs$path[i])), silent = TRUE)
  if (!inherits(h, "try-error")) inputs[i, md5 := h]
}
f_in <- file.path(OUTPUT_DIR, "14_manifest", "STEP08_INPUT_FILE_MANIFEST.csv")
fwrite(inputs, f_in); register_out(f_in)

outdt2 <- if (length(.OUTFILES)) unique(rbindlist(.OUTFILES), by = "file") else
  data.table(file = character(0), exists = logical(0), size_kb = numeric(0))
if (nrow(outdt2)) {
  outdt2[, exists := file.exists(file)]
  outdt2[, size_kb := ifelse(exists, round(file.size(file) / 1024, 1), NA_real_)]
}
f_out <- file.path(OUTPUT_DIR, "14_manifest", "STEP08_OUTPUT_FILE_MANIFEST.csv")
fwrite(outdt2, f_out); register_out(f_out)

pkgs <- c("data.table", "arrow", "ggplot2",
          if (HAS_XLSX) "writexl", if (HAS_SF) "sf")
pkgtxt <- vapply(pkgs, function(p) paste0(p, " ", as.character(utils::packageVersion(p))),
                 character(1))

man <- c(
  "STEP 08 REPRODUCIBILITY MANIFEST", strrep("=", 78),
  paste0("Script            : ", SCRIPT_NAME),
  paste0("Script version    : ", SCRIPT_VERSION),
  paste0("Run start (UTC)   : ", format(t_start, format = "%Y-%m-%d %H:%M:%S")),
  paste0("Run end   (UTC)   : ", format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")),
  paste0("Elapsed           : ",
         sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs")))),
  paste0("Software          : ", R.version.string, " on ", R.version$platform),
  paste0("Packages          : ", paste(pkgtxt, collapse = " | ")),
  "", "LOCKED PARAMETERS",
  paste0("  Primary period        : ", PRIMARY_START_YEAR, "-", PRIMARY_END_YEAR),
  paste0("  Context period        : ", CONTEXT_START_YEAR, "-", CONTEXT_END_YEAR),
  paste0("  Primary variant       : ", PRIMARY_VARIANT),
  paste0("  Sensitivity variants  : ", paste(SENSITIVITY_VARIANTS, collapse = ", ")),
  paste0("  FDR threshold q       : ", FDR_Q),
  paste0("  Raw alpha (reported)  : ", RAW_ALPHA),
  paste0("  Completeness rule     : >= ", COMPLETENESS_THRESHOLD * 100, "%"),
  paste0("  Grids                 : ", EXP_GRIDS),
  paste0("  Indices               : ", paste(INDICES, collapse = ", ")),
  paste0("  Months (climatology)  : ", paste(unname(MONTH_LABELS), collapse = ", ")),
  paste0("  CF <= 0 policy        : test_undefined; CF never floored"),
  "", "CARTOGRAPHY",
  paste0("  State boundary file   : ", ifelse(is.na(SHAPEFILE), "<not supplied>", SHAPEFILE)),
  paste0("  Boundary required     : ", SHAPEFILE_REQUIRED),
  paste0("  Boundary loaded       : ", .SHP_INFO$loaded),
  paste0("  Boundary features     : ", .SHP_INFO$n_features),
  paste0("  CRS as read           : ", .SHP_INFO$crs_in),
  paste0("  CRS used for plotting : ", .SHP_INFO$crs_out,
         "  (transformed = ", .SHP_INFO$transformed, ")"),
  paste0("  Invalid geoms repaired: ", .SHP_INFO$invalid_before),
  paste0("  Maps with boundary    : ", .SHP_INFO$layers_added),
  paste0("  Draw method           : geom_path over sf::st_coordinates, so the locked",
         " coord_fixed() layout is preserved"),
  "", "EXECUTION",
  paste0("  RUN_SECTION           : ", RUN_SECTION),
  paste0("  Sections executed     : ", paste(sort(.SECTIONS_RUN), collapse = ", ")),
  paste0("  Output directory      : ", OUTPUT_DIR),
  paste0("  Output files written  : ", nrow(outdt2)),
  paste0("  MAKE_CSV/XLSX/PNG/PDF : ", MAKE_CSV, "/", MAKE_XLSX, "/", MAKE_PNG,
         "/", MAKE_PDF, "  DPI=", DPI, "  OVERWRITE=", OVERWRITE),
  "", "INPUT FILES (md5 where available)",
  paste0("  ", inputs$role, " : ", inputs$path,
         "  [", ifelse(inputs$exists, "present", "ABSENT"), "]",
         ifelse(is.na(inputs$md5), "", paste0("  md5=", inputs$md5)),
         ifelse(is.na(inputs$modified), "", paste0("  modified=", inputs$modified))),
  "", paste0("WARNINGS (", length(.WARNINGS), ")"),
  if (length(.WARNINGS)) paste0("  - ", .WARNINGS) else "  none",
  "", paste0("QC STATUS         : ", STEP08_STATUS,
             "   (FATAL ", n_fatal, " | WARN ", n_warn, ")"),
  "", "STATEMENT",
  "  STEP 08 performed no statistical re-estimation. All trend statistics were",
  "  read from the locked STEP 07 / 07A / 07B outputs. No percentile threshold,",
  "  bootstrap correction, p-value or q-value was recomputed in this script.",
  "  The v1.1.0 revision added only the official state boundary overlay to the",
  "  map-rendering block; no scale, axis, label, stipple, layout, dimension,",
  "  font, title or panel arrangement was altered.",
  strrep("=", 78))
f_man <- file.path(OUTPUT_DIR, "14_manifest", "STEP08_REPRODUCIBILITY_MANIFEST.txt")
writeLines(man, f_man)

.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
f_log <- file.path(OUTPUT_DIR, "14_manifest", "STEP08_RUN_LOG.txt")
writeLines(.LOG, f_log)

log_head("STEP 08 COMPLETE")
log_msg("Sections executed : ", paste(sort(.SECTIONS_RUN), collapse = ", "))
log_msg("Outputs written   : ", nrow(outdt2))
log_msg("Warnings          : ", length(.WARNINGS))
log_msg("STEP 08 STATUS    : ", STEP08_STATUS)
log_msg("QC report         : ", f_qc)
log_msg("Manifest          : ", f_man)
log_msg("Run log           : ", f_log)
cat("\nSTEP 08 STATUS:", STEP08_STATUS, "\n")

if (identical(STEP08_STATUS, "FAIL"))
  warning("STEP 08 finished with FATAL QC failures. See STEP08_FINAL_QC.txt.", call. = FALSE)

invisible(list(status = STEP08_STATUS, qc = qc, inputs = inputs, outputs = outdt2))
###############################################################################
## END OF SCRIPT
###############################################################################