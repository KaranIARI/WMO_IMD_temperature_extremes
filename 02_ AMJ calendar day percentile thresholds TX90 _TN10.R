# ============================================================
# CHECK VALID-OBSERVATION QC FILE STATUS
# ============================================================

library(arrow)

QC_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_valid_observation_QC_1981_2010.parquet"

cat("\n========================================\n")
cat("VALID-OBSERVATION QC FILE CHECK\n")
cat("========================================\n")

if (file.exists(QC_FILE)) {
  
  cat("FILE EXISTS: YES\n\n")
  
  info <- file.info(QC_FILE)
  
  cat("File size:",
      round(info$size / 1024^2, 3),
      "MB\n")
  
  cat("Created:",
      as.character(info$ctime),
      "\n")
  
  cat("Modified:",
      as.character(info$mtime),
      "\n")
  
  # Read file
  qc <- read_parquet(QC_FILE)
  
  cat("\nRows:", nrow(qc), "\n")
  cat("Columns:", ncol(qc), "\n")
  
  cat("\nColumn names:\n")
  print(names(qc))
  
  cat("\nFirst 5 rows:\n")
  print(head(qc, 5))
  
  cat("\n========================================\n")
  
  if (nrow(qc) == 91 * 283) {
    cat("STATUS: VALID-OBSERVATION QC COMPLETED\n")
    cat("Expected rows = 25,753\n")
    cat("Actual rows   =", nrow(qc), "\n")
  } else {
    cat("STATUS: FILE EXISTS BUT ROW COUNT IS WRONG\n")
    cat("Expected rows = 25,753\n")
    cat("Actual rows   =", nrow(qc), "\n")
  }
  
} else {
  
  cat("FILE EXISTS: NO\n")
  cat("\nThe last QC script did NOT reach its final save step,\n")
  cat("or the output path/name was different.\n")
}

cat("========================================\n")
# ============================================================
# FINAL VALID-COUNT SUMMARY
# ============================================================

library(arrow)
library(data.table)

QC_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_valid_observation_QC_1981_2010.parquet"

qc <- as.data.table(read_parquet(QC_FILE))

cat("\n========================================\n")
cat("VALID OBSERVATION SUMMARY\n")
cat("========================================\n")

cat("\nTmax valid observations:\n")
print(
  quantile(
    qc$Tmax_valid,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
)

cat("\nTmin valid observations:\n")
print(
  quantile(
    qc$Tmin_valid,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
)

cat("\nMinimum Tmax valid:", min(qc$Tmax_valid), "/150\n")
cat("Minimum Tmin valid:", min(qc$Tmin_valid), "/150\n")

cat("\nTmax grid-days with <90% valid:",
    sum(qc$Tmax_valid_pct < 90), "\n")

cat("Tmin grid-days with <90% valid:",
    sum(qc$Tmin_valid_pct < 90), "\n")

cat("\nTmax grid-days with <80% valid:",
    sum(qc$Tmax_valid_pct < 80), "\n")

cat("Tmin grid-days with <80% valid:",
    sum(qc$Tmin_valid_pct < 80), "\n")

cat("\nTmax grid-days with <75% valid:",
    sum(qc$Tmax_valid_pct < 75), "\n")

cat("Tmin grid-days with <75% valid:",
    sum(qc$Tmin_valid_pct < 75), "\n")

cat("========================================\n")
###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 02 : AMJ calendar-day percentile thresholds (TX90 / TN10)
##  Baseline: 1981-2010   |   Target season: 01 Apr - 30 Jun (91 days)
##
##  INPUT  : F:/WMO_IMD_R/WMO_IMD/data/AMJ_buffer_1981_2010_30Mar_02Jul.parquet
##  OUTPUT : F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet
##
##  METHOD (locked - do not alter without documenting):
##    * For every target calendar day d (01Apr..30Jun) and every grid g,
##      pool all daily values falling in the 5-day window [d-2, d+2]
##      across the 30 baseline years 1981-2010.
##    * Pool size (possible) = 30 years x 5 days = 150 observations.
##    * NA = source missing value. NEVER interpolated, imputed or filled.
##    * Minimum validity = 80% of 150 = 120 valid (non-NA) observations.
##      If n_valid < 120  ->  threshold = NA (no partial estimate).
##    * TX90 = 90th percentile of Tmax pool ; TN10 = 10th percentile of Tmin pool.
##    * Expected output = 91 target days x 283 grids = 25,753 rows.
##
##  DESIGN NOTES:
##    (a) Calendar-day matching uses MONTH-DAY (MMDD) keys, NOT day-of-year.
##        Day-of-year shifts by 1 after 29-Feb in leap years; MMDD does not.
##        The whole window range (30 Mar - 02 Jul) sits after 29-Feb, so MMDD
##        matching gives exactly 5 window days x 30 years = 150 for every day.
##    (b) The previous `invalid 'trim' argument` error came from calling
##        format(<IDate>, "%m%d") -- the format string was consumed positionally
##        as the `trim` argument by the dispatched method. FIXED HERE by never
##        calling format() on a date column: month/day are extracted with the
##        integer accessors data.table::month() / data.table::mday(), and every
##        remaining format() call passes `format=` as a NAMED argument.
##    (c) Performance: the data are NOT re-filtered 91 times. A tiny 455-row
##        window map (91 targets x 5 offsets) is key-joined once to the buffer,
##        expanding it to 3,862,950 rows, then aggregated in a single grouped
##        data.table pass. Result is numerically identical to looping.
##    (d) Zhang et al. (2005) out-of-base bootstrap is deliberately NOT applied:
##        it is required when computing exceedance *indices* inside the base
##        period, not when publishing the threshold table itself.
##
##  This script is non-interactive and runs start to finish under Rscript.
###############################################################################

## ---------------------------------------------------------------------------
## 0. ENVIRONMENT
## ---------------------------------------------------------------------------
rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(TZ = "UTC")

AUTO_INSTALL <- TRUE
CRAN_REPO    <- "https://cloud.r-project.org"

.need <- c("data.table", "arrow")
for (p in .need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (AUTO_INSTALL) {
      message("Installing missing package: ", p)
      install.packages(p, repos = CRAN_REPO)
    }
    if (!requireNamespace(p, quietly = TRUE))
      stop("Required package '", p, "' is not available. Install it and re-run.",
           call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

setDTthreads(0L)   # use all available cores for data.table

## ---------------------------------------------------------------------------
## 1. LOCKED PARAMETERS  (single source of truth)
## ---------------------------------------------------------------------------
IN_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_buffer_1981_2010_30Mar_02Jul.parquet"
OUT_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet"
QC_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010_QC.csv"
LOG_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010_LOG.txt"

BASE_YR_START <- 1981L
BASE_YR_END   <- 2010L
N_YEARS       <- BASE_YR_END - BASE_YR_START + 1L      # 30

WIN_HALF      <- 2L                                    # +/- 2 days
WIN_WIDTH     <- 2L * WIN_HALF + 1L                    # 5

N_POSSIBLE    <- N_YEARS * WIN_WIDTH                   # 150
MIN_FRAC      <- 0.80
MIN_VALID     <- ceiling(MIN_FRAC * N_POSSIBLE)        # 120

P_TX          <- 0.90                                  # TX90
P_TN          <- 0.10                                  # TN10

## Percentile estimator -------------------------------------------------------
## quantile(type = 8) = Hyndman & Fan (1996) definition 8, the median-unbiased
## plotting-position estimator used by climdex.pcic / RClimDex / ETCCDI.
## Changing this constant changes the published threshold values.
PCTL_TYPE     <- 8L

EXP_IN_ROWS   <- 806550L                               # 283 x 95 x 30
EXP_GRIDS     <- 283L
EXP_TARGETS   <- 91L
EXP_OUT_ROWS  <- 25753L                                # 283 x 91
EXP_BUF_DAYS  <- 95L                                   # 30Mar..02Jul inclusive

## Physical plausibility envelope (QC warning only, never alters values)
TX90_LO <- 15; TX90_HI <- 60
TN10_LO <- -15; TN10_HI <- 40

## ---------------------------------------------------------------------------
## 2. LOGGING + CHECK HELPERS
## ---------------------------------------------------------------------------
.LOG <- character(0)

log_msg <- function(..., blank_before = FALSE) {
  txt <- paste0(..., collapse = "")
  stamp <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")   # NAMED `format=`
  line <- paste0("[", stamp, "] ", txt)
  if (blank_before) cat("\n")
  cat(line, "\n", sep = "")
  .LOG <<- c(.LOG, if (blank_before) "" else NULL, line)
  invisible(line)
}

log_head <- function(txt) {
  bar <- strrep("-", 76)
  cat("\n", bar, "\n", sep = "")
  cat(txt, "\n", sep = "")
  cat(bar, "\n", sep = "")
  .LOG <<- c(.LOG, "", bar, txt, bar)
  invisible(NULL)
}

CHECKS <- list()

record_check <- function(name, passed, detail, severity = c("FATAL", "WARN")) {
  severity <- match.arg(severity)
  status <- if (isTRUE(passed)) "PASS" else severity
  CHECKS[[length(CHECKS) + 1L]] <<- data.table(
    check = name, status = status, detail = detail
  )
  log_msg(sprintf("  %-5s | %-38s | %s", status, name, detail))
  if (!isTRUE(passed) && severity == "FATAL")
    stop("FATAL QC failure: ", name, " -- ", detail, call. = FALSE)
  invisible(isTRUE(passed))
}

t_start <- Sys.time()
log_head("STEP 02 | AMJ TX90 / TN10 baseline percentile thresholds (1981-2010)")
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("Baseline         : ", BASE_YR_START, "-", BASE_YR_END, " (", N_YEARS, " years)")
log_msg("Window           : target day +/- ", WIN_HALF, " days (", WIN_WIDTH, "-day window)")
log_msg("Possible obs     : ", N_YEARS, " x ", WIN_WIDTH, " = ", N_POSSIBLE, " per grid-day")
log_msg("Validity rule    : n_valid >= ", MIN_VALID, " (", MIN_FRAC * 100, "% of ", N_POSSIBLE, ") else NA")
log_msg("Percentiles      : Tmax p", P_TX * 100, " / Tmin p", P_TN * 100,
        " via stats::quantile(type = ", PCTL_TYPE, ")")
log_msg("NA policy        : source missing values preserved; no interpolation/imputation")
log_msg("Input            : ", IN_FILE)
log_msg("Output           : ", OUT_FILE)

## ---------------------------------------------------------------------------
## 3. READ VALIDATED BUFFER (read-only; never modified or recomputed)
## ---------------------------------------------------------------------------
log_head("3. READ INPUT")

if (!file.exists(IN_FILE))
  stop("Input parquet not found: ", IN_FILE, call. = FALSE)

buf <- as.data.table(arrow::read_parquet(IN_FILE))
log_msg("Rows read        : ", format(nrow(buf), big.mark = ","))
log_msg("Columns present  : ", paste(names(buf), collapse = ", "))

## --- resolve column names case-insensitively --------------------------------
pick_col <- function(dt, candidates, label) {
  hit <- names(dt)[tolower(names(dt)) %in% tolower(candidates)]
  if (length(hit) < 1L)
    stop("Could not find a '", label, "' column. Looked for: ",
         paste(candidates, collapse = ", "), ". Present: ",
         paste(names(dt), collapse = ", "), call. = FALSE)
  hit[1L]
}

C_GRID <- pick_col(buf, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id")
C_LAT  <- pick_col(buf, c("lat", "latitude", "y"), "lat")
C_LON  <- pick_col(buf, c("lon", "long", "longitude", "x"), "lon")
C_DATE <- pick_col(buf, c("date", "dates", "obs_date", "time", "day"), "date")
C_TMAX <- pick_col(buf, c("tmax", "tx", "t_max", "max_temp", "temp_max", "tmax_c"), "tmax")
C_TMIN <- pick_col(buf, c("tmin", "tn", "t_min", "min_temp", "temp_min", "tmin_c"), "tmin")

log_msg("Column mapping   : grid_id<-", C_GRID, " lat<-", C_LAT, " lon<-", C_LON,
        " date<-", C_DATE, " tmax<-", C_TMAX, " tmin<-", C_TMIN)

buf <- buf[, c(C_GRID, C_LAT, C_LON, C_DATE, C_TMAX, C_TMIN), with = FALSE]
setnames(buf, c("grid_id", "lat", "lon", "date", "tmax", "tmin"))

## --- robust date coercion to IDate (no format() on date columns) ------------
coerce_idate <- function(x) {
  if (inherits(x, "IDate"))   return(x)
  if (inherits(x, "Date"))    return(as.IDate(x))
  if (inherits(x, "POSIXct")) return(as.IDate(as.Date(x, tz = "UTC")))
  if (is.character(x))        return(as.IDate(as.Date(x, format = "%Y-%m-%d")))
  if (is.numeric(x)) {
    if (all(is.na(x) | (x > 1e7 & x < 1e8)))   # YYYYMMDD integer
      return(as.IDate(as.Date(as.character(as.integer(x)), format = "%Y%m%d")))
    return(as.IDate(as.Date(x, origin = "1970-01-01")))
  }
  stop("Unsupported date column class: ", paste(class(x), collapse = "/"), call. = FALSE)
}
buf[, date := coerce_idate(date)]
if (anyNA(buf$date)) stop("Date coercion produced NA values -- inspect input.", call. = FALSE)

buf[, `:=`(
  yr = as.integer(data.table::year(date)),
  md = as.integer(data.table::month(date)) * 100L + as.integer(data.table::mday(date))
)]

buf[, `:=`(tmax = as.numeric(tmax), tmin = as.numeric(tmin))]
if (!is.integer(buf$grid_id) && !is.character(buf$grid_id))
  buf[, grid_id := as.character(grid_id)]

## ---------------------------------------------------------------------------
## 4. PRE-FLIGHT INPUT VALIDATION (hard gates)
## ---------------------------------------------------------------------------
log_head("4. INPUT VALIDATION")

n_in     <- nrow(buf)
n_grid   <- uniqueN(buf$grid_id)
yrs      <- sort(unique(buf$yr))
n_bufday <- uniqueN(buf$md)

record_check("input_row_count", n_in == EXP_IN_ROWS,
             sprintf("observed %s vs expected %s",
                     format(n_in, big.mark = ","), format(EXP_IN_ROWS, big.mark = ",")))

record_check("grid_count", n_grid == EXP_GRIDS,
             sprintf("observed %d vs expected %d", n_grid, EXP_GRIDS))

record_check("baseline_years", identical(yrs, BASE_YR_START:BASE_YR_END),
             sprintf("observed %d years [%d-%d]", length(yrs), min(yrs), max(yrs)))

record_check("buffer_calendar_days", n_bufday == EXP_BUF_DAYS,
             sprintf("observed %d unique MMDD vs expected %d (30Mar-02Jul)",
                     n_bufday, EXP_BUF_DAYS))

record_check("no_duplicate_grid_date",
             !anyDuplicated(buf, by = c("grid_id", "date")),
             sprintf("%d duplicate grid-date rows",
                     sum(duplicated(buf, by = c("grid_id", "date")))))

meta <- unique(buf[, .(grid_id, lat, lon)])
record_check("unique_coords_per_grid", nrow(meta) == EXP_GRIDS,
             sprintf("%d unique grid/lat/lon triplets vs %d grids", nrow(meta), EXP_GRIDS))

record_check("coords_not_missing", !anyNA(meta$lat) && !anyNA(meta$lon),
             sprintf("%d NA lat, %d NA lon", sum(is.na(meta$lat)), sum(is.na(meta$lon))))

log_msg("Source NA (Tmax) : ", format(sum(is.na(buf$tmax)), big.mark = ","),
        sprintf("  (%.3f%% of rows)", 100 * mean(is.na(buf$tmax))))
log_msg("Source NA (Tmin) : ", format(sum(is.na(buf$tmin)), big.mark = ","),
        sprintf("  (%.3f%% of rows)", 100 * mean(is.na(buf$tmin))))

## ---------------------------------------------------------------------------
## 5. BUILD TARGET-DAY x WINDOW MAP  (91 targets x 5 offsets = 455 rows)
## ---------------------------------------------------------------------------
log_head("5. BUILD 5-DAY WINDOW MAP")

## A non-leap reference year is used purely to enumerate calendar labels.
REF_YR   <- 2001L
tgt_dates <- seq(as.Date(paste0(REF_YR, "-04-01")),
                 as.Date(paste0(REF_YR, "-06-30")), by = "day")

record_check("target_day_count", length(tgt_dates) == EXP_TARGETS,
             sprintf("observed %d vs expected %d (01Apr-30Jun)",
                     length(tgt_dates), EXP_TARGETS))

targets <- data.table(
  target_index = seq_along(tgt_dates),
  target_md    = as.integer(data.table::month(tgt_dates)) * 100L +
    as.integer(data.table::mday(tgt_dates))
)
targets[, target_day := sprintf("%02d-%02d", target_md %/% 100L, target_md %% 100L)]

win <- CJ(target_index = targets$target_index, offset = -WIN_HALF:WIN_HALF)
win <- merge(win, targets, by = "target_index", sort = FALSE)
win[, win_date := tgt_dates[target_index] + offset]
win[, md := as.integer(data.table::month(win_date)) * 100L +
      as.integer(data.table::mday(win_date))]
win[, win_date := NULL]

record_check("window_map_size", nrow(win) == EXP_TARGETS * WIN_WIDTH,
             sprintf("observed %d vs expected %d", nrow(win), EXP_TARGETS * WIN_WIDTH))

record_check("window_days_within_buffer",
             all(win$md %in% unique(buf$md)),
             sprintf("%d window MMDD not present in buffer",
                     sum(!win$md %in% unique(buf$md))))

log_msg("Window map built : ", nrow(win), " (target_day, window_day) pairs")
log_msg("First target     : ", targets$target_day[1], "  window MMDD = ",
        paste(sort(win[target_index == 1L, md]), collapse = ", "))
log_msg("Last target      : ", targets$target_day[EXP_TARGETS], "  window MMDD = ",
        paste(sort(win[target_index == EXP_TARGETS, md]), collapse = ", "))

## ---------------------------------------------------------------------------
## 6. SINGLE KEY-JOIN EXPANSION + ONE GROUPED PASS
##    (replaces 91 repeated full-data filters; identical scientific result)
## ---------------------------------------------------------------------------
log_head("6. THRESHOLD COMPUTATION")

setkey(buf, md)
setkey(win, md)

t0 <- Sys.time()
pool <- buf[win, .(grid_id, target_index, tmax, tmin),
            allow.cartesian = TRUE, nomatch = 0L]
log_msg("Expanded pool    : ", format(nrow(pool), big.mark = ","), " rows in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("Expected pool    : ", format(EXP_OUT_ROWS * N_POSSIBLE, big.mark = ","),
        " rows (25,753 grid-days x 150)")

record_check("expanded_pool_rows", nrow(pool) == EXP_OUT_ROWS * N_POSSIBLE,
             sprintf("observed %s vs expected %s",
                     format(nrow(pool), big.mark = ","),
                     format(EXP_OUT_ROWS * N_POSSIBLE, big.mark = ",")))

t0 <- Sys.time()
res <- pool[, {
  vx <- tmax[!is.na(tmax)]
  vn <- tmin[!is.na(tmin)]
  nx <- length(vx)
  nn <- length(vn)
  .(n_obs                = .N,
    n_valid_tmax         = nx,
    n_valid_tmin         = nn,
    completeness_tmax_pct = round(100 * nx / N_POSSIBLE, 2),
    completeness_tmin_pct = round(100 * nn / N_POSSIBLE, 2),
    tx90 = if (nx >= MIN_VALID)
      as.numeric(stats::quantile(vx, probs = P_TX, type = PCTL_TYPE, names = FALSE))
    else NA_real_,
    tn10 = if (nn >= MIN_VALID)
      as.numeric(stats::quantile(vn, probs = P_TN, type = PCTL_TYPE, names = FALSE))
    else NA_real_)
}, by = .(grid_id, target_index)]
log_msg("Aggregated       : ", format(nrow(res), big.mark = ","), " grid-days in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

rm(pool); invisible(gc(verbose = FALSE))

## --- attach labels + coordinates -------------------------------------------
res <- merge(res, targets, by = "target_index", all.x = TRUE, sort = FALSE)
res <- merge(res, meta,    by = "grid_id",      all.x = TRUE, sort = FALSE)

res[, `:=`(
  baseline_start = BASE_YR_START,
  baseline_end   = BASE_YR_END,
  n_possible     = N_POSSIBLE,
  min_valid_req  = MIN_VALID,
  window_days    = WIN_WIDTH,
  pctl_type      = PCTL_TYPE,
  tx90_prob      = P_TX,
  tn10_prob      = P_TN,
  tx90_valid     = n_valid_tmax >= MIN_VALID,
  tn10_valid     = n_valid_tmin >= MIN_VALID
)]

setcolorder(res, c(
  "grid_id", "lat", "lon",
  "target_index", "target_md", "target_day",
  "n_obs", "n_possible", "min_valid_req",
  "n_valid_tmax", "completeness_tmax_pct", "tx90", "tx90_valid",
  "n_valid_tmin", "completeness_tmin_pct", "tn10", "tn10_valid",
  "window_days", "baseline_start", "baseline_end",
  "pctl_type", "tx90_prob", "tn10_prob"
))
setorder(res, grid_id, target_index)

## ---------------------------------------------------------------------------
## 7. FINAL QUALITY CONTROL
## ---------------------------------------------------------------------------
log_head("7. FINAL QC")

record_check("output_row_count", nrow(res) == EXP_OUT_ROWS,
             sprintf("observed %s vs expected %s",
                     format(nrow(res), big.mark = ","),
                     format(EXP_OUT_ROWS, big.mark = ",")))

record_check("output_grid_count", uniqueN(res$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(res$grid_id), EXP_GRIDS))

record_check("output_target_day_count", uniqueN(res$target_day) == EXP_TARGETS,
             sprintf("observed %d vs expected %d", uniqueN(res$target_day), EXP_TARGETS))

n_dup <- sum(duplicated(res, by = c("grid_id", "target_day")))
record_check("no_duplicate_grid_day", n_dup == 0L,
             sprintf("%d duplicated (grid_id, target_day) rows", n_dup))

record_check("complete_grid_x_day_cross",
             nrow(res) == uniqueN(res$grid_id) * uniqueN(res$target_day),
             sprintf("%d rows vs %d x %d full cross",
                     nrow(res), uniqueN(res$grid_id), uniqueN(res$target_day)))

bad_nobs <- res[n_obs != N_POSSIBLE, .N]
record_check("pool_size_is_150", bad_nobs == 0L,
             sprintf("%d grid-days with n_obs != %d", bad_nobs, N_POSSIBLE))

record_check("n_valid_within_bounds",
             res[, all(n_valid_tmax >= 0 & n_valid_tmax <= N_POSSIBLE &
                         n_valid_tmin >= 0 & n_valid_tmin <= N_POSSIBLE)],
             sprintf("Tmax range %d-%d ; Tmin range %d-%d",
                     min(res$n_valid_tmax), max(res$n_valid_tmax),
                     min(res$n_valid_tmin), max(res$n_valid_tmin)))

## validity rule enforced in both directions
v1 <- res[n_valid_tmax <  MIN_VALID & !is.na(tx90), .N]
v2 <- res[n_valid_tmax >= MIN_VALID &  is.na(tx90), .N]
v3 <- res[n_valid_tmin <  MIN_VALID & !is.na(tn10), .N]
v4 <- res[n_valid_tmin >= MIN_VALID &  is.na(tn10), .N]
record_check("validity_rule_tx90", v1 == 0L && v2 == 0L,
             sprintf("%d thresholds below n=%d, %d NA above n=%d", v1, MIN_VALID, v2, MIN_VALID))
record_check("validity_rule_tn10", v3 == 0L && v4 == 0L,
             sprintf("%d thresholds below n=%d, %d NA above n=%d", v3, MIN_VALID, v4, MIN_VALID))

na_tx <- res[is.na(tx90), .N]
na_tn <- res[is.na(tn10), .N]
log_msg("TX90 NA thresholds : ", format(na_tx, big.mark = ","),
        sprintf("  (%.2f%%)", 100 * na_tx / nrow(res)))
log_msg("TN10 NA thresholds : ", format(na_tn, big.mark = ","),
        sprintf("  (%.2f%%)", 100 * na_tn / nrow(res)))
record_check("na_thresholds_documented", TRUE,
             sprintf("TX90 NA = %d ; TN10 NA = %d (insufficient valid observations)",
                     na_tx, na_tn))

## plausibility envelope -- WARN only, values are never altered
oor_tx <- res[!is.na(tx90) & (tx90 < TX90_LO | tx90 > TX90_HI), .N]
oor_tn <- res[!is.na(tn10) & (tn10 < TN10_LO | tn10 > TN10_HI), .N]
record_check("tx90_plausible_range", oor_tx == 0L,
             sprintf("%d values outside [%g, %g] degC ; observed [%.2f, %.2f]",
                     oor_tx, TX90_LO, TX90_HI,
                     suppressWarnings(min(res$tx90, na.rm = TRUE)),
                     suppressWarnings(max(res$tx90, na.rm = TRUE))),
             severity = "WARN")
record_check("tn10_plausible_range", oor_tn == 0L,
             sprintf("%d values outside [%g, %g] degC ; observed [%.2f, %.2f]",
                     oor_tn, TN10_LO, TN10_HI,
                     suppressWarnings(min(res$tn10, na.rm = TRUE)),
                     suppressWarnings(max(res$tn10, na.rm = TRUE))),
             severity = "WARN")

n_inv <- res[!is.na(tx90) & !is.na(tn10) & tx90 <= tn10, .N]
record_check("tx90_greater_than_tn10", n_inv == 0L,
             sprintf("%d grid-days with TX90 <= TN10", n_inv), severity = "WARN")

record_check("coords_carried_through", !anyNA(res$lat) && !anyNA(res$lon),
             sprintf("%d NA lat, %d NA lon", sum(is.na(res$lat)), sum(is.na(res$lon))))

## descriptive summary
log_msg("TX90 summary     : ", paste(sprintf("%.2f", quantile(res$tx90, c(0, .25, .5, .75, 1),
                                                              na.rm = TRUE, names = FALSE)), collapse = " | "), "   (min|q1|med|q3|max degC)")
log_msg("TN10 summary     : ", paste(sprintf("%.2f", quantile(res$tn10, c(0, .25, .5, .75, 1),
                                                              na.rm = TRUE, names = FALSE)), collapse = " | "), "   (min|q1|med|q3|max degC)")
log_msg("Grids with any NA TX90 : ", uniqueN(res[is.na(tx90), grid_id]), " / ", EXP_GRIDS)
log_msg("Grids with any NA TN10 : ", uniqueN(res[is.na(tn10), grid_id]), " / ", EXP_GRIDS)

## ---------------------------------------------------------------------------
## 8. WRITE OUTPUTS
## ---------------------------------------------------------------------------
log_head("8. WRITE OUTPUTS")

out_dir <- dirname(OUT_FILE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

arrow::write_parquet(res, OUT_FILE, compression = "snappy")
log_msg("Thresholds saved : ", OUT_FILE, "  (",
        format(nrow(res), big.mark = ","), " rows x ", ncol(res), " cols, ",
        sprintf("%.2f MB", file.size(OUT_FILE) / 1024^2), ")")

qc_tab <- rbindlist(CHECKS)
qc_tab[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
              script   = "02_AMJ_TX90_TN10_thresholds_1981_2010.R")]
fwrite(qc_tab, QC_FILE)
log_msg("QC summary saved : ", QC_FILE, "  (", nrow(qc_tab), " checks)")

## read-back verification of the written file
chk <- as.data.table(arrow::read_parquet(OUT_FILE))
record_check("readback_row_count", nrow(chk) == EXP_OUT_ROWS,
             sprintf("re-read %s rows", format(nrow(chk), big.mark = ",")))
record_check("readback_values_identical",
             isTRUE(all.equal(chk$tx90, res$tx90)) && isTRUE(all.equal(chk$tn10, res$tn10)),
             "TX90/TN10 identical after round-trip")
rm(chk)

n_fail <- qc_tab[status == "FATAL", .N]
n_warn <- rbindlist(CHECKS)[status == "WARN", .N]

log_head("RUN COMPLETE")
log_msg("Checks passed    : ", rbindlist(CHECKS)[status == "PASS", .N], " / ", length(CHECKS))
log_msg("Warnings         : ", n_warn)
log_msg("Fatal            : ", n_fail)
log_msg("Elapsed          : ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

.LOG <- c(.LOG, "", "--- sessionInfo() ---",
          capture.output(utils::sessionInfo()))
writeLines(.LOG, LOG_FILE)
cat("\nConsole log written to:", LOG_FILE, "\n")

## final object left in the workspace for interactive inspection: `res`
invisible(res)
###############################################################################
## END OF SCRIPT
###############################################################################
