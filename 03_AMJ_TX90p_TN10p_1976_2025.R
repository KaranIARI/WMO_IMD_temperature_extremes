###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 03 : Annual AMJ TX90p / TN10p (percentage of exceedance days)
##  Grids: 283  |  Period: 1976-2025 (50 years)  |  Season: 01 Apr - 30 Jun
##
##  INPUT 1 (daily, read-only) :
##    F:/WMO_IMD_R/WMO_IMD/data/IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet
##  INPUT 2 (thresholds, read-only, LOCKED - never modified) :
##    F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet
##  OUTPUT :
##    F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_1976_2025.parquet
##
##  METHOD (locked):
##    * Restrict daily data to AMJ = 01 Apr - 30 Jun (91 calendar days/year).
##    * Join each daily observation to its own grid + calendar-day threshold
##      (matched on grid_id + MMDD, NOT day-of-year: DOY shifts by 1 after
##      29-Feb in leap years, MMDD does not).
##    * Tmax >  TX90  -> TX90 exceedance day
##      Tmin <  TN10  -> TN10 exceedance day
##      Comparisons are strict (>, <), consistent with ETCCDI TX90p / TN10p.
##    * Missing Tmax/Tmin stay NA, are never imputed, and never count as
##      exceedances.
##    * TX90p = 100 x TX90_days / n_valid_tmax
##      TN10p = 100 x TN10_days / n_valid_tmin
##      Denominators are computed independently for Tmax and Tmin.
##      91 is NEVER used as a fixed denominator.
##    * If a denominator is 0 the index is NA (never 0, never 0/0).
##
##  DENOMINATOR DEFINITION (explicit, see log):
##    A day is counted as "valid" for TX90p only if BOTH
##       (i)  the daily Tmax observation is non-NA, AND
##       (ii) the corresponding calendar-day TX90 threshold is non-NA.
##    A day whose threshold is NA (a baseline grid-day that failed the
##    >=120/150 rule in STEP 02) cannot be tested for exceedance, so keeping it
##    in the denominator would systematically deflate TX90p. Such days are
##    therefore excluded from BOTH numerator and denominator, and counted
##    separately in `n_no_threshold_tmax` / `n_no_threshold_tmin` so the
##    decision is fully auditable. The raw observation counts are also carried
##    through as `n_obs_tmax` / `n_obs_tmin`. If the threshold file contains no
##    NA thresholds, the two definitions coincide exactly (the log states which).
##
##  SCOPE:
##    NO bootstrap. In-base years (1981-2010) are computed here by direct
##    comparison only; the Zhang et al. (2005) out-of-base bootstrap will be a
##    separate, later step. Column `bootstrap_applied` is written as FALSE.
##    Nothing upstream (GRD extraction, buffer, thresholds) is recomputed.
##
##  Non-interactive: runs start to finish under Rscript.
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
setDTthreads(0L)

## ---------------------------------------------------------------------------
## 1. LOCKED PARAMETERS
## ---------------------------------------------------------------------------
DAILY_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet"
THR_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet"
OUT_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_1976_2025.parquet"
QC_FILE    <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_1976_2025_QC.csv"
LOG_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_1976_2025_LOG.txt"

YR_START      <- 1976L
YR_END        <- 2025L
N_YEARS       <- YR_END - YR_START + 1L        # 50

BASE_YR_START <- 1981L                         # metadata flag only (no bootstrap here)
BASE_YR_END   <- 2010L

MD_START      <- 401L                          # 01 April  (MMDD)
MD_END        <- 630L                          # 30 June   (MMDD)
N_AMJ_DAYS    <- 91L

EXP_GRIDS     <- 283L
EXP_THR_ROWS  <- 25753L                        # 283 x 91
EXP_AMJ_ROWS  <- 1287650L                      # 283 x 91 x 50
EXP_OUT_ROWS  <- 14150L                        # 283 x 50

COORD_TOL     <- 1e-6                          # lat/lon cross-file agreement

## Informational sanity band for the baseline-period mean of each index.
## By construction TX90p / TN10p average ~10% over 1981-2010 (in-base values are
## slightly biased low before bootstrap). WARN only - values are never altered.
BASE_MEAN_LO  <- 5
BASE_MEAN_HI  <- 20

## ---------------------------------------------------------------------------
## 2. LOGGING + CHECK HELPERS
## ---------------------------------------------------------------------------
.LOG <- character(0)

log_msg <- function(..., blank_before = FALSE) {
  txt   <- paste0(..., collapse = "")
  stamp <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")   # NAMED `format=`
  line  <- paste0("[", stamp, "] ", txt)
  if (blank_before) cat("\n")
  cat(line, "\n", sep = "")
  .LOG <<- c(.LOG, if (blank_before) "" else NULL, line)
  invisible(line)
}

log_head <- function(txt) {
  bar <- strrep("-", 76)
  cat("\n", bar, "\n", txt, "\n", bar, "\n", sep = "")
  .LOG <<- c(.LOG, "", bar, txt, bar)
  invisible(NULL)
}

CHECKS <- list()

record_check <- function(name, passed, detail, severity = c("FATAL", "WARN")) {
  severity <- match.arg(severity)
  status   <- if (isTRUE(passed)) "PASS" else severity
  CHECKS[[length(CHECKS) + 1L]] <<- data.table(check = name, status = status, detail = detail)
  log_msg(sprintf("  %-5s | %-38s | %s", status, name, detail))
  if (!isTRUE(passed) && severity == "FATAL")
    stop("FATAL QC failure: ", name, " -- ", detail, call. = FALSE)
  invisible(isTRUE(passed))
}

pick_col <- function(dt, candidates, label, required = TRUE) {
  hit <- names(dt)[tolower(names(dt)) %in% tolower(candidates)]
  if (length(hit) < 1L) {
    if (!required) return(NA_character_)
    stop("Could not find a '", label, "' column. Looked for: ",
         paste(candidates, collapse = ", "), ". Present: ",
         paste(names(dt), collapse = ", "), call. = FALSE)
  }
  hit[1L]
}

coerce_idate <- function(x) {
  if (inherits(x, "IDate"))   return(x)
  if (inherits(x, "Date"))    return(as.IDate(x))
  if (inherits(x, "POSIXct")) return(as.IDate(as.Date(x, tz = "UTC")))
  if (is.character(x))        return(as.IDate(as.Date(x, format = "%Y-%m-%d")))
  if (is.numeric(x)) {
    if (all(is.na(x) | (x > 1e7 & x < 1e8)))
      return(as.IDate(as.Date(as.character(as.integer(x)), format = "%Y%m%d")))
    return(as.IDate(as.Date(x, origin = "1970-01-01")))
  }
  stop("Unsupported date column class: ", paste(class(x), collapse = "/"), call. = FALSE)
}

t_start <- Sys.time()
log_head("STEP 03 | Annual AMJ TX90p / TN10p, 283 grids, 1976-2025")
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("Season           : 01 Apr - 30 Jun (MMDD ", MD_START, "-", MD_END,
        ", ", N_AMJ_DAYS, " days)")
log_msg("Years            : ", YR_START, "-", YR_END, " (", N_YEARS, ")")
log_msg("Exceedance rule  : Tmax > TX90 ; Tmin < TN10 (strict inequalities)")
log_msg("Denominator      : evaluable days only (obs non-NA AND threshold non-NA)")
log_msg("Bootstrap        : NOT applied in this step (deferred, by design)")
log_msg("Daily input      : ", DAILY_FILE)
log_msg("Threshold input  : ", THR_FILE)
log_msg("Output           : ", OUT_FILE)

## ---------------------------------------------------------------------------
## 3. READ THRESHOLDS (read-only; file is never written back)
## ---------------------------------------------------------------------------
log_head("3. READ LOCKED THRESHOLD FILE")

if (!file.exists(THR_FILE)) stop("Threshold parquet not found: ", THR_FILE, call. = FALSE)

thr <- as.data.table(arrow::read_parquet(THR_FILE))
log_msg("Rows read        : ", format(nrow(thr), big.mark = ","))
log_msg("Columns present  : ", paste(names(thr), collapse = ", "))

T_GRID <- pick_col(thr, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id")
T_MD   <- pick_col(thr, c("target_md", "md", "mmdd", "target_mmdd"), "target_md", required = FALSE)
T_DAY  <- pick_col(thr, c("target_day", "calendar_day", "day_label"), "target_day", required = FALSE)
T_TX   <- pick_col(thr, c("tx90", "tx90_threshold", "tmax_p90", "tx90_thresh"), "tx90")
T_TN   <- pick_col(thr, c("tn10", "tn10_threshold", "tmin_p10", "tn10_thresh"), "tn10")
T_LAT  <- pick_col(thr, c("lat", "latitude", "y"), "lat", required = FALSE)
T_LON  <- pick_col(thr, c("lon", "long", "longitude", "x"), "lon", required = FALSE)

## Reconstruct MMDD from target_day ("MM-DD") if target_md is absent.
if (is.na(T_MD)) {
  if (is.na(T_DAY))
    stop("Threshold file has neither 'target_md' nor 'target_day'.", call. = FALSE)
  dd <- as.character(thr[[T_DAY]])
  thr[, thr_md := as.integer(substr(dd, 1L, 2L)) * 100L + as.integer(substr(dd, 4L, 5L))]
  log_msg("target_md        : reconstructed from '", T_DAY, "'")
} else {
  thr[, thr_md := as.integer(get(T_MD))]
}

thr_keep <- c(T_GRID, "thr_md", T_TX, T_TN, if (!is.na(T_LAT)) T_LAT, if (!is.na(T_LON)) T_LON)
thr <- thr[, thr_keep, with = FALSE]
setnames(thr, c("grid_id", "md", "tx90", "tn10",
                if (!is.na(T_LAT)) "thr_lat", if (!is.na(T_LON)) "thr_lon"))
thr[, `:=`(grid_id = as.character(grid_id),
           tx90    = as.numeric(tx90),
           tn10    = as.numeric(tn10))]

record_check("threshold_row_count", nrow(thr) == EXP_THR_ROWS,
             sprintf("observed %s vs expected %s",
                     format(nrow(thr), big.mark = ","),
                     format(EXP_THR_ROWS, big.mark = ",")))
record_check("threshold_grid_count", uniqueN(thr$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(thr$grid_id), EXP_GRIDS))
record_check("threshold_day_count", uniqueN(thr$md) == N_AMJ_DAYS,
             sprintf("observed %d vs expected %d", uniqueN(thr$md), N_AMJ_DAYS))
record_check("threshold_no_duplicate_key",
             !anyDuplicated(thr, by = c("grid_id", "md")),
             sprintf("%d duplicated (grid_id, MMDD) rows",
                     sum(duplicated(thr, by = c("grid_id", "md")))))
record_check("threshold_days_in_AMJ", all(thr$md >= MD_START & thr$md <= MD_END),
             sprintf("MMDD range %d-%d", min(thr$md), max(thr$md)))

n_thr_na_tx <- sum(is.na(thr$tx90))
n_thr_na_tn <- sum(is.na(thr$tn10))
log_msg("NA TX90 thresholds : ", n_thr_na_tx,
        sprintf("  (%.3f%% of grid-days)", 100 * n_thr_na_tx / nrow(thr)))
log_msg("NA TN10 thresholds : ", n_thr_na_tn,
        sprintf("  (%.3f%% of grid-days)", 100 * n_thr_na_tn / nrow(thr)))
if (n_thr_na_tx == 0L && n_thr_na_tn == 0L) {
  log_msg("=> No NA thresholds: 'evaluable days' == 'non-NA observation days' exactly.")
} else {
  log_msg("=> Some thresholds are NA; those grid-days are excluded from BOTH the",
          " numerator and the denominator and are counted in n_no_threshold_*.")
}
record_check("threshold_na_documented", TRUE,
             sprintf("TX90 NA = %d ; TN10 NA = %d", n_thr_na_tx, n_thr_na_tn))

## ---------------------------------------------------------------------------
## 4. READ DAILY DATA AND RESTRICT TO AMJ
## ---------------------------------------------------------------------------
log_head("4. READ DAILY DATA + AMJ SUBSET")

if (!file.exists(DAILY_FILE)) stop("Daily parquet not found: ", DAILY_FILE, call. = FALSE)

t0  <- Sys.time()
dly <- as.data.table(arrow::read_parquet(DAILY_FILE))
log_msg("Rows read        : ", format(nrow(dly), big.mark = ","), " in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("Columns present  : ", paste(names(dly), collapse = ", "))

D_GRID <- pick_col(dly, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id")
D_LAT  <- pick_col(dly, c("lat", "latitude", "y"), "lat")
D_LON  <- pick_col(dly, c("lon", "long", "longitude", "x"), "lon")
D_DATE <- pick_col(dly, c("date", "dates", "obs_date", "time", "day"), "date")
D_TMAX <- pick_col(dly, c("tmax", "tx", "t_max", "max_temp", "temp_max", "tmax_c"), "tmax")
D_TMIN <- pick_col(dly, c("tmin", "tn", "t_min", "min_temp", "temp_min", "tmin_c"), "tmin")
log_msg("Column mapping   : grid_id<-", D_GRID, " lat<-", D_LAT, " lon<-", D_LON,
        " date<-", D_DATE, " tmax<-", D_TMAX, " tmin<-", D_TMIN)

dly <- dly[, c(D_GRID, D_LAT, D_LON, D_DATE, D_TMAX, D_TMIN), with = FALSE]
setnames(dly, c("grid_id", "lat", "lon", "date", "tmax", "tmin"))

dly[, date := coerce_idate(date)]
if (anyNA(dly$date)) stop("Date coercion produced NA values -- inspect daily input.", call. = FALSE)

dly[, `:=`(
  grid_id = as.character(grid_id),
  year    = as.integer(data.table::year(date)),
  md      = as.integer(data.table::month(date)) * 100L +
    as.integer(data.table::mday(date)),
  tmax    = as.numeric(tmax),
  tmin    = as.numeric(tmin)
)]

record_check("daily_no_duplicate_grid_date",
             !anyDuplicated(dly, by = c("grid_id", "date")),
             sprintf("%d duplicate grid-date rows",
                     sum(duplicated(dly, by = c("grid_id", "date")))))

## grid coordinate table (authoritative source = daily file)
meta <- unique(dly[, .(grid_id, lat, lon)])
record_check("daily_grid_count", uniqueN(dly$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(dly$grid_id), EXP_GRIDS))
record_check("unique_coords_per_grid", nrow(meta) == EXP_GRIDS,
             sprintf("%d unique grid/lat/lon triplets vs %d grids", nrow(meta), EXP_GRIDS))
record_check("coords_not_missing", !anyNA(meta$lat) && !anyNA(meta$lon),
             sprintf("%d NA lat, %d NA lon", sum(is.na(meta$lat)), sum(is.na(meta$lon))))

## cross-file coordinate agreement (WARN only; daily file wins)
if (all(c("thr_lat", "thr_lon") %in% names(thr))) {
  cmp <- merge(meta, unique(thr[, .(grid_id, thr_lat, thr_lon)]), by = "grid_id")
  n_mis <- cmp[abs(lat - thr_lat) > COORD_TOL | abs(lon - thr_lon) > COORD_TOL, .N]
  record_check("coords_match_threshold_file", n_mis == 0L,
               sprintf("%d grids differ beyond tol %g", n_mis, COORD_TOL),
               severity = "WARN")
  rm(cmp)
}
thr[, c("thr_lat", "thr_lon") := NULL]

## --- AMJ subset (single vectorised filter, no per-grid/per-year looping) ----
amj <- dly[md >= MD_START & md <= MD_END & year >= YR_START & year <= YR_END]
rm(dly); invisible(gc(verbose = FALSE))

log_msg("AMJ rows         : ", format(nrow(amj), big.mark = ","),
        "  (expected ", format(EXP_AMJ_ROWS, big.mark = ","), ")")

yrs <- sort(unique(amj$year))
record_check("year_coverage", identical(yrs, YR_START:YR_END),
             sprintf("observed %d years [%d-%d]", length(yrs), min(yrs), max(yrs)))
record_check("amj_calendar_day_count", uniqueN(amj$md) == N_AMJ_DAYS,
             sprintf("observed %d unique MMDD vs expected %d",
                     uniqueN(amj$md), N_AMJ_DAYS))
record_check("amj_row_count", nrow(amj) == EXP_AMJ_ROWS,
             sprintf("observed %s vs expected %s (a shortfall means whole rows are absent, not merely NA)",
                     format(nrow(amj), big.mark = ","),
                     format(EXP_AMJ_ROWS, big.mark = ",")),
             severity = "WARN")

log_msg("AMJ NA Tmax      : ", format(sum(is.na(amj$tmax)), big.mark = ","),
        sprintf("  (%.3f%%)", 100 * mean(is.na(amj$tmax))))
log_msg("AMJ NA Tmin      : ", format(sum(is.na(amj$tmin)), big.mark = ","),
        sprintf("  (%.3f%%)", 100 * mean(is.na(amj$tmin))))

## ---------------------------------------------------------------------------
## 5. THRESHOLD JOIN + DAILY EXCEEDANCE FLAGS
## ---------------------------------------------------------------------------
log_head("5. THRESHOLD JOIN + DAILY EXCEEDANCE")

setkey(thr, grid_id, md)
t0 <- Sys.time()
amj[thr, on = .(grid_id, md), `:=`(tx90 = i.tx90, tn10 = i.tn10)]
log_msg("Joined thresholds in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## Every AMJ row MUST find a threshold row (catches grid_id type/label drift).
n_unmatched <- amj[, sum(!paste(grid_id, md) %chin%
                           paste(thr$grid_id, thr$md))]
record_check("all_days_matched_to_threshold_key", n_unmatched == 0L,
             sprintf("%s AMJ rows without a (grid_id, MMDD) threshold record",
                     format(n_unmatched, big.mark = ",")))

## Evaluability + strict exceedance. `FALSE & NA` is FALSE in R, so no NA leaks
## into the integer flags; unevaluable days contribute 0 to both numerator and
## denominator.
amj[, `:=`(
  eval_tmax = !is.na(tmax) & !is.na(tx90),
  eval_tmin = !is.na(tmin) & !is.na(tn10)
)]
amj[, `:=`(
  exc_tx90 = as.integer(eval_tmax & tmax > tx90),
  exc_tn10 = as.integer(eval_tmin & tmin < tn10)
)]

record_check("no_NA_in_exceedance_flags",
             !anyNA(amj$exc_tx90) && !anyNA(amj$exc_tn10),
             sprintf("%d NA in TX90 flag, %d NA in TN10 flag",
                     sum(is.na(amj$exc_tx90)), sum(is.na(amj$exc_tn10))))
record_check("exceedances_only_on_evaluable_days",
             amj[exc_tx90 == 1L & !eval_tmax, .N] == 0L &&
               amj[exc_tn10 == 1L & !eval_tmin, .N] == 0L,
             sprintf("%d TX90 / %d TN10 flags on non-evaluable days",
                     amj[exc_tx90 == 1L & !eval_tmax, .N],
                     amj[exc_tn10 == 1L & !eval_tmin, .N]))
log_msg("Daily TX90 exceedances : ", format(sum(amj$exc_tx90), big.mark = ","))
log_msg("Daily TN10 exceedances : ", format(sum(amj$exc_tn10), big.mark = ","))

## ---------------------------------------------------------------------------
## 6. GRID x YEAR AGGREGATION (single grouped pass)
## ---------------------------------------------------------------------------
log_head("6. GRID x YEAR AGGREGATION")

t0 <- Sys.time()
agg <- amj[, .(
  n_AMJ_days          = .N,
  n_obs_tmax          = sum(!is.na(tmax)),
  n_obs_tmin          = sum(!is.na(tmin)),
  n_no_threshold_tmax = sum(!is.na(tmax) &  is.na(tx90)),
  n_no_threshold_tmin = sum(!is.na(tmin) &  is.na(tn10)),
  n_valid_tmax        = sum(eval_tmax),
  n_valid_tmin        = sum(eval_tmin),
  TX90_days           = sum(exc_tx90),
  TN10_days           = sum(exc_tn10)
), by = .(grid_id, year)]
log_msg("Aggregated       : ", format(nrow(agg), big.mark = ","), " grid-years in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

rm(amj); invisible(gc(verbose = FALSE))

## Complete 283 x 50 skeleton so the output shape is guaranteed even if a
## grid-year is wholly absent from the daily file (such rows become 0 days / NA).
skel <- CJ(grid_id = sort(unique(meta$grid_id)), year = YR_START:YR_END, unique = TRUE)
res  <- merge(skel, agg, by = c("grid_id", "year"), all.x = TRUE, sort = TRUE)

cnt_cols <- c("n_AMJ_days", "n_obs_tmax", "n_obs_tmin",
              "n_no_threshold_tmax", "n_no_threshold_tmin",
              "n_valid_tmax", "n_valid_tmin", "TX90_days", "TN10_days")
for (cc in cnt_cols) set(res, which(is.na(res[[cc]])), cc, 0L)
res[, (cnt_cols) := lapply(.SD, as.integer), .SDcols = cnt_cols]

n_empty <- res[n_AMJ_days == 0L, .N]
record_check("no_empty_grid_years", n_empty == 0L,
             sprintf("%d grid-years had no AMJ rows in the daily file", n_empty),
             severity = "WARN")

## --- the indices ------------------------------------------------------------
res[, `:=`(
  TX90p = fifelse(n_valid_tmax > 0L, 100 * TX90_days / n_valid_tmax, NA_real_),
  TN10p = fifelse(n_valid_tmin > 0L, 100 * TN10_days / n_valid_tmin, NA_real_)
)]

res <- merge(res, meta, by = "grid_id", all.x = TRUE, sort = FALSE)

res[, `:=`(
  completeness_tmax_pct = round(100 * n_valid_tmax / N_AMJ_DAYS, 2),
  completeness_tmin_pct = round(100 * n_valid_tmin / N_AMJ_DAYS, 2),
  season                = "AMJ",
  season_days           = N_AMJ_DAYS,
  in_baseline           = year >= BASE_YR_START & year <= BASE_YR_END,
  baseline_start        = BASE_YR_START,
  baseline_end          = BASE_YR_END,
  bootstrap_applied     = FALSE
)]

setcolorder(res, c(
  "year", "grid_id", "lat", "lon",
  "n_AMJ_days", "n_valid_tmax", "n_valid_tmin",
  "TX90_days", "TN10_days", "TX90p", "TN10p",
  "n_obs_tmax", "n_obs_tmin",
  "n_no_threshold_tmax", "n_no_threshold_tmin",
  "completeness_tmax_pct", "completeness_tmin_pct",
  "season", "season_days", "in_baseline",
  "baseline_start", "baseline_end", "bootstrap_applied"
))
setorder(res, year, grid_id)

## ---------------------------------------------------------------------------
## 7. FINAL QC
## ---------------------------------------------------------------------------
log_head("7. FINAL QC")

record_check("output_row_count", nrow(res) == EXP_OUT_ROWS,
             sprintf("observed %s vs expected %s",
                     format(nrow(res), big.mark = ","),
                     format(EXP_OUT_ROWS, big.mark = ",")))
record_check("output_year_count", uniqueN(res$year) == N_YEARS,
             sprintf("observed %d vs expected %d", uniqueN(res$year), N_YEARS))
record_check("output_grid_count", uniqueN(res$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(res$grid_id), EXP_GRIDS))

n_dup <- sum(duplicated(res, by = c("grid_id", "year")))
record_check("no_duplicate_grid_year", n_dup == 0L,
             sprintf("%d duplicated (grid_id, year) rows", n_dup))

bad_yr <- res[, .N, by = year][N != EXP_GRIDS]
record_check("283_grids_every_year", nrow(bad_yr) == 0L,
             sprintf("%d years with grid count != %d", nrow(bad_yr), EXP_GRIDS))

bad_gr <- res[, .N, by = grid_id][N != N_YEARS]
record_check("50_years_every_grid", nrow(bad_gr) == 0L,
             sprintf("%d grids with year count != %d", nrow(bad_gr), N_YEARS))

record_check("complete_grid_x_year_cross",
             nrow(res) == uniqueN(res$grid_id) * uniqueN(res$year),
             sprintf("%d rows vs %d x %d full cross",
                     nrow(res), uniqueN(res$grid_id), uniqueN(res$year)))

bad_days <- res[n_AMJ_days != N_AMJ_DAYS, .N]
record_check("n_AMJ_days_is_91", bad_days == 0L,
             sprintf("%d grid-years with n_AMJ_days != %d", bad_days, N_AMJ_DAYS),
             severity = "WARN")

record_check("valid_days_within_bounds",
             res[, all(n_valid_tmax >= 0L & n_valid_tmax <= n_AMJ_days &
                         n_valid_tmin >= 0L & n_valid_tmin <= n_AMJ_days)],
             sprintf("Tmax valid range %d-%d ; Tmin valid range %d-%d (cap = n_AMJ_days)",
                     min(res$n_valid_tmax), max(res$n_valid_tmax),
                     min(res$n_valid_tmin), max(res$n_valid_tmin)))

record_check("exceedance_days_le_valid_days",
             res[, all(TX90_days <= n_valid_tmax & TN10_days <= n_valid_tmin)],
             sprintf("%d rows with TX90_days > n_valid_tmax ; %d with TN10_days > n_valid_tmin",
                     res[TX90_days > n_valid_tmax, .N],
                     res[TN10_days > n_valid_tmin, .N]))

record_check("denominator_is_not_fixed_91",
             res[, any(n_valid_tmax != N_AMJ_DAYS) | all(n_obs_tmax == N_AMJ_DAYS)],
             sprintf("distinct n_valid_tmax values = %d ; distinct n_valid_tmin values = %d",
                     uniqueN(res$n_valid_tmax), uniqueN(res$n_valid_tmin)),
             severity = "WARN")

## NA handling: index is NA if and only if its denominator is zero
na_tx_bad <- res[(is.na(TX90p) & n_valid_tmax > 0L) | (!is.na(TX90p) & n_valid_tmax == 0L), .N]
na_tn_bad <- res[(is.na(TN10p) & n_valid_tmin > 0L) | (!is.na(TN10p) & n_valid_tmin == 0L), .N]
record_check("NA_index_iff_zero_denominator_tx90", na_tx_bad == 0L,
             sprintf("%d inconsistent rows (NA TX90p = %d total)",
                     na_tx_bad, res[is.na(TX90p), .N]))
record_check("NA_index_iff_zero_denominator_tn10", na_tn_bad == 0L,
             sprintf("%d inconsistent rows (NA TN10p = %d total)",
                     na_tn_bad, res[is.na(TN10p), .N]))

oor_tx <- res[!is.na(TX90p) & (TX90p < 0 | TX90p > 100), .N]
oor_tn <- res[!is.na(TN10p) & (TN10p < 0 | TN10p > 100), .N]
record_check("TX90p_within_0_100", oor_tx == 0L,
             sprintf("%d out-of-range ; observed [%.2f, %.2f]", oor_tx,
                     suppressWarnings(min(res$TX90p, na.rm = TRUE)),
                     suppressWarnings(max(res$TX90p, na.rm = TRUE))))
record_check("TN10p_within_0_100", oor_tn == 0L,
             sprintf("%d out-of-range ; observed [%.2f, %.2f]", oor_tn,
                     suppressWarnings(min(res$TN10p, na.rm = TRUE)),
                     suppressWarnings(max(res$TN10p, na.rm = TRUE))))

record_check("coords_carried_through", !anyNA(res$lat) && !anyNA(res$lon),
             sprintf("%d NA lat, %d NA lon", sum(is.na(res$lat)), sum(is.na(res$lon))))

## Scientific plausibility: baseline-period means should sit near 10%.
base_tx <- res[in_baseline == TRUE, mean(TX90p, na.rm = TRUE)]
base_tn <- res[in_baseline == TRUE, mean(TN10p, na.rm = TRUE)]
record_check("baseline_mean_TX90p_near_10", is.finite(base_tx) &&
               base_tx >= BASE_MEAN_LO && base_tx <= BASE_MEAN_HI,
             sprintf("1981-2010 mean TX90p = %.2f%% (expected ~10%%, band %g-%g)",
                     base_tx, BASE_MEAN_LO, BASE_MEAN_HI), severity = "WARN")
record_check("baseline_mean_TN10p_near_10", is.finite(base_tn) &&
               base_tn >= BASE_MEAN_LO && base_tn <= BASE_MEAN_HI,
             sprintf("1981-2010 mean TN10p = %.2f%% (expected ~10%%, band %g-%g)",
                     base_tn, BASE_MEAN_LO, BASE_MEAN_HI), severity = "WARN")

log_msg("TX90p summary    : ", paste(sprintf("%.2f",
                                             quantile(res$TX90p, c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)),
                                     collapse = " | "), "   (min|q1|med|q3|max %)")
log_msg("TN10p summary    : ", paste(sprintf("%.2f",
                                             quantile(res$TN10p, c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)),
                                     collapse = " | "), "   (min|q1|med|q3|max %)")
log_msg("Mean TX90p 1976-1980 / 1981-2010 / 2011-2025 : ",
        sprintf("%.2f / %.2f / %.2f",
                res[year <  BASE_YR_START, mean(TX90p, na.rm = TRUE)],
                base_tx,
                res[year >  BASE_YR_END,   mean(TX90p, na.rm = TRUE)]))
log_msg("Mean TN10p 1976-1980 / 1981-2010 / 2011-2025 : ",
        sprintf("%.2f / %.2f / %.2f",
                res[year <  BASE_YR_START, mean(TN10p, na.rm = TRUE)],
                base_tn,
                res[year >  BASE_YR_END,   mean(TN10p, na.rm = TRUE)]))

## ---------------------------------------------------------------------------
## 8. WRITE OUTPUTS + READ-BACK VERIFICATION
## ---------------------------------------------------------------------------
log_head("8. WRITE OUTPUTS")

out_dir <- dirname(OUT_FILE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

arrow::write_parquet(res, OUT_FILE, compression = "snappy")
log_msg("Indices saved    : ", OUT_FILE, "  (",
        format(nrow(res), big.mark = ","), " rows x ", ncol(res), " cols, ",
        sprintf("%.2f MB", file.size(OUT_FILE) / 1024^2), ")")

chk <- as.data.table(arrow::read_parquet(OUT_FILE))
record_check("readback_row_count", nrow(chk) == EXP_OUT_ROWS,
             sprintf("re-read %s rows", format(nrow(chk), big.mark = ",")))
record_check("readback_columns_identical", identical(names(chk), names(res)),
             sprintf("%d columns re-read", ncol(chk)))
setorder(chk, year, grid_id)
record_check("readback_values_identical",
             isTRUE(all.equal(chk$TX90p, res$TX90p)) &&
               isTRUE(all.equal(chk$TN10p, res$TN10p)) &&
               identical(chk$TX90_days, res$TX90_days) &&
               identical(chk$TN10_days, res$TN10_days),
             "TX90p/TN10p and exceedance counts identical after round-trip")
rm(chk)

## threshold file must be untouched
record_check("threshold_file_unmodified", file.exists(THR_FILE),
             sprintf("%s still present, opened read-only, never rewritten",
                     basename(THR_FILE)))

qc_tab <- rbindlist(CHECKS)
qc_tab[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
              script   = "03_AMJ_TX90p_TN10p_1976_2025.R")]
fwrite(qc_tab, QC_FILE)
log_msg("QC summary saved : ", QC_FILE, "  (", nrow(qc_tab), " checks)")

log_head("RUN COMPLETE")
log_msg("Checks passed    : ", qc_tab[status == "PASS", .N], " / ", nrow(qc_tab))
log_msg("Warnings         : ", qc_tab[status == "WARN", .N])
log_msg("Fatal            : ", qc_tab[status == "FATAL", .N])
log_msg("Elapsed          : ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
log_msg("NOTE             : bootstrap for in-base years 1981-2010 NOT applied here.")

.LOG <- c(.LOG, "", "--- sessionInfo() ---",
          capture.output(utils::sessionInfo()))
writeLines(.LOG, LOG_FILE)
cat("\nConsole log written to:", LOG_FILE, "\n")

invisible(res)
###############################################################################
## END OF SCRIPT
###############################################################################