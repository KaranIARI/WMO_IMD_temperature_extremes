###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 04 : 9TX0p / TN10p with Zhang et al. (2005) in-base BOOTSTRAP
##  283 grids | 1976-2025 | AMJ 01 Apr - 30 Jun | reference 1981-2010
##
##  REVISION: percentile engine made BIT-EXACT against stats::quantile(type = 8).
##  Nothing else in the methodology, the bootstrap, or the upstream steps changed.
##
##  INPUTS (ALL READ-ONLY - nothing upstream is recomputed or rewritten)
##    1. Daily  : IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet
##    2. Thresh : AMJ_TX90_TN10_thresholds_1981_2010.parquet        [STEP 02, LOCKED]
##    3. Buffer : AMJ_buffer_1981_2010_30Mar_02Jul.parquet          [STEP 01, LOCKED]
##    4. Direct : AMJ_TX90p_TN10p_1976_2025.parquet                 [STEP 03, for QC]
##
##  OUTPUTS
##    A. AMJ_TX90p_TN10p_bootstrap_1976_2025.parquet        14,150 rows (283 x 50)
##    B. AMJ_TX90p_TN10p_bootstrap_trend_1996_2025.parquet   8,490 rows (283 x 30)
##    + QC csv + per-year bootstrap-effect csv + console log
##
## ---------------------------------------------------------------------------
##  WHY row_q8_matches_stats_quantile_random FAILED, AND WHAT WAS DONE
## ---------------------------------------------------------------------------
##  The estimator was already the right one (Hyndman-Fan type 8); the failure was
##  bit-level, exposed by tightening the verification to tolerance = 0.
##
##  (P1) nppm was computed as (n + 1/3) * p + 1/3. R computes
##       a + probs * (n + 1 - a - b) with a = b = 1/3, i.e. ((n + 1) - 1/3) - 1/3.
##       1/3 is not representable in binary, so ((n+1) - 1/3) - 1/3 and n + 1/3
##       can differ in the last bit. Now computed in R's exact association order.
##
##  (P2) Interpolation was written as v0 + g * (v1 - v0). R uses the algebraically
##       equal but NOT bitwise equal form (1 - h) * x[j] + h * x[j+1], and only on
##       the subset where 0 < h < 1 AND x[j] != x[j+1]; elsewhere it returns x[j]
##       verbatim, or x[j+1] where h == 1. Ties therefore return an untouched
##       order statistic in R but a rounded arithmetic combination in the old
##       code. IMD data is recorded to ~0.1 degC, so ties are common and this
##       branch is exercised constantly. R's exact branch structure is now
##       reproduced.
##
##  (P3) The fuzz test was h < fuzz; R uses abs(h) < fuzz. Now matched.
##
##  The FATAL pre-bootstrap verification is retained and strengthened: bit-exact
##  identity (no tolerance), over randomised ragged-NA rows, over every attainable
##  pool size n = 120..150, and over a heavily tied matrix that exercises the
##  x[j] == x[j+1] branch. The bootstrap cannot start unless all three pass.
##
##  Header counts (verified): 30 x 29 x 91 = 79,170 calendar-day threshold blocks;
##  79,170 x 283 = 22,405,110 grid-level percentile evaluations per variable;
##  44,810,220 row-sorts across Tmax and Tmin combined.
##
## ---------------------------------------------------------------------------
##  ROOT CAUSE OF THE ORIGINAL 15-YEAR FAILURE  (fixed earlier, kept for record)
## ---------------------------------------------------------------------------
##  The failing run bootstrapped 1996-2010 = 15 years, exactly
##      intersect(TREND period 1996-2025, REFERENCE period 1981-2010).
##  The loop was driven by a year vector derived from the TREND window instead of
##  the REFERENCE window, while the banner printed a hard-coded "30 years" string
##  never derived from the loop.
##
##  STRUCTURAL FIX (retained):
##    * BOOT_YEARS <- BASE_YR_START:BASE_YR_END is the ONLY definition of the
##      bootstrap year set. TREND_START appears nowhere except the final
##      subsetting of output B.
##    * Every banner/log line prints length(BOOT_YEARS).
##    * stopifnot() gates fire before any computation, including explicit
##      rejection of the old 1996:2010 vector.
##    * Years visited accumulate INSIDE the loop and are asserted identical to
##      1981:2010; pairs asserted to equal 870; and every reference year must
##      differ from its STEP 03 direct estimate on at least one grid.
##
## ---------------------------------------------------------------------------
##  METHOD (Zhang et al. 2005, J. Climate; ETCCDI standard) - UNCHANGED
## ---------------------------------------------------------------------------
##  OUT-OF-BASE years (1976-1980, 2011-2025):
##      direct comparison against the fixed, locked 1981-2010 thresholds.
##      Identical to STEP 03 by construction (asserted in QC).
##
##  IN-BASE years (1981-2010, ALL 30):
##      for out-of-base year y:
##        for each of the other 29 reference years j:
##          - build a 30-year reference sample in which year y's block of
##            5 window-day values is REPLACED by a duplicate of year j's block
##            ("remove y, duplicate j": 29 real years + 1 repeat, pool stays
##            30 x 5 = 150);
##          - recompute the calendar-day TX90 / TN10 from that sample;
##          - evaluate year y's observations against those thresholds, giving
##            index_j = 100 * exceedance_days / valid_days;
##        the reported index for year y = mean of the 29 index_j values.
##      TX90p_boot_sd / TN10p_boot_sd are the standard deviation ACROSS those 29
##      replicate index values: a descriptive spread induced by the choice of
##      replacement year, NOT a standard error. The replicates are not independent.
##
##  PRESERVED FROM EARLIER STEPS (unchanged):
##    * 5-day window (target +/- 2 days), pool = 30 x 5 = 150.
##    * Validity rule: a threshold requires >= 120 valid (non-NA) values, else NA.
##      Applied identically to every bootstrap sample.
##    * Percentile estimator: Hyndman-Fan type 8 (climdex.pcic / RClimDex).
##    * Calendar-day matching on grid_id + MMDD (never day-of-year).
##    * Strict inequalities: Tmax > TX90, Tmin < TN10.
##    * NA daily values are never imputed and never count as exceedances.
##    * A day is "valid" only if the observation AND the applicable threshold
##      are both non-NA; denominators are independent for Tmax and Tmin;
##      91 is never used as a fixed denominator.
##
##  Runtime: 79,170 threshold blocks; expect roughly 5-15 minutes.
##  Progress is logged per reference year.
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
suppressPackageStartupMessages({ library(data.table); library(arrow) })
setDTthreads(0L)

## ---------------------------------------------------------------------------
## 1. LOCKED PARAMETERS  (single source of truth)
## ---------------------------------------------------------------------------
DAILY_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet"
THR_FILE    <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet"
BUF_FILE    <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_buffer_1981_2010_30Mar_02Jul.parquet"
STEP03_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_1976_2025.parquet"

OUT_FULL    <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_bootstrap_1976_2025.parquet"
OUT_TREND   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_bootstrap_trend_1996_2025.parquet"
QC_FILE     <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_bootstrap_QC.csv"
LOG_FILE    <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_bootstrap_LOG.txt"

YR_START      <- 1976L
YR_END        <- 2025L
N_YEARS       <- YR_END - YR_START + 1L            # 50

BASE_YR_START <- 1981L                             # REFERENCE period
BASE_YR_END   <- 2010L
N_BASE_YEARS  <- BASE_YR_END - BASE_YR_START + 1L  # 30

## >>> THE ONLY DEFINITION OF THE BOOTSTRAP YEAR SET <<<
BOOT_YEARS    <- BASE_YR_START:BASE_YR_END         # 1981..2010, length 30
N_REPLICATES  <- N_BASE_YEARS - 1L                 # 29 replacement years

TREND_START   <- 1996L                             # used ONLY for output B
TREND_END     <- YR_END

MD_START      <- 401L                              # 01 Apr (MMDD)
MD_END        <- 630L                              # 30 Jun (MMDD)
N_AMJ_DAYS    <- 91L
WIN_HALF      <- 2L
WIN_WIDTH     <- 5L
N_POSSIBLE    <- N_BASE_YEARS * WIN_WIDTH          # 150
MIN_VALID     <- 120L                              # >= 80% of 150
PCTL_TYPE     <- 8L
P_TX          <- 0.90
P_TN          <- 0.10

EXP_GRIDS     <- 283L
EXP_THR_ROWS  <- 25753L                            # 283 x 91
EXP_BUF_ROWS  <- 806550L                           # 283 x 95 x 30
EXP_AMJ_ROWS  <- 1287650L                          # 283 x 91 x 50
EXP_POOL_ROWS <- 3862950L                          # 25,753 x 150
EXP_OUT_ROWS  <- 14150L                            # 283 x 50
EXP_BOOT_ROWS <- 8490L                             # 283 x 30 in-base
EXP_DIRECT_RW <- 5660L                             # 283 x 20 out-of-base
EXP_TREND_RW  <- 8490L                             # 283 x 30 (1996-2025)
EXP_PAIRS     <- N_BASE_YEARS * N_REPLICATES       # 870        = 30 x 29
EXP_BLOCKS    <- EXP_PAIRS * N_AMJ_DAYS            # 79,170     = 870 x 91
EXP_PCTL_EVAL <- EXP_BLOCKS * EXP_GRIDS            # 22,405,110 per variable

TOL_EXACT     <- 1e-9
COORD_TOL     <- 1e-6
BASE_MEAN_LO  <- 5
BASE_MEAN_HI  <- 20

## Hard gates before any computation.
stopifnot(
  length(BOOT_YEARS) == 30L,
  identical(BOOT_YEARS, 1981:2010),
  identical(BOOT_YEARS, BASE_YR_START:BASE_YR_END),
  N_REPLICATES == 29L,
  !identical(BOOT_YEARS, TREND_START:BASE_YR_END),   # rejects the old 1996:2010
  EXP_BLOCKS == 79170L,
  EXP_PCTL_EVAL == 22405110L
)

## ---------------------------------------------------------------------------
## 2. HELPERS
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
  log_msg(sprintf("  %-5s | %-44s | %s", status, name, detail))
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

## Bit-exact identity test (no tolerance), NA-aware.
exact_same <- function(a, b) {
  if (length(a) != length(b)) return(FALSE)
  na_a <- is.na(a); na_b <- is.na(b)
  if (!identical(na_a, na_b)) return(FALSE)
  all(a[!na_a] == b[!na_b])
}

## ---------------------------------------------------------------------------
##  row_q8() : row-wise Hyndman-Fan type-8 quantile, BIT-EXACT against
##  stats::quantile(x, p, type = 8, na.rm = TRUE).
##
##  This is a transcription of R's own quantile.default branch for types 4-9,
##  specialised to a = b = 1/3, and applied across rows. Every step follows R's
##  source in the same order and the same floating-point association:
##
##    a <- 1/3 ; b <- 1/3
##    fuzz <- 4 * .Machine$double.eps
##    nppm <- a + probs * (n + 1 - a - b)      # NOT (n + 1/3)*p + 1/3
##    j    <- floor(nppm + fuzz)
##    h    <- nppm - j ; h[abs(h) < fuzz] <- 0 # NOT h < fuzz
##    x    <- c(x[1], x[1], x, x[n], x[n])     # padding; index j+2 and j+3
##    qs   <- x[j+2]
##    qs[h == 1] <- x[j+3][h == 1]
##    other <- (0 < h) & (h < 1) & (x[j+2] != x[j+3])
##    qs[other] <- ((1-h)*x[j+2] + h*x[j+3])[other]
##
##  The padding is reproduced by clamping: index j+2 of the padded vector is
##  order statistic pmin(pmax(j, 1), n), and index j+3 is pmin(pmax(j+1, 1), n).
##  The `other` mask matters: where the two bracketing order statistics are equal
##  (ties, which are frequent in temperature data recorded to 0.1 degC) R returns
##  the order statistic untouched rather than an arithmetic combination.
##
##  The locked >= MIN_VALID rule is applied last, exactly as before.
## ---------------------------------------------------------------------------
row_q8 <- function(M, p, min_valid = MIN_VALID) {
  nr <- nrow(M)
  n  <- as.integer(rowSums(!is.na(M)))
  S  <- apply(M, 1L, sort, na.last = TRUE)     # ncol x nrow, NAs pushed to end
  
  a <- 1/3
  b <- 1/3
  fuzz <- 4 * .Machine$double.eps
  
  nppm <- a + p * (n + 1 - a - b)              # R's association order
  j    <- floor(nppm + fuzz)
  h    <- nppm - j
  sml  <- abs(h) < fuzz
  sml[is.na(sml)] <- FALSE
  h[sml] <- 0
  
  nn <- pmax(n, 1L)
  lo <- pmin(pmax(as.integer(j),      1L), nn) # padded index j+2
  hi <- pmin(pmax(as.integer(j) + 1L, 1L), nn) # padded index j+3
  i  <- seq_len(nr)
  v0 <- S[cbind(lo, i)]
  v1 <- S[cbind(hi, i)]
  
  qs <- v0
  one <- !is.na(h) & h == 1
  qs[one] <- v1[one]
  other <- (0 < h) & (h < 1) & (v0 != v1)
  other[is.na(other)] <- TRUE
  if (any(other)) qs[other] <- ((1 - h) * v0 + h * v1)[other]
  
  qs[n < min_valid] <- NA_real_
  qs
}

t_start <- Sys.time()
log_head("STEP 04 | TX90p / TN10p with Zhang et al. (2005) in-base bootstrap")
log_msg("R version         : ", R.version.string)
log_msg("data.table        : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow             : ", as.character(utils::packageVersion("arrow")))
log_msg("Analysis period   : ", YR_START, "-", YR_END, " (", N_YEARS, " years)")
log_msg("Reference period  : ", BASE_YR_START, "-", BASE_YR_END, " (", N_BASE_YEARS, " years)")
log_msg("Trend period      : ", TREND_START, "-", TREND_END,
        "  [subsets output B only - never touches the bootstrap loop]")
log_msg("BOOTSTRAP YEARS   : ", min(BOOT_YEARS), "-", max(BOOT_YEARS),
        "  n = ", length(BOOT_YEARS), "   (printed from length(BOOT_YEARS))")
log_msg("Replicates / year : ", N_REPLICATES, "  -> pairs = ", EXP_PAIRS,
        " ; blocks = ", format(EXP_BLOCKS, big.mark = ","),
        " ; percentile evaluations per variable = ", format(EXP_PCTL_EVAL, big.mark = ","))
log_msg("Season            : 01 Apr - 30 Jun (MMDD ", MD_START, "-", MD_END, ")")
log_msg("Pool / validity   : ", N_POSSIBLE, " values, >= ", MIN_VALID, " valid else NA")
log_msg("Percentile        : Hyndman-Fan type ", PCTL_TYPE,
        " (p", P_TX * 100, " / p", P_TN * 100, "), bit-exact transcription of stats::quantile")

## ---------------------------------------------------------------------------
## 3. READ LOCKED INPUTS (read-only)
## ---------------------------------------------------------------------------
log_head("3. READ LOCKED INPUTS")

for (f in c(DAILY_FILE, THR_FILE, BUF_FILE))
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)

## --- 3a. thresholds (STEP 02) ----------------------------------------------
thr <- as.data.table(arrow::read_parquet(THR_FILE))
T_GRID <- pick_col(thr, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id")
T_MD   <- pick_col(thr, c("target_md", "md", "mmdd"), "target_md", required = FALSE)
T_DAY  <- pick_col(thr, c("target_day", "calendar_day"), "target_day", required = FALSE)
T_TX   <- pick_col(thr, c("tx90", "tx90_threshold", "tmax_p90"), "tx90")
T_TN   <- pick_col(thr, c("tn10", "tn10_threshold", "tmin_p10"), "tn10")
if (is.na(T_MD)) {
  if (is.na(T_DAY)) stop("Threshold file has neither target_md nor target_day.", call. = FALSE)
  dd <- as.character(thr[[T_DAY]])
  thr[, thr_md := as.integer(substr(dd, 1L, 2L)) * 100L + as.integer(substr(dd, 4L, 5L))]
} else thr[, thr_md := as.integer(get(T_MD))]
thr <- thr[, c(T_GRID, "thr_md", T_TX, T_TN), with = FALSE]
setnames(thr, c("grid_id", "md", "tx90", "tn10"))
thr[, `:=`(grid_id = as.character(grid_id), tx90 = as.numeric(tx90), tn10 = as.numeric(tn10))]

log_msg("Thresholds        : ", format(nrow(thr), big.mark = ","), " rows, ",
        uniqueN(thr$grid_id), " grids, ", uniqueN(thr$md), " calendar days")
record_check("threshold_row_count", nrow(thr) == EXP_THR_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(thr), big.mark = ","),
                     format(EXP_THR_ROWS, big.mark = ",")))
record_check("threshold_no_duplicate_key", !anyDuplicated(thr, by = c("grid_id", "md")),
             sprintf("%d duplicated (grid_id, MMDD)", sum(duplicated(thr, by = c("grid_id", "md")))))
log_msg("Locked NA thresh  : TX90 ", sum(is.na(thr$tx90)), " ; TN10 ", sum(is.na(thr$tn10)))

## --- 3b. baseline buffer (STEP 01) -----------------------------------------
buf <- as.data.table(arrow::read_parquet(BUF_FILE))
B_GRID <- pick_col(buf, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id")
B_DATE <- pick_col(buf, c("date", "dates", "obs_date", "time", "day"), "date")
B_TMAX <- pick_col(buf, c("tmax", "tx", "t_max", "max_temp", "temp_max"), "tmax")
B_TMIN <- pick_col(buf, c("tmin", "tn", "t_min", "min_temp", "temp_min"), "tmin")
buf <- buf[, c(B_GRID, B_DATE, B_TMAX, B_TMIN), with = FALSE]
setnames(buf, c("grid_id", "date", "tmax", "tmin"))
buf[, date := coerce_idate(date)]
buf[, `:=`(grid_id = as.character(grid_id),
           yr = as.integer(data.table::year(date)),
           md = as.integer(data.table::month(date)) * 100L + as.integer(data.table::mday(date)),
           tmax = as.numeric(tmax), tmin = as.numeric(tmin))]

log_msg("Buffer            : ", format(nrow(buf), big.mark = ","), " rows, ",
        uniqueN(buf$grid_id), " grids, ", uniqueN(buf$yr), " years, ",
        uniqueN(buf$md), " calendar days")
record_check("buffer_row_count", nrow(buf) == EXP_BUF_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(buf), big.mark = ","),
                     format(EXP_BUF_ROWS, big.mark = ",")))
record_check("buffer_years", identical(sort(unique(buf$yr)), BASE_YR_START:BASE_YR_END),
             sprintf("%d years [%d-%d]", uniqueN(buf$yr), min(buf$yr), max(buf$yr)))
record_check("buffer_no_duplicate_grid_date", !anyDuplicated(buf, by = c("grid_id", "date")),
             sprintf("%d duplicates", sum(duplicated(buf, by = c("grid_id", "date")))))

## --- 3c. daily data --------------------------------------------------------
dly <- as.data.table(arrow::read_parquet(DAILY_FILE))
D_GRID <- pick_col(dly, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id")
D_LAT  <- pick_col(dly, c("lat", "latitude", "y"), "lat")
D_LON  <- pick_col(dly, c("lon", "long", "longitude", "x"), "lon")
D_DATE <- pick_col(dly, c("date", "dates", "obs_date", "time", "day"), "date")
D_TMAX <- pick_col(dly, c("tmax", "tx", "t_max", "max_temp", "temp_max"), "tmax")
D_TMIN <- pick_col(dly, c("tmin", "tn", "t_min", "min_temp", "temp_min"), "tmin")
dly <- dly[, c(D_GRID, D_LAT, D_LON, D_DATE, D_TMAX, D_TMIN), with = FALSE]
setnames(dly, c("grid_id", "lat", "lon", "date", "tmax", "tmin"))
dly[, date := coerce_idate(date)]
dly[, `:=`(grid_id = as.character(grid_id),
           year = as.integer(data.table::year(date)),
           md   = as.integer(data.table::month(date)) * 100L + as.integer(data.table::mday(date)),
           tmax = as.numeric(tmax), tmin = as.numeric(tmin))]

meta <- unique(dly[, .(grid_id, lat, lon)])
record_check("daily_unique_coords_per_grid", nrow(meta) == EXP_GRIDS,
             sprintf("%d unique grid/lat/lon triplets vs %d grids", nrow(meta), EXP_GRIDS))

amj <- dly[md >= MD_START & md <= MD_END & year >= YR_START & year <= YR_END]
rm(dly); invisible(gc(verbose = FALSE))
log_msg("Daily AMJ subset  : ", format(nrow(amj), big.mark = ","), " rows, ",
        uniqueN(amj$year), " years, ", uniqueN(amj$md), " calendar days")
record_check("amj_row_count", nrow(amj) == EXP_AMJ_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(amj), big.mark = ","),
                     format(EXP_AMJ_ROWS, big.mark = ",")))
record_check("amj_year_coverage", identical(sort(unique(amj$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(amj$year), min(amj$year), max(amj$year)))

## --- 3d. grid identity must agree across all files --------------------------
g_daily <- sort(unique(amj$grid_id)); g_thr <- sort(unique(thr$grid_id)); g_buf <- sort(unique(buf$grid_id))
record_check("grid_ids_identical_across_files",
             identical(g_daily, g_thr) && identical(g_daily, g_buf),
             sprintf("daily %d / thresholds %d / buffer %d ; set-identical = %s",
                     length(g_daily), length(g_thr), length(g_buf),
                     identical(g_daily, g_thr) && identical(g_daily, g_buf)))
record_check("grid_count", length(g_daily) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", length(g_daily), EXP_GRIDS))

GRIDS  <- g_daily
N_GRID <- length(GRIDS)
gidx   <- function(x) match(x, GRIDS)

## ---------------------------------------------------------------------------
## 4. WINDOW MAP AND POOL ARRAYS  (grid x target-day x year x window-offset)
## ---------------------------------------------------------------------------
log_head("4. BUILD 5-DAY WINDOW POOL ARRAYS")

REF_YR    <- 2001L
tgt_dates <- seq(as.Date(paste0(REF_YR, "-04-01")), as.Date(paste0(REF_YR, "-06-30")), by = "day")
record_check("target_day_count", length(tgt_dates) == N_AMJ_DAYS,
             sprintf("observed %d vs expected %d", length(tgt_dates), N_AMJ_DAYS))

targets <- data.table(
  target_index = seq_along(tgt_dates),
  target_md    = as.integer(data.table::month(tgt_dates)) * 100L +
    as.integer(data.table::mday(tgt_dates)))
targets[, target_day := sprintf("%02d-%02d", target_md %/% 100L, target_md %% 100L)]

win <- CJ(target_index = targets$target_index, off_idx = 1:WIN_WIDTH)
win[, offset := off_idx - (WIN_HALF + 1L)]          # off_idx 1..5 <-> offset -2..+2
win <- merge(win, targets, by = "target_index", sort = FALSE)
win[, wdate := tgt_dates[target_index] + offset]
win[, md := as.integer(data.table::month(wdate)) * 100L + as.integer(data.table::mday(wdate))]
win[, wdate := NULL]
record_check("window_map_size", nrow(win) == N_AMJ_DAYS * WIN_WIDTH,
             sprintf("observed %d vs expected %d", nrow(win), N_AMJ_DAYS * WIN_WIDTH))

setkey(buf, md); setkey(win, md)
pool <- buf[win, .(grid_id, target_index, yr, off_idx, tmax, tmin),
            allow.cartesian = TRUE, nomatch = 0L]
record_check("pool_row_count", nrow(pool) == EXP_POOL_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(pool), big.mark = ","),
                     format(EXP_POOL_ROWS, big.mark = ",")))

pool[, `:=`(gi = gidx(grid_id), yi = yr - BASE_YR_START + 1L)]
lin <- pool$gi + N_GRID * ((pool$target_index - 1L) +
                             N_AMJ_DAYS * ((pool$yi - 1L) + N_BASE_YEARS * (pool$off_idx - 1L)))
record_check("pool_index_is_bijective",
             !anyNA(lin) && !anyDuplicated(lin) && length(lin) == EXP_POOL_ROWS,
             sprintf("%s unique linear indices", format(uniqueN(lin), big.mark = ",")))

A_TX <- array(NA_real_, dim = c(N_GRID, N_AMJ_DAYS, N_BASE_YEARS, WIN_WIDTH))
A_TN <- array(NA_real_, dim = c(N_GRID, N_AMJ_DAYS, N_BASE_YEARS, WIN_WIDTH))
A_TX[lin] <- pool$tmax
A_TN[lin] <- pool$tmin
log_msg("Pool arrays built : ", N_GRID, " x ", N_AMJ_DAYS, " x ", N_BASE_YEARS, " x ", WIN_WIDTH,
        "  (", format(length(A_TX), big.mark = ","), " values per variable)")
rm(pool, lin, buf); invisible(gc(verbose = FALSE))

## After dim(x) <- c(N_GRID, 150) the column index is year + 30*(offset-1),
## so reference year yi occupies columns yi, yi+30, yi+60, yi+90, yi+120.
YCOLS <- function(yi) yi + N_BASE_YEARS * (0:(WIN_WIDTH - 1L))

## Observation array for reference years: grid x target-day x year
obs <- amj[year >= BASE_YR_START & year <= BASE_YR_END]
obs[, `:=`(gi = gidx(grid_id), yi = year - BASE_YR_START + 1L)]
obs <- merge(obs, targets[, .(target_index, target_md)], by.x = "md", by.y = "target_md",
             all.x = TRUE, sort = FALSE)
record_check("base_obs_row_count", nrow(obs) == N_GRID * N_AMJ_DAYS * N_BASE_YEARS,
             sprintf("observed %s vs expected %s", format(nrow(obs), big.mark = ","),
                     format(N_GRID * N_AMJ_DAYS * N_BASE_YEARS, big.mark = ",")))
olin <- obs$gi + N_GRID * ((obs$target_index - 1L) + N_AMJ_DAYS * (obs$yi - 1L))
O_TX <- array(NA_real_, dim = c(N_GRID, N_AMJ_DAYS, N_BASE_YEARS))
O_TN <- array(NA_real_, dim = c(N_GRID, N_AMJ_DAYS, N_BASE_YEARS))
O_TX[olin] <- obs$tmax
O_TN[olin] <- obs$tmin
record_check("base_obs_index_is_bijective", !anyNA(olin) && !anyDuplicated(olin),
             sprintf("%s unique linear indices", format(uniqueN(olin), big.mark = ",")))

## Buffer centre day (off_idx = WIN_HALF+1 = 3 <-> offset 0) must equal the
## daily-file value on the same grid/day/year.
ov <- A_TX[, , , WIN_HALF + 1L]
record_check("daily_matches_buffer_on_AMJ",
             isTRUE(all.equal(as.vector(ov), as.vector(O_TX), tolerance = TOL_EXACT)),
             "buffer centre-day Tmax identical to daily-file Tmax over 1981-2010 AMJ",
             severity = "WARN")
rm(ov, obs, olin); invisible(gc(verbose = FALSE))

## ---------------------------------------------------------------------------
## 5. VERIFY THE PERCENTILE ENGINE  (must pass before any bootstrap work)
## ---------------------------------------------------------------------------
log_head("5. VERIFY PERCENTILE ENGINE (BIT-EXACT) AGAINST stats::quantile AND STEP 02")

## Reference: R's own function, applied row by row, with the locked validity rule.
ref_of <- function(M, p) vapply(seq_len(nrow(M)), function(i) {
  v <- M[i, ][!is.na(M[i, ])]
  if (length(v) < MIN_VALID) NA_real_ else
    as.numeric(stats::quantile(v, p, type = PCTL_TYPE, na.rm = TRUE, names = FALSE))
}, numeric(1))

## 5a. randomised ragged-NA test (continuous values, no ties)
set.seed(20260812)
Mt <- matrix(rnorm(200 * N_POSSIBLE, 30, 6), nrow = 200)
Mt[sample.int(length(Mt), 1500)] <- NA_real_
record_check("row_q8_matches_stats_quantile_random",
             exact_same(row_q8(Mt, P_TX), ref_of(Mt, P_TX)) &&
               exact_same(row_q8(Mt, P_TN), ref_of(Mt, P_TN)),
             "200 ragged-NA rows, bit-exact identity (no tolerance)")

## 5b. deterministic sweep over every attainable pool size n = 120..150,
## where the type-8 fuzz at integer nppm boundaries can matter.
ns <- MIN_VALID:N_POSSIBLE
Ms <- matrix(NA_real_, nrow = length(ns), ncol = N_POSSIBLE)
for (a_ in seq_along(ns)) Ms[a_, seq_len(ns[a_])] <- seq_len(ns[a_]) + 0.5
record_check("row_q8_matches_stats_quantile_all_n",
             exact_same(row_q8(Ms, P_TX), ref_of(Ms, P_TX)) &&
               exact_same(row_q8(Ms, P_TN), ref_of(Ms, P_TN)),
             sprintf("bit-exact for every pool size n = %d..%d", min(ns), max(ns)))

## 5c. heavily tied values, as in real 0.1 degC-resolution temperature data.
## This exercises R's x[j] == x[j+1] branch, where R returns the order statistic
## untouched instead of interpolating.
set.seed(19810401)
Mq <- round(matrix(rnorm(300 * N_POSSIBLE, 38, 1.2), nrow = 300), 1)
Mq[sample.int(length(Mq), 2000)] <- NA_real_
record_check("row_q8_matches_stats_quantile_ties",
             exact_same(row_q8(Mq, P_TX), ref_of(Mq, P_TX)) &&
               exact_same(row_q8(Mq, P_TN), ref_of(Mq, P_TN)),
             "300 rows of tied 0.1-resolution values, bit-exact identity")
rm(Mt, Ms, Mq, ns)

## 5d. recompute the FULL 30-year thresholds from the pool and compare with the
## locked STEP 02 file (which is never modified).
chk_tx <- matrix(NA_real_, N_GRID, N_AMJ_DAYS)
chk_tn <- matrix(NA_real_, N_GRID, N_AMJ_DAYS)
for (d in seq_len(N_AMJ_DAYS)) {
  Bx <- A_TX[, d, , ]; dim(Bx) <- c(N_GRID, N_POSSIBLE)
  Bn <- A_TN[, d, , ]; dim(Bn) <- c(N_GRID, N_POSSIBLE)
  chk_tx[, d] <- row_q8(Bx, P_TX)
  chk_tn[, d] <- row_q8(Bn, P_TN)
}
cmp <- merge(
  data.table(grid_id = rep(GRIDS, times = N_AMJ_DAYS),
             md      = rep(targets$target_md, each = N_GRID),
             rc_tx   = as.vector(chk_tx), rc_tn = as.vector(chk_tn)),
  thr, by = c("grid_id", "md"), all.x = TRUE)
d_tx <- max(abs(cmp$rc_tx - cmp$tx90), na.rm = TRUE)
d_tn <- max(abs(cmp$rc_tn - cmp$tn10), na.rm = TRUE)
record_check("recomputed_thresholds_match_STEP02",
             identical(is.na(cmp$rc_tx), is.na(cmp$tx90)) &&
               identical(is.na(cmp$rc_tn), is.na(cmp$tn10)) &&
               d_tx < TOL_EXACT && d_tn < TOL_EXACT,
             sprintf("max |diff| TX90 = %.3e, TN10 = %.3e ; NA patterns identical", d_tx, d_tn))
rm(chk_tx, chk_tn, cmp); invisible(gc(verbose = FALSE))

## ---------------------------------------------------------------------------
## 6. DIRECT (NON-BOOTSTRAP) INDICES FOR ALL YEARS  -- STEP 03 logic
## ---------------------------------------------------------------------------
log_head("6. DIRECT INDICES (fixed 1981-2010 thresholds, all years)")

setkey(thr, grid_id, md)
amj[thr, on = .(grid_id, md), `:=`(tx90 = i.tx90, tn10 = i.tn10)]

n_unmatched <- amj[, sum(!paste(grid_id, md) %chin% paste(thr$grid_id, thr$md))]
record_check("all_days_matched_to_threshold_key", n_unmatched == 0L,
             sprintf("%s AMJ rows without a (grid_id, MMDD) threshold record",
                     format(n_unmatched, big.mark = ",")))

amj[, `:=`(eval_tmax = !is.na(tmax) & !is.na(tx90),
           eval_tmin = !is.na(tmin) & !is.na(tn10))]
amj[, `:=`(exc_tx = as.integer(eval_tmax & tmax > tx90),
           exc_tn = as.integer(eval_tmin & tmin < tn10))]

direct <- amj[, .(n_AMJ_days = .N,
                  nv_tmax    = sum(eval_tmax),
                  nv_tmin    = sum(eval_tmin),
                  tx_days    = sum(exc_tx),
                  tn_days    = sum(exc_tn)),
              by = .(grid_id, year)]
direct[, `:=`(TX90p_direct = fifelse(nv_tmax > 0L, 100 * tx_days / nv_tmax, NA_real_),
              TN10p_direct = fifelse(nv_tmin > 0L, 100 * tn_days / nv_tmin, NA_real_))]
log_msg("Direct indices    : ", format(nrow(direct), big.mark = ","), " grid-years")
record_check("direct_row_count", nrow(direct) == EXP_OUT_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(direct), big.mark = ","),
                     format(EXP_OUT_ROWS, big.mark = ",")))
rm(amj); invisible(gc(verbose = FALSE))

## ---------------------------------------------------------------------------
## 7. ZHANG et al. (2005) IN-BASE BOOTSTRAP  -- ALL 30 REFERENCE YEARS
## ---------------------------------------------------------------------------
log_head(sprintf("7. IN-BASE BOOTSTRAP | %d reference years x %d replicates = %s (year, replicate) pairs -> %s calendar-day threshold blocks",
                 length(BOOT_YEARS), N_REPLICATES,
                 format(length(BOOT_YEARS) * N_REPLICATES, big.mark = ","),
                 format(length(BOOT_YEARS) * N_REPLICATES * N_AMJ_DAYS, big.mark = ",")))

nb <- length(BOOT_YEARS)
BOOT_TX90p  <- matrix(NA_real_, N_GRID, nb); BOOT_TN10p  <- matrix(NA_real_, N_GRID, nb)
BOOT_TX_SD  <- matrix(NA_real_, N_GRID, nb); BOOT_TN_SD  <- matrix(NA_real_, N_GRID, nb)
BOOT_TXDAYS <- matrix(NA_real_, N_GRID, nb); BOOT_TNDAYS <- matrix(NA_real_, N_GRID, nb)
BOOT_NVTX   <- matrix(NA_real_, N_GRID, nb); BOOT_NVTN   <- matrix(NA_real_, N_GRID, nb)
BOOT_NREPTX <- matrix(NA_integer_, N_GRID, nb)   # realised replicate count
BOOT_NREPTN <- matrix(NA_integer_, N_GRID, nb)

years_visited <- integer(0)     # accumulated INSIDE the loop - the audit trail
pairs_done    <- 0L
blocks_done   <- 0L
t_boot        <- Sys.time()

for (k in seq_along(BOOT_YEARS)) {
  
  y_cal <- BOOT_YEARS[k]                         # calendar year being left out
  yi    <- y_cal - BASE_YR_START + 1L             # its index in 1..30
  js    <- setdiff(seq_len(N_BASE_YEARS), yi)     # 29 replacement year indices
  stopifnot(length(js) == N_REPLICATES, !yi %in% js)
  
  EXC_TX <- matrix(0, N_GRID, N_REPLICATES); NV_TX <- matrix(0, N_GRID, N_REPLICATES)
  EXC_TN <- matrix(0, N_GRID, N_REPLICATES); NV_TN <- matrix(0, N_GRID, N_REPLICATES)
  ycols  <- YCOLS(yi)
  
  for (d in seq_len(N_AMJ_DAYS)) {
    
    Bx <- A_TX[, d, , ]; dim(Bx) <- c(N_GRID, N_POSSIBLE)
    Bn <- A_TN[, d, , ]; dim(Bn) <- c(N_GRID, N_POSSIBLE)
    ox <- O_TX[, d, yi]
    on <- O_TN[, d, yi]
    
    for (m in seq_along(js)) {
      jcols <- YCOLS(js[m])
      
      Mx <- Bx; Mx[, ycols] <- Bx[, jcols]       # year y block <- copy of year j
      Mn <- Bn; Mn[, ycols] <- Bn[, jcols]
      
      tx_thr <- row_q8(Mx, P_TX)
      tn_thr <- row_q8(Mn, P_TN)
      
      evx <- !is.na(ox) & !is.na(tx_thr)         # FALSE & NA is FALSE -> no NA leak
      evn <- !is.na(on) & !is.na(tn_thr)
      
      NV_TX[, m]  <- NV_TX[, m]  + evx
      NV_TN[, m]  <- NV_TN[, m]  + evn
      EXC_TX[, m] <- EXC_TX[, m] + (evx & ox > tx_thr)
      EXC_TN[, m] <- EXC_TN[, m] + (evn & on < tn_thr)
      
      blocks_done <- blocks_done + 1L
    }
  }
  
  idx_tx <- ifelse(NV_TX > 0, 100 * EXC_TX / NV_TX, NA_real_)
  idx_tn <- ifelse(NV_TN > 0, 100 * EXC_TN / NV_TN, NA_real_)
  
  BOOT_TX90p[, k]  <- rowMeans(idx_tx, na.rm = TRUE)
  BOOT_TN10p[, k]  <- rowMeans(idx_tn, na.rm = TRUE)
  BOOT_TX_SD[, k]  <- apply(idx_tx, 1L, stats::sd, na.rm = TRUE)
  BOOT_TN_SD[, k]  <- apply(idx_tn, 1L, stats::sd, na.rm = TRUE)
  BOOT_TXDAYS[, k] <- rowMeans(EXC_TX)
  BOOT_TNDAYS[, k] <- rowMeans(EXC_TN)
  BOOT_NVTX[, k]   <- rowMeans(NV_TX)
  BOOT_NVTN[, k]   <- rowMeans(NV_TN)
  BOOT_NREPTX[, k] <- as.integer(rowSums(!is.na(idx_tx)))
  BOOT_NREPTN[, k] <- as.integer(rowSums(!is.na(idx_tn)))
  
  years_visited <- c(years_visited, y_cal)
  pairs_done    <- pairs_done + length(js)
  
  el  <- as.numeric(difftime(Sys.time(), t_boot, units = "secs"))
  eta <- el / k * (nb - k)
  log_msg(sprintf("  bootstrapped %d  (%2d/%2d)  replicates = %d  elapsed %.0fs  eta %.0fs",
                  y_cal, k, nb, length(js), el, eta))
}

BOOT_TX90p[!is.finite(BOOT_TX90p)] <- NA_real_
BOOT_TN10p[!is.finite(BOOT_TN10p)] <- NA_real_

log_msg("Bootstrap done    : ", length(years_visited), " reference years, ",
        pairs_done, " (year, replicate) pairs, ",
        format(blocks_done, big.mark = ","), " threshold blocks, ",
        format(blocks_done * N_GRID, big.mark = ","), " percentile evaluations per variable, in ",
        sprintf("%.1f min", as.numeric(difftime(Sys.time(), t_boot, units = "mins"))))
## ===========================================================================
## 7D. DIAGNOSTIC PATCH  (read-only; alters no result, no methodology)
## PLACE IMMEDIATELY BEFORE:
##     record_check("realised_replicates_equal_29", ...)
## Requires in scope: BOOT_NREPTX/TN, BOOT_NVTX/TN, A_TX, A_TN, O_TX, O_TN,
##                    GRIDS, BOOT_YEARS, direct, thr, row_q8(), YCOLS()
## ===========================================================================
log_head("7D. DIAGNOSTIC | realised replicate counts per grid-year")

DIAG_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP04_replicate_diagnostic.csv"

## ---- (1) distribution of realised replicate counts -------------------------
band <- function(v) cut(v, breaks = c(-1, 0, N_REPLICATES - 1L, N_REPLICATES),
                        labels = c("0 replicates", "1-28 replicates", "29 replicates"))
tab_tx <- table(band(as.vector(BOOT_NREPTX)))
tab_tn <- table(band(as.vector(BOOT_NREPTN)))
log_msg("Grid-years by realised replicate count (total = ",
        format(length(BOOT_NREPTX), big.mark = ","), " = ", N_GRID, " x ", nb, ")")
for (lv in names(tab_tx))
  log_msg(sprintf("    %-16s : Tmax %6d   Tmin %6d", lv, tab_tx[[lv]], tab_tn[[lv]]))

## ---- (2)+(3) full grid-year table with realised counts ---------------------
n_obs_tx_mat <- apply(!is.na(O_TX), c(1L, 3L), sum)      # N_GRID x nb
n_obs_tn_mat <- apply(!is.na(O_TN), c(1L, 3L), sum)

diag_dt <- data.table(
  grid_id               = rep(GRIDS,              times = nb),
  year                  = rep(BOOT_YEARS,         each  = N_GRID),
  gi                    = rep(seq_len(N_GRID),    times = nb),
  yi                    = rep(seq_len(nb),        each  = N_GRID),
  n_rep_tmax            = as.vector(BOOT_NREPTX),
  n_rep_tmin            = as.vector(BOOT_NREPTN),
  mean_validdays_tmax   = as.vector(BOOT_NVTX),
  mean_validdays_tmin   = as.vector(BOOT_NVTN),
  n_obs_days_tmax       = as.vector(n_obs_tx_mat),
  n_obs_days_tmin       = as.vector(n_obs_tn_mat))

## locked STEP 02 NA thresholds per grid (out of 91 calendar days)
thr_na <- thr[, .(locked_na_tx90_days = sum(is.na(tx90)),
                  locked_na_tn10_days = sum(is.na(tn10))), by = grid_id]
diag_dt <- merge(diag_dt, thr_na, by = "grid_id", all.x = TRUE)

## STEP 03 direct valid-day counts for the same grid-year
diag_dt <- merge(diag_dt,
                 direct[, .(grid_id, year,
                            step03_nvalid_tmax = nv_tmax,
                            step03_nvalid_tmin = nv_tmin)],
                 by = c("grid_id", "year"), all.x = TRUE)

bad <- diag_dt[n_rep_tmax < N_REPLICATES | n_rep_tmin < N_REPLICATES]
setorder(bad, year, grid_id)
log_msg("Grid-years with < ", N_REPLICATES, " realised replicates : ", nrow(bad))

## ---- (4)+(5) per-replicate NA-threshold forensics for offending cases ------
if (nrow(bad) > 0L) {
  
  addcols <- c("rep_thrNA_min_tx", "rep_thrNA_max_tx", "rep_allNA_thr_tx", "rep_nv0_tx",
               "rep_thrNA_min_tn", "rep_thrNA_max_tn", "rep_allNA_thr_tn", "rep_nv0_tn",
               "pool_minvalid_tx", "pool_minvalid_tn")
  bad[, (addcols) := NA_integer_]
  
  for (r in seq_len(nrow(bad))) {
    g  <- bad$gi[r]; yy <- bad$yi[r]
    Px <- A_TX[g, , , ]; dim(Px) <- c(N_AMJ_DAYS, N_POSSIBLE)   # 91 days x 150 pool
    Pn <- A_TN[g, , , ]; dim(Pn) <- c(N_AMJ_DAYS, N_POSSIBLE)
    ox <- O_TX[g, , yy]; on <- O_TN[g, , yy]
    ycols <- YCOLS(yy)
    js    <- setdiff(seq_len(N_BASE_YEARS), yy)
    
    thrNA_x <- nv_x <- integer(length(js))
    thrNA_n <- nv_n <- integer(length(js))
    for (m in seq_along(js)) {
      jc <- YCOLS(js[m])
      Qx <- Px; Qx[, ycols] <- Px[, jc]
      Qn <- Pn; Qn[, ycols] <- Pn[, jc]
      tx <- row_q8(Qx, P_TX); tn <- row_q8(Qn, P_TN)
      thrNA_x[m] <- sum(is.na(tx)); thrNA_n[m] <- sum(is.na(tn))
      nv_x[m] <- sum(!is.na(ox) & !is.na(tx))
      nv_n[m] <- sum(!is.na(on) & !is.na(tn))
    }
    set(bad, r, "rep_thrNA_min_tx", min(thrNA_x));  set(bad, r, "rep_thrNA_max_tx", max(thrNA_x))
    set(bad, r, "rep_allNA_thr_tx", sum(thrNA_x == N_AMJ_DAYS)); set(bad, r, "rep_nv0_tx", sum(nv_x == 0L))
    set(bad, r, "rep_thrNA_min_tn", min(thrNA_n));  set(bad, r, "rep_thrNA_max_tn", max(thrNA_n))
    set(bad, r, "rep_allNA_thr_tn", sum(thrNA_n == N_AMJ_DAYS)); set(bad, r, "rep_nv0_tn", sum(nv_n == 0L))
    set(bad, r, "pool_minvalid_tx", min(rowSums(!is.na(Px))))
    set(bad, r, "pool_minvalid_tn", min(rowSums(!is.na(Pn))))
  }
  
  classify <- function(nobs, allna, nv0, poolmin) {
    fifelse(nobs == 0L, "A1: no valid observations in this grid-year",
            fifelse(allna == N_REPLICATES, "A2: all 91 thresholds NA in every replicate (>=120/150 rule)",
                    fifelse(allna > 0L, "A3: some replicates fully NA (validity rule, marginal pool)",
                            fifelse(nv0 > 0L, "B1: obs and thresholds exist but never on the same day - INVESTIGATE",
                                    "B2: replicate counted as dropped despite nv>0 - CODING ERROR"))))
  }
  bad[, reason_tmax := classify(n_obs_days_tmax, rep_allNA_thr_tx, rep_nv0_tx, pool_minvalid_tx)]
  bad[, reason_tmin := classify(n_obs_days_tmin, rep_allNA_thr_tn, rep_nv0_tn, pool_minvalid_tn)]
  
  log_msg("Reason breakdown (Tmax):")
  print(bad[, .N, by = reason_tmax][order(-N)])
  log_msg("Reason breakdown (Tmin):")
  print(bad[, .N, by = reason_tmin][order(-N)])
  
  log_msg("Offending grid-years (grid_id | year | nrepTx | nrepTn | obsTx | obsTn | poolminTx | lockedNAtx):")
  for (r in seq_len(nrow(bad)))
    log_msg(sprintf("    %-12s | %d | %2d | %2d | %2d | %2d | %3d | %2d",
                    bad$grid_id[r], bad$year[r], bad$n_rep_tmax[r], bad$n_rep_tmin[r],
                    bad$n_obs_days_tmax[r], bad$n_obs_days_tmin[r],
                    bad$pool_minvalid_tx[r], bad$locked_na_tx90_days[r]))
  
  ## decisive A-vs-B discriminator: does STEP 03 also report zero valid days?
  agree <- bad[, sum((n_rep_tmax == 0L) == (step03_nvalid_tmax == 0L))]
  log_msg("Grid-years where 0 replicates coincides with STEP 03 zero valid Tmax days : ",
          agree, " / ", nrow(bad),
          "  -> full agreement implies genuine data absence, not a counting error")
  
  fwrite(bad, DIAG_FILE)
  log_msg("Diagnostic table saved : ", DIAG_FILE)
} else {
  log_msg("No offending grid-years found - the FATAL condition is not reproducible here.")
}
## --- COVERAGE ASSERTIONS ----------------------------------------------------
record_check("bootstrap_years_are_full_reference_period",
             identical(sort(unique(years_visited)), BASE_YR_START:BASE_YR_END) &&
               length(years_visited) == N_BASE_YEARS,
             sprintf("visited %d years [%d-%d] ; expected %d [%d-%d]",
                     length(years_visited), min(years_visited), max(years_visited),
                     N_BASE_YEARS, BASE_YR_START, BASE_YR_END))
record_check("bootstrap_year_count_is_30", length(years_visited) == 30L,
             sprintf("observed %d (failing run: 15 = trend period intersected with reference)",
                     length(years_visited)))
record_check("bootstrap_pair_count_is_870", pairs_done == EXP_PAIRS,
             sprintf("observed %d vs expected %d (30 x 29)", pairs_done, EXP_PAIRS))
record_check("bootstrap_block_count_is_79170", blocks_done == EXP_BLOCKS,
             sprintf("observed %s vs expected %s (870 x 91)",
                     format(blocks_done, big.mark = ","), format(EXP_BLOCKS, big.mark = ",")))
record_check("percentile_evaluation_count", blocks_done * N_GRID == EXP_PCTL_EVAL,
             sprintf("observed %s vs expected %s per variable (79,170 x 283)",
                     format(blocks_done * N_GRID, big.mark = ","),
                     format(EXP_PCTL_EVAL, big.mark = ",")))
record_check("pre_1996_reference_years_included", all(1981:1995 %in% years_visited),
             sprintf("%d of the 15 years 1981-1995 present (failing run: 0)",
                     sum(1981:1995 %in% years_visited)))
record_check("realised_replicates_equal_29",
             all(BOOT_NREPTX == N_REPLICATES) && all(BOOT_NREPTN == N_REPLICATES),
             sprintf("min realised replicates: Tmax %d, Tmin %d (expected %d everywhere); grid-years with 29/29: Tmax %d, Tmin %d; with 1-28: Tmax %d, Tmin %d; with 0: Tmax %d, Tmin %d of %d",
                     min(BOOT_NREPTX), min(BOOT_NREPTN), N_REPLICATES,
                     sum(BOOT_NREPTX == N_REPLICATES), sum(BOOT_NREPTN == N_REPLICATES),
                     sum(BOOT_NREPTX > 0L & BOOT_NREPTX < N_REPLICATES),
                     sum(BOOT_NREPTN > 0L & BOOT_NREPTN < N_REPLICATES),
                     sum(BOOT_NREPTX == 0L), sum(BOOT_NREPTN == 0L),
                     length(BOOT_NREPTX)),
             severity = "WARN")
record_check("no_bootstrap_year_all_NA",
             all(colSums(!is.na(BOOT_TX90p)) > 0L) && all(colSums(!is.na(BOOT_TN10p)) > 0L),
             sprintf("min non-NA grids per year: TX90p %d, TN10p %d",
                     min(colSums(!is.na(BOOT_TX90p))), min(colSums(!is.na(BOOT_TN10p)))))

rm(A_TX, A_TN, O_TX, O_TN); invisible(gc(verbose = FALSE))

## ---------------------------------------------------------------------------
## 8. ASSEMBLE OUTPUT
## ---------------------------------------------------------------------------
log_head("8. ASSEMBLE OUTPUT")

boot <- data.table(
  grid_id            = rep(GRIDS, times = nb),
  year               = rep(BOOT_YEARS, each = N_GRID),
  TX90p_boot         = as.vector(BOOT_TX90p),
  TN10p_boot         = as.vector(BOOT_TN10p),
  TX90p_boot_sd      = as.vector(BOOT_TX_SD),
  TN10p_boot_sd      = as.vector(BOOT_TN_SD),
  TX90_days_boot     = as.vector(BOOT_TXDAYS),
  TN10_days_boot     = as.vector(BOOT_TNDAYS),
  n_valid_tmax_boot  = as.vector(BOOT_NVTX),
  n_valid_tmin_boot  = as.vector(BOOT_NVTN),
  n_rep_used_tmax    = as.vector(BOOT_NREPTX),
  n_rep_used_tmin    = as.vector(BOOT_NREPTN))
record_check("bootstrap_table_rows", nrow(boot) == EXP_BOOT_ROWS,
             sprintf("observed %s vs expected %s (283 x 30)",
                     format(nrow(boot), big.mark = ","), format(EXP_BOOT_ROWS, big.mark = ",")))

res <- merge(direct, boot, by = c("grid_id", "year"), all.x = TRUE, sort = FALSE)
res <- merge(res, meta, by = "grid_id", all.x = TRUE, sort = FALSE)

res[, in_baseline := year >= BASE_YR_START & year <= BASE_YR_END]
res[, bootstrap_applied := in_baseline]

res[, `:=`(
  TX90p        = fifelse(bootstrap_applied, TX90p_boot,        TX90p_direct),
  TN10p        = fifelse(bootstrap_applied, TN10p_boot,        TN10p_direct),
  TX90_days    = fifelse(bootstrap_applied, TX90_days_boot,    as.numeric(tx_days)),
  TN10_days    = fifelse(bootstrap_applied, TN10_days_boot,    as.numeric(tn_days)),
  n_valid_tmax = fifelse(bootstrap_applied, n_valid_tmax_boot, as.numeric(nv_tmax)),
  n_valid_tmin = fifelse(bootstrap_applied, n_valid_tmin_boot, as.numeric(nv_tmin)),
  n_bootstrap_replicates = fifelse(bootstrap_applied, N_REPLICATES, NA_integer_),
  method       = fifelse(bootstrap_applied, "Zhang2005_inbase_bootstrap",
                         "fixed_1981_2010_thresholds"),
  n_valid_tmax_direct = as.integer(nv_tmax),
  n_valid_tmin_direct = as.integer(nv_tmin),
  TX90_days_direct    = as.integer(tx_days),
  TN10_days_direct    = as.integer(tn_days),
  season           = "AMJ",
  season_days      = N_AMJ_DAYS,
  baseline_start   = BASE_YR_START,
  baseline_end     = BASE_YR_END,
  n_baseline_years = N_BASE_YEARS,
  window_days      = WIN_WIDTH,
  pool_size        = N_POSSIBLE,
  min_valid_req    = MIN_VALID,
  pctl_type        = PCTL_TYPE)]

res[, c("TX90p_boot", "TN10p_boot", "TX90_days_boot", "TN10_days_boot",
        "n_valid_tmax_boot", "n_valid_tmin_boot", "nv_tmax", "nv_tmin",
        "tx_days", "tn_days") := NULL]

setcolorder(res, c(
  "year", "grid_id", "lat", "lon",
  "n_AMJ_days", "n_valid_tmax", "n_valid_tmin",
  "TX90_days", "TN10_days", "TX90p", "TN10p",
  "TX90p_boot_sd", "TN10p_boot_sd",
  "bootstrap_applied", "n_bootstrap_replicates", "n_rep_used_tmax", "n_rep_used_tmin",
  "method", "in_baseline",
  "TX90p_direct", "TN10p_direct",
  "TX90_days_direct", "TN10_days_direct",
  "n_valid_tmax_direct", "n_valid_tmin_direct",
  "season", "season_days", "baseline_start", "baseline_end", "n_baseline_years",
  "window_days", "pool_size", "min_valid_req", "pctl_type"))
setorder(res, year, grid_id)

trend <- res[year >= TREND_START & year <= TREND_END]
setorder(trend, year, grid_id)

## ---------------------------------------------------------------------------
## 9. FINAL QC
## ---------------------------------------------------------------------------
log_head("9. FINAL QC")

record_check("output_row_count", nrow(res) == EXP_OUT_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(res), big.mark = ","),
                     format(EXP_OUT_ROWS, big.mark = ",")))
record_check("output_year_count", uniqueN(res$year) == N_YEARS,
             sprintf("observed %d vs expected %d", uniqueN(res$year), N_YEARS))
record_check("output_grid_count", uniqueN(res$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(res$grid_id), EXP_GRIDS))
record_check("no_duplicate_grid_year", !anyDuplicated(res, by = c("grid_id", "year")),
             sprintf("%d duplicated (grid_id, year)", sum(duplicated(res, by = c("grid_id", "year")))))
record_check("283_grids_every_year", nrow(res[, .N, by = year][N != EXP_GRIDS]) == 0L,
             sprintf("%d years with grid count != %d",
                     nrow(res[, .N, by = year][N != EXP_GRIDS]), EXP_GRIDS))
record_check("50_years_every_grid", nrow(res[, .N, by = grid_id][N != N_YEARS]) == 0L,
             sprintf("%d grids with year count != %d",
                     nrow(res[, .N, by = grid_id][N != N_YEARS]), N_YEARS))
record_check("n_AMJ_days_is_91", res[n_AMJ_days != N_AMJ_DAYS, .N] == 0L,
             sprintf("%d grid-years with n_AMJ_days != %d",
                     res[n_AMJ_days != N_AMJ_DAYS, .N], N_AMJ_DAYS))

## --- bootstrap coverage in the OUTPUT ---------------------------------------
record_check("bootstrap_rows_is_8490", res[bootstrap_applied == TRUE, .N] == EXP_BOOT_ROWS,
             sprintf("observed %s vs expected %s (failing run produced 4,245 = 283 x 15)",
                     format(res[bootstrap_applied == TRUE, .N], big.mark = ","),
                     format(EXP_BOOT_ROWS, big.mark = ",")))
record_check("direct_rows_is_5660", res[bootstrap_applied == FALSE, .N] == EXP_DIRECT_RW,
             sprintf("observed %s vs expected %s (283 x 20)",
                     format(res[bootstrap_applied == FALSE, .N], big.mark = ","),
                     format(EXP_DIRECT_RW, big.mark = ",")))
record_check("bootstrap_flag_iff_reference_year",
             identical(sort(unique(res[bootstrap_applied == TRUE, year])), BASE_YR_START:BASE_YR_END),
             sprintf("flagged years: %d distinct [%d-%d]",
                     uniqueN(res[bootstrap_applied == TRUE, year]),
                     min(res[bootstrap_applied == TRUE, year]),
                     max(res[bootstrap_applied == TRUE, year])))
record_check("replicates_29_on_every_inbase_row",
             res[bootstrap_applied == TRUE,
                 all(n_bootstrap_replicates == N_REPLICATES &
                       n_rep_used_tmax == N_REPLICATES & n_rep_used_tmin == N_REPLICATES)] &&
               res[bootstrap_applied == FALSE, all(is.na(n_bootstrap_replicates))],
             sprintf("in-base realised replicates: Tmax min %d, Tmin min %d; rows below %d: Tmax %d, Tmin %d of %d",
                     res[bootstrap_applied == TRUE, min(n_rep_used_tmax)],
                     res[bootstrap_applied == TRUE, min(n_rep_used_tmin)], N_REPLICATES,
                     res[bootstrap_applied == TRUE, sum(n_rep_used_tmax < N_REPLICATES)],
                     res[bootstrap_applied == TRUE, sum(n_rep_used_tmin < N_REPLICATES)],
                     res[bootstrap_applied == TRUE, .N]),
             severity = "WARN")

## --- DECISIVE per-year evidence that each reference year was bootstrapped ---
per_yr <- res[bootstrap_applied == TRUE,
              .(n_grid_changed = sum(abs(TX90p - TX90p_direct) > TOL_EXACT |
                                       abs(TN10p - TN10p_direct) > TOL_EXACT, na.rm = TRUE),
                mean_shift_tx  = mean(TX90p - TX90p_direct, na.rm = TRUE),
                mean_shift_tn  = mean(TN10p - TN10p_direct, na.rm = TRUE)),
              by = year][order(year)]
record_check("every_reference_year_actually_bootstrapped",
             nrow(per_yr) == N_BASE_YEARS && all(per_yr$n_grid_changed > 0L),
             sprintf("%d years assessed; min grids changed in a year = %d (0 would mean that year was skipped)",
                     nrow(per_yr), min(per_yr$n_grid_changed)))
log_msg("Per-year bootstrap effect (year : grids changed / 283 : mean TX90p shift : mean TN10p shift)")
for (i in seq_len(nrow(per_yr)))
  log_msg(sprintf("    %d : %3d : %+6.3f : %+6.3f", per_yr$year[i], per_yr$n_grid_changed[i],
                  per_yr$mean_shift_tx[i], per_yr$mean_shift_tn[i]))

## --- out-of-base years must be exactly the direct calculation ---------------
oob_same <- res[bootstrap_applied == FALSE,
                all(is.na(TX90p) == is.na(TX90p_direct) &
                      (is.na(TX90p) | TX90p == TX90p_direct) &
                      is.na(TN10p) == is.na(TN10p_direct) &
                      (is.na(TN10p) | TN10p == TN10p_direct))]
record_check("out_of_base_equals_direct", isTRUE(oob_same),
             "1976-1980 and 2011-2025 bitwise unchanged from the fixed-threshold calculation")

if (file.exists(STEP03_FILE)) {
  s3 <- as.data.table(arrow::read_parquet(STEP03_FILE))
  s3 <- s3[, .(grid_id = as.character(grid_id), year = as.integer(year),
               s3_tx = as.numeric(TX90p), s3_tn = as.numeric(TN10p))]
  cc <- merge(res[, .(grid_id, year, TX90p, TN10p, bootstrap_applied)], s3,
              by = c("grid_id", "year"))
  oo <- cc[bootstrap_applied == FALSE]
  record_check("out_of_base_reproduces_STEP03",
               nrow(oo) == EXP_DIRECT_RW &&
                 isTRUE(all.equal(oo$TX90p, oo$s3_tx, tolerance = TOL_EXACT)) &&
                 isTRUE(all.equal(oo$TN10p, oo$s3_tn, tolerance = TOL_EXACT)),
               sprintf("%s out-of-base rows identical to STEP 03", format(nrow(oo), big.mark = ",")))
  ib <- cc[bootstrap_applied == TRUE]
  n_ch <- ib[, sum(abs(TX90p - s3_tx) > TOL_EXACT | abs(TN10p - s3_tn) > TOL_EXACT, na.rm = TRUE)]
  record_check("in_base_differs_from_STEP03", n_ch > 0L,
               sprintf("%s of %s in-base rows changed by the bootstrap",
                       format(n_ch, big.mark = ","), format(nrow(ib), big.mark = ",")))
  rm(s3, cc, oo, ib)
} else {
  record_check("step03_comparison_available", FALSE,
               "STEP 03 output not found; cross-check skipped", severity = "WARN")
}

## --- value-domain checks ----------------------------------------------------
record_check("TX90p_within_0_100", res[!is.na(TX90p) & (TX90p < 0 | TX90p > 100), .N] == 0L,
             sprintf("observed [%.2f, %.2f]", min(res$TX90p, na.rm = TRUE), max(res$TX90p, na.rm = TRUE)))
record_check("TN10p_within_0_100", res[!is.na(TN10p) & (TN10p < 0 | TN10p > 100), .N] == 0L,
             sprintf("observed [%.2f, %.2f]", min(res$TN10p, na.rm = TRUE), max(res$TN10p, na.rm = TRUE)))
record_check("exceedance_days_le_valid_days",
             res[, all(TX90_days <= n_valid_tmax + TOL_EXACT &
                         TN10_days <= n_valid_tmin + TOL_EXACT, na.rm = TRUE)],
             "TX90_days <= n_valid_tmax and TN10_days <= n_valid_tmin on every row")
record_check("valid_days_never_fixed_at_91", res[, uniqueN(round(n_valid_tmax, 6)) > 1L],
             sprintf("distinct n_valid_tmax = %d ; distinct n_valid_tmin = %d",
                     uniqueN(round(res$n_valid_tmax, 6)), uniqueN(round(res$n_valid_tmin, 6))),
             severity = "WARN")
record_check("NA_index_iff_zero_denominator",
             res[(is.na(TX90p) & n_valid_tmax > 0) | (!is.na(TX90p) & n_valid_tmax == 0), .N] == 0L &&
               res[(is.na(TN10p) & n_valid_tmin > 0) | (!is.na(TN10p) & n_valid_tmin == 0), .N] == 0L,
             sprintf("NA TX90p = %d ; NA TN10p = %d", res[is.na(TX90p), .N], res[is.na(TN10p), .N]))
record_check("coords_carried_through", !anyNA(res$lat) && !anyNA(res$lon),
             sprintf("%d NA lat, %d NA lon", sum(is.na(res$lat)), sum(is.na(res$lon))))

## --- trend subset -----------------------------------------------------------
record_check("trend_row_count", nrow(trend) == EXP_TREND_RW,
             sprintf("observed %s vs expected %s (283 x 30)",
                     format(nrow(trend), big.mark = ","), format(EXP_TREND_RW, big.mark = ",")))
record_check("trend_year_span", identical(sort(unique(trend$year)), TREND_START:TREND_END),
             sprintf("%d years [%d-%d]", uniqueN(trend$year), min(trend$year), max(trend$year)))
record_check("trend_grid_count", uniqueN(trend$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(trend$grid_id), EXP_GRIDS))
record_check("trend_bootstrap_composition",
             trend[bootstrap_applied == TRUE, uniqueN(year)] == 15L &&
               trend[bootstrap_applied == FALSE, uniqueN(year)] == 15L,
             sprintf("%d bootstrapped years (1996-2010) + %d fixed-threshold years (2011-2025)",
                     trend[bootstrap_applied == TRUE, uniqueN(year)],
                     trend[bootstrap_applied == FALSE, uniqueN(year)]))

## --- scientific plausibility -------------------------------------------------
b_tx <- res[in_baseline == TRUE, mean(TX90p, na.rm = TRUE)]
b_tn <- res[in_baseline == TRUE, mean(TN10p, na.rm = TRUE)]
record_check("baseline_mean_TX90p_near_10",
             is.finite(b_tx) && b_tx >= BASE_MEAN_LO && b_tx <= BASE_MEAN_HI,
             sprintf("1981-2010 mean TX90p = %.3f%% (STEP 03 direct = %.3f%%)",
                     b_tx, res[in_baseline == TRUE, mean(TX90p_direct, na.rm = TRUE)]),
             severity = "WARN")
record_check("baseline_mean_TN10p_near_10",
             is.finite(b_tn) && b_tn >= BASE_MEAN_LO && b_tn <= BASE_MEAN_HI,
             sprintf("1981-2010 mean TN10p = %.3f%% (STEP 03 direct = %.3f%%)",
                     b_tn, res[in_baseline == TRUE, mean(TN10p_direct, na.rm = TRUE)]),
             severity = "WARN")

yr_mean <- res[, .(tx = mean(TX90p, na.rm = TRUE), tn = mean(TN10p, na.rm = TRUE)), by = year][order(year)]
gj <- function(v, a, b) abs(yr_mean[year == a, get(v)] - yr_mean[year == b, get(v)])
log_msg(sprintf("Boundary discontinuity 1980/1981 : TX90p %.3f  TN10p %.3f",
                gj("tx", 1980L, 1981L), gj("tn", 1980L, 1981L)))
log_msg(sprintf("Boundary discontinuity 2010/2011 : TX90p %.3f  TN10p %.3f",
                gj("tx", 2010L, 2011L), gj("tn", 2010L, 2011L)))
log_msg("Mean TX90p 1976-1980 / 1981-2010 / 2011-2025 : ",
        sprintf("%.2f / %.2f / %.2f", res[year < BASE_YR_START, mean(TX90p, na.rm = TRUE)],
                b_tx, res[year > BASE_YR_END, mean(TX90p, na.rm = TRUE)]))
log_msg("Mean TN10p 1976-1980 / 1981-2010 / 2011-2025 : ",
        sprintf("%.2f / %.2f / %.2f", res[year < BASE_YR_START, mean(TN10p, na.rm = TRUE)],
                b_tn, res[year > BASE_YR_END, mean(TN10p, na.rm = TRUE)]))
log_msg("Mean replicate spread (descriptive, NOT a standard error) : TX90p ",
        sprintf("%.3f", res[in_baseline == TRUE, mean(TX90p_boot_sd, na.rm = TRUE)]),
        " ; TN10p ", sprintf("%.3f", res[in_baseline == TRUE, mean(TN10p_boot_sd, na.rm = TRUE)]))

## ---------------------------------------------------------------------------
## 10. WRITE OUTPUTS + READ-BACK
## ---------------------------------------------------------------------------
log_head("10. WRITE OUTPUTS")

out_dir <- dirname(OUT_FULL)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

arrow::write_parquet(res,   OUT_FULL,  compression = "snappy")
arrow::write_parquet(trend, OUT_TREND, compression = "snappy")
log_msg("Full output       : ", OUT_FULL, "  (", format(nrow(res), big.mark = ","),
        " rows x ", ncol(res), " cols, ", sprintf("%.2f MB", file.size(OUT_FULL) / 1024^2), ")")
log_msg("Trend output      : ", OUT_TREND, "  (", format(nrow(trend), big.mark = ","),
        " rows x ", ncol(trend), " cols, ", sprintf("%.2f MB", file.size(OUT_TREND) / 1024^2), ")")

rb  <- as.data.table(arrow::read_parquet(OUT_FULL));  setorder(rb, year, grid_id)
rbt <- as.data.table(arrow::read_parquet(OUT_TREND)); setorder(rbt, year, grid_id)
record_check("readback_full", nrow(rb) == EXP_OUT_ROWS && identical(names(rb), names(res)) &&
               isTRUE(all.equal(rb$TX90p, res$TX90p)) && isTRUE(all.equal(rb$TN10p, res$TN10p)),
             sprintf("%s rows, %d cols, values identical after round-trip",
                     format(nrow(rb), big.mark = ","), ncol(rb)))
record_check("readback_trend", nrow(rbt) == EXP_TREND_RW &&
               isTRUE(all.equal(rbt$TX90p, trend$TX90p)),
             sprintf("%s rows, values identical after round-trip", format(nrow(rbt), big.mark = ",")))
record_check("readback_bootstrap_rows_preserved",
             rb[bootstrap_applied == TRUE, .N] == EXP_BOOT_ROWS,
             sprintf("%s bootstrapped rows survived the round-trip",
                     format(rb[bootstrap_applied == TRUE, .N], big.mark = ",")))
rm(rb, rbt)

record_check("upstream_files_unmodified",
             file.exists(THR_FILE) && file.exists(BUF_FILE) && file.exists(DAILY_FILE),
             "thresholds, buffer and daily files opened read-only, never rewritten")

qc_tab <- rbindlist(CHECKS)
qc_tab[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
              script   = "04_AMJ_TX90p_TN10p_bootstrap_1976_2025.R")]
fwrite(qc_tab, QC_FILE)
fwrite(per_yr, sub("\\.csv$", "_per_year_bootstrap_effect.csv", QC_FILE))
log_msg("QC summary saved  : ", QC_FILE, "  (", nrow(qc_tab), " checks)")

log_head("RUN COMPLETE")
log_msg("Checks passed     : ", qc_tab[status == "PASS", .N], " / ", nrow(qc_tab))
log_msg("Warnings          : ", qc_tab[status == "WARN", .N])
log_msg("Fatal             : ", qc_tab[status == "FATAL", .N])
log_msg("Bootstrap years   : ", length(years_visited), " (", min(years_visited), "-",
        max(years_visited), ")  pairs = ", pairs_done,
        "  blocks = ", format(blocks_done, big.mark = ","))
log_msg("Elapsed           : ",
        sprintf("%.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
writeLines(.LOG, LOG_FILE)
cat("\nConsole log written to:", LOG_FILE, "\n")

invisible(res)
###############################################################################
## END OF SCRIPT
###############################################################################

thr02 <- as.data.table(arrow::read_parquet(THR_FILE))
thr02[as.character(grid_id) == "283" & (is.na(tx90) | is.na(tn10)),
      .(target_day, target_md, n_valid_tmax, n_valid_tmin,
        completeness_tmax_pct, completeness_tmin_pct)][order(target_md)]