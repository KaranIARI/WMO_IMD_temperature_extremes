###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 05 : Validation, flagging, and LOCKED analysis-ready dataset
##  283 grids | 1976-2025 | AMJ 01 Apr - 30 Jun | reference 1981-2010
##
##  SCOPE
##    * Reads STEP 04 outputs READ-ONLY. Nothing upstream is modified.
##    * Recalculates NOTHING: TX90p / TN10p are carried through unchanged.
##    * No change to the 80% / >=120-of-150 rule, the 150-value pool, NA
##      handling, percentile type 8, or the Zhang (2005) bootstrap.
##    * All 14,150 grid-year rows are PRESERVED. Problem rows are FLAGGED,
##      never deleted. Grid 283 is NOT removed.
##    * The two STEP 04 WARNs are treated as documented data limitations
##      unless QC finds a NEW error.
##    * NO Mann-Kendall / Sen slope here. Trend estimation is STEP 06.
##
##  INPUTS (read-only)
##    1. AMJ_TX90p_TN10p_bootstrap_1976_2025.parquet        [STEP 04, 14,150]
##    2. AMJ_TX90p_TN10p_bootstrap_trend_1996_2025.parquet  [STEP 04,  8,490]
##    3. AMJ_TX90_TN10_thresholds_1981_2010.parquet         [STEP 02, locked]
##    4. AMJ_TX90p_TN10p_1976_2025.parquet                  [STEP 03, cross-check]
##
##  OUTPUTS
##    A. AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet        14,150 rows
##    B. AMJ_TX90p_TN10p_ANALYSIS_READY_trend_1996_2025.parquet   8,490 rows
##    C. AMJ_STEP05_QC.csv               - every check, PASS/WARN/FATAL
##    D. AMJ_STEP05_grid_summary.csv     - 283 rows, per-grid usability
##    E. AMJ_STEP05_data_limitations.csv - every flagged grid-year, with reason
##    F. AMJ_STEP05_LOG.txt              - full console log + sessionInfo
##
##  FLAGS ADDED (advisory only - nothing is filtered or altered by them)
##    flag_no_valid_tmax / flag_no_valid_tmin  denominator == 0, index is NA
##    flag_reduced_replicates                  in-base row with < 29 realised
##    flag_low_completeness_tmax / _tmin       valid days < ANNUAL_MIN_VALID
##    flag_locked_na_thresholds                grid has NA days in STEP 02
##    data_quality_tmax / data_quality_tmin    NO_DATA > LOW_COMPLETENESS >
##                                             REDUCED_REPLICATES > OK
##    usable_tmax / usable_tmin                index non-NA (use for STEP 06)
##
##  ANNUAL_MIN_VALID is a NEW, ADVISORY reporting threshold introduced in this
##  step alone. It applies 80% to the 91 AMJ days (ceiling(0.80 * 91) = 73) by
##  analogy with the locked baseline rule. It filters nothing, changes no index,
##  and can be set to 0 to disable flagging entirely.
##
##  Schema policy: required STEP 04 columns must be present by exact name.
##  If any is missing the script STOPS with an explicit schema error and never
##  guesses a substitute.
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
      stop("Required package '", p, "' is not available. Install it and re-run.", call. = FALSE)
  }
}
suppressPackageStartupMessages({ library(data.table); library(arrow) })
setDTthreads(0L)

## ---------------------------------------------------------------------------
## 1. PATHS AND LOCKED CONSTANTS
## ---------------------------------------------------------------------------
IN_FULL   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_bootstrap_1976_2025.parquet"
IN_TREND  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_bootstrap_trend_1996_2025.parquet"
THR_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet"
S03_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_1976_2025.parquet"

OUT_FULL  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet"
OUT_TREND <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_ANALYSIS_READY_trend_1996_2025.parquet"
QC_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP05_QC.csv"
GRID_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP05_grid_summary.csv"
LIM_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP05_data_limitations.csv"
LOG_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP05_LOG.txt"

YR_START <- 1976L; YR_END <- 2025L; N_YEARS <- 50L
BASE_START <- 1981L; BASE_END <- 2010L; N_BASE_YEARS <- 30L
TREND_START <- 1996L; TREND_END <- 2025L; N_TREND_YEARS <- 30L
N_AMJ_DAYS <- 91L
EXP_GRIDS <- 283L
EXP_FULL_ROWS  <- 14150L      # 283 x 50
EXP_TREND_ROWS <- 8490L       # 283 x 30
EXP_BOOT_ROWS  <- 8490L       # 283 x 30 in-base
EXP_DIRECT_ROWS<- 5660L       # 283 x 20 out-of-base
N_REPLICATES   <- 29L
POOL_SIZE      <- 150L
MIN_VALID_POOL <- 120L
PCTL_TYPE      <- 8L

## Advisory annual-completeness flag (this step only; filters nothing).
ANNUAL_MIN_FRAC  <- 0.80
ANNUAL_MIN_VALID <- as.integer(ceiling(ANNUAL_MIN_FRAC * N_AMJ_DAYS))   # 73
TOL <- 1e-9

REQUIRED_COLS <- c(
  "year", "grid_id", "lat", "lon",
  "n_AMJ_days", "n_valid_tmax", "n_valid_tmin",
  "TX90_days", "TN10_days", "TX90p", "TN10p",
  "TX90p_boot_sd", "TN10p_boot_sd",
  "bootstrap_applied", "n_bootstrap_replicates",
  "n_rep_used_tmax", "n_rep_used_tmin",
  "method", "in_baseline",
  "TX90p_direct", "TN10p_direct",
  "TX90_days_direct", "TN10_days_direct",
  "n_valid_tmax_direct", "n_valid_tmin_direct",
  "season", "season_days", "baseline_start", "baseline_end",
  "n_baseline_years", "window_days", "pool_size", "min_valid_req", "pctl_type")

## ---------------------------------------------------------------------------
## 2. HELPERS
## ---------------------------------------------------------------------------
.LOG <- character(0)
log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "] ",
                 paste0(..., collapse = ""))
  cat(line, "\n", sep = ""); .LOG <<- c(.LOG, line); invisible(line)
}
log_head <- function(txt) {
  bar <- strrep("-", 76)
  cat("\n", bar, "\n", txt, "\n", bar, "\n", sep = "")
  .LOG <<- c(.LOG, "", bar, txt, bar); invisible(NULL)
}
CHECKS <- list()
record_check <- function(name, passed, detail, severity = c("FATAL", "WARN")) {
  severity <- match.arg(severity)
  status <- if (isTRUE(passed)) "PASS" else severity
  CHECKS[[length(CHECKS) + 1L]] <<- data.table(check = name, status = status, detail = detail)
  log_msg(sprintf("  %-5s | %-42s | %s", status, name, detail))
  if (!isTRUE(passed) && severity == "FATAL")
    stop("FATAL QC failure: ", name, " -- ", detail, call. = FALSE)
  invisible(isTRUE(passed))
}
assert_schema <- function(dt, required, label) {
  miss <- setdiff(required, names(dt))
  if (length(miss))
    stop("SCHEMA ERROR in ", label, ".\n  Missing required column(s): ",
         paste(miss, collapse = ", "),
         "\n  Columns present: ", paste(names(dt), collapse = ", "),
         "\n  STEP 05 will not guess substitutes. Re-run STEP 04 or correct the file.",
         call. = FALSE)
  invisible(TRUE)
}

t_start <- Sys.time()
log_head("STEP 05 | Validation, flagging, and locked analysis-ready dataset")
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("Full period      : ", YR_START, "-", YR_END, " (", N_YEARS, " years)")
log_msg("Reference period : ", BASE_START, "-", BASE_END, " (", N_BASE_YEARS, " years)")
log_msg("Trend period     : ", TREND_START, "-", TREND_END, " (", N_TREND_YEARS, " years)")
log_msg("Policy           : read-only inputs, no recalculation, no row deletion, flags only")
log_msg("Advisory flag    : annual completeness < ", ANNUAL_MIN_VALID, "/", N_AMJ_DAYS,
        " days (", ANNUAL_MIN_FRAC * 100, "%) - reporting only, filters nothing")

## ---------------------------------------------------------------------------
## 3. READ STEP 04 OUTPUTS (READ-ONLY) AND ENFORCE SCHEMA
## ---------------------------------------------------------------------------
log_head("3. READ STEP 04 OUTPUTS (READ-ONLY)")

for (f in c(IN_FULL, IN_TREND))
  if (!file.exists(f)) stop("Required STEP 04 input not found: ", f, call. = FALSE)

full  <- as.data.table(arrow::read_parquet(IN_FULL))
trend <- as.data.table(arrow::read_parquet(IN_TREND))
log_msg("Full  file       : ", format(nrow(full),  big.mark = ","), " rows x ", ncol(full),  " cols")
log_msg("Trend file       : ", format(nrow(trend), big.mark = ","), " rows x ", ncol(trend), " cols")

assert_schema(full,  REQUIRED_COLS, basename(IN_FULL))
assert_schema(trend, REQUIRED_COLS, basename(IN_TREND))
record_check("schema_full_file", TRUE,
             sprintf("all %d required columns present", length(REQUIRED_COLS)))
record_check("schema_trend_file", TRUE,
             sprintf("all %d required columns present", length(REQUIRED_COLS)))

full[,  `:=`(grid_id = as.character(grid_id), year = as.integer(year))]
trend[, `:=`(grid_id = as.character(grid_id), year = as.integer(year))]
setorder(full, year, grid_id); setorder(trend, year, grid_id)

## metadata columns must be constant and match the locked settings
record_check("locked_settings_constant",
             full[, uniqueN(pool_size) == 1L && uniqueN(min_valid_req) == 1L &&
                    uniqueN(pctl_type) == 1L && uniqueN(window_days) == 1L &&
                    uniqueN(season_days) == 1L &&
                    pool_size[1] == POOL_SIZE && min_valid_req[1] == MIN_VALID_POOL &&
                    pctl_type[1] == PCTL_TYPE && season_days[1] == N_AMJ_DAYS],
             sprintf("pool=%s min_valid=%s pctl_type=%s window=%s season_days=%s",
                     full$pool_size[1], full$min_valid_req[1], full$pctl_type[1],
                     full$window_days[1], full$season_days[1]))
record_check("locked_baseline_constant",
             full[, all(baseline_start == BASE_START) && all(baseline_end == BASE_END) &&
                    all(n_baseline_years == N_BASE_YEARS)],
             sprintf("baseline %d-%d, %d years", full$baseline_start[1],
                     full$baseline_end[1], full$n_baseline_years[1]))
record_check("season_is_AMJ", full[, all(season == "AMJ")],
             sprintf("distinct season labels = %d", uniqueN(full$season)))

## ---------------------------------------------------------------------------
## 4. STRUCTURAL QC
## ---------------------------------------------------------------------------
log_head("4. STRUCTURAL QC")

record_check("full_row_count", nrow(full) == EXP_FULL_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(full), big.mark = ","),
                     format(EXP_FULL_ROWS, big.mark = ",")))
record_check("full_grid_count", uniqueN(full$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(full$grid_id), EXP_GRIDS))
record_check("full_year_coverage", identical(sort(unique(full$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(full$year), min(full$year), max(full$year)))
record_check("no_duplicate_grid_year", !anyDuplicated(full, by = c("grid_id", "year")),
             sprintf("%d duplicated (grid_id, year)", sum(duplicated(full, by = c("grid_id", "year")))))
record_check("283_grids_every_year", nrow(full[, .N, by = year][N != EXP_GRIDS]) == 0L,
             sprintf("%d years with grid count != %d",
                     nrow(full[, .N, by = year][N != EXP_GRIDS]), EXP_GRIDS))
record_check("50_years_every_grid", nrow(full[, .N, by = grid_id][N != N_YEARS]) == 0L,
             sprintf("%d grids with year count != %d",
                     nrow(full[, .N, by = grid_id][N != N_YEARS]), N_YEARS))
record_check("complete_grid_x_year_cross",
             nrow(full) == uniqueN(full$grid_id) * uniqueN(full$year),
             sprintf("%d rows vs %d x %d", nrow(full), uniqueN(full$grid_id), uniqueN(full$year)))
record_check("coords_present_and_unique",
             !anyNA(full$lat) && !anyNA(full$lon) &&
               nrow(unique(full[, .(grid_id, lat, lon)])) == EXP_GRIDS,
             sprintf("%d unique grid/lat/lon triplets, %d NA coords",
                     nrow(unique(full[, .(grid_id, lat, lon)])),
                     sum(is.na(full$lat)) + sum(is.na(full$lon))))
record_check("n_AMJ_days_is_91", full[n_AMJ_days != N_AMJ_DAYS, .N] == 0L,
             sprintf("%d grid-years with n_AMJ_days != %d",
                     full[n_AMJ_days != N_AMJ_DAYS, .N], N_AMJ_DAYS))

## ---------------------------------------------------------------------------
## 5. BOOTSTRAP-FLAG AND METHOD QC
## ---------------------------------------------------------------------------
log_head("5. BOOTSTRAP FLAG AND METHOD QC")

record_check("bootstrap_flag_iff_reference_year",
             identical(sort(unique(full[bootstrap_applied == TRUE, year])), BASE_START:BASE_END) &&
               full[bootstrap_applied != in_baseline, .N] == 0L,
             sprintf("%d flagged years [%d-%d]; flag == in_baseline on all rows",
                     uniqueN(full[bootstrap_applied == TRUE, year]),
                     min(full[bootstrap_applied == TRUE, year]),
                     max(full[bootstrap_applied == TRUE, year])))
record_check("bootstrap_row_count", full[bootstrap_applied == TRUE, .N] == EXP_BOOT_ROWS,
             sprintf("observed %s vs expected %s",
                     format(full[bootstrap_applied == TRUE, .N], big.mark = ","),
                     format(EXP_BOOT_ROWS, big.mark = ",")))
record_check("direct_row_count", full[bootstrap_applied == FALSE, .N] == EXP_DIRECT_ROWS,
             sprintf("observed %s vs expected %s",
                     format(full[bootstrap_applied == FALSE, .N], big.mark = ","),
                     format(EXP_DIRECT_ROWS, big.mark = ",")))
record_check("method_label_matches_flag",
             full[bootstrap_applied == TRUE,  all(method == "Zhang2005_inbase_bootstrap")] &&
               full[bootstrap_applied == FALSE, all(method == "fixed_1981_2010_thresholds")],
             sprintf("distinct method labels = %d", uniqueN(full$method)))
record_check("design_replicates_29_in_base",
             full[bootstrap_applied == TRUE, all(n_bootstrap_replicates == N_REPLICATES)] &&
               full[bootstrap_applied == FALSE, all(is.na(n_bootstrap_replicates))],
             sprintf("in-base design replicate value(s): %s",
                     paste(unique(full[bootstrap_applied == TRUE, n_bootstrap_replicates]), collapse = ",")))
record_check("realised_replicates_le_design",
             full[bootstrap_applied == TRUE,
                  all(n_rep_used_tmax <= N_REPLICATES & n_rep_used_tmin <= N_REPLICATES &
                        n_rep_used_tmax >= 0L & n_rep_used_tmin >= 0L)],
             sprintf("realised range Tmax %d-%d, Tmin %d-%d",
                     full[bootstrap_applied == TRUE, min(n_rep_used_tmax)],
                     full[bootstrap_applied == TRUE, max(n_rep_used_tmax)],
                     full[bootstrap_applied == TRUE, min(n_rep_used_tmin)],
                     full[bootstrap_applied == TRUE, max(n_rep_used_tmin)]))
n_red <- full[bootstrap_applied == TRUE &
                (n_rep_used_tmax < N_REPLICATES | n_rep_used_tmin < N_REPLICATES), .N]
record_check("reduced_replicate_rows_documented", TRUE,
             sprintf("%d in-base grid-years below %d realised replicates (STEP 04 WARN, carried as a flag)",
                     n_red, N_REPLICATES))
record_check("boot_sd_only_in_base",
             full[bootstrap_applied == FALSE, all(is.na(TX90p_boot_sd) & is.na(TN10p_boot_sd))],
             sprintf("%d out-of-base rows carry a bootstrap SD",
                     full[bootstrap_applied == FALSE & (!is.na(TX90p_boot_sd) | !is.na(TN10p_boot_sd)), .N]),
             severity = "WARN")

## ---------------------------------------------------------------------------
## 6. DENOMINATOR, NA-CONSISTENCY AND VALUE-RANGE QC
## ---------------------------------------------------------------------------
log_head("6. DENOMINATOR, NA CONSISTENCY AND VALUE RANGES")

record_check("valid_days_within_bounds",
             full[, all(n_valid_tmax >= 0 & n_valid_tmax <= n_AMJ_days + TOL &
                          n_valid_tmin >= 0 & n_valid_tmin <= n_AMJ_days + TOL)],
             sprintf("Tmax %.2f-%.2f ; Tmin %.2f-%.2f (cap %d)",
                     min(full$n_valid_tmax), max(full$n_valid_tmax),
                     min(full$n_valid_tmin), max(full$n_valid_tmin), N_AMJ_DAYS))
record_check("exceedance_days_le_valid_days",
             full[, all(TX90_days <= n_valid_tmax + TOL & TN10_days <= n_valid_tmin + TOL &
                          TX90_days >= -TOL & TN10_days >= -TOL)],
             sprintf("%d rows with TX90_days > n_valid_tmax ; %d with TN10_days > n_valid_tmin",
                     full[TX90_days > n_valid_tmax + TOL, .N],
                     full[TN10_days > n_valid_tmin + TOL, .N]))
record_check("denominator_not_fixed_at_91",
             full[, uniqueN(round(n_valid_tmax, 6)) > 1L],
             sprintf("distinct n_valid_tmax = %d ; distinct n_valid_tmin = %d",
                     uniqueN(round(full$n_valid_tmax, 6)), uniqueN(round(full$n_valid_tmin, 6))))
record_check("NA_index_iff_zero_denominator",
             full[(is.na(TX90p) & n_valid_tmax > TOL) | (!is.na(TX90p) & n_valid_tmax <= TOL), .N] == 0L &&
               full[(is.na(TN10p) & n_valid_tmin > TOL) | (!is.na(TN10p) & n_valid_tmin <= TOL), .N] == 0L,
             sprintf("NA TX90p = %d ; NA TN10p = %d ; zero-denominator rows Tmax = %d, Tmin = %d",
                     full[is.na(TX90p), .N], full[is.na(TN10p), .N],
                     full[n_valid_tmax <= TOL, .N], full[n_valid_tmin <= TOL, .N]))
record_check("TX90p_within_0_100", full[!is.na(TX90p) & (TX90p < -TOL | TX90p > 100 + TOL), .N] == 0L,
             sprintf("observed [%.3f, %.3f]", min(full$TX90p, na.rm = TRUE), max(full$TX90p, na.rm = TRUE)))
record_check("TN10p_within_0_100", full[!is.na(TN10p) & (TN10p < -TOL | TN10p > 100 + TOL), .N] == 0L,
             sprintf("observed [%.3f, %.3f]", min(full$TN10p, na.rm = TRUE), max(full$TN10p, na.rm = TRUE)))
record_check("boot_sd_non_negative",
             full[!is.na(TX90p_boot_sd) & TX90p_boot_sd < -TOL, .N] == 0L &&
               full[!is.na(TN10p_boot_sd) & TN10p_boot_sd < -TOL, .N] == 0L,
             sprintf("max SD: TX90p %.3f ; TN10p %.3f",
                     max(full$TX90p_boot_sd, na.rm = TRUE), max(full$TN10p_boot_sd, na.rm = TRUE)))
record_check("out_of_base_equals_direct_columns",
             full[bootstrap_applied == FALSE,
                  all(is.na(TX90p) == is.na(TX90p_direct) & (is.na(TX90p) | TX90p == TX90p_direct) &
                        is.na(TN10p) == is.na(TN10p_direct) & (is.na(TN10p) | TN10p == TN10p_direct))],
             "1976-1980 and 2011-2025 identical to the fixed-threshold columns")
record_check("in_base_differs_from_direct",
             full[bootstrap_applied == TRUE,
                  sum(abs(TX90p - TX90p_direct) > TOL | abs(TN10p - TN10p_direct) > TOL,
                      na.rm = TRUE)] > 0L,
             sprintf("%s in-base rows changed by the bootstrap",
                     format(full[bootstrap_applied == TRUE,
                                 sum(abs(TX90p - TX90p_direct) > TOL |
                                       abs(TN10p - TN10p_direct) > TOL, na.rm = TRUE)], big.mark = ",")))

## ---------------------------------------------------------------------------
## 7. TREND-FILE CONSISTENCY WITH THE FULL FILE
## ---------------------------------------------------------------------------
log_head("7. TREND FILE CONSISTENCY")

record_check("trend_row_count", nrow(trend) == EXP_TREND_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(trend), big.mark = ","),
                     format(EXP_TREND_ROWS, big.mark = ",")))
record_check("trend_year_span", identical(sort(unique(trend$year)), TREND_START:TREND_END),
             sprintf("%d years [%d-%d]", uniqueN(trend$year), min(trend$year), max(trend$year)))
record_check("trend_grid_count", uniqueN(trend$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(trend$grid_id), EXP_GRIDS))
record_check("trend_no_duplicate_grid_year", !anyDuplicated(trend, by = c("grid_id", "year")),
             sprintf("%d duplicates", sum(duplicated(trend, by = c("grid_id", "year")))))

sub_full <- full[year >= TREND_START & year <= TREND_END]
cmp <- merge(sub_full[, .(grid_id, year, f_tx = TX90p, f_tn = TN10p,
                          f_nvx = n_valid_tmax, f_nvn = n_valid_tmin)],
             trend[, .(grid_id, year, t_tx = TX90p, t_tn = TN10p,
                       t_nvx = n_valid_tmax, t_nvn = n_valid_tmin)],
             by = c("grid_id", "year"))
eq <- function(a, b) all(is.na(a) == is.na(b) & (is.na(a) | a == b))
record_check("trend_is_exact_subset_of_full",
             nrow(cmp) == EXP_TREND_ROWS && eq(cmp$f_tx, cmp$t_tx) && eq(cmp$f_tn, cmp$t_tn) &&
               eq(cmp$f_nvx, cmp$t_nvx) && eq(cmp$f_nvn, cmp$t_nvn),
             sprintf("%s rows matched; indices and denominators bitwise identical",
                     format(nrow(cmp), big.mark = ",")))
record_check("trend_bootstrap_composition",
             trend[bootstrap_applied == TRUE,  uniqueN(year)] == 15L &&
               trend[bootstrap_applied == FALSE, uniqueN(year)] == 15L,
             sprintf("%d bootstrapped years (1996-2010) + %d fixed-threshold years (2011-2025)",
                     trend[bootstrap_applied == TRUE,  uniqueN(year)],
                     trend[bootstrap_applied == FALSE, uniqueN(year)]))
rm(cmp, sub_full)

## ---------------------------------------------------------------------------
## 8. CROSS-CHECK AGAINST STEP 02 AND STEP 03 (read-only, optional)
## ---------------------------------------------------------------------------
log_head("8. CROSS-CHECK AGAINST STEP 02 AND STEP 03")

locked_na <- NULL
if (file.exists(THR_FILE)) {
  t02 <- as.data.table(arrow::read_parquet(THR_FILE))
  if (all(c("grid_id", "tx90", "tn10") %in% names(t02))) {
    t02[, grid_id := as.character(grid_id)]
    locked_na <- t02[, .(locked_na_tx90_days = sum(is.na(tx90)),
                         locked_na_tn10_days = sum(is.na(tn10))), by = grid_id]
    tot_tx <- sum(locked_na$locked_na_tx90_days); tot_tn <- sum(locked_na$locked_na_tn10_days)
    aff <- locked_na[locked_na_tx90_days > 0L | locked_na_tn10_days > 0L][order(-locked_na_tx90_days)]
    record_check("step02_threshold_rows", nrow(t02) == EXP_GRIDS * N_AMJ_DAYS,
                 sprintf("%s rows (283 x 91)", format(nrow(t02), big.mark = ",")))
    log_msg("Locked NA thresholds : TX90 ", tot_tx, " ; TN10 ", tot_tn,
            " ; grids affected = ", nrow(aff))
    for (i in seq_len(nrow(aff)))
      log_msg(sprintf("    grid %-6s : %2d NA TX90 days, %2d NA TN10 days",
                      aff$grid_id[i], aff$locked_na_tx90_days[i], aff$locked_na_tn10_days[i]))
    record_check("locked_na_confined_to_known_grids", nrow(aff) <= 5L,
                 sprintf("%d grids carry NA thresholds: %s", nrow(aff),
                         paste(aff$grid_id, collapse = ", ")), severity = "WARN")
  } else {
    record_check("step02_schema_readable", FALSE,
                 "threshold file lacks grid_id/tx90/tn10; cross-check skipped", severity = "WARN")
  }
  rm(t02)
} else {
  record_check("step02_file_available", FALSE,
               "STEP 02 threshold file not found; cross-check skipped", severity = "WARN")
}

if (file.exists(S03_FILE)) {
  s03 <- as.data.table(arrow::read_parquet(S03_FILE))
  if (all(c("grid_id", "year", "TX90p", "TN10p") %in% names(s03))) {
    s03 <- s03[, .(grid_id = as.character(grid_id), year = as.integer(year),
                   s3_tx = as.numeric(TX90p), s3_tn = as.numeric(TN10p))]
    cc <- merge(full[, .(grid_id, year, TX90p, TN10p, bootstrap_applied)], s03,
                by = c("grid_id", "year"))
    oo <- cc[bootstrap_applied == FALSE]
    record_check("out_of_base_reproduces_STEP03",
                 nrow(oo) == EXP_DIRECT_ROWS &&
                   isTRUE(all.equal(oo$TX90p, oo$s3_tx, tolerance = TOL)) &&
                   isTRUE(all.equal(oo$TN10p, oo$s3_tn, tolerance = TOL)),
                 sprintf("%s out-of-base rows identical to STEP 03", format(nrow(oo), big.mark = ",")))
    rm(s03, cc, oo)
  } else {
    record_check("step03_schema_readable", FALSE,
                 "STEP 03 file lacks expected columns; cross-check skipped", severity = "WARN")
  }
} else {
  record_check("step03_file_available", FALSE,
               "STEP 03 output not found; cross-check skipped", severity = "WARN")
}

## ---------------------------------------------------------------------------
## 9. BUILD THE ANALYSIS-READY TABLE  (flags only; no row is removed)
## ---------------------------------------------------------------------------
log_head("9. BUILD ANALYSIS-READY TABLE (FLAGS ONLY, NO DELETION)")

ar <- copy(full)

if (!is.null(locked_na)) {
  ar <- merge(ar, locked_na, by = "grid_id", all.x = TRUE, sort = FALSE)
  ar[is.na(locked_na_tx90_days), locked_na_tx90_days := 0L]
  ar[is.na(locked_na_tn10_days), locked_na_tn10_days := 0L]
} else {
  ar[, `:=`(locked_na_tx90_days = NA_integer_, locked_na_tn10_days = NA_integer_)]
}

ar[, `:=`(
  completeness_tmax_pct = round(100 * n_valid_tmax / N_AMJ_DAYS, 2),
  completeness_tmin_pct = round(100 * n_valid_tmin / N_AMJ_DAYS, 2),
  usable_tmax = !is.na(TX90p),
  usable_tmin = !is.na(TN10p),
  flag_no_valid_tmax = n_valid_tmax <= TOL,
  flag_no_valid_tmin = n_valid_tmin <= TOL,
  flag_reduced_replicates = bootstrap_applied == TRUE &
    (n_rep_used_tmax < N_REPLICATES | n_rep_used_tmin < N_REPLICATES),
  flag_low_completeness_tmax = n_valid_tmax > TOL & n_valid_tmax < ANNUAL_MIN_VALID,
  flag_low_completeness_tmin = n_valid_tmin > TOL & n_valid_tmin < ANNUAL_MIN_VALID,
  flag_locked_na_thresholds = !is.na(locked_na_tx90_days) &
    (locked_na_tx90_days > 0L | locked_na_tn10_days > 0L))]

ar[, data_quality_tmax := fifelse(flag_no_valid_tmax, "NO_DATA",
                                  fifelse(flag_low_completeness_tmax, "LOW_COMPLETENESS",
                                          fifelse(flag_reduced_replicates, "REDUCED_REPLICATES", "OK")))]
ar[, data_quality_tmin := fifelse(flag_no_valid_tmin, "NO_DATA",
                                  fifelse(flag_low_completeness_tmin, "LOW_COMPLETENESS",
                                          fifelse(flag_reduced_replicates, "REDUCED_REPLICATES", "OK")))]
ar[, any_flag := flag_no_valid_tmax | flag_no_valid_tmin | flag_reduced_replicates |
     flag_low_completeness_tmax | flag_low_completeness_tmin]

ar[, `:=`(in_trend_period = year >= TREND_START & year <= TREND_END,
          step05_annual_min_valid = ANNUAL_MIN_VALID,
          step05_stamp = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
          step05_status = "ANALYSIS_READY_LOCKED")]

setcolorder(ar, c(
  "year", "grid_id", "lat", "lon",
  "TX90p", "TN10p", "usable_tmax", "usable_tmin",
  "n_AMJ_days", "n_valid_tmax", "n_valid_tmin",
  "completeness_tmax_pct", "completeness_tmin_pct",
  "TX90_days", "TN10_days",
  "data_quality_tmax", "data_quality_tmin", "any_flag",
  "flag_no_valid_tmax", "flag_no_valid_tmin", "flag_reduced_replicates",
  "flag_low_completeness_tmax", "flag_low_completeness_tmin", "flag_locked_na_thresholds",
  "locked_na_tx90_days", "locked_na_tn10_days",
  "bootstrap_applied", "in_baseline", "in_trend_period",
  "n_bootstrap_replicates", "n_rep_used_tmax", "n_rep_used_tmin",
  "TX90p_boot_sd", "TN10p_boot_sd", "method",
  "TX90p_direct", "TN10p_direct", "TX90_days_direct", "TN10_days_direct",
  "n_valid_tmax_direct", "n_valid_tmin_direct",
  "season", "season_days", "baseline_start", "baseline_end", "n_baseline_years",
  "window_days", "pool_size", "min_valid_req", "pctl_type",
  "step05_annual_min_valid", "step05_stamp", "step05_status"))
setorder(ar, year, grid_id)

ar_trend <- ar[in_trend_period == TRUE]
setorder(ar_trend, year, grid_id)

record_check("no_rows_lost_in_flagging", nrow(ar) == EXP_FULL_ROWS,
             sprintf("%s rows preserved (deletion is not performed at this step)",
                     format(nrow(ar), big.mark = ",")))
record_check("indices_unchanged_by_step05",
             eq(ar$TX90p, full$TX90p) && eq(ar$TN10p, full$TN10p) &&
               eq(ar$n_valid_tmax, full$n_valid_tmax) && eq(ar$n_valid_tmin, full$n_valid_tmin),
             "TX90p / TN10p / denominators bitwise identical to STEP 04")
record_check("grid_283_retained", ar[grid_id == "283", .N] == N_YEARS,
             sprintf("grid 283 present with %d of %d years (not removed)",
                     ar[grid_id == "283", .N], N_YEARS))

## ---------------------------------------------------------------------------
## 10. FLAG SUMMARY, GRID 283 REPORT, PER-GRID SUMMARY
## ---------------------------------------------------------------------------
log_head("10. FLAG SUMMARY AND DATA LIMITATIONS")

log_msg("Flagged grid-years (of ", format(nrow(ar), big.mark = ","), "):")
log_msg("    no valid Tmax days        : ", ar[flag_no_valid_tmax == TRUE, .N])
log_msg("    no valid Tmin days        : ", ar[flag_no_valid_tmin == TRUE, .N])
log_msg("    reduced replicates (<29)  : ", ar[flag_reduced_replicates == TRUE, .N])
log_msg("    low completeness Tmax     : ", ar[flag_low_completeness_tmax == TRUE, .N])
log_msg("    low completeness Tmin     : ", ar[flag_low_completeness_tmin == TRUE, .N])
log_msg("    any flag                  : ", ar[any_flag == TRUE, .N])
log_msg("    usable TX90p / TN10p      : ", ar[usable_tmax == TRUE, .N], " / ", ar[usable_tmin == TRUE, .N])
log_msg("Trend period 1996-2025 (of ", format(nrow(ar_trend), big.mark = ","), "):")
log_msg("    any flag                  : ", ar_trend[any_flag == TRUE, .N])
log_msg("    usable TX90p / TN10p      : ", ar_trend[usable_tmax == TRUE, .N], " / ",
        ar_trend[usable_tmin == TRUE, .N])

log_msg("Data-quality class counts (Tmax):")
print(ar[, .N, by = data_quality_tmax][order(-N)])
log_msg("Data-quality class counts (Tmin):")
print(ar[, .N, by = data_quality_tmin][order(-N)])

## --- grid 283 dedicated report ----------------------------------------------
g283 <- ar[grid_id == "283"]
if (nrow(g283) > 0L) {
  log_msg("GRID 283 report ------------------------------------------------------")
  log_msg("    years total / usable TX90p / usable TN10p : ", nrow(g283), " / ",
          g283[usable_tmax == TRUE, .N], " / ", g283[usable_tmin == TRUE, .N])
  log_msg("    locked NA threshold days                  : TX90 ",
          g283$locked_na_tx90_days[1], " ; TN10 ", g283$locked_na_tn10_days[1], " of ", N_AMJ_DAYS)
  log_msg("    years with zero valid Tmax days           : ",
          paste(sort(g283[flag_no_valid_tmax == TRUE, year]), collapse = ", "))
  log_msg("    in-base years with < 29 replicates        : ",
          g283[flag_reduced_replicates == TRUE, .N])
  log_msg("    trend-period usable years (1996-2025)     : ",
          g283[in_trend_period == TRUE & usable_tmax == TRUE, .N], " (TX90p) / ",
          g283[in_trend_period == TRUE & usable_tmin == TRUE, .N], " (TN10p) of ", N_TREND_YEARS)
  record_check("grid_283_trend_coverage_reported", TRUE,
               sprintf("%d of %d trend years usable for TX90p; inclusion is a STEP 06 decision",
                       g283[in_trend_period == TRUE & usable_tmax == TRUE, .N], N_TREND_YEARS),
               severity = "WARN")
}

## --- per-grid summary --------------------------------------------------------
grid_sum <- ar[, .(
  lat = lat[1], lon = lon[1],
  n_years                 = .N,
  n_usable_tmax_full      = sum(usable_tmax),
  n_usable_tmin_full      = sum(usable_tmin),
  n_usable_tmax_trend     = sum(usable_tmax & in_trend_period),
  n_usable_tmin_trend     = sum(usable_tmin & in_trend_period),
  n_flagged_years         = sum(any_flag),
  n_reduced_replicates    = sum(flag_reduced_replicates),
  n_no_data_years         = sum(flag_no_valid_tmax | flag_no_valid_tmin),
  locked_na_tx90_days     = locked_na_tx90_days[1],
  locked_na_tn10_days     = locked_na_tn10_days[1],
  mean_TX90p_full         = mean(TX90p, na.rm = TRUE),
  mean_TN10p_full         = mean(TN10p, na.rm = TRUE),
  mean_TX90p_trend        = mean(TX90p[in_trend_period], na.rm = TRUE),
  mean_TN10p_trend        = mean(TN10p[in_trend_period], na.rm = TRUE)
), by = grid_id][order(as.integer(grid_id))]
grid_sum[, `:=`(complete_series_tmax_trend = n_usable_tmax_trend == N_TREND_YEARS,
                complete_series_tmin_trend = n_usable_tmin_trend == N_TREND_YEARS)]

record_check("grid_summary_row_count", nrow(grid_sum) == EXP_GRIDS,
             sprintf("%d grids summarised", nrow(grid_sum)))
log_msg("Grids with a complete 30-year trend series : TX90p ",
        grid_sum[complete_series_tmax_trend == TRUE, .N], " / ", EXP_GRIDS,
        " ; TN10p ", grid_sum[complete_series_tmin_trend == TRUE, .N], " / ", EXP_GRIDS)

## --- data-limitation table (one row per flagged grid-year) -------------------
lims <- ar[any_flag == TRUE, .(
  grid_id, year, lat, lon, TX90p, TN10p,
  n_valid_tmax, n_valid_tmin, completeness_tmax_pct, completeness_tmin_pct,
  bootstrap_applied, n_rep_used_tmax, n_rep_used_tmin,
  locked_na_tx90_days, locked_na_tn10_days,
  data_quality_tmax, data_quality_tmin, in_trend_period)]
lims[, limitation := fifelse(data_quality_tmax == "NO_DATA" | data_quality_tmin == "NO_DATA",
                             "no valid AMJ observations in this grid-year",
                             fifelse(data_quality_tmax == "LOW_COMPLETENESS" |
                                       data_quality_tmin == "LOW_COMPLETENESS",
                                     sprintf("annual valid days below advisory %d/%d", ANNUAL_MIN_VALID, N_AMJ_DAYS),
                                     "fewer than 29 bootstrap replicates (baseline pool below 120/150)"))]
setorder(lims, grid_id, year)
log_msg("Data-limitation rows : ", nrow(lims))

## --- scientific plausibility (reporting only) --------------------------------
log_msg("Mean TX90p 1976-1980 / 1981-2010 / 2011-2025 : ",
        sprintf("%.2f / %.2f / %.2f",
                ar[year < BASE_START, mean(TX90p, na.rm = TRUE)],
                ar[in_baseline == TRUE, mean(TX90p, na.rm = TRUE)],
                ar[year > BASE_END, mean(TX90p, na.rm = TRUE)]))
log_msg("Mean TN10p 1976-1980 / 1981-2010 / 2011-2025 : ",
        sprintf("%.2f / %.2f / %.2f",
                ar[year < BASE_START, mean(TN10p, na.rm = TRUE)],
                ar[in_baseline == TRUE, mean(TN10p, na.rm = TRUE)],
                ar[year > BASE_END, mean(TN10p, na.rm = TRUE)]))
log_msg("Mean TX90p / TN10p over trend period 1996-2025 : ",
        sprintf("%.2f / %.2f", ar_trend[, mean(TX90p, na.rm = TRUE)],
                ar_trend[, mean(TN10p, na.rm = TRUE)]))
bt <- ar[in_baseline == TRUE, mean(TX90p, na.rm = TRUE)]
bn <- ar[in_baseline == TRUE, mean(TN10p, na.rm = TRUE)]
record_check("baseline_mean_TX90p_plausible", is.finite(bt) && bt >= 5 && bt <= 20,
             sprintf("1981-2010 mean TX90p = %.3f%%", bt), severity = "WARN")
record_check("baseline_mean_TN10p_plausible", is.finite(bn) && bn >= 5 && bn <= 20,
             sprintf("1981-2010 mean TN10p = %.3f%%", bn), severity = "WARN")

## ---------------------------------------------------------------------------
## 11. WRITE LOCKED OUTPUTS + READ-BACK
## ---------------------------------------------------------------------------
log_head("11. WRITE LOCKED ANALYSIS-READY OUTPUTS")

out_dir <- dirname(OUT_FULL)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

arrow::write_parquet(ar,       OUT_FULL,  compression = "snappy")
arrow::write_parquet(ar_trend, OUT_TREND, compression = "snappy")
fwrite(grid_sum, GRID_FILE)
fwrite(lims,     LIM_FILE)
log_msg("Analysis-ready full  : ", OUT_FULL, "  (", format(nrow(ar), big.mark = ","),
        " rows x ", ncol(ar), " cols, ", sprintf("%.2f MB", file.size(OUT_FULL) / 1024^2), ")")
log_msg("Analysis-ready trend : ", OUT_TREND, "  (", format(nrow(ar_trend), big.mark = ","),
        " rows x ", ncol(ar_trend), " cols, ", sprintf("%.2f MB", file.size(OUT_TREND) / 1024^2), ")")
log_msg("Grid summary         : ", GRID_FILE, "  (", nrow(grid_sum), " rows)")
log_msg("Data limitations     : ", LIM_FILE,  "  (", nrow(lims), " rows)")

rb  <- as.data.table(arrow::read_parquet(OUT_FULL));  setorder(rb, year, grid_id)
rbt <- as.data.table(arrow::read_parquet(OUT_TREND)); setorder(rbt, year, grid_id)
record_check("readback_full", nrow(rb) == EXP_FULL_ROWS && identical(names(rb), names(ar)) &&
               eq(rb$TX90p, ar$TX90p) && eq(rb$TN10p, ar$TN10p),
             sprintf("%s rows, %d cols, indices identical after round-trip",
                     format(nrow(rb), big.mark = ","), ncol(rb)))
record_check("readback_trend", nrow(rbt) == EXP_TREND_ROWS && eq(rbt$TX90p, ar_trend$TX90p),
             sprintf("%s rows, indices identical after round-trip", format(nrow(rbt), big.mark = ",")))
record_check("readback_flags_preserved",
             rb[any_flag == TRUE, .N] == ar[any_flag == TRUE, .N] &&
               rb[usable_tmax == TRUE, .N] == ar[usable_tmax == TRUE, .N],
             sprintf("%d flagged and %d usable-TX90p rows survived the round-trip",
                     rb[any_flag == TRUE, .N], rb[usable_tmax == TRUE, .N]))
rm(rb, rbt)

record_check("step04_inputs_unmodified",
             file.exists(IN_FULL) && file.exists(IN_TREND),
             "STEP 04 outputs opened read-only, never rewritten")

qc_tab <- rbindlist(CHECKS)
qc_tab[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
              script   = "05_AMJ_TX90p_TN10p_analysis_ready.R")]
fwrite(qc_tab, QC_FILE)
log_msg("QC summary saved     : ", QC_FILE, "  (", nrow(qc_tab), " checks)")

log_head("RUN COMPLETE")
log_msg("Checks passed        : ", qc_tab[status == "PASS", .N], " / ", nrow(qc_tab))
log_msg("Warnings             : ", qc_tab[status == "WARN", .N])
log_msg("Fatal                : ", qc_tab[status == "FATAL", .N])
log_msg("Rows preserved       : ", format(nrow(ar), big.mark = ","), " full / ",
        format(nrow(ar_trend), big.mark = ","), " trend")
log_msg("Elapsed              : ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
log_msg("NOTE                 : no MK / Sen trend estimation performed; that is STEP 06.")

.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
writeLines(.LOG, LOG_FILE)
cat("\nConsole log written to:", LOG_FILE, "\n")

invisible(ar)
###############################################################################
## END OF SCRIPT
###############################################################################