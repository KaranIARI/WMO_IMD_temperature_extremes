###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 06 : AMJ absolute / duration indices
##            WSDI, TXx, TXn, TNx, TNn, DTR
##  283 grids | 1976-2025 | AMJ 01 Apr - 30 Jun (91 days)
##
##  REVISION: full debug pass. Scientific methodology and assumptions UNCHANGED.
##
## ---------------------------------------------------------------------------
##  BUGS FOUND AND FIXED IN THIS REVISION  (none alter the science)
## ---------------------------------------------------------------------------
##  B1  add_test() recycling. T11 passes length-2 got/want, so data.table()
##      recycled the scalar columns and emitted TWO rows for one test (hence
##      "12/12" for 11 tests). got/want are now collapsed to single strings, so
##      every test occupies exactly one row and `ok` is always scalar.
##  B2  Length-2 subscript in three record_check() conditions. `tst[test == nm,
##      ok]` returned c(TRUE,TRUE) for T11 and `&&` rejects length > 1 under
##      R >= 4.3. Replaced by ok_of(), which is scalar by construction and also
##      fails closed if a test name is misspelled (logical(0) would otherwise
##      make all() return TRUE and silently pass a check that never ran).
##  B3  Row adjacency was assumed to mean calendar adjacency. If a grid-year
##      were missing whole rows (absent, not NA), two non-consecutive dates
##      would sit next to each other and a spell could bridge the gap. This
##      contradicts the locked rule that unevaluable days must never be
##      bridged - an absent row is more unevaluable than an NA row. compute_wsdi
##      now breaks a run whenever ord jumps by more than one day. With a
##      complete 91-day season this changes nothing; it only removes a silent
##      failure mode. Two new unit tests (T12, T13) cover it.
##  B4  ifelse() evaluates both branches, so the "all NA" guards in several
##      record_check details still called min()/max() on empty input, emitting
##      warnings and Inf. Replaced by scalar-safe helpers.
##  B5  setcolorder() was given a fixed column list that omitted the optional
##      s05_TX90_days_direct column and would error outright if any listed
##      column were absent. Now intersect/setdiff based.
##  B6  Column name `.brk` used a dot prefix inside a data.table `by`. Renamed
##      spell_break to avoid collision with data.table's reserved symbols.
##  B7  options(error = NULL) is now set explicitly, so a lingering
##      options(error = recover) from an earlier session cannot drop the run
##      into a browser() prompt that looks like a hang.
##  B8  QC-CSV assembly re-derived got/want via vapply on what are now plain
##      character scalars. Simplified to direct column use.
##
## ---------------------------------------------------------------------------
##  INPUTS (ALL READ-ONLY - nothing upstream is recomputed or rewritten)
##    1. STEP 05 : AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet
##                 -> authoritative grid roster (283), lat/lon, and
##                    TX90_days_direct used as an independent cross-check.
##    2. STEP 02 : AMJ_TX90_TN10_thresholds_1981_2010.parquet   [LOCKED]
##                 -> calendar-day TX90 thresholds for the WSDI warm-day test.
##    3. DAILY   : IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet
##                 -> the daily Tmax / Tmin observations.
##
##  WHY THE DAILY FILE IS REQUIRED
##    The STEP 05 parquet is aggregated to grid-year (14,150 rows) and holds no
##    daily temperatures. TXx, TXn, TNx, TNn, DTR and WSDI are all daily
##    constructs. STEP 05 remains the primary reference for the grid roster,
##    coordinates, row shape and cross-checks; TX90p / TN10p are never touched.
##
##  OUTPUTS
##    A. AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet   14,150 rows
##    B. AMJ_STEP06_QC.csv        + AMJ_STEP06_QC_incomplete_grid_years.csv
##    C. AMJ_STEP06_LOG.txt       - full console log + sessionInfo
##
## ---------------------------------------------------------------------------
##  INDEX DEFINITIONS  (unchanged)
## ---------------------------------------------------------------------------
##  WSDI (Warm Spell Duration Index, AMJ-restricted)
##    * The AMJ season is ONE continuous 91-day sequence per grid-year. The
##      spell counter is NEVER reset at 30 Apr / 01 May or 31 May / 01 Jun; a
##      spell may cross month boundaries freely.
##    * 31 Mar -> 01 Apr and 30 Jun -> 01 Jul are NEVER connected: only the 91
##      AMJ days enter the sequence.
##    * Warm day = Tmax > the corresponding calendar-day TX90 threshold
##      (strict inequality, matched on grid_id + MMDD).
##    * A day with NA Tmax OR NA TX90 threshold is UNEVALUABLE and BREAKS the
##      spell, exactly as a cold day does. Never bridged, never imputed.
##    * WSDI = total days belonging to spells of >= 6 consecutive warm days.
##    * WSDI is NA only when the grid-year has zero evaluable Tmax days;
##      otherwise a genuine absence of spells gives WSDI = 0.
##
##  TXx / TXn = max / min valid AMJ Tmax     TNx / TNn = max / min valid AMJ Tmin
##    Over non-NA observations only. No threshold involved.
##  DTR = mean(Tmax - Tmin) over days where BOTH are valid; NA if no paired day.
##
##  Any index returns NA when it has no valid input. Nothing is imputed. All
##  14,150 grid-year rows are preserved; incomplete grid-years are FLAGGED.
##
##  Report as AMJ-RESTRICTED WSDI: spells straddling 01 Apr or 30 Jun are
##  truncated by design, so values are not comparable to annual ETCCDI WSDI.
##
##  OUT OF SCOPE: Mann-Kendall, Hamed-Rao, Sen slope, FDR, trend
##  classification, mapping.
##
##  Non-interactive, deterministic: runs start to finish under Rscript.
###############################################################################

## ---------------------------------------------------------------------------
## 0. ENVIRONMENT
## ---------------------------------------------------------------------------
rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE, scipen = 999)
options(error = NULL)          # [B7] neutralise any inherited error handler
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
S05_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90p_TN10p_ANALYSIS_READY_1976_2025.parquet"
THR_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_TX90_TN10_thresholds_1981_2010.parquet"
DAILY_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/IMD_283grids_Tmax_Tmin_MarJul_1976_2025.parquet"

OUT_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet"
QC_FILE    <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_QC.csv"
INC_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_QC_incomplete_grid_years.csv"
LOG_FILE   <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_LOG.txt"

YR_START <- 1976L; YR_END <- 2025L; N_YEARS <- 50L
BASE_START <- 1981L; BASE_END <- 2010L
TREND_START <- 1996L; TREND_END <- 2025L
MD_START <- 401L; MD_END <- 630L
N_AMJ_DAYS <- 91L
EXP_GRIDS <- 283L
EXP_ROWS  <- 14150L                     # 283 x 50
EXP_DAILY_AMJ_ROWS <- 1287650L          # 283 x 91 x 50
WSDI_MIN_LEN <- 6L                      # >= 6 consecutive warm days
TOL <- 1e-9

## Plausibility envelopes (reporting / WARN only; nothing is altered)
TX_LO <- 0;   TX_HI <- 60
TN_LO <- -25; TN_HI <- 45
DTR_LO <- 0;  DTR_HI <- 30

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
  ok <- isTRUE(all(passed)) && length(passed) > 0L
  status <- if (ok) "PASS" else severity
  CHECKS[[length(CHECKS) + 1L]] <<- data.table(check = name, status = status, detail = detail)
  log_msg(sprintf("  %-5s | %-44s | %s", status, name, detail))
  if (!ok && severity == "FATAL")
    stop("FATAL QC failure: ", name, " -- ", detail, call. = FALSE)
  invisible(ok)
}

## [B4] scalar-safe summaries: never call min/max on empty input
n_fin   <- function(v) sum(is.finite(v))
sc_min  <- function(v, f = "%.2f") { v <- v[is.finite(v)]; if (!length(v)) "NA" else sprintf(f, min(v)) }
sc_max  <- function(v, f = "%.2f") { v <- v[is.finite(v)]; if (!length(v)) "NA" else sprintf(f, max(v)) }
sc_mean <- function(v, f = "%.2f") { v <- v[is.finite(v)]; if (!length(v)) "NA" else sprintf(f, mean(v)) }
sc_med  <- function(v, f = "%.1f") { v <- v[is.finite(v)]; if (!length(v)) "NA" else sprintf(f, stats::median(v)) }
rng_txt <- function(v, f = "%.2f") { v <- v[is.finite(v)]
if (!length(v)) "all NA" else sprintf(paste0("[", f, ", ", f, "]"), min(v), max(v)) }

pick_col <- function(dt, candidates, label, file_label, required = TRUE) {
  hit <- names(dt)[tolower(names(dt)) %in% tolower(candidates)]
  if (length(hit) < 1L) {
    if (!required) return(NA_character_)
    stop("SCHEMA ERROR in ", file_label, ".\n  Could not resolve a '", label,
         "' column. Looked for: ", paste(candidates, collapse = ", "),
         "\n  Columns present: ", paste(names(dt), collapse = ", "),
         "\n  STEP 06 will not guess a substitute.", call. = FALSE)
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
safe_max <- function(v) { v <- v[!is.na(v)]; if (!length(v)) NA_real_ else max(v) }
safe_min <- function(v) { v <- v[!is.na(v)]; if (!length(v)) NA_real_ else min(v) }

## ---------------------------------------------------------------------------
##  WARM-DAY DERIVATION AND WSDI ENGINE
##  The deterministic unit tests and the production run call these exact
##  functions, so the tests validate the code that produces the output.
## ---------------------------------------------------------------------------

## Warm only if BOTH the observation and the threshold exist and Tmax > TX90.
## Unevaluable days become FALSE, which is what breaks a spell. `FALSE & NA` is
## FALSE in R, so no NA can leak into the run logic.
derive_warm <- function(tmax, thr) {
  ev <- !is.na(tmax) & !is.na(thr)
  ev & tmax > thr
}

## WSDI over one continuous ordered sequence per (grid_id, year).
##   `ord` must be a calendar day number (consecutive days differ by exactly 1).
##   A run is broken by (a) a non-warm day, or (b) [B3] a jump in `ord`, which
##   means a calendar day is absent from the table entirely. Both are treated
##   identically: the spell terminates and never bridges.
##   Month boundaries are adjacent rows in the same group and pass through
##   untouched; group boundaries stop any spell at 01 Apr / 30 Jun and between
##   grids and years.
compute_wsdi <- function(d, min_len = WSDI_MIN_LEN) {
  stopifnot(all(c("grid_id", "year", "ord", "warm") %in% names(d)))
  if (anyNA(d$warm)) stop("compute_wsdi(): `warm` must be logical with no NA.", call. = FALSE)
  if (anyNA(d$ord))  stop("compute_wsdi(): `ord` must not contain NA.", call. = FALSE)
  d <- copy(d)
  setorder(d, grid_id, year, ord)
  d[, gap_before := c(TRUE, diff(ord) != 1L), by = .(grid_id, year)]     # [B3]
  d[, spell_break := cumsum(!warm | gap_before), by = .(grid_id, year)]  # [B6]
  runs <- d[warm == TRUE, .(len = .N), by = .(grid_id, year, spell_break)]
  if (!nrow(runs))
    return(data.table(grid_id = character(0), year = integer(0),
                      WSDI = integer(0), longest_warm_spell = integer(0),
                      n_spells_ge_min = integer(0), n_warm_days = integer(0)))
  runs[, .(WSDI               = as.integer(sum(len[len >= min_len])),
           longest_warm_spell = as.integer(max(len)),
           n_spells_ge_min    = as.integer(sum(len >= min_len)),
           n_warm_days        = as.integer(sum(len))),
       by = .(grid_id, year)]
}

t_start <- Sys.time()
log_head("STEP 06 | AMJ WSDI, TXx, TXn, TNx, TNn, DTR (1976-2025, 283 grids)")
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("Season           : 01 Apr - 30 Jun, MMDD ", MD_START, "-", MD_END,
        " (", N_AMJ_DAYS, " days), ONE continuous sequence per grid-year")
log_msg("WSDI rule        : >= ", WSDI_MIN_LEN,
        " consecutive warm days; unevaluable days and absent days break spells")
log_msg("Warm day         : Tmax > calendar-day TX90 (strict), matched on grid_id + MMDD")
log_msg("Out of scope     : MK / Hamed-Rao / Sen / FDR / classification / mapping")

## ---------------------------------------------------------------------------
## 3. DETERMINISTIC UNIT TESTS  (must pass before any real data is touched)
## ---------------------------------------------------------------------------
log_head("3. DETERMINISTIC WSDI UNIT TESTS")

.ref_dates <- seq(as.Date("2001-04-01"), as.Date("2001-06-30"), by = "day")
stopifnot(length(.ref_dates) == N_AMJ_DAYS)
.ref_md <- as.integer(data.table::month(.ref_dates)) * 100L +
  as.integer(data.table::mday(.ref_dates))

## Synthetic grid-year, then the SAME derive_warm() + compute_wsdi() path.
## `ord` is a true day number so the [B3] gap rule is exercised realistically.
make_case <- function(id, warm_md, na_md = integer(0), thr_na_md = integer(0),
                      drop_md = integer(0)) {
  d <- data.table(grid_id = id, year = 2001L,
                  ord = as.integer(as.IDate(.ref_dates)), md = .ref_md)
  d[, thr  := 35.0]
  d[, tmax := 30.0]                       # below threshold -> not warm
  d[md %in% warm_md,   tmax := 40.0]      # above threshold -> warm
  d[md %in% na_md,     tmax := NA_real_]  # missing observation -> unevaluable
  d[md %in% thr_na_md, thr  := NA_real_]  # missing threshold   -> unevaluable
  if (length(drop_md)) d <- d[!md %in% drop_md]   # row absent entirely
  d[, warm := derive_warm(tmax, thr)]
  d[]
}
run_case <- function(d) {
  r <- compute_wsdi(d)
  if (!nrow(r)) return(list(WSDI = 0L, longest = 0L))
  list(WSDI = r$WSDI[1], longest = r$longest_warm_spell[1])
}

TESTS <- list()
## [B1] got/want collapsed to single strings -> exactly one row per test.
add_test <- function(name, got, want, note) {
  ok <- identical(as.integer(got), as.integer(want))
  TESTS[[length(TESTS) + 1L]] <<- data.table(
    test = name,
    got  = paste(as.integer(got),  collapse = ","),
    want = paste(as.integer(want), collapse = ","),
    ok   = ok,
    note = note)
}

## T1 - THE MANDATED CASE: 27 Apr - 02 May, six consecutive days across the
## April/May boundary. Must give WSDI = 6.
t1 <- run_case(make_case("T1", c(427L, 428L, 429L, 430L, 501L, 502L)))
add_test("T1_apr27_to_may02_six_day_spell", t1$WSDI, 6L,
         "spell crosses 30 Apr / 01 May without reset")

## T2 - same spell, 30 Apr has NA Tmax -> 3 + 2 fragments -> 0.
t2 <- run_case(make_case("T2", c(427L, 428L, 429L, 430L, 501L, 502L), na_md = 430L))
add_test("T2_NA_tmax_breaks_spell", t2$WSDI, 0L,
         "NA Tmax on 30 Apr breaks continuity; fragments discarded")

## T3 - same spell, 30 Apr THRESHOLD is NA -> 0.
t3 <- run_case(make_case("T3", c(427L, 428L, 429L, 430L, 501L, 502L), thr_na_md = 430L))
add_test("T3_NA_threshold_breaks_spell", t3$WSDI, 0L,
         "NA TX90 threshold is unevaluable and breaks continuity")

## T4 - 29 May - 04 Jun, seven days across the May/June boundary.
t4 <- run_case(make_case("T4", c(529L, 530L, 531L, 601L, 602L, 603L, 604L)))
add_test("T4_may29_to_jun04_seven_day_spell", t4$WSDI, 7L,
         "spell crosses 31 May / 01 Jun without reset")

## T5 - a five-day spell contributes nothing.
t5 <- run_case(make_case("T5", c(510L, 511L, 512L, 513L, 514L)))
add_test("T5_five_day_spell_excluded", t5$WSDI, 0L, "below the 6-day minimum")

## T6 - two qualifying spells (6 and 7 days) sum to 13.
t6 <- run_case(make_case("T6", c(405L:410L, 620L:626L)))
add_test("T6_two_spells_sum", t6$WSDI, 13L, "6 + 7 days, both >= minimum")

## T7 - first six AMJ days warm; nothing inherited from March.
t7 <- run_case(make_case("T7", c(401L:406L)))
add_test("T7_season_start_no_march_bridge", t7$WSDI, 6L,
         "01-06 Apr counted; 31 Mar -> 01 Apr never connected")

## T8 - last five AMJ days warm; nothing continues into July.
t8 <- run_case(make_case("T8", c(626L:630L)))
add_test("T8_season_end_no_july_bridge", t8$WSDI, 0L,
         "5 days at season end; 30 Jun -> 01 Jul never connected")

## T9 - every AMJ day warm.
t9 <- run_case(make_case("T9", .ref_md))
add_test("T9_all_warm", t9$WSDI, N_AMJ_DAYS, "entire 91-day season is one spell")

## T10 - no warm day at all.
t10 <- run_case(make_case("T10", integer(0)))
add_test("T10_no_warm_days", t10$WSDI, 0L, "no spell present")

## T11 - two grid-years processed together must not bleed into each other.
mixA <- make_case("A", c(625L:630L)); mixB <- make_case("B", c(401L:405L))
mix  <- compute_wsdi(rbind(mixA, mixB))
add_test("T11_no_bleed_between_grid_years",
         c(mix[grid_id == "A", WSDI], mix[grid_id == "B", WSDI]), c(6L, 0L),
         "grid A 6-day spell counted; grid B 5-day spell not extended by grid A")

## T12 - [B3] the 30 Apr ROW IS ABSENT (not NA). Row adjacency would wrongly
## bridge 29 Apr to 01 May and yield 6; the gap rule must give 0.
t12 <- run_case(make_case("T12", c(427L, 428L, 429L, 430L, 501L, 502L), drop_md = 430L))
add_test("T12_absent_row_breaks_spell", t12$WSDI, 0L,
         "missing calendar day is a gap, not an adjacency; never bridged")

## T13 - scrambled input order must give the same answer as T1 (setorder works).
set.seed(1L)
scr <- make_case("T13", c(427L, 428L, 429L, 430L, 501L, 502L))
scr <- scr[sample.int(nrow(scr))]
add_test("T13_row_order_independent", run_case(scr)$WSDI, 6L,
         "compute_wsdi sorts internally; input row order is irrelevant")

tst <- rbindlist(TESTS)

## [B2] scalar-by-construction accessor; fails closed on an unknown test name.
ok_of <- function(nm) {
  v <- tst[test == nm, ok]
  length(v) == 1L && isTRUE(v)
}

for (i in seq_len(nrow(tst)))
  log_msg(sprintf("  %-4s | %-38s | got %-8s want %-8s | %s",
                  if (tst$ok[i]) "PASS" else "FAIL", tst$test[i],
                  tst$got[i], tst$want[i], tst$note[i]))

record_check("unit_test_table_one_row_per_test",
             nrow(tst) == length(TESTS) && !anyDuplicated(tst$test),
             sprintf("%d tests, %d rows, %d duplicated test names",
                     length(TESTS), nrow(tst), sum(duplicated(tst$test))))
record_check("wsdi_deterministic_unit_tests", all(tst$ok),
             sprintf("%d of %d passed; T1 (27 Apr - 02 May) = %s, required 6",
                     sum(tst$ok), nrow(tst), tst[test == "T1_apr27_to_may02_six_day_spell", got]))
record_check("wsdi_month_boundary_continuity",
             ok_of("T1_apr27_to_may02_six_day_spell") &&
               ok_of("T4_may29_to_jun04_seven_day_spell"),
             "spells cross 30 Apr/01 May and 31 May/01 Jun without the counter resetting")
record_check("wsdi_NA_days_break_spells",
             ok_of("T2_NA_tmax_breaks_spell") &&
               ok_of("T3_NA_threshold_breaks_spell") &&
               ok_of("T12_absent_row_breaks_spell"),
             "NA Tmax, NA TX90 and absent rows all break continuity; never bridged")
record_check("wsdi_season_edges_not_bridged",
             ok_of("T7_season_start_no_march_bridge") &&
               ok_of("T8_season_end_no_july_bridge") &&
               ok_of("T11_no_bleed_between_grid_years"),
             "31 Mar -> 01 Apr and 30 Jun -> 01 Jul never connected; no cross-group bleed")
record_check("wsdi_order_independent", ok_of("T13_row_order_independent"),
             "result invariant to input row order")
rm(mixA, mixB, mix, scr)

## ---------------------------------------------------------------------------
## 4. READ INPUTS (READ-ONLY) AND RESOLVE SCHEMAS
## ---------------------------------------------------------------------------
log_head("4. READ INPUTS AND RESOLVE SCHEMAS")

for (f in c(S05_FILE, THR_FILE, DAILY_FILE))
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)
.mtime_before <- file.mtime(c(S05_FILE, THR_FILE, DAILY_FILE))

## --- 4a. STEP 05 analysis-ready ---------------------------------------------
s05 <- as.data.table(arrow::read_parquet(S05_FILE))
log_msg("STEP 05 file     : ", format(nrow(s05), big.mark = ","), " rows x ", ncol(s05), " cols")
log_msg("STEP 05 columns  : ", paste(names(s05), collapse = ", "))
S_GRID <- pick_col(s05, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id", "STEP 05")
S_YEAR <- pick_col(s05, c("year", "yr"), "year", "STEP 05")
S_LAT  <- pick_col(s05, c("lat", "latitude", "y"), "lat", "STEP 05")
S_LON  <- pick_col(s05, c("lon", "long", "longitude", "x"), "lon", "STEP 05")
S_TXD  <- pick_col(s05, "TX90_days_direct", "TX90_days_direct", "STEP 05", required = FALSE)
log_msg("STEP 05 mapping  : grid_id<-", S_GRID, " year<-", S_YEAR, " lat<-", S_LAT,
        " lon<-", S_LON, " TX90_days_direct<-", ifelse(is.na(S_TXD), "(absent)", S_TXD))

keep05 <- c(S_GRID, S_YEAR, S_LAT, S_LON, if (!is.na(S_TXD)) S_TXD)
s05 <- s05[, keep05, with = FALSE]
setnames(s05, c("grid_id", "year", "lat", "lon",
                if (!is.na(S_TXD)) "s05_TX90_days_direct"))
s05[, `:=`(grid_id = as.character(grid_id), year = as.integer(year))]

record_check("step05_row_count", nrow(s05) == EXP_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(s05), big.mark = ","),
                     format(EXP_ROWS, big.mark = ",")))
record_check("step05_grid_count", uniqueN(s05$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(s05$grid_id), EXP_GRIDS))
record_check("step05_year_coverage", identical(sort(unique(s05$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(s05$year), min(s05$year), max(s05$year)))
record_check("step05_no_duplicate_grid_year", !anyDuplicated(s05, by = c("grid_id", "year")),
             sprintf("%d duplicates", sum(duplicated(s05, by = c("grid_id", "year")))))

GRIDS <- sort(unique(s05$grid_id))
meta  <- unique(s05[, .(grid_id, lat, lon)])
record_check("step05_unique_coords_per_grid", nrow(meta) == EXP_GRIDS,
             sprintf("%d unique grid/lat/lon triplets", nrow(meta)))

## --- 4b. STEP 02 locked thresholds ------------------------------------------
thr <- as.data.table(arrow::read_parquet(THR_FILE))
log_msg("STEP 02 file     : ", format(nrow(thr), big.mark = ","), " rows x ", ncol(thr), " cols")
T_GRID <- pick_col(thr, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id", "STEP 02")
T_MD   <- pick_col(thr, c("target_md", "md", "mmdd"), "target_md", "STEP 02", required = FALSE)
T_DAY  <- pick_col(thr, c("target_day", "calendar_day"), "target_day", "STEP 02", required = FALSE)
T_TX   <- pick_col(thr, c("tx90", "tx90_threshold", "tmax_p90"), "tx90", "STEP 02")
if (is.na(T_MD)) {
  if (is.na(T_DAY))
    stop("SCHEMA ERROR in STEP 02: neither 'target_md' nor 'target_day' present.", call. = FALSE)
  dd <- as.character(thr[[T_DAY]])
  thr[, thr_md := as.integer(substr(dd, 1L, 2L)) * 100L + as.integer(substr(dd, 4L, 5L))]
} else thr[, thr_md := as.integer(get(T_MD))]
thr <- thr[, c(T_GRID, "thr_md", T_TX), with = FALSE]
setnames(thr, c("grid_id", "md", "tx90"))
thr[, `:=`(grid_id = as.character(grid_id), tx90 = as.numeric(tx90))]
log_msg("STEP 02 mapping  : grid_id<-", T_GRID, " MMDD<-",
        ifelse(is.na(T_MD), paste0("derived from ", T_DAY), T_MD), " tx90<-", T_TX)

record_check("threshold_row_count", nrow(thr) == EXP_GRIDS * N_AMJ_DAYS,
             sprintf("observed %s vs expected %s (283 x 91)",
                     format(nrow(thr), big.mark = ","),
                     format(EXP_GRIDS * N_AMJ_DAYS, big.mark = ",")))
record_check("threshold_no_duplicate_key", !anyDuplicated(thr, by = c("grid_id", "md")),
             sprintf("%d duplicated (grid_id, MMDD)", sum(duplicated(thr, by = c("grid_id", "md")))))
record_check("threshold_days_in_AMJ", all(thr$md >= MD_START & thr$md <= MD_END),
             sprintf("MMDD range %d-%d", min(thr$md), max(thr$md)))
record_check("threshold_grids_match_step05", identical(sort(unique(thr$grid_id)), GRIDS),
             sprintf("%d threshold grids vs %d STEP 05 grids, set-identical",
                     uniqueN(thr$grid_id), length(GRIDS)))
log_msg("Locked NA TX90   : ", sum(is.na(thr$tx90)),
        " calendar grid-days (unevaluable for WSDI, by design)")

## --- 4c. daily observations --------------------------------------------------
dly <- as.data.table(arrow::read_parquet(DAILY_FILE))
log_msg("Daily file       : ", format(nrow(dly), big.mark = ","), " rows x ", ncol(dly), " cols")
log_msg("Daily columns    : ", paste(names(dly), collapse = ", "))
D_GRID <- pick_col(dly, c("grid_id", "gridid", "grid", "cell_id", "id"), "grid_id", "daily")
D_DATE <- pick_col(dly, c("date", "dates", "obs_date", "time", "day"), "date", "daily")
D_TMAX <- pick_col(dly, c("tmax", "tx", "t_max", "max_temp", "temp_max"), "tmax", "daily")
D_TMIN <- pick_col(dly, c("tmin", "tn", "t_min", "min_temp", "temp_min"), "tmin", "daily")
dly <- dly[, c(D_GRID, D_DATE, D_TMAX, D_TMIN), with = FALSE]
setnames(dly, c("grid_id", "date", "tmax", "tmin"))
dly[, date := coerce_idate(date)]
if (anyNA(dly$date)) stop("Date coercion produced NA values in the daily file.", call. = FALSE)
dly[, `:=`(grid_id = as.character(grid_id),
           year = as.integer(data.table::year(date)),
           md   = as.integer(data.table::month(date)) * 100L + as.integer(data.table::mday(date)),
           tmax = as.numeric(tmax), tmin = as.numeric(tmin))]
log_msg("Daily mapping    : grid_id<-", D_GRID, " date<-", D_DATE,
        " tmax<-", D_TMAX, " tmin<-", D_TMIN)

amj <- dly[md >= MD_START & md <= MD_END & year >= YR_START & year <= YR_END]
rm(dly); invisible(gc(verbose = FALSE))
log_msg("AMJ daily subset : ", format(nrow(amj), big.mark = ","), " rows")

record_check("daily_no_duplicate_grid_date", !anyDuplicated(amj, by = c("grid_id", "date")),
             sprintf("%d duplicate (grid_id, date) rows",
                     sum(duplicated(amj, by = c("grid_id", "date")))))
record_check("daily_grid_count", uniqueN(amj$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(amj$grid_id), EXP_GRIDS))
record_check("daily_year_coverage", identical(sort(unique(amj$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(amj$year), min(amj$year), max(amj$year)))
record_check("daily_amj_day_count", uniqueN(amj$md) == N_AMJ_DAYS,
             sprintf("%d unique MMDD vs expected %d", uniqueN(amj$md), N_AMJ_DAYS))
record_check("daily_amj_row_count", nrow(amj) == EXP_DAILY_AMJ_ROWS,
             sprintf("observed %s vs expected %s (a shortfall means absent rows, not NA values)",
                     format(nrow(amj), big.mark = ","),
                     format(EXP_DAILY_AMJ_ROWS, big.mark = ",")), severity = "WARN")
record_check("daily_grids_match_step05", identical(sort(unique(amj$grid_id)), GRIDS),
             sprintf("daily %d grids vs STEP 05 %d grids, set-identical",
                     uniqueN(amj$grid_id), length(GRIDS)))
log_msg("NA Tmax in AMJ   : ", format(sum(is.na(amj$tmax)), big.mark = ","),
        sprintf("  (%.3f%%)", 100 * mean(is.na(amj$tmax))))
log_msg("NA Tmin in AMJ   : ", format(sum(is.na(amj$tmin)), big.mark = ","),
        sprintf("  (%.3f%%)", 100 * mean(is.na(amj$tmin))))

## ---------------------------------------------------------------------------
## 5. THRESHOLD JOIN AND WARM-DAY DERIVATION
## ---------------------------------------------------------------------------
log_head("5. THRESHOLD JOIN AND WARM-DAY DERIVATION")

setkey(thr, grid_id, md)
amj[thr, on = .(grid_id, md), tx90 := i.tx90]

n_unmatched <- amj[, sum(!paste(grid_id, md) %chin% paste(thr$grid_id, thr$md))]
record_check("zero_unmatched_threshold_keys", n_unmatched == 0L,
             sprintf("%s AMJ rows without a (grid_id, MMDD) threshold record",
                     format(n_unmatched, big.mark = ",")))

amj[, `:=`(evaluable_tmax = !is.na(tmax) & !is.na(tx90),
           paired         = !is.na(tmax) & !is.na(tmin))]
amj[, warm := derive_warm(tmax, tx90)]
record_check("warm_flag_has_no_NA", !anyNA(amj$warm),
             sprintf("%d NA in warm flag", sum(is.na(amj$warm))))
record_check("warm_implies_evaluable", amj[warm == TRUE & evaluable_tmax == FALSE, .N] == 0L,
             sprintf("%d warm days on unevaluable rows",
                     amj[warm == TRUE & evaluable_tmax == FALSE, .N]))
log_msg("Evaluable Tmax   : ", format(sum(amj$evaluable_tmax), big.mark = ","),
        " ; warm days = ", format(sum(amj$warm), big.mark = ","))
log_msg("Paired Tmax+Tmin : ", format(sum(amj$paired), big.mark = ","))

## ---------------------------------------------------------------------------
## 6. COMPUTE INDICES
## ---------------------------------------------------------------------------
log_head("6. COMPUTE WSDI / TXx / TXn / TNx / TNn / DTR")

amj[, ord := as.integer(date)]          # calendar day number: consecutive = +1
setorder(amj, grid_id, year, ord)

## Report any within-season calendar gaps the [B3] rule will act on.
gapchk <- amj[, .(n_gaps = sum(c(FALSE, diff(ord) != 1L))), by = .(grid_id, year)]
n_gap_gy <- gapchk[n_gaps > 0L, .N]
record_check("no_within_season_calendar_gaps", n_gap_gy == 0L,
             sprintf("%d grid-years contain absent calendar days; the gap rule breaks spells there",
                     n_gap_gy), severity = "WARN")
rm(gapchk)

## --- 6a. WSDI (the engine the unit tests validated) --------------------------
t0 <- Sys.time()
wsdi <- compute_wsdi(amj[, .(grid_id, year, ord, warm)])
log_msg("WSDI computed    : ", format(nrow(wsdi), big.mark = ","),
        " grid-years with >= 1 warm day, in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## --- 6b. absolute extremes, DTR and all counts -------------------------------
t0 <- Sys.time()
idx <- amj[, {
  dv <- tmax - tmin
  dv <- dv[!is.na(dv)]
  .(n_AMJ_days        = .N,
    n_obs_tmax        = sum(!is.na(tmax)),
    n_obs_tmin        = sum(!is.na(tmin)),
    n_evaluable_tmax  = sum(evaluable_tmax),
    n_no_threshold_tx = sum(!is.na(tmax) & is.na(tx90)),
    n_paired_days     = length(dv),
    TXx = safe_max(tmax), TXn = safe_min(tmax),
    TNx = safe_max(tmin), TNn = safe_min(tmin),
    DTR     = if (length(dv)) mean(dv) else NA_real_,
    DTR_min = if (length(dv)) min(dv)  else NA_real_,
    DTR_max = if (length(dv)) max(dv)  else NA_real_)
}, by = .(grid_id, year)]
log_msg("Extremes / DTR   : ", format(nrow(idx), big.mark = ","), " grid-years in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
rm(amj); invisible(gc(verbose = FALSE))

## --- 6c. assemble on the full 283 x 50 skeleton ------------------------------
res <- CJ(grid_id = GRIDS, year = YR_START:YR_END, unique = TRUE)
res <- merge(res, idx,  by = c("grid_id", "year"), all.x = TRUE, sort = FALSE)
res <- merge(res, wsdi, by = c("grid_id", "year"), all.x = TRUE, sort = FALSE)
res <- merge(res, meta, by = "grid_id", all.x = TRUE, sort = FALSE)
rm(idx, wsdi); invisible(gc(verbose = FALSE))

cnt0 <- c("n_AMJ_days", "n_obs_tmax", "n_obs_tmin", "n_evaluable_tmax",
          "n_no_threshold_tx", "n_paired_days", "n_warm_days",
          "longest_warm_spell", "n_spells_ge_min")
for (cc in cnt0) {
  if (!cc %in% names(res)) res[, (cc) := 0L]
  set(res, which(is.na(res[[cc]])), cc, 0L)
  res[, (cc) := as.integer(get(cc))]
}

## WSDI: 0 is a real value when evaluable days exist; NA only when none do.
res[, WSDI := as.integer(WSDI)]
res[is.na(WSDI) & n_evaluable_tmax > 0L, WSDI := 0L]
res[n_evaluable_tmax == 0L, `:=`(WSDI = NA_integer_, longest_warm_spell = NA_integer_,
                                 n_spells_ge_min = NA_integer_, n_warm_days = NA_integer_)]

## --- 6d. completeness and flags (advisory; nothing is deleted) ---------------
res[, `:=`(
  completeness_tmax_pct      = round(100 * n_obs_tmax       / N_AMJ_DAYS, 2),
  completeness_tmin_pct      = round(100 * n_obs_tmin       / N_AMJ_DAYS, 2),
  completeness_paired_pct    = round(100 * n_paired_days    / N_AMJ_DAYS, 2),
  completeness_evaluable_pct = round(100 * n_evaluable_tmax / N_AMJ_DAYS, 2),
  flag_incomplete_season     = n_AMJ_days != N_AMJ_DAYS,
  flag_no_obs_tmax           = n_obs_tmax == 0L,
  flag_no_obs_tmin           = n_obs_tmin == 0L,
  flag_no_paired_days        = n_paired_days == 0L,
  flag_no_evaluable_tmax     = n_evaluable_tmax == 0L,
  usable_TXx = !is.na(TXx), usable_TXn = !is.na(TXn),
  usable_TNx = !is.na(TNx), usable_TNn = !is.na(TNn),
  usable_DTR = !is.na(DTR), usable_WSDI = !is.na(WSDI),
  in_baseline     = year >= BASE_START & year <= BASE_END,
  in_trend_period = year >= TREND_START & year <= TREND_END,
  season = "AMJ", season_days = N_AMJ_DAYS,
  wsdi_min_spell_len = WSDI_MIN_LEN,
  wsdi_sequence = "continuous_91day_AMJ_no_month_reset",
  step06_stamp = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"))]
res[, any_flag := flag_incomplete_season | flag_no_obs_tmax | flag_no_obs_tmin |
      flag_no_paired_days | flag_no_evaluable_tmax]

if ("s05_TX90_days_direct" %in% names(s05))
  res <- merge(res, s05[, .(grid_id, year, s05_TX90_days_direct)],
               by = c("grid_id", "year"), all.x = TRUE, sort = FALSE)

## [B5] partial, order-safe column arrangement
ord_cols <- c(
  "year", "grid_id", "lat", "lon",
  "WSDI", "TXx", "TXn", "TNx", "TNn", "DTR",
  "usable_WSDI", "usable_TXx", "usable_TXn", "usable_TNx", "usable_TNn", "usable_DTR",
  "n_AMJ_days", "n_obs_tmax", "n_obs_tmin", "n_evaluable_tmax", "n_paired_days",
  "n_no_threshold_tx", "n_warm_days", "longest_warm_spell", "n_spells_ge_min",
  "completeness_tmax_pct", "completeness_tmin_pct",
  "completeness_paired_pct", "completeness_evaluable_pct",
  "DTR_min", "DTR_max",
  "any_flag", "flag_incomplete_season", "flag_no_obs_tmax", "flag_no_obs_tmin",
  "flag_no_paired_days", "flag_no_evaluable_tmax",
  "in_baseline", "in_trend_period",
  "season", "season_days", "wsdi_min_spell_len", "wsdi_sequence", "step06_stamp")
setcolorder(res, c(intersect(ord_cols, names(res)), setdiff(names(res), ord_cols)))
setorder(res, year, grid_id)

## ---------------------------------------------------------------------------
## 7. FINAL QC
## ---------------------------------------------------------------------------
log_head("7. FINAL QC")

record_check("output_row_count", nrow(res) == EXP_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(res), big.mark = ","),
                     format(EXP_ROWS, big.mark = ",")))
record_check("output_grid_count", uniqueN(res$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(res$grid_id), EXP_GRIDS))
record_check("output_year_coverage", identical(sort(unique(res$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(res$year), min(res$year), max(res$year)))
record_check("no_duplicate_grid_year", !anyDuplicated(res, by = c("grid_id", "year")),
             sprintf("%d duplicated (grid_id, year)", sum(duplicated(res, by = c("grid_id", "year")))))
record_check("283_grids_every_year", nrow(res[, .N, by = year][N != EXP_GRIDS]) == 0L,
             sprintf("%d years with grid count != %d",
                     nrow(res[, .N, by = year][N != EXP_GRIDS]), EXP_GRIDS))
record_check("50_years_every_grid", nrow(res[, .N, by = grid_id][N != N_YEARS]) == 0L,
             sprintf("%d grids with year count != %d",
                     nrow(res[, .N, by = grid_id][N != N_YEARS]), N_YEARS))
record_check("complete_grid_x_year_cross",
             nrow(res) == uniqueN(res$grid_id) * uniqueN(res$year),
             sprintf("%d rows vs %d x %d", nrow(res), uniqueN(res$grid_id), uniqueN(res$year)))
record_check("n_AMJ_days_is_91", res[n_AMJ_days != N_AMJ_DAYS, .N] == 0L,
             sprintf("%d grid-years with n_AMJ_days != %d",
                     res[n_AMJ_days != N_AMJ_DAYS, .N], N_AMJ_DAYS), severity = "WARN")
record_check("coords_carried_through", !anyNA(res$lat) && !anyNA(res$lon),
             sprintf("%d NA lat, %d NA lon", sum(is.na(res$lat)), sum(is.na(res$lon))))

## --- WSDI domain -------------------------------------------------------------
record_check("WSDI_non_negative", res[!is.na(WSDI) & WSDI < 0L, .N] == 0L,
             sprintf("min WSDI = %s", sc_min(res$WSDI, "%.0f")))
record_check("WSDI_le_evaluable_tmax_days",
             res[!is.na(WSDI) & WSDI > n_evaluable_tmax, .N] == 0L,
             sprintf("%d violations ; max WSDI = %s",
                     res[!is.na(WSDI) & WSDI > n_evaluable_tmax, .N], sc_max(res$WSDI, "%.0f")))
record_check("WSDI_le_warm_days", res[!is.na(WSDI) & WSDI > n_warm_days, .N] == 0L,
             sprintf("%d rows with WSDI > n_warm_days",
                     res[!is.na(WSDI) & WSDI > n_warm_days, .N]))
record_check("WSDI_le_91", res[!is.na(WSDI) & WSDI > N_AMJ_DAYS, .N] == 0L,
             sprintf("%d rows with WSDI > %d", res[!is.na(WSDI) & WSDI > N_AMJ_DAYS, .N], N_AMJ_DAYS))
record_check("WSDI_NA_iff_no_evaluable_days",
             res[(is.na(WSDI) & n_evaluable_tmax > 0L) |
                   (!is.na(WSDI) & n_evaluable_tmax == 0L), .N] == 0L,
             sprintf("NA WSDI = %d ; zero-evaluable grid-years = %d",
                     res[is.na(WSDI), .N], res[n_evaluable_tmax == 0L, .N]))
record_check("WSDI_zero_iff_no_qualifying_spell",
             res[!is.na(WSDI) & ((WSDI == 0L) != (n_spells_ge_min == 0L)), .N] == 0L,
             sprintf("%d rows inconsistent between WSDI and n_spells_ge_min",
                     res[!is.na(WSDI) & ((WSDI == 0L) != (n_spells_ge_min == 0L)), .N]))
record_check("longest_spell_consistent",
             res[!is.na(WSDI) & n_spells_ge_min > 0L & longest_warm_spell < WSDI_MIN_LEN, .N] == 0L,
             sprintf("max longest_warm_spell = %s", sc_max(res$longest_warm_spell, "%.0f")))
record_check("warm_days_le_evaluable",
             res[!is.na(n_warm_days) & n_warm_days > n_evaluable_tmax, .N] == 0L,
             sprintf("%d rows with n_warm_days > n_evaluable_tmax",
                     res[!is.na(n_warm_days) & n_warm_days > n_evaluable_tmax, .N]))

## --- independent cross-check: warm days vs STEP 05 TX90_days_direct ----------
if ("s05_TX90_days_direct" %in% names(res)) {
  cmpw <- res[!is.na(s05_TX90_days_direct) & !is.na(n_warm_days)]
  nbad <- cmpw[abs(n_warm_days - s05_TX90_days_direct) > TOL, .N]
  record_check("warm_days_match_STEP05_TX90_days_direct", nbad == 0L,
               sprintf("%d of %s comparable grid-years disagree (same locked thresholds, same strict >)",
                       nbad, format(nrow(cmpw), big.mark = ",")))
  rm(cmpw)
} else {
  record_check("step05_warm_day_crosscheck_available", FALSE,
               "TX90_days_direct absent from STEP 05; cross-check skipped", severity = "WARN")
}

## --- absolute extremes and DTR ------------------------------------------------
record_check("TXx_ge_TXn", res[!is.na(TXx) & !is.na(TXn) & TXx < TXn - TOL, .N] == 0L,
             sprintf("%d violations ; TXx %s, TXn %s",
                     res[!is.na(TXx) & !is.na(TXn) & TXx < TXn - TOL, .N],
                     rng_txt(res$TXx), rng_txt(res$TXn)))
record_check("TNx_ge_TNn", res[!is.na(TNx) & !is.na(TNn) & TNx < TNn - TOL, .N] == 0L,
             sprintf("%d violations ; TNx %s, TNn %s",
                     res[!is.na(TNx) & !is.na(TNn) & TNx < TNn - TOL, .N],
                     rng_txt(res$TNx), rng_txt(res$TNn)))
record_check("extremes_NA_iff_no_observations",
             res[(is.na(TXx) & n_obs_tmax > 0L) | (!is.na(TXx) & n_obs_tmax == 0L), .N] == 0L &&
               res[(is.na(TXn) & n_obs_tmax > 0L) | (!is.na(TXn) & n_obs_tmax == 0L), .N] == 0L &&
               res[(is.na(TNx) & n_obs_tmin > 0L) | (!is.na(TNx) & n_obs_tmin == 0L), .N] == 0L &&
               res[(is.na(TNn) & n_obs_tmin > 0L) | (!is.na(TNn) & n_obs_tmin == 0L), .N] == 0L,
             sprintf("NA: TXx %d, TXn %d, TNx %d, TNn %d ; zero-obs grid-years Tmax %d, Tmin %d",
                     res[is.na(TXx), .N], res[is.na(TXn), .N], res[is.na(TNx), .N],
                     res[is.na(TNn), .N], res[n_obs_tmax == 0L, .N], res[n_obs_tmin == 0L, .N]))
record_check("no_infinite_values",
             res[, sum(is.infinite(TXx) | is.infinite(TXn) | is.infinite(TNx) |
                         is.infinite(TNn) | is.infinite(DTR))] == 0L,
             "no +/-Inf produced by max()/min() on all-NA groups")
record_check("DTR_finite_iff_paired_days_positive",
             res[n_paired_days > 0L & (is.na(DTR) | !is.finite(DTR)), .N] == 0L &&
               res[n_paired_days == 0L & !is.na(DTR), .N] == 0L,
             sprintf("%d rows with paired days but non-finite DTR ; %d rows with zero paired days but non-NA DTR",
                     res[n_paired_days > 0L & (is.na(DTR) | !is.finite(DTR)), .N],
                     res[n_paired_days == 0L & !is.na(DTR), .N]))
record_check("DTR_bracketed_by_daily_extremes",
             res[!is.na(DTR) & (DTR < DTR_min - TOL | DTR > DTR_max + TOL), .N] == 0L,
             sprintf("mean DTR always within [daily min, daily max]; DTR %s", rng_txt(res$DTR)))
record_check("DTR_daily_minimum_non_negative",
             res[!is.na(DTR_min) & DTR_min < -TOL, .N] == 0L,
             sprintf("%d grid-years contain a day with Tmax < Tmin",
                     res[!is.na(DTR_min) & DTR_min < -TOL, .N]), severity = "WARN")

## --- plausibility envelopes (WARN only) --------------------------------------
record_check("TXx_TXn_plausible_range",
             res[!is.na(TXx) & (TXx < TX_LO | TXx > TX_HI), .N] == 0L &&
               res[!is.na(TXn) & (TXn < TX_LO | TXn > TX_HI), .N] == 0L,
             sprintf("outside [%g, %g] degC: TXx %d, TXn %d", TX_LO, TX_HI,
                     res[!is.na(TXx) & (TXx < TX_LO | TXx > TX_HI), .N],
                     res[!is.na(TXn) & (TXn < TX_LO | TXn > TX_HI), .N]), severity = "WARN")
record_check("TNx_TNn_plausible_range",
             res[!is.na(TNx) & (TNx < TN_LO | TNx > TN_HI), .N] == 0L &&
               res[!is.na(TNn) & (TNn < TN_LO | TNn > TN_HI), .N] == 0L,
             sprintf("outside [%g, %g] degC: TNx %d, TNn %d", TN_LO, TN_HI,
                     res[!is.na(TNx) & (TNx < TN_LO | TNx > TN_HI), .N],
                     res[!is.na(TNn) & (TNn < TN_LO | TNn > TN_HI), .N]), severity = "WARN")
record_check("DTR_plausible_range",
             res[!is.na(DTR) & (DTR < DTR_LO | DTR > DTR_HI), .N] == 0L,
             sprintf("outside [%g, %g] degC: %d", DTR_LO, DTR_HI,
                     res[!is.na(DTR) & (DTR < DTR_LO | DTR > DTR_HI), .N]), severity = "WARN")

## --- incomplete grid-years explicitly reported -------------------------------
log_head("7b. INCOMPLETE GRID-YEARS (flagged, never deleted)")
log_msg("Flagged grid-years of ", format(nrow(res), big.mark = ","), ":")
log_msg("    incomplete season (< 91 rows) : ", res[flag_incomplete_season == TRUE, .N])
log_msg("    zero Tmax observations        : ", res[flag_no_obs_tmax == TRUE, .N])
log_msg("    zero Tmin observations        : ", res[flag_no_obs_tmin == TRUE, .N])
log_msg("    zero paired Tmax+Tmin days    : ", res[flag_no_paired_days == TRUE, .N])
log_msg("    zero evaluable Tmax days      : ", res[flag_no_evaluable_tmax == TRUE, .N])
log_msg("    any flag                      : ", res[any_flag == TRUE, .N])
log_msg("Usable: WSDI ", res[usable_WSDI == TRUE, .N], " ; TXx ", res[usable_TXx == TRUE, .N],
        " ; TXn ", res[usable_TXn == TRUE, .N], " ; TNx ", res[usable_TNx == TRUE, .N],
        " ; TNn ", res[usable_TNn == TRUE, .N], " ; DTR ", res[usable_DTR == TRUE, .N])

bad_rows <- res[any_flag == TRUE, .(grid_id, year, n_AMJ_days, n_obs_tmax, n_obs_tmin,
                                    n_evaluable_tmax, n_paired_days, n_no_threshold_tx,
                                    WSDI, TXx, TXn, TNx, TNn, DTR)]
setorder(bad_rows, grid_id, year)
if (nrow(bad_rows) > 0L) {
  log_msg("Flagged detail (grid_id | year | nAMJ | obsTx | obsTn | evalTx | paired):")
  for (i in seq_len(min(nrow(bad_rows), 200L)))
    log_msg(sprintf("    %-8s | %d | %2d | %2d | %2d | %2d | %2d",
                    bad_rows$grid_id[i], bad_rows$year[i], bad_rows$n_AMJ_days[i],
                    bad_rows$n_obs_tmax[i], bad_rows$n_obs_tmin[i],
                    bad_rows$n_evaluable_tmax[i], bad_rows$n_paired_days[i]))
  if (nrow(bad_rows) > 200L)
    log_msg("    ... ", nrow(bad_rows) - 200L, " further rows in ", basename(INC_FILE))
}
record_check("incomplete_grid_years_reported", TRUE,
             sprintf("%d flagged grid-years enumerated; all retained in the output", nrow(bad_rows)))
record_check("all_rows_retained", nrow(res) == EXP_ROWS,
             sprintf("%s rows preserved; no grid-year deleted", format(nrow(res), big.mark = ",")))

## --- descriptive summary ------------------------------------------------------
log_msg("WSDI  mean / median / max : ", sc_mean(res$WSDI), " / ", sc_med(res$WSDI),
        " / ", sc_max(res$WSDI, "%.0f"), "   (n usable = ", n_fin(res$WSDI), ")")
log_msg("TXx   mean / max          : ", sc_mean(res$TXx), " / ", sc_max(res$TXx))
log_msg("TXn   mean / min          : ", sc_mean(res$TXn), " / ", sc_min(res$TXn))
log_msg("TNx   mean / max          : ", sc_mean(res$TNx), " / ", sc_max(res$TNx))
log_msg("TNn   mean / min          : ", sc_mean(res$TNn), " / ", sc_min(res$TNn))
log_msg("DTR   mean / range        : ", sc_mean(res$DTR), " / ", rng_txt(res$DTR))
log_msg("Mean WSDI 1976-1980 / 1981-2010 / 2011-2025 : ",
        sc_mean(res[year < BASE_START, WSDI]), " / ",
        sc_mean(res[in_baseline == TRUE, WSDI]), " / ",
        sc_mean(res[year > BASE_END, WSDI]))
log_msg("Mean TXx 1976-1980 / 1981-2010 / 2011-2025  : ",
        sc_mean(res[year < BASE_START, TXx]), " / ",
        sc_mean(res[in_baseline == TRUE, TXx]), " / ",
        sc_mean(res[year > BASE_END, TXx]))

## ---------------------------------------------------------------------------
## 8. WRITE OUTPUTS + READ-BACK + UPSTREAM INTEGRITY
## ---------------------------------------------------------------------------
log_head("8. WRITE OUTPUTS")

out_dir <- dirname(OUT_FILE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
arrow::write_parquet(res, OUT_FILE, compression = "snappy")
log_msg("STEP 06 output   : ", OUT_FILE, "  (", format(nrow(res), big.mark = ","),
        " rows x ", ncol(res), " cols, ", sprintf("%.2f MB", file.size(OUT_FILE) / 1024^2), ")")

eqv <- function(a, b) length(a) == length(b) && all(is.na(a) == is.na(b) & (is.na(a) | a == b))
rb <- as.data.table(arrow::read_parquet(OUT_FILE)); setorder(rb, year, grid_id)
record_check("readback_output",
             nrow(rb) == EXP_ROWS && identical(names(rb), names(res)) &&
               eqv(rb$WSDI, res$WSDI) && eqv(rb$TXx, res$TXx) && eqv(rb$TXn, res$TXn) &&
               eqv(rb$TNx, res$TNx) && eqv(rb$TNn, res$TNn) &&
               isTRUE(all.equal(rb$DTR, res$DTR)),
             sprintf("%s rows, %d cols, all six indices identical after round-trip",
                     format(nrow(rb), big.mark = ","), ncol(rb)))
rm(rb)

.mtime_after <- file.mtime(c(S05_FILE, THR_FILE, DAILY_FILE))
record_check("upstream_files_unchanged", identical(.mtime_before, .mtime_after),
             "STEP 05, STEP 02 and daily files: modification times identical before and after")

qc_tab <- rbindlist(CHECKS)
qc_tab[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
              script   = "06_AMJ_WSDI_TXx_TXn_TNx_TNn_DTR.R")]
## [B8] tst$got / tst$want are already character scalars
qc_tab <- rbind(qc_tab,
                data.table(check  = paste0("unit_", tst$test),
                           status = fifelse(tst$ok, "PASS", "FATAL"),
                           detail = sprintf("got %s, want %s | %s", tst$got, tst$want, tst$note),
                           run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
                           script = "06_AMJ_WSDI_TXx_TXn_TNx_TNn_DTR.R"))
fwrite(qc_tab, QC_FILE)
fwrite(bad_rows, INC_FILE)
log_msg("QC summary saved : ", QC_FILE, "  (", nrow(qc_tab), " entries)")
log_msg("Incomplete rows  : ", INC_FILE, "  (", nrow(bad_rows), " rows)")

log_head("RUN COMPLETE")
log_msg("Unit tests       : ", sum(tst$ok), " / ", nrow(tst), " passed")
log_msg("Checks passed    : ", qc_tab[status == "PASS", .N], " / ", nrow(qc_tab))
log_msg("Warnings         : ", qc_tab[status == "WARN", .N])
log_msg("Fatal            : ", qc_tab[status == "FATAL", .N])
log_msg("Rows written     : ", format(nrow(res), big.mark = ","), " (283 x 50)")
log_msg("Elapsed          : ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
log_msg("NOTE             : no MK / Hamed-Rao / Sen / FDR / classification / mapping performed.")
log_msg("NOTE             : report as AMJ-restricted WSDI (01 Apr / 30 Jun truncation).")

.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
writeLines(.LOG, LOG_FILE)
cat("\nConsole log written to:", LOG_FILE, "\n")

invisible(res)
###############################################################################
## END OF SCRIPT
###############################################################################