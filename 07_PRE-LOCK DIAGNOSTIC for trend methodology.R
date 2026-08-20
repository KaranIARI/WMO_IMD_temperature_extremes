###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 07A : PRE-LOCK DIAGNOSTIC for trend methodology
##             (no final trend results are produced or published here)
##
##  283 grids | AMJ | indices: WSDI, TXx, TXn, TNx, TNn, DTR
##  Periods  : FULL 1976-2025 (n = 50)   TREND 1996-2025 (n = 30)
##
## ---------------------------------------------------------------------------
##  PURPOSE
##    This script exists to INFORM four methodological decisions that are NOT
##    yet locked. It deliberately computes every candidate and reports how the
##    answers differ, rather than assuming any of them:
##      D1  minimum valid-year rule        (sweep 50% .. 100%, not just 80%)
##      D2  autocorrelation structure      (detrended-rank ACF by lag)
##      D3  Hamed-Rao lag truncation       (none / lag 1 / lag <=3 / all lags)
##      D4  multiplicity control           (raw p vs BH-FDR at q = .05/.10/.20)
##    Nothing here should be quoted as a result. STEP 07 will fix one option
##    per decision on the evidence below and recompute cleanly.
##
##  INPUT  (the ONLY input, opened READ-ONLY, never modified)
##    F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet
##
##  OUTPUTS (all diagnostic; none is an analysis product)
##    1. AMJ_STEP07A_validity_by_grid.csv     grid x index x period completeness
##    2. AMJ_STEP07A_completeness_sweep.csv   eligible grids under each rule
##    3. AMJ_STEP07A_acf_summary.csv          ACF by lag, pooled over grids
##    4. AMJ_STEP07A_acf_by_grid.csv          ACF lags 1-10 per grid
##    5. AMJ_STEP07A_mk_variants.csv          per-series S, Var, CF, Z, p, Sen
##    6. AMJ_STEP07A_variant_comparison.csv   discordance between lag choices
##    7. AMJ_STEP07A_fdr_comparison.csv       raw vs BH counts
##    8. AMJ_STEP07A_wsdi_distribution.csv    tie / zero structure of WSDI
##    9. AMJ_STEP07A_QC.csv  +  AMJ_STEP07A_LOG.txt
##
## ---------------------------------------------------------------------------
##  MODIFIED MANN-KENDALL (Hamed & Rao 1998) AS IMPLEMENTED HERE
## ---------------------------------------------------------------------------
##   1. S    = sum_{i<j} sign(x_j - x_i) over non-missing pairs.
##   2. Var(S) with the TIE correction:
##        Var = [ n(n-1)(2n+5) - sum_g t_g(t_g-1)(2t_g+5) ] / 18
##      The tie term is mandatory here: WSDI is integer-valued and heavily
##      tied, and the untied formula would be simply wrong for it.
##   3. Theil-Sen slope beta = median over all pairs of (x_j-x_i)/(t_j-t_i).
##   4. Detrend:  x'_t = x_t - beta * t.
##   5. RANK the detrended series (ranks, not raw values - a common
##      implementation error), keeping NA in place.
##   6. Autocorrelation rho(i) of those ranks, computed PAIRWISE-COMPLETE at
##      true calendar lag, so a missing year does not silently shift the lag.
##   7. Screen rho(i) against the white-noise bound +/- z_{1-alpha/2} / sqrt(n),
##      as in modifiedmk::mmkh() and pymannkendall. Non-significant rho(i) are
##      set to zero. (Anderson (1942) bounds are retained in the code only as a
##      non-default branch for reproducing the earlier run.)
##
##   8. Eligible lags are 1 .. (n-1). The Hamed-Rao correction uses the
##      weights (n-i)(n-i-1)(n-i-2); therefore the terms at i = n-2 and
##      i = n-1 have zero weight and do not affect the correction factor.
##   9. Var*(S) = Var_ties(S) * CF.
##  10. Z = (S - sign(S)) / sqrt(Var*(S))   [continuity correction ON]
##      p = 2 * Phi(-|Z|).
##
##  CF < 1 (variance DEFLATION, arising when significant rho are negative) is
##  permitted here because it is faithful to the original method, but BOTH the
##  faithful and the floored (CF := max(CF,1)) p-values are reported so the
##  consequence of that policy can be seen before it is locked.
##
##  Hamed-Rao is applied UNIFORMLY to every series. No conditional
##  pre-test-then-choose procedure is used: selecting the estimator on the
##  basis of an autocorrelation test distorts the very type-I error it is
##  meant to protect.
##
##  CAVEAT recorded per series: the CF formula assumes evenly spaced data.
##  Series containing gaps are counted and flagged (`has_gaps`) so their
##  influence can be judged.
##
##  Non-interactive, deterministic, read-only. Runs under Rscript.
###############################################################################

## ---------------------------------------------------------------------------
## 0. ENVIRONMENT
## ---------------------------------------------------------------------------
rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE, scipen = 999)
options(error = NULL)
Sys.setenv(TZ = "UTC")
set.seed(20260814L)          # only used by the self-tests; nothing is resampled

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
## 1. PATHS, CANDIDATE SETTINGS (nothing here is locked)
## ---------------------------------------------------------------------------
IN_FILE  <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet"
OUT_DIR  <- "F:/WMO_IMD_R/WMO_IMD/data"

F_VALID  <- file.path(OUT_DIR, "AMJ_STEP07A_validity_by_grid.csv")
F_SWEEP  <- file.path(OUT_DIR, "AMJ_STEP07A_completeness_sweep.csv")
F_ACFS   <- file.path(OUT_DIR, "AMJ_STEP07A_acf_summary.csv")
F_ACFG   <- file.path(OUT_DIR, "AMJ_STEP07A_acf_by_grid.csv")
F_MK     <- file.path(OUT_DIR, "AMJ_STEP07A_mk_variants.csv")
F_VCMP   <- file.path(OUT_DIR, "AMJ_STEP07A_variant_comparison.csv")
F_FDR    <- file.path(OUT_DIR, "AMJ_STEP07A_fdr_comparison.csv")
F_WSDI   <- file.path(OUT_DIR, "AMJ_STEP07A_wsdi_distribution.csv")
F_QC     <- file.path(OUT_DIR, "AMJ_STEP07A_QC.csv")
F_LOG    <- file.path(OUT_DIR, "AMJ_STEP07A_LOG.txt")

INDICES  <- c("WSDI", "TXx", "TXn", "TNx", "TNn", "DTR")
USABLE   <- c(WSDI = "usable_WSDI", TXx = "usable_TXx", TXn = "usable_TXn",
              TNx = "usable_TNx",  TNn = "usable_TNn", DTR = "usable_DTR")
UNITS    <- c(WSDI = "days/decade", TXx = "degC/decade", TXn = "degC/decade",
              TNx = "degC/decade",  TNn = "degC/decade", DTR = "degC/decade")

PERIODS <- list(
  FULL_1976_2025  = c(1976L, 2025L),
  TREND_1996_2025 = c(1996L, 2025L))

EXP_GRIDS <- 283L
EXP_ROWS  <- 14150L
YR_START  <- 1976L; YR_END <- 2025L

## D1 candidates - completeness rules to sweep (NOT locked)
COMPLETENESS_RULES <- c(0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
RULES_FOR_FDR      <- c(0.70, 0.80, 0.90, 1.00)
MIN_N_ABSOLUTE     <- 10L      # below this a trend test is meaningless at all

## D3 candidates - lag truncation (NOT locked). NA = all eligible lags 1..(n-3)
LAG_VARIANTS <- list(MK_original = 0L, HR_lag1 = 1L, HR_lag3 = 3L,
                     HR_lag5 = 5L, HR_all = NA_integer_)

ACF_REPORT_LAGS <- 10L         # per-grid ACF written out to this lag

## D4 candidates - multiplicity control (NOT locked)
ALPHA_RAW <- 0.05
Q_LEVELS  <- c(0.05, 0.10, 0.20)

## ---------------------------------------------------------------------------
## 2. HELPERS AND LOGGING
## ---------------------------------------------------------------------------
.LOG <- character(0)
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
CHECKS <- list()
record_check <- function(name, passed, detail, severity = c("FATAL", "WARN")) {
  severity <- match.arg(severity)
  ok <- length(passed) == 1L && isTRUE(passed)
  status <- if (ok) "PASS" else severity
  CHECKS[[length(CHECKS) + 1L]] <<- data.table(check = name, status = status, detail = detail)
  log_msg(sprintf("  %-5s | %-44s | %s", status, name, detail))
  if (!ok && severity == "FATAL")
    stop("FATAL QC failure: ", name, " -- ", detail, call. = FALSE)
  invisible(ok)
}
sc <- function(v, f = "%.3f") { v <- v[is.finite(v)]; if (!length(v)) "NA" else sprintf(f, mean(v)) }
sc_med <- function(v, f = "%.3f") { v <- v[is.finite(v)]; if (!length(v)) "NA" else sprintf(f, stats::median(v)) }
rng <- function(v, f = "%.3f") { v <- v[is.finite(v)]
if (!length(v)) "all NA" else sprintf(paste0("[", f, ", ", f, "]"), min(v), max(v)) }


## ===========================================================================
## SECTION 2 (ENGINE) - REPLACE FROM "MANN-KENDALL / SEN / HAMED-RAO ENGINE"
## DOWN TO (AND INCLUDING) THE OLD mk_series() DEFINITION
## ===========================================================================

## Screening rule for retaining autocorrelations. "standard" reproduces
## modifiedmk::mmkh() and pymannkendall exactly. "anderson" is retained only so
## the earlier behaviour can be reproduced for the record; it is NOT the default
## and must not be used for the published result without justification.
HR_SCREEN <- "standard"        # one of: "standard", "anderson"
ACF_ALPHA <- 0.05              # two-sided level for the screening bound

## S and the Theil-Sen slope share one pairwise-difference matrix. Unchanged.
mk_S_sen <- function(x, t) {
  ok <- !is.na(x)
  xv <- x[ok]; tv <- t[ok]
  n <- length(xv)
  if (n < 2L) return(list(S = NA_real_, sen = NA_real_, n = n))
  dx <- outer(xv, xv, "-")            # dx[i,j] = x_i - x_j
  dt <- outer(tv, tv, "-")
  lt <- lower.tri(dx)                 # i > j  =>  t_i > t_j
  list(S   = as.numeric(sum(sign(dx[lt]))),
       sen = as.numeric(stats::median(dx[lt] / dt[lt])),
       n   = n)
}

## Tie-corrected variance of S. Identical to mmkh's construction. Unchanged.
mk_var_ties <- function(x) {
  xv <- x[!is.na(x)]
  n <- length(xv)
  if (n < 3L) return(NA_real_)
  v <- n * (n - 1) * (2 * n + 5)
  tg <- as.numeric(table(xv)); tg <- tg[tg > 1]
  if (length(tg)) v <- v - sum(tg * (tg - 1) * (2 * tg + 5))
  v / 18
}

## [BUG 5] Pairwise-complete rank autocorrelation at TRUE calendar lag.
## Common mean and common denominator, matching stats::acf()'s normalisation,
## so a complete series reproduces acf(rank(xd))$acf[-1] exactly.
## Undefined lags return 0, not NA.
acf_hr <- function(r, max_lag) {
  n <- length(r)
  out <- rep(0, max(max_lag, 0L))
  if (max_lag < 1L) return(out)
  mu <- mean(r, na.rm = TRUE)
  den <- sum((r - mu)^2, na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(out)
  for (k in seq_len(min(max_lag, n - 1L))) {
    a <- r[1:(n - k)]; b <- r[(k + 1):n]
    ok <- !is.na(a) & !is.na(b)
    if (any(ok)) out[k] <- sum((a[ok] - mu) * (b[ok] - mu)) / den
  }
  out[!is.finite(out)] <- 0
  out
}

## [BUG 1] Screening bound. "standard" == modifiedmk / pymannkendall.
hr_significant_rho <- function(rho, n, alpha = ACF_ALPHA, rule = HR_SCREEN) {
  k <- seq_along(rho)
  if (identical(rule, "standard")) {
    b <- stats::qnorm(1 - alpha / 2) / sqrt(n)
    return(is.finite(rho) & abs(rho) > b)
  }
  ## legacy Anderson (1942), kept only for reproducing the earlier run
  z <- stats::qnorm(1 - alpha / 2); d <- n - k
  lo <- (-1 - z * sqrt(pmax(d - 1, 0))) / d
  hi <- (-1 + z * sqrt(pmax(d - 1, 0))) / d
  is.finite(rho) & d > 1 & (rho < lo | rho > hi)
}

## [BUG 3] Correction factor, summing i = 1..max_lag as mmkh does.
## The weight is zero at i = n-2 and n-1, so extending to n-1 is harmless and
## matches the reference exactly.
hr_cf <- function(rho, sig, n, max_lag) {
  if (max_lag < 1L) return(1)
  k <- seq_len(max_lag)
  r <- rho[k]; r[!sig[k]] <- 0; r[!is.finite(r)] <- 0
  w <- (n - k) * (n - k - 1) * (n - k - 2)
  w[w < 0] <- 0
  1 + (2 / (n * (n - 1) * (n - 2))) * sum(w * r)
}

## [BUG 2] Full analysis of one series.
##   x : full CALENDAR-length vector over the period, NA where the year is absent
##   t : 1..n_series
##   n_series drives the ACF, the screening bound and the CF weights
##   n_obs    drives S and Var(S)
## For a complete series n_series == n_obs and this reduces exactly to mmkh().
mk_series <- function(x, t, lag_variants = LAG_VARIANTS, alpha = ACF_ALPHA) {
  n_series <- length(x)
  base <- mk_S_sen(x, t)
  n_obs <- base$n
  varS  <- mk_var_ties(x)
  empty_acf <- rep(NA_real_, ACF_REPORT_LAGS)
  
  if (is.na(base$S) || is.na(varS) || n_obs < MIN_N_ABSOLUTE || n_series < 4L)
    return(list(ok = FALSE, n_obs = n_obs, n_series = n_series,
                S = base$S, varS = varS, sen = base$sen,
                acf = empty_acf, rows = NULL))
  
  ## Sen-detrend on the calendar index, then RANK the detrended values
  xd <- x - base$sen * t
  rk <- rep(NA_real_, n_series)
  rk[!is.na(xd)] <- rank(xd[!is.na(xd)], ties.method = "average")
  
  max_lag_all <- n_series - 1L                  # mmkh convention
  rho <- acf_hr(rk, max_lag_all)
  sig <- hr_significant_rho(rho, n_series, alpha, HR_SCREEN)
  
  out_acf <- empty_acf
  nrep <- min(ACF_REPORT_LAGS, max_lag_all)
  out_acf[seq_len(nrep)] <- rho[seq_len(nrep)]
  
  rows <- vector("list", length(lag_variants))
  for (v in seq_along(lag_variants)) {
    nm <- names(lag_variants)[v]; L <- lag_variants[[v]]
    if (is.na(L)) L <- max_lag_all
    L <- min(L, max_lag_all)
    
    if (L < 1L) {
      cf <- 1; nsig <- 0L; ridx <- integer(0)
    } else {
      cf   <- hr_cf(rho, sig, n_series, L)
      ridx <- which(sig[seq_len(L)])
      nsig <- length(ridx)
    }
    
    cf_status <- if (!is.finite(cf)) "non_finite" else
      if (cf <= 0)        "non_positive" else
        if (cf < 1)         "deflated" else "ok"
    
    zp <- function(cfx) {
      vx <- varS * cfx
      if (!is.finite(vx) || vx <= 0) return(c(NA_real_, NA_real_))
      z <- if (base$S > 0) (base$S - 1) / sqrt(vx) else
        if (base$S < 0) (base$S + 1) / sqrt(vx) else 0
      c(z, 2 * stats::pnorm(-abs(z)))
    }
    a  <- zp(cf)
    bf <- zp(max(cf, 1))                       # sensitivity only, never primary
    
    ## diagnostics that explain a non-positive CF
    neg <- if (length(ridx)) ridx[rho[ridx] < 0] else integer(0)
    wneg <- if (length(neg))
      sum((n_series - neg) * (n_series - neg - 1) * (n_series - neg - 2) * rho[neg]) *
      (2 / (n_series * (n_series - 1) * (n_series - 2))) else 0
    
    rows[[v]] <- data.table(
      variant = nm, lags_used = L, n_lags_sig = as.integer(nsig),
      CF = cf, cf_status = cf_status, CF_floored = max(cf, 1),
      varS_corrected = varS * cf,
      Z = a[1], p = a[2], Z_cf_floored = bf[1], p_cf_floored = bf[2],
      rho_lag1 = if (max_lag_all >= 1L) rho[1] else NA_real_,
      rho_lag1_significant = if (max_lag_all >= 1L) sig[1] else NA,
      min_sig_rho = if (nsig) min(rho[ridx]) else NA_real_,
      lag_of_min_sig_rho = if (nsig) ridx[which.min(rho[ridx])] else NA_integer_,
      n_sig_rho_negative = length(neg),
      cf_negative_contribution = wneg,
      screen_rule = HR_SCREEN,
      screen_bound = stats::qnorm(1 - alpha / 2) / sqrt(n_series))
  }
  
  list(ok = TRUE, n_obs = n_obs, n_series = n_series,
       S = base$S, varS = varS, sen = base$sen,
       acf = out_acf, rho_all = rho, sig = sig, rows = rbindlist(rows))
}

t_start <- Sys.time()
log_head("STEP 07A | PRE-LOCK TREND DIAGNOSTIC (no final results produced)")
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("Input (readonly) : ", IN_FILE)
log_msg("Indices          : ", paste(INDICES, collapse = ", "))
log_msg("Periods          : ", paste(names(PERIODS), collapse = " ; "))
log_msg("D1 sweep         : completeness rules ", paste(COMPLETENESS_RULES * 100, collapse = "/"), "%")
log_msg("D3 sweep         : ", paste(names(LAG_VARIANTS), collapse = ", "),
        "  (HR_all = lags 1..n-1)")
log_msg("D4 sweep         : raw p < ", ALPHA_RAW, " vs BH q <= ",
        paste(Q_LEVELS, collapse = "/"))
log_msg("NOTE             : nothing in this script is locked; STEP 07 decides.")

## ---------------------------------------------------------------------------
## 3. READ STEP 06 (READ-ONLY) AND VALIDATE
## ---------------------------------------------------------------------------
log_head("3. READ STEP 06 OUTPUT (READ-ONLY)")

if (!file.exists(IN_FILE)) stop("STEP 06 parquet not found: ", IN_FILE, call. = FALSE)
.mtime_before <- file.mtime(IN_FILE)

d6 <- as.data.table(arrow::read_parquet(IN_FILE))
log_msg("Rows read        : ", format(nrow(d6), big.mark = ","), " x ", ncol(d6), " cols")

need <- c("grid_id", "year", "lat", "lon", INDICES, unname(USABLE))
miss <- setdiff(need, names(d6))
if (length(miss))
  stop("SCHEMA ERROR in STEP 06 parquet.\n  Missing: ", paste(miss, collapse = ", "),
       "\n  Present: ", paste(names(d6), collapse = ", "), call. = FALSE)
record_check("schema_step06", TRUE, sprintf("all %d required columns present", length(need)))

d6[, `:=`(grid_id = as.character(grid_id), year = as.integer(year))]
setorder(d6, grid_id, year)

record_check("input_row_count", nrow(d6) == EXP_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(d6), big.mark = ","),
                     format(EXP_ROWS, big.mark = ",")))
record_check("input_grid_count", uniqueN(d6$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(d6$grid_id), EXP_GRIDS))
record_check("input_year_coverage", identical(sort(unique(d6$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(d6$year), min(d6$year), max(d6$year)))
record_check("input_no_duplicate_grid_year", !anyDuplicated(d6, by = c("grid_id", "year")),
             sprintf("%d duplicates", sum(duplicated(d6, by = c("grid_id", "year")))))
record_check("usable_flags_match_NA_pattern",
             all(vapply(INDICES, function(ix)
               all(d6[[USABLE[[ix]]]] == !is.na(d6[[ix]])), logical(1))),
             "usable_* flags agree with the NA pattern of every index")

GRIDS <- sort(unique(d6$grid_id))
meta  <- unique(d6[, .(grid_id, lat, lon)])


## ===========================================================================
## SECTION 4 - ENGINE SELF-TESTS  (replaces the old Section 4 entirely)
## ===========================================================================
# ===========================================================================
log_head("4. ENGINE SELF-TESTS (must pass before any diagnostic is computed)")

TT <- list()
addt <- function(nm, got, want, tol = 1e-9, note = "") {
  ok <- is.finite(got) && is.finite(want) && abs(got - want) <= tol
  TT[[length(TT) + 1L]] <<- data.table(test = nm, got = sprintf("%.10g", got),
                                       want = sprintf("%.10g", want), ok = ok, note = note)
}

## Warning-free, deterministic AR(1) generator (replaces stats::arima.sim).
## phi = 0 is safe here and yields plain white noise.
sim_ar1 <- function(n, phi, seed, burn = 200L) {
  set.seed(seed)
  e <- stats::rnorm(n + burn)
  as.numeric(stats::filter(e, filter = phi, method = "recursive"))[(burn + 1L):(burn + n)]
}

## --- S, Sen, Var(S) ---------------------------------------------------------
set.seed(1001L)
xr <- round(stats::rnorm(25, 30, 5), 1); tr <- seq_along(xr)
bf <- 0; for (i in 1:24) for (j in (i + 1):25) bf <- bf + sign(xr[j] - xr[i])
addt("S_matches_brute_force", mk_S_sen(xr, tr)$S, bf, 0, "vectorised S == double loop")
addt("S_perfect_increase", mk_S_sen(1:10, 1:10)$S, 45, 0, "n(n-1)/2 for n=10")
addt("S_perfect_decrease", mk_S_sen(10:1, 1:10)$S, -45, 0, "sign symmetry")
addt("S_constant_series",  mk_S_sen(rep(5, 10), 1:10)$S, 0, 0, "all ties -> S = 0")
addt("sen_exact_linear",   mk_S_sen(2 * (1:20) + 7, 1:20)$sen, 2, 1e-12, "y = 2t + 7")
addt("sen_exact_negative", mk_S_sen(-1.5 * (1:20), 1:20)$sen, -1.5, 1e-12, "y = -1.5t")

xt <- c(1, 1, 1, 2, 2, 3, 4, 5, 5, 5); nt <- 10
vt <- (nt * (nt - 1) * (2 * nt + 5) - (3 * 2 * 11 + 2 * 1 * 9 + 3 * 2 * 11)) / 18
addt("var_ties_formula", mk_var_ties(xt), vt, 1e-12, "3 tie groups: 3,2,3")
addt("var_no_ties", mk_var_ties(1:10), 10 * 9 * 25 / 18, 1e-12, "untied reduces correctly")

## --- ACF reproduces stats::acf() on a complete series -----------------------
set.seed(1002L)
rr <- rank(stats::rnorm(40))
addt("acf_hr_matches_stats_acf",
     max(abs(acf_hr(rr, 10) -
               as.numeric(stats::acf(rr, lag.max = 10, plot = FALSE)$acf)[-1])),
     0, 1e-12, "same normalisation as stats::acf on complete data")

## AR(1) detection: 1000 points gives a tight lag-1 estimate. No arima.sim.
ar1 <- sim_ar1(1000L, 0.7, seed = 101L)
addt("acf_detects_AR1", acf_hr(rank(ar1), 1)[1], 0.7, 0.10,
     "lag-1 rank ACF near the true AR coefficient (n = 1000, warning-free)")

ar1s <- ar1[1:100]
ag <- ar1s; ag[c(20, 55, 80)] <- NA
rg <- rep(NA_real_, 100); rg[!is.na(ag)] <- rank(ag[!is.na(ag)])
addt("acf_pairwise_handles_gaps",
     sign(acf_hr(rg, 1)[1]), sign(acf_hr(rank(ar1s), 1)[1]), 0,
     "pairwise-complete ACF keeps the true calendar lag")

## --- screening bound is the standard one ------------------------------------
xb <- sim_ar1(30L, 0.3, seed = 7L) + 0.05 * (1:30)
addt("screen_bound_is_z_over_sqrt_n",
     stats::qnorm(0.975) / sqrt(30),
     unique(mk_series(xb, 1:30, list(HR_all = NA_integer_))$rows$screen_bound),
     1e-12, "matches mmkh sig <- qnorm((1+ci)/2)/sqrt(n)")

## --- 4a  DEGENERATE-EXACT: CF == 1 with no randomness whatsoever -------------
## x = 1:30 has Sen slope exactly 1, so xd = x - 1*t is identically 0. All ranks
## tie, the ACF denominator is 0, acf_hr returns zeros, nothing can be retained,
## and CF is 1 exactly. Var(S) is untied and positive, so Z is well defined.
lin <- as.numeric(1:30)
r_lin <- mk_series(lin, 1:30, list(MK_original = 0L, HR_all = NA_integer_))
addt("CF_equals_one_exact_linear_series",
     r_lin$rows[variant == "HR_all", CF], 1, 0,
     sprintf("detrended residuals identically zero; n_lags_sig = %d",
             r_lin$rows[variant == "HR_all", n_lags_sig]))
addt("HR_reduces_to_MK_exact_linear_series",
     r_lin$rows[variant == "HR_all", Z], r_lin$rows[variant == "MK_original", Z], 0,
     "identical Z when CF is exactly 1 (not a self-comparison)")

## --- 4b  NON-DEGENERATE END-TO-END: deterministic search for a clean series --
## Fixed seeds in fixed order, so the same series is found on every run and
## every machine. P(a given seed retains zero of 29 lags) = 0.95^29 = 0.226,
## so a hit is expected within the first handful of seeds.
find_clean_series <- function(n = 30L, slope = 0.5, max_try = 500L) {
  for (s in seq_len(max_try)) {
    xx <- slope * (1:n) + sim_ar1(n, 0, seed = 5000L + s)   # phi = 0: white noise
    rr <- mk_series(xx, 1:n, list(MK_original = 0L, HR_all = NA_integer_))
    if (isTRUE(rr$ok) && rr$rows[variant == "HR_all", n_lags_sig] == 0L)
      return(list(seed = 5000L + s, try = s, x = xx, res = rr))
  }
  NULL
}
clean <- find_clean_series()
if (is.null(clean)) {
  addt("CF_equals_one_when_no_rho_retained", NA_real_, 1, 0,
       "no clean seed found in 500 tries - investigate the screening rule")
  addt("HR_reduces_to_MK_when_CF_is_one", NA_real_, 1, 0, "precondition unavailable")
} else {
  addt("CF_equals_one_when_no_rho_retained",
       clean$res$rows[variant == "HR_all", CF], 1, 1e-12,
       sprintf("seed %d (try %d), non-degenerate series, n_lags_sig = %d",
               clean$seed, clean$try, clean$res$rows[variant == "HR_all", n_lags_sig]))
  addt("HR_reduces_to_MK_when_CF_is_one",
       clean$res$rows[variant == "HR_all", Z],
       clean$res$rows[variant == "MK_original", Z], 1e-12,
       "HR collapses to uncorrected MK when nothing is retained")
}

## Guard the guard: the search must have verified its own precondition.
addt("clean_series_precondition_verified",
     if (is.null(clean)) NA_real_ else
       as.numeric(clean$res$rows[variant == "HR_all", n_lags_sig] == 0L),
     1, 0, "the tested series provably retains zero lags")

## --- CF formula -------------------------------------------------------------
rho_t <- c(0.4, 0.2, 0.1); n_t <- 20; k_t <- 1:3
cf_t <- 1 + (2 / (n_t * (n_t - 1) * (n_t - 2))) *
  sum((n_t - k_t) * (n_t - k_t - 1) * (n_t - k_t - 2) * rho_t)
addt("hr_cf_formula", hr_cf(rho_t, rep(TRUE, 3), n_t, 3L), cf_t, 1e-12,
     "weights (n-i)(n-i-1)(n-i-2)")
addt("hr_cf_ignores_unscreened_rho", hr_cf(rho_t, rep(FALSE, 3), n_t, 3L), 1, 1e-12,
     "zeroed rho contribute nothing")
addt("lag1_weight_is_2n_minus_3_over_n",
     hr_cf(c(-1, rep(0, 28)), c(TRUE, rep(FALSE, 28)), 30L, 29L),
     1 - 2 * (30 - 3) / 30, 1e-12,
     "a single rho(1) = -1 drives CF to 1 - 1.80 = -0.80 at n = 30")

## --- DECISIVE: exact reproduction of modifiedmk::mmkh() ---------------------
## Deterministic phi sequence and warning-free generator; phi = 0 is now safe.
if (requireNamespace("modifiedmk", quietly = TRUE)) {
  phis <- rep(c(-0.6, -0.3, 0.0, 0.3, 0.7), each = 8L)
  dmax <- 0; ncmp <- 0L
  for (i in seq_along(phis)) {
    xs <- sim_ar1(30L, phis[i], seed = 90000L + i) + 0.03 * (1:30)
    mine <- mk_series(xs, 1:30, list(HR_all = NA_integer_))$rows[1]
    ref  <- tryCatch(modifiedmk::mmkh(xs), error = function(e) NULL)
    if (is.null(ref)) next
    pref <- unname(ref["new P-value"])
    ## mmkh does not guard CF <= 0 either; compare only where both are defined
    if (is.finite(pref) && is.finite(mine$p)) {
      dmax <- max(dmax, abs(mine$p - pref)); ncmp <- ncmp + 1L
    }
  }
  addt("reproduces_modifiedmk_mmkh", dmax, 0, 1e-8,
       sprintf("identical p over %d complete AR(1) series, phi in [-0.6, 0.7]", ncmp))
} else {
  log_msg("  NOTE  | modifiedmk not installed; the decisive cross-check is skipped.")
  log_msg("        | install.packages('modifiedmk') before the publication run.")
}

tt <- rbindlist(TT)
for (i in seq_len(nrow(tt)))
  log_msg(sprintf("  %-4s | %-40s | got %-14s want %-14s | %s",
                  if (tt$ok[i]) "PASS" else "FAIL", tt$test[i], tt$got[i], tt$want[i], tt$note[i]))
record_check("engine_self_tests", all(tt$ok),
             sprintf("%d of %d passed", sum(tt$ok), nrow(tt)))
record_check("crosscheck_modifiedmk_available",
             "reproduces_modifiedmk_mmkh" %in% tt$test,
             if ("reproduces_modifiedmk_mmkh" %in% tt$test)
               sprintf("max |p difference| vs mmkh = %s",
                       tt[test == "reproduces_modifiedmk_mmkh", got])
             else "modifiedmk not installed; install before the publication run",
             severity = "WARN")
record_check("screening_rule_documented", identical(HR_SCREEN, "standard"),
             sprintf("HR_SCREEN = '%s' (standard = z/sqrt(n), as in mmkh and pymannkendall)",
                     HR_SCREEN), severity = "WARN")



## ===========================================================================
## SECTION 5 - D1 VALIDITY AND COMPLETENESS  (unchanged; included for drop-in)
## ===========================================================================
log_head("5. D1 | VALID-YEAR COUNTS AND COMPLETENESS SWEEP")

valid <- rbindlist(lapply(names(PERIODS), function(pn) {
  yr <- PERIODS[[pn]]
  sub <- d6[year >= yr[1] & year <= yr[2]]
  n_period <- yr[2] - yr[1] + 1L
  rbindlist(lapply(INDICES, function(ix) {
    s <- sub[, .(n_period_years = n_period,
                 n_valid = sum(!is.na(get(ix))),
                 first_valid_year = suppressWarnings(min(year[!is.na(get(ix))])),
                 last_valid_year  = suppressWarnings(max(year[!is.na(get(ix))])),
                 n_gaps = { yy <- sort(year[!is.na(get(ix))])
                 if (length(yy) < 2L) 0L else sum(diff(yy) != 1L) }),
             by = grid_id]
    s[, `:=`(period = pn, index = ix)]
    s[!is.finite(first_valid_year), first_valid_year := NA_integer_]
    s[!is.finite(last_valid_year),  last_valid_year  := NA_integer_]
    s[, completeness_pct := round(100 * n_valid / n_period_years, 2)]
    s[, has_gaps := n_gaps > 0L]
    s[]
  }))
}))
valid <- merge(valid, meta, by = "grid_id", all.x = TRUE, sort = FALSE)
for (r in COMPLETENESS_RULES)
  valid[, (sprintf("eligible_%02d", round(r * 100))) :=
          n_valid >= ceiling(r * n_period_years) & n_valid >= MIN_N_ABSOLUTE]
setcolorder(valid, c("period", "index", "grid_id", "lat", "lon",
                     "n_period_years", "n_valid", "completeness_pct",
                     "n_gaps", "has_gaps", "first_valid_year", "last_valid_year"))
setorder(valid, period, index, grid_id)

record_check("validity_row_count",
             nrow(valid) == EXP_GRIDS * length(INDICES) * length(PERIODS),
             sprintf("%s rows = %d grids x %d indices x %d periods",
                     format(nrow(valid), big.mark = ","), EXP_GRIDS,
                     length(INDICES), length(PERIODS)))

sweep <- rbindlist(lapply(COMPLETENESS_RULES, function(r) {
  cn <- sprintf("eligible_%02d", round(r * 100))
  valid[, .(rule_frac = r,
            min_years_required = ceiling(r * n_period_years[1]),
            n_grids_eligible = sum(get(cn)),
            n_grids_excluded = sum(!get(cn)),
            pct_eligible = round(100 * mean(get(cn)), 2),
            n_with_gaps_among_eligible = sum(get(cn) & has_gaps)),
        by = .(period, index)]
}))
setorder(sweep, period, index, rule_frac)

log_msg("Eligible grids of ", EXP_GRIDS, " by completeness rule:")
for (pn in names(PERIODS)) for (ix in INDICES) {
  z <- sweep[period == pn & index == ix]
  log_msg(sprintf("  %-16s %-5s : %s", pn, ix,
                  paste(sprintf("%d%%=%d", round(z$rule_frac * 100), z$n_grids_eligible),
                        collapse = "  ")))
}
log_msg("Grid-index-period series containing internal gaps : ",
        valid[has_gaps == TRUE, .N], " of ", nrow(valid))



## ===========================================================================
## SECTION 6 - D2 + D3 AUTOCORRELATION AND HAMED-RAO VARIANTS
## (replaces the old Section 6 entirely)
## ===========================================================================
log_head("6. D2 + D3 | AUTOCORRELATION AND HAMED-RAO VARIANTS")

MKROWS <- vector("list", length(PERIODS) * length(INDICES) * EXP_GRIDS)
ACFROW <- vector("list", length(PERIODS) * length(INDICES) * EXP_GRIDS)
kk <- 0L; t0 <- Sys.time()

for (pn in names(PERIODS)) {
  yr  <- PERIODS[[pn]]
  yrs <- yr[1]:yr[2]
  sub <- d6[year >= yr[1] & year <= yr[2]]
  setkey(sub, grid_id, year)
  for (ix in INDICES) {
    w <- dcast(sub[, .(grid_id, year, v = get(ix))], grid_id ~ year, value.var = "v")
    setkey(w, grid_id); w <- w[J(GRIDS)]
    M <- as.matrix(w[, -1L, with = FALSE])
    stopifnot(ncol(M) == length(yrs))          # calendar completeness of the frame
    tvec <- as.numeric(seq_along(yrs))
    for (g in seq_len(nrow(M))) {
      kk <- kk + 1L
      xg <- as.numeric(M[g, ])
      r  <- mk_series(xg, tvec)
      ACFROW[[kk]] <- data.table(period = pn, index = ix, grid_id = GRIDS[g],
                                 n_valid = r$n_obs, n_series = r$n_series,
                                 lag = seq_len(ACF_REPORT_LAGS), rho = r$acf)
      if (!isTRUE(r$ok)) next
      rr <- copy(r$rows)
      rr[, `:=`(period = pn, index = ix, grid_id = GRIDS[g],
                n_valid = r$n_obs, n_series = r$n_series,
                S = r$S, varS_uncorrected = r$varS,
                sen_per_year = r$sen, sen_per_decade = r$sen * 10)]
      MKROWS[[kk]] <- rr
    }
  }
  log_msg("Processed period ", pn, " in ",
          sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

mk    <- rbindlist(MKROWS, use.names = TRUE, fill = TRUE)
acf_g <- rbindlist(ACFROW, use.names = TRUE, fill = TRUE)
rm(MKROWS, ACFROW); invisible(gc(verbose = FALSE))

mk <- merge(mk, valid[, .(period, index, grid_id, completeness_pct, n_gaps, has_gaps,
                          n_period_years)],
            by = c("period", "index", "grid_id"), all.x = TRUE, sort = FALSE)
mk <- merge(mk, meta, by = "grid_id", all.x = TRUE, sort = FALSE)
mk[, units := UNITS[index]]
setcolorder(mk, c("period", "index", "grid_id", "lat", "lon", "variant",
                  "n_valid", "n_series", "n_period_years", "completeness_pct",
                  "n_gaps", "has_gaps",
                  "S", "varS_uncorrected", "lags_used", "n_lags_sig",
                  "CF", "cf_status", "CF_floored", "varS_corrected",
                  "Z", "p", "Z_cf_floored", "p_cf_floored",
                  "rho_lag1", "rho_lag1_significant", "min_sig_rho",
                  "lag_of_min_sig_rho", "n_sig_rho_negative",
                  "cf_negative_contribution", "screen_rule", "screen_bound",
                  "sen_per_year", "sen_per_decade", "units"))
setorder(mk, period, index, grid_id, variant)

## --- QC on the engine output -------------------------------------------------
record_check("mk_p_in_unit_interval",
             mk[!is.na(p) & (p < 0 | p > 1), .N] == 0L,
             sprintf("%d p-values outside [0,1]; range %s",
                     mk[!is.na(p) & (p < 0 | p > 1), .N], rng(mk$p)))

## [RETARGETED] Non-positive corrected variance is permitted ONLY where the CF
## itself is non-positive, which is a documented property of Hamed-Rao. Any
## other route to a non-positive variance is a bug and still halts the run.
n_np   <- mk[cf_status == "non_positive", .N]
n_unex <- mk[!is.na(varS_corrected) & varS_corrected <= 0 & cf_status != "non_positive", .N]
record_check("mk_variance_non_positive_only_via_CF", n_unex == 0L,
             sprintf("%d unexplained non-positive variances (must be 0); %d explained by CF <= 0",
                     n_unex, n_np))
record_check("mk_p_is_NA_exactly_where_CF_non_positive",
             mk[cf_status == "non_positive" & !is.na(p), .N] == 0L &&
               mk[cf_status != "non_positive" & is.finite(varS_corrected) & is.na(p), .N] == 0L,
             sprintf("%d CF<=0 rows carry a p-value; %d finite-variance rows lack one",
                     mk[cf_status == "non_positive" & !is.na(p), .N],
                     mk[cf_status != "non_positive" & is.finite(varS_corrected) & is.na(p), .N]))
record_check("mk_cf_non_positive_reported", TRUE,
             sprintf("CF <= 0 in %d of %s rows (%.3f%%); no flooring applied",
                     n_np, format(nrow(mk), big.mark = ","),
                     100 * n_np / max(nrow(mk), 1L)),
             severity = "WARN")
record_check("mk_original_has_CF_one",
             mk[variant == "MK_original", all(abs(CF - 1) < 1e-12)],
             "the uncorrected variant carries CF = 1 exactly")
record_check("mk_variant_coverage", mk[, uniqueN(variant)] == length(LAG_VARIANTS),
             sprintf("%d variants: %s", mk[, uniqueN(variant)],
                     paste(sort(unique(mk$variant)), collapse = ", ")))
record_check("mk_n_series_equals_period_length",
             mk[n_series != n_period_years, .N] == 0L,
             sprintf("%d rows where the ACF time base differs from the period length",
                     mk[n_series != n_period_years, .N]))

log_msg("MK rows computed : ", format(nrow(mk), big.mark = ","), " in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("CF status counts : ",
        paste(sprintf("%s=%d", names(table(mk$cf_status)), as.integer(table(mk$cf_status))),
              collapse = "  "))

## --- where do the CF <= 0 cases live? ----------------------------------------
if (n_np > 0L) {
  npx <- mk[cf_status == "non_positive"]
  log_msg("CF <= 0 breakdown (index / period / variant : count):")
  for (i in seq_len(nrow(npx[, .N, by = .(index, period, variant)])))
  { z <- npx[, .N, by = .(index, period, variant)][i]
  log_msg(sprintf("    %-5s %-16s %-11s : %d", z$index, z$period, z$variant, z$N)) }
  log_msg(sprintf("    driving rho(1): median %.3f, min %.3f ; median n_sig_rho_negative %.1f",
                  stats::median(npx$rho_lag1, na.rm = TRUE),
                  min(npx$rho_lag1, na.rm = TRUE),
                  stats::median(npx$n_sig_rho_negative, na.rm = TRUE)))
  fwrite(npx, file.path(OUT_DIR, "AMJ_STEP07A_cf_non_positive_cases.csv"))
  log_msg("Saved            : AMJ_STEP07A_cf_non_positive_cases.csv (", nrow(npx), " rows)")
}

## --- ACF summary pooled over grids -------------------------------------------
acf_g <- acf_g[!is.na(rho)]
acf_sum <- acf_g[, .(n_grids = .N,
                     mean_rho    = mean(rho),
                     median_rho  = stats::median(rho),
                     q10_rho     = stats::quantile(rho, 0.10, names = FALSE),
                     q90_rho     = stats::quantile(rho, 0.90, names = FALSE),
                     max_abs_rho = max(abs(rho)),
                     pct_beyond_bound = round(100 * mean(
                       abs(rho) > stats::qnorm(1 - ACF_ALPHA / 2) / sqrt(n_series)), 2)),
                 by = .(period, index, lag)]
sig_by_lag <- mk[variant == "HR_all",
                 .(mean_n_lags_sig = mean(n_lags_sig),
                   pct_any_sig_lag = round(100 * mean(n_lags_sig > 0L), 2)),
                 by = .(period, index)]
acf_sum <- merge(acf_sum, sig_by_lag, by = c("period", "index"), all.x = TRUE)
setorder(acf_sum, period, index, lag)

log_msg("Lag-1 detrended-rank ACF (mean / median / max|rho| / %% beyond bound):")
for (pn in names(PERIODS)) for (ix in INDICES) {
  z <- acf_sum[period == pn & index == ix & lag == 1L]
  if (nrow(z))
    log_msg(sprintf("  %-16s %-5s : %6.3f / %6.3f / %5.3f / %5.1f%%   any sig lag: %5.1f%%",
                    pn, ix, z$mean_rho, z$median_rho, z$max_abs_rho,
                    z$pct_beyond_bound, z$pct_any_sig_lag))
}

log_msg("Grids where a SIGNIFICANT rho occurs beyond lag 3 (HR_all vs HR_lag3):")
beyond3 <- merge(mk[variant == "HR_lag3", .(period, index, grid_id, s3 = n_lags_sig)],
                 mk[variant == "HR_all",  .(period, index, grid_id, sa = n_lags_sig)],
                 by = c("period", "index", "grid_id"))
b3 <- beyond3[, .(n_grids_extra_sig_lags = sum(sa > s3),
                  pct = round(100 * mean(sa > s3), 2)), by = .(period, index)]
for (i in seq_len(nrow(b3)))
  log_msg(sprintf("  %-16s %-5s : %3d grids (%.1f%%)", b3$period[i], b3$index[i],
                  b3$n_grids_extra_sig_lags[i], b3$pct[i]))

## ---------------------------------------------------------------------------
## 7. D3 - VARIANT DISCORDANCE
## ---------------------------------------------------------------------------

log_head("7. D3 | HOW MANY GRIDS CHANGE WITH THE LAG CHOICE")

wide <- dcast(mk[, .(period, index, grid_id, variant, p, CF, Z)],
              period + index + grid_id ~ variant, value.var = c("p", "CF", "Z"))

vcmp <- rbindlist(lapply(names(LAG_VARIANTS), function(a) {
  rbindlist(lapply(names(LAG_VARIANTS), function(b) {
    if (a >= b) return(NULL)
    pa <- paste0("p_", a); pb <- paste0("p_", b)
    if (!all(c(pa, pb) %in% names(wide))) return(NULL)
    wide[, .(variant_a = a, variant_b = b,
             n_grids = sum(!is.na(get(pa)) & !is.na(get(pb))),
             n_sig_a = sum(get(pa) < ALPHA_RAW, na.rm = TRUE),
             n_sig_b = sum(get(pb) < ALPHA_RAW, na.rm = TRUE),
             n_discordant = sum((get(pa) < ALPHA_RAW) != (get(pb) < ALPHA_RAW), na.rm = TRUE),
             median_abs_dp = stats::median(abs(get(pa) - get(pb)), na.rm = TRUE),
             max_abs_dp = suppressWarnings(max(abs(get(pa) - get(pb)), na.rm = TRUE))),
         by = .(period, index)]
  }))
}))
vcmp[!is.finite(max_abs_dp), max_abs_dp := NA_real_]
vcmp[, pct_discordant := round(100 * n_discordant / pmax(n_grids, 1L), 2)]
setorder(vcmp, period, index, variant_a, variant_b)

log_msg("Raw-significance discordance between lag choices (grids of 283):")
for (i in seq_len(nrow(vcmp[variant_a == "HR_lag3" | variant_b == "HR_lag3"]))) {
  z <- vcmp[variant_a == "HR_lag3" | variant_b == "HR_lag3"][i]
  log_msg(sprintf("  %-16s %-5s | %-11s vs %-11s : sig %3d vs %3d, discordant %3d (%.1f%%)",
                  z$period, z$index, z$variant_a, z$variant_b,
                  z$n_sig_a, z$n_sig_b, z$n_discordant, z$pct_discordant))
}

cf_stats <- mk[variant != "MK_original",
               .(n = .N,
                 mean_CF = mean(CF, na.rm = TRUE),
                 median_CF = stats::median(CF, na.rm = TRUE),
                 min_CF = min(CF, na.rm = TRUE),
                 max_CF = max(CF, na.rm = TRUE),
                 n_deflated = sum(cf_status == "deflated"),
                 n_non_positive = sum(cf_status == "non_positive"),
                 pct_below_1 = round(100 * mean(CF < 1, na.rm = TRUE), 2),
                 n_sig_change_if_floored =
                   sum((p < ALPHA_RAW) != (p_cf_floored < ALPHA_RAW), na.rm = TRUE)),
               by = .(period, index, variant)]
setorder(cf_stats, period, index, variant)

log_msg("CF summary (variance deflation policy evidence):")
for (i in seq_len(nrow(cf_stats))) {
  z <- cf_stats[i]
  log_msg(sprintf(
    "  %-16s %-5s %-11s : CF med %.3f range [%.3f, %.3f], CF<1 in %3d (%.1f%%), sig changes if floored: %d",
    z$period, z$index, z$variant,
    z$median_CF, z$min_CF, z$max_CF,
    z$n_deflated + z$n_non_positive,
    z$pct_below_1,
    z$n_sig_change_if_floored
  ))
}

## ---------------------------------------------------------------------------
## 8. D4 - RAW vs BENJAMINI-HOCHBERG
## ---------------------------------------------------------------------------
log_head("8. D4 | RAW p < 0.05 vs BENJAMINI-HOCHBERG FDR")

fdr <- rbindlist(lapply(RULES_FOR_FDR, function(rule) {
  cn <- sprintf("eligible_%02d", round(rule * 100))
  elig <- valid[get(cn) == TRUE, .(period, index, grid_id)]
  mkx  <- merge(mk, elig, by = c("period", "index", "grid_id"))
  rbindlist(lapply(Q_LEVELS, function(q) {
    mkx[!is.na(p), {
      padj <- stats::p.adjust(p, method = "BH")
      .(completeness_rule = rule, q_level = q,
        n_tested = .N,
        n_expected_false_at_raw = round(ALPHA_RAW * .N, 1),
        n_sig_raw = sum(p < ALPHA_RAW),
        n_sig_bh  = sum(padj <= q),
        n_lost_to_bh = sum(p < ALPHA_RAW & padj > q),
        n_gained_by_bh = sum(p >= ALPHA_RAW & padj <= q),
        min_p = min(p), min_q = min(padj),
        n_sig_raw_positive = sum(p < ALPHA_RAW & sen_per_decade > 0),
        n_sig_raw_negative = sum(p < ALPHA_RAW & sen_per_decade < 0),
        n_sig_bh_positive  = sum(padj <= q & sen_per_decade > 0),
        n_sig_bh_negative  = sum(padj <= q & sen_per_decade < 0))
    }, by = .(period, index, variant)]
  }))
}))
setorder(fdr, completeness_rule, period, index, variant, q_level)

log_msg("Raw vs BH at the 80% completeness rule, HR_lag3 (grids of 283):")
zz <- fdr[completeness_rule == 0.80 & variant == "HR_lag3"]
for (i in seq_len(nrow(zz))) {
  z <- zz[i]
  log_msg(sprintf("  %-16s %-5s q=%.2f : tested %3d, expected-false %4.1f, raw %3d, BH %3d (lost %3d, gained %3d)",
                  z$period, z$index, z$q_level, z$n_tested, z$n_expected_false_at_raw,
                  z$n_sig_raw, z$n_sig_bh, z$n_lost_to_bh, z$n_gained_by_bh))
}

record_check("bh_never_exceeds_raw_at_same_alpha",
             fdr[q_level == ALPHA_RAW & n_sig_bh > n_sig_raw, .N] == 0L,
             sprintf("%d rows where BH at q=0.05 exceeded raw p<0.05 (would indicate a bug)",
                     fdr[q_level == ALPHA_RAW & n_sig_bh > n_sig_raw, .N]))
record_check("bh_monotone_in_q",
             {
               chk <- dcast(fdr[, .(period, index, variant, completeness_rule, q_level, n_sig_bh)],
                            period + index + variant + completeness_rule ~ q_level,
                            value.var = "n_sig_bh")
               qn <- as.character(Q_LEVELS)
               all(chk[[qn[1]]] <= chk[[qn[2]]]) && all(chk[[qn[2]]] <= chk[[qn[3]]])
             },
             "BH counts increase monotonically with q, as they must")

## ---------------------------------------------------------------------------
## 9. WSDI DISTRIBUTIONAL DIAGNOSTIC
## ---------------------------------------------------------------------------
log_head("9. WSDI TIE AND ZERO STRUCTURE")

wsdi_dist <- rbindlist(lapply(names(PERIODS), function(pn) {
  yr <- PERIODS[[pn]]
  sub <- d6[year >= yr[1] & year <= yr[2]]
  s <- sub[, .(n_valid = sum(!is.na(WSDI)),
               n_zero  = sum(WSDI == 0L, na.rm = TRUE),
               pct_zero = round(100 * mean(WSDI == 0L, na.rm = TRUE), 2),
               max_tie_group = { v <- WSDI[!is.na(WSDI)]
               if (!length(v)) NA_integer_ else max(table(v)) },
               n_distinct = uniqueN(WSDI[!is.na(WSDI)]),
               max_WSDI = suppressWarnings(max(WSDI, na.rm = TRUE))),
           by = grid_id]
  s[!is.finite(max_WSDI), max_WSDI := NA_integer_]
  s[, period := pn][]
}))
wsdi_dist <- merge(wsdi_dist,
                   mk[index == "WSDI" & variant == "HR_lag3",
                      .(period, grid_id, sen_per_decade, p)],
                   by = c("period", "grid_id"), all.x = TRUE)
wsdi_dist[, sen_exactly_zero := is.finite(sen_per_decade) & sen_per_decade == 0]
setorder(wsdi_dist, period, grid_id)

for (pn in names(PERIODS)) {
  z <- wsdi_dist[period == pn]
  log_msg(sprintf("  %-16s : mean %% zero WSDI years = %.1f%% ; grids with Sen slope exactly 0 = %d of %d ; median distinct values per series = %.0f",
                  pn, mean(z$pct_zero, na.rm = TRUE), sum(z$sen_exactly_zero, na.rm = TRUE),
                  nrow(z), stats::median(z$n_distinct, na.rm = TRUE)))
}
record_check("wsdi_tie_structure_reported", TRUE,
             sprintf("zero-inflation and tie structure recorded; tie-corrected Var(S) is therefore essential for WSDI"))

## ---------------------------------------------------------------------------
## 10. WRITE DIAGNOSTIC OUTPUTS
## ---------------------------------------------------------------------------
log_head("10. WRITE DIAGNOSTIC OUTPUTS")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
fwrite(valid,      F_VALID)
fwrite(sweep,      F_SWEEP)
fwrite(acf_sum,    F_ACFS)
fwrite(acf_g,      F_ACFG)
fwrite(mk,         F_MK)
fwrite(vcmp,       F_VCMP)
fwrite(fdr,        F_FDR)
fwrite(wsdi_dist,  F_WSDI)
fwrite(cf_stats,   file.path(OUT_DIR, "AMJ_STEP07A_cf_summary.csv"))

for (f in c(F_VALID, F_SWEEP, F_ACFS, F_ACFG, F_MK, F_VCMP, F_FDR, F_WSDI))
  log_msg("Saved            : ", basename(f), sprintf("  (%.2f MB)", file.size(f) / 1024^2))

.mtime_after <- file.mtime(IN_FILE)
record_check("input_unchanged", identical(.mtime_before, .mtime_after),
             "STEP 06 parquet opened read-only; modification time unchanged")
record_check("no_trend_results_published", TRUE,
             "diagnostic CSVs only; no analysis-ready trend product written")

qc <- rbindlist(CHECKS)
qc[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
          script = "07A_trend_diagnostic.R")]
qc <- rbind(qc, data.table(check = paste0("selftest_", tt$test),
                           status = fifelse(tt$ok, "PASS", "FATAL"),
                           detail = sprintf("got %s, want %s | %s", tt$got, tt$want, tt$note),
                           run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
                           script = "07A_trend_diagnostic.R"))
fwrite(qc, F_QC)
log_msg("QC summary       : ", basename(F_QC), "  (", nrow(qc), " entries)")

## ---------------------------------------------------------------------------
## 11. DECISION BRIEF - what the numbers above imply for STEP 07
## ---------------------------------------------------------------------------
log_head("11. DECISION BRIEF (read these numbers, then lock STEP 07)")

log_msg("D1 minimum-n rule")
log_msg("    Compare AMJ_STEP07A_completeness_sweep.csv across rules. Choose the")
log_msg("    loosest rule at which the eligible-grid count stops changing")
log_msg("    materially; if 80% and 90% give the same set, the rule is not doing")
log_msg("    any work and the looser one is easier to defend.")
for (pn in names(PERIODS)) {
  z <- sweep[period == pn & rule_frac %in% c(0.70, 0.80, 0.90, 1.00)]
  log_msg(sprintf("    %-16s : %s", pn,
                  paste(sprintf("%s %d%%->%d", z$index, round(z$rule_frac * 100),
                                z$n_grids_eligible), collapse = "  ")))
}
log_msg("D2/D3 lag truncation")
log_msg("    If the significant-rho count beyond lag 3 is near zero and the")
log_msg("    HR_lag3-vs-HR_all discordance is a handful of grids, lag<=3 and")
log_msg("    all-lags are equivalent in practice: pick lag<=3 and report the")
log_msg("    other as a sensitivity. If discordance is large, prefer HR_all,")
log_msg("    since the (n-i)(n-i-1)(n-i-2) weight already suppresses high lags.")
log_msg("    Also inspect the CF<1 column: if variance deflation is common, the")
log_msg("    floor-at-1 policy must be stated explicitly in the Methods.")
log_msg("D4 multiplicity")
log_msg("    Compare n_sig_raw against n_expected_false_at_raw. Where raw counts")
log_msg("    barely exceed the ~14-per-283 expected by chance, no spatial cluster")
log_msg("    should be interpreted without FDR. Wilks (2016) recommends")
log_msg("    q = 2 x alpha_global, i.e. q = 0.10 for a 0.05 global level.")
log_msg("Additional points already evidenced above")
log_msg("    - WSDI is heavily tied and zero-inflated; the tie-corrected Var(S)")
log_msg("      is mandatory and many Sen slopes will be exactly zero.")
log_msg("    - Series with internal gaps violate the even-spacing assumption of")
log_msg("      the CF formula; the has_gaps column identifies them.")
log_msg("    - Declare 1996-2025 primary and 1976-2025 contextual BEFORE seeing")
log_msg("      the final results, and FDR-adjust each period separately.")

log_head("RUN COMPLETE")
log_msg("Self-tests       : ", sum(tt$ok), " / ", nrow(tt), " passed")
log_msg("Checks passed    : ", qc[status == "PASS", .N], " / ", nrow(qc))
log_msg("Warnings         : ", qc[status == "WARN", .N])
log_msg("Fatal            : ", qc[status == "FATAL", .N])
log_msg("Series analysed  : ", format(uniqueN(mk[, .(period, index, grid_id)]), big.mark = ","),
        " ; MK rows = ", format(nrow(mk), big.mark = ","))
log_msg("Elapsed          : ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
log_msg("REMINDER         : these are diagnostics, not results. Lock D1-D4, then run STEP 07.")

.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
writeLines(.LOG, F_LOG)
cat("\nConsole log written to:", F_LOG, "\n")

invisible(list(valid = valid, sweep = sweep, acf = acf_sum, mk = mk,
               variants = vcmp, fdr = fdr, wsdi = wsdi_dist))
###############################################################################
## END OF SCRIPT
###############################################################################