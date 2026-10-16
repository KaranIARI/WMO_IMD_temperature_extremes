###############################################################################
##  IMD 1-degree GRD -> Temperature Extremes Analysis
##  STEP 07 : FINAL TREND ANALYSIS
##            Modified Mann-Kendall (Hamed & Rao 1998) + Theil-Sen + BH-FDR
##
##  283 grids | AMJ (01 Apr - 30 Jun) | indices: WSDI, TXx, TXn, TNx, TNn, DTR
##
## ===========================================================================
##  LOCKED PARAMETERS - fixed by the STEP 07A diagnostic, not revisited here
## ===========================================================================
##  PRIMARY PERIOD        1996-2025 (n = 30). Pre-registered study design.
##  CONTEXT PERIOD        1976-2025 (n = 50). Sensitivity/supplement ONLY.
##                        It does NOT replace the primary period anywhere.
##  COMPLETENESS          >= 80% of period years (>= 24 of 30; >= 40 of 50).
##                        07A showed all 283 grids satisfy this at 1996-2025,
##                        so the rule excludes no grid; it is retained because
##                        it was pre-specified and matches the baseline logic.
##  PRIMARY TEST          Hamed-Rao Modified Mann-Kendall, lag truncation 3.
##                        Chosen over HR_all because at n = 30 the all-lag
##                        variant screens 29 lags at alpha = 0.05, retains ~1.5
##                        spurious lags per series, and - since the mean lag-1
##                        rank ACF is NEGATIVE for every index at 1996-2025
##                        (-0.04 to -0.14) - systematically DEFLATES the
##                        variance. HR_all is therefore the anti-conservative
##                        choice: CF < 1 in 17-26% of grids under HR_all versus
##                        2.8-20% under HR_lag3, and all CF <= 0 cases fall in
##                        HR_all / HR_lag5.
##  SENSITIVITY TESTS     HR_all and MK_original, reported alongside, never
##                        substituted for the primary result.
##  SCREENING BOUND       +/- z_{1-alpha/2} / sqrt(n), alpha = 0.05, matching
##                        modifiedmk::mmkh() and pymannkendall. Verified to
##                        1e-16 against mmkh() in STEP 07A.
##  VARIANCE              Tie-corrected Var(S). Mandatory: WSDI is integer and
##                        heavily tied (68% zero years at 1996-2025).
##  CF POLICY             Faithful to Hamed-Rao. CF is NOT floored at 1.
##                        CF <= 0 yields an undefined test: p, q and CI are NA
##                        and the grid is classified "test_undefined", NOT
##                        "no significant trend".
##  MULTIPLE TESTING      Benjamini-Hochberg FDR, q <= 0.10 (Wilks 2016:
##                        q = 2 x alpha_global for a 0.05 global level).
##                        One family = one index x one period x one variant,
##                        over testable grids only. Periods are NEVER pooled.
##  SLOPE                 Theil-Sen, reported per DECADE with explicit units.
##                        CI from Gilbert (1987) using the SAME corrected
##                        variance that produced the p-value.
##  CLASSIFICATION        Significance (categorical, from q) and magnitude
##                        (continuous, with units) are kept strictly separate.
##                        No "strong"/"weak" labels and no universal magnitude
##                        threshold anywhere in the pipeline.
##
##  UNITS   WSDI days/decade ; TXx TXn TNx TNn DTR degC/decade
##
## ===========================================================================
##  INPUT (READ-ONLY - never modified)
##    F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet
##
##  OUTPUTS
##    A. AMJ_STEP07_trends_all_variants.parquet   grid x index x period x variant
##    B. AMJ_STEP07_trends_PRIMARY_1996_2025.csv  primary result, HR_lag3
##    C. AMJ_STEP07_maps_long.csv                 tidy, one row per grid-index
##    D. AMJ_STEP07_maps_wide_1996_2025.csv       one row per grid, all indices
##    E. AMJ_STEP07_summary_by_index.csv          raw vs BH counts, class counts
##    F. AMJ_STEP07_sensitivity_variants.csv      HR_lag3 vs HR_all vs MK
##    G. AMJ_STEP07_wsdi_zero_inflation.csv       WSDI tie / zero structure
##    H. AMJ_STEP07_QC.csv  +  AMJ_STEP07_LOG.txt
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
## 1. LOCKED CONSTANTS
## ---------------------------------------------------------------------------
IN_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet"
OUT_DIR <- "F:/WMO_IMD_R/WMO_IMD/data"

F_ALL   <- file.path(OUT_DIR, "AMJ_STEP07_trends_all_variants.parquet")
F_PRIM  <- file.path(OUT_DIR, "AMJ_STEP07_trends_PRIMARY_1996_2025.csv")
F_MAPL  <- file.path(OUT_DIR, "AMJ_STEP07_maps_long.csv")
F_MAPW  <- file.path(OUT_DIR, "AMJ_STEP07_maps_wide_1996_2025.csv")
F_SUMM  <- file.path(OUT_DIR, "AMJ_STEP07_summary_by_index.csv")
F_SENS  <- file.path(OUT_DIR, "AMJ_STEP07_sensitivity_variants.csv")
F_WSDI  <- file.path(OUT_DIR, "AMJ_STEP07_wsdi_zero_inflation.csv")
F_QC    <- file.path(OUT_DIR, "AMJ_STEP07_QC.csv")
F_LOG   <- file.path(OUT_DIR, "AMJ_STEP07_LOG.txt")

INDICES <- c("WSDI", "TXx", "TXn", "TNx", "TNn", "DTR")
USABLE  <- c(WSDI = "usable_WSDI", TXx = "usable_TXx", TXn = "usable_TXn",
             TNx = "usable_TNx",  TNn = "usable_TNn", DTR = "usable_DTR")
UNITS   <- c(WSDI = "days/decade", TXx = "degC/decade", TXn = "degC/decade",
             TNx = "degC/decade",  TNn = "degC/decade", DTR = "degC/decade")

PRIMARY_PERIOD <- "TREND_1996_2025"
PERIODS <- list(TREND_1996_2025 = c(1996L, 2025L),   # PRIMARY
                FULL_1976_2025  = c(1976L, 2025L))   # context / sensitivity
PERIOD_ROLE <- c(TREND_1996_2025 = "PRIMARY", FULL_1976_2025 = "CONTEXT_SENSITIVITY")

COMPLETENESS_FRAC <- 0.80
MIN_N_ABSOLUTE    <- 10L

PRIMARY_VARIANT <- "HR_lag3"
LAG_VARIANTS    <- list(HR_lag3 = 3L, HR_all = NA_integer_, MK_original = 0L)

HR_SCREEN <- "standard"        # z/sqrt(n), as in modifiedmk::mmkh
ACF_ALPHA <- 0.05
CI_ALPHA  <- 0.05              # 95% Theil-Sen confidence interval
Q_PRIMARY <- 0.10              # BH-FDR threshold (Wilks 2016)
ALPHA_RAW <- 0.05              # reported only, for the raw-vs-BH comparison
ACF_REPORT_LAGS <- 10L

EXP_GRIDS <- 283L
EXP_ROWS  <- 14150L
YR_START  <- 1976L; YR_END <- 2025L

CLASS_LEVELS <- c("significant_increase", "significant_decrease",
                  "significant_zero_slope", "no_significant_trend",
                  "test_undefined", "insufficient_data")

## Hard gates on the locked settings.
stopifnot(identical(PRIMARY_PERIOD, "TREND_1996_2025"),
          identical(PRIMARY_VARIANT, "HR_lag3"),
          identical(HR_SCREEN, "standard"),
          COMPLETENESS_FRAC == 0.80, Q_PRIMARY == 0.10,
          PRIMARY_VARIANT %in% names(LAG_VARIANTS),
          identical(PERIODS[[PRIMARY_PERIOD]], c(1996L, 2025L)))

## ---------------------------------------------------------------------------
## 2. LOGGING, CHECKS, AND THE STEP 07A ENGINE (verbatim, unchanged)
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
  log_msg(sprintf("  %-5s | %-46s | %s", status, name, detail))
  if (!ok && severity == "FATAL")
    stop("FATAL QC failure: ", name, " -- ", detail, call. = FALSE)
  invisible(ok)
}
rng <- function(v, f = "%.3f") { v <- v[is.finite(v)]
if (!length(v)) "all NA" else sprintf(paste0("[", f, ", ", f, "]"), min(v), max(v)) }

## ---- ENGINE (identical to the verified STEP 07A implementation) -------------
mk_S_sen <- function(x, t) {
  ok <- !is.na(x); xv <- x[ok]; tv <- t[ok]; n <- length(xv)
  if (n < 2L) return(list(S = NA_real_, sen = NA_real_, n = n))
  dx <- outer(xv, xv, "-"); dt <- outer(tv, tv, "-"); lt <- lower.tri(dx)
  list(S = as.numeric(sum(sign(dx[lt]))),
       sen = as.numeric(stats::median(dx[lt] / dt[lt])), n = n)
}
mk_var_ties <- function(x) {
  xv <- x[!is.na(x)]; n <- length(xv)
  if (n < 3L) return(NA_real_)
  v <- n * (n - 1) * (2 * n + 5)
  tg <- as.numeric(table(xv)); tg <- tg[tg > 1]
  if (length(tg)) v <- v - sum(tg * (tg - 1) * (2 * tg + 5))
  v / 18
}
acf_hr <- function(r, max_lag) {
  n <- length(r); out <- rep(0, max(max_lag, 0L))
  if (max_lag < 1L) return(out)
  mu <- mean(r, na.rm = TRUE); den <- sum((r - mu)^2, na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(out)
  for (k in seq_len(min(max_lag, n - 1L))) {
    a <- r[1:(n - k)]; b <- r[(k + 1):n]; ok <- !is.na(a) & !is.na(b)
    if (any(ok)) out[k] <- sum((a[ok] - mu) * (b[ok] - mu)) / den
  }
  out[!is.finite(out)] <- 0
  out
}
hr_significant_rho <- function(rho, n, alpha = ACF_ALPHA, rule = HR_SCREEN) {
  k <- seq_along(rho)
  if (identical(rule, "standard")) {
    b <- stats::qnorm(1 - alpha / 2) / sqrt(n)
    return(is.finite(rho) & abs(rho) > b)
  }
  z <- stats::qnorm(1 - alpha / 2); d <- n - k
  lo <- (-1 - z * sqrt(pmax(d - 1, 0))) / d
  hi <- (-1 + z * sqrt(pmax(d - 1, 0))) / d
  is.finite(rho) & d > 1 & (rho < lo | rho > hi)
}
hr_cf <- function(rho, sig, n, max_lag) {
  if (max_lag < 1L) return(1)
  k <- seq_len(max_lag)
  r <- rho[k]; r[!sig[k]] <- 0; r[!is.finite(r)] <- 0
  w <- (n - k) * (n - k - 1) * (n - k - 2); w[w < 0] <- 0
  1 + (2 / (n * (n - 1) * (n - 2))) * sum(w * r)
}
mk_series <- function(x, t, lag_variants = LAG_VARIANTS, alpha = ACF_ALPHA) {
  n_series <- length(x); base <- mk_S_sen(x, t); n_obs <- base$n
  varS <- mk_var_ties(x); empty_acf <- rep(NA_real_, ACF_REPORT_LAGS)
  if (is.na(base$S) || is.na(varS) || n_obs < MIN_N_ABSOLUTE || n_series < 4L)
    return(list(ok = FALSE, n_obs = n_obs, n_series = n_series, S = base$S,
                varS = varS, sen = base$sen, acf = empty_acf, rows = NULL))
  xd <- x - base$sen * t
  rk <- rep(NA_real_, n_series)
  rk[!is.na(xd)] <- rank(xd[!is.na(xd)], ties.method = "average")
  max_lag_all <- n_series - 1L
  rho <- acf_hr(rk, max_lag_all)
  sig <- hr_significant_rho(rho, n_series, alpha, HR_SCREEN)
  out_acf <- empty_acf; nrep <- min(ACF_REPORT_LAGS, max_lag_all)
  out_acf[seq_len(nrep)] <- rho[seq_len(nrep)]
  rows <- vector("list", length(lag_variants))
  for (v in seq_along(lag_variants)) {
    nm <- names(lag_variants)[v]; L <- lag_variants[[v]]
    if (is.na(L)) L <- max_lag_all
    L <- min(L, max_lag_all)
    if (L < 1L) { cf <- 1; nsig <- 0L; ridx <- integer(0) }
    else { cf <- hr_cf(rho, sig, n_series, L); ridx <- which(sig[seq_len(L)]); nsig <- length(ridx) }
    cf_status <- if (!is.finite(cf)) "non_finite" else
      if (cf <= 0) "non_positive" else if (cf < 1) "deflated" else "ok"
    vx <- varS * cf
    if (!is.finite(vx) || vx <= 0) { z <- NA_real_; pv <- NA_real_ } else {
      z <- if (base$S > 0) (base$S - 1) / sqrt(vx) else
        if (base$S < 0) (base$S + 1) / sqrt(vx) else 0
      pv <- 2 * stats::pnorm(-abs(z))
    }
    rows[[v]] <- data.table(
      variant = nm, lags_used = L, n_lags_sig = as.integer(nsig),
      CF = cf, cf_status = cf_status, varS_corrected = vx, Z = z, p = pv,
      rho_lag1 = if (max_lag_all >= 1L) rho[1] else NA_real_,
      n_sig_rho_negative = if (nsig) sum(rho[ridx] < 0) else 0L,
      screen_rule = HR_SCREEN,
      screen_bound = stats::qnorm(1 - alpha / 2) / sqrt(n_series))
  }
  list(ok = TRUE, n_obs = n_obs, n_series = n_series, S = base$S, varS = varS,
       sen = base$sen, acf = out_acf, rows = rbindlist(rows))
}

## ---- Theil-Sen confidence interval, Gilbert (1987) --------------------------
## C = z_{1-alpha/2} * sqrt(Var*(S)) where Var*(S) is the SAME Hamed-Rao
## corrected variance used for the p-value, so the interval and the test are
## mutually consistent. With N ordered pairwise slopes,
##   M1 = (N - C)/2 , M2 = (N + C)/2
##   lower = slope[floor(M1)] , upper = slope[ceiling(M2) + 1]
## (conservative rounding; indices clamped to [1, N]). NA when the variance is
## not positive, i.e. exactly where the test itself is undefined.
sen_ci <- function(x, t, varS_corr, alpha = CI_ALPHA) {
  ok <- !is.na(x); xv <- x[ok]; tv <- t[ok]; n <- length(xv)
  if (n < 3L || !is.finite(varS_corr) || varS_corr <= 0) return(c(NA_real_, NA_real_))
  dx <- outer(xv, xv, "-"); dt <- outer(tv, tv, "-"); lt <- lower.tri(dx)
  sl <- sort(dx[lt] / dt[lt]); N <- length(sl)
  C  <- stats::qnorm(1 - alpha / 2) * sqrt(varS_corr)
  k1 <- max(1L, min(N, as.integer(floor((N - C) / 2))))
  k2 <- max(1L, min(N, as.integer(ceiling((N + C) / 2)) + 1L))
  c(sl[k1], sl[k2])
}

t_start <- Sys.time()
log_head("STEP 07 | FINAL TREND ANALYSIS (locked methodology)")
log_msg("R version        : ", R.version.string)
log_msg("data.table       : ", as.character(utils::packageVersion("data.table")))
log_msg("arrow            : ", as.character(utils::packageVersion("arrow")))
log_msg("PRIMARY period   : ", PRIMARY_PERIOD, " (", PERIODS[[PRIMARY_PERIOD]][1], "-",
        PERIODS[[PRIMARY_PERIOD]][2], ", n = 30)")
log_msg("Context period   : FULL_1976_2025 (sensitivity / supplement only)")
log_msg("Completeness     : >= ", COMPLETENESS_FRAC * 100, "% of period years")
log_msg("PRIMARY test     : Hamed-Rao MMK, ", PRIMARY_VARIANT,
        " ; screening +/- z/sqrt(n), alpha = ", ACF_ALPHA)
log_msg("Sensitivity      : ", paste(setdiff(names(LAG_VARIANTS), PRIMARY_VARIANT), collapse = ", "))
log_msg("CF policy        : faithful (NOT floored); CF <= 0 -> test_undefined")
log_msg("Variance         : tie-corrected Var(S)")
log_msg("Multiple testing : Benjamini-Hochberg, q <= ", Q_PRIMARY,
        " (family = index x period x variant)")
log_msg("Slope            : Theil-Sen per decade, Gilbert (1987) CI at ",
        (1 - CI_ALPHA) * 100, "% from the corrected variance")

## ---------------------------------------------------------------------------
## 3. READ STEP 06 (READ-ONLY)
## ---------------------------------------------------------------------------
log_head("3. READ STEP 06 OUTPUT (READ-ONLY)")

if (!file.exists(IN_FILE)) stop("STEP 06 parquet not found: ", IN_FILE, call. = FALSE)
.mtime_before <- file.mtime(IN_FILE)

d6 <- as.data.table(arrow::read_parquet(IN_FILE))
need <- c("grid_id", "year", "lat", "lon", INDICES, unname(USABLE))
miss <- setdiff(need, names(d6))
if (length(miss))
  stop("SCHEMA ERROR in STEP 06 parquet.\n  Missing: ", paste(miss, collapse = ", "),
       "\n  Present: ", paste(names(d6), collapse = ", "), call. = FALSE)
d6[, `:=`(grid_id = as.character(grid_id), year = as.integer(year))]
setorder(d6, grid_id, year)

record_check("schema_step06", TRUE, sprintf("all %d required columns present", length(need)))
record_check("input_row_count", nrow(d6) == EXP_ROWS,
             sprintf("observed %s vs expected %s", format(nrow(d6), big.mark = ","),
                     format(EXP_ROWS, big.mark = ",")))
record_check("input_grid_count", uniqueN(d6$grid_id) == EXP_GRIDS,
             sprintf("observed %d vs expected %d", uniqueN(d6$grid_id), EXP_GRIDS))
record_check("input_year_coverage", identical(sort(unique(d6$year)), YR_START:YR_END),
             sprintf("%d years [%d-%d]", uniqueN(d6$year), min(d6$year), max(d6$year)))
record_check("input_no_duplicate_grid_year", !anyDuplicated(d6, by = c("grid_id", "year")),
             sprintf("%d duplicates", sum(duplicated(d6, by = c("grid_id", "year")))))

GRIDS <- sort(unique(d6$grid_id))
meta  <- unique(d6[, .(grid_id, lat, lon)])
record_check("unique_coords_per_grid", nrow(meta) == EXP_GRIDS,
             sprintf("%d unique grid/lat/lon triplets", nrow(meta)))

## ---------------------------------------------------------------------------
## 4. ENGINE VERIFICATION (must pass before any trend is produced)
## ---------------------------------------------------------------------------
log_head("4. ENGINE VERIFICATION")

TT <- list()
addt <- function(nm, got, want, tol = 1e-9, note = "") {
  ok <- is.finite(got) && is.finite(want) && abs(got - want) <= tol
  TT[[length(TT) + 1L]] <<- data.table(test = nm, got = sprintf("%.10g", got),
                                       want = sprintf("%.10g", want), ok = ok, note = note)
}
sim_ar1 <- function(n, phi, seed, burn = 200L) {
  set.seed(seed); e <- stats::rnorm(n + burn)
  as.numeric(stats::filter(e, filter = phi, method = "recursive"))[(burn + 1L):(burn + n)]
}

addt("S_perfect_increase", mk_S_sen(1:10, 1:10)$S, 45, 0, "n(n-1)/2")
addt("sen_exact_linear", mk_S_sen(2 * (1:20) + 7, 1:20)$sen, 2, 1e-12, "y = 2t + 7")
xt <- c(1, 1, 1, 2, 2, 3, 4, 5, 5, 5); nt <- 10
addt("var_ties_formula", mk_var_ties(xt),
     (nt * (nt - 1) * (2 * nt + 5) - (3 * 2 * 11 + 2 * 1 * 9 + 3 * 2 * 11)) / 18,
     1e-12, "tie correction retained for WSDI")
lin <- as.numeric(1:30)
r_lin <- mk_series(lin, 1:30, list(MK_original = 0L, HR_all = NA_integer_))
addt("CF_one_exact_linear", r_lin$rows[variant == "HR_all", CF], 1, 0,
     "detrended residuals identically zero")
addt("HR_reduces_to_MK", r_lin$rows[variant == "HR_all", Z],
     r_lin$rows[variant == "MK_original", Z], 0, "identical Z when CF = 1")
addt("screen_bound_is_z_over_sqrt_n", stats::qnorm(0.975) / sqrt(30),
     unique(mk_series(sim_ar1(30L, 0.3, 7L) + 0.05 * (1:30), 1:30,
                      list(HR_all = NA_integer_))$rows$screen_bound),
     1e-12, "matches modifiedmk::mmkh")
rho_t <- c(0.4, 0.2, 0.1); n_t <- 20; k_t <- 1:3
addt("hr_cf_formula", hr_cf(rho_t, rep(TRUE, 3), n_t, 3L),
     1 + (2 / (n_t * (n_t - 1) * (n_t - 2))) *
       sum((n_t - k_t) * (n_t - k_t - 1) * (n_t - k_t - 2) * rho_t), 1e-12,
     "weights (n-i)(n-i-1)(n-i-2)")
## Sen CI must bracket the slope on a clean linear-plus-noise series
xs <- 0.5 * (1:30) + sim_ar1(30L, 0, 11L)
rs <- mk_series(xs, 1:30, list(HR_lag3 = 3L))
ci <- sen_ci(xs, 1:30, rs$rows[1, varS_corrected])
addt("sen_ci_brackets_slope",
     as.numeric(ci[1] <= rs$sen && rs$sen <= ci[2]), 1, 0,
     sprintf("slope %.4f in [%.4f, %.4f]", rs$sen, ci[1], ci[2]))

if (requireNamespace("modifiedmk", quietly = TRUE)) {
  phis <- rep(c(-0.6, -0.3, 0.0, 0.3, 0.7), each = 6L); dmax <- 0; ncmp <- 0L
  for (i in seq_along(phis)) {
    xx <- sim_ar1(30L, phis[i], seed = 90000L + i) + 0.03 * (1:30)
    mine <- mk_series(xx, 1:30, list(HR_all = NA_integer_))$rows[1]
    ref <- tryCatch(modifiedmk::mmkh(xx), error = function(e) NULL)
    if (is.null(ref)) next
    pr <- unname(ref["new P-value"])
    if (is.finite(pr) && is.finite(mine$p)) { dmax <- max(dmax, abs(mine$p - pr)); ncmp <- ncmp + 1L }
  }
  addt("reproduces_modifiedmk_mmkh", dmax, 0, 1e-8,
       sprintf("identical p over %d AR(1) series", ncmp))
} else {
  log_msg("  NOTE  | modifiedmk not installed; the external cross-check is skipped.")
}

tt <- rbindlist(TT)
for (i in seq_len(nrow(tt)))
  log_msg(sprintf("  %-4s | %-32s | got %-14s want %-14s | %s",
                  if (tt$ok[i]) "PASS" else "FAIL", tt$test[i], tt$got[i], tt$want[i], tt$note[i]))
record_check("engine_verification", all(tt$ok),
             sprintf("%d of %d passed", sum(tt$ok), nrow(tt)))

## ---------------------------------------------------------------------------
## 5. TREND COMPUTATION
## ---------------------------------------------------------------------------
log_head("5. TREND COMPUTATION (Mann-Kendall + Hamed-Rao + Theil-Sen)")

ROWS <- vector("list", length(PERIODS) * length(INDICES) * EXP_GRIDS)
kk <- 0L; t0 <- Sys.time()

for (pn in names(PERIODS)) {
  yr <- PERIODS[[pn]]; yrs <- yr[1]:yr[2]
  n_period <- length(yrs)
  min_years <- as.integer(ceiling(COMPLETENESS_FRAC * n_period))
  sub <- d6[year >= yr[1] & year <= yr[2]]
  for (ix in INDICES) {
    w <- dcast(sub[, .(grid_id, year, v = get(ix))], grid_id ~ year, value.var = "v")
    setkey(w, grid_id); w <- w[J(GRIDS)]
    M <- as.matrix(w[, -1L, with = FALSE])
    stopifnot(ncol(M) == n_period)
    tvec <- as.numeric(seq_along(yrs))
    for (g in seq_len(nrow(M))) {
      kk <- kk + 1L
      xg <- as.numeric(M[g, ])
      n_valid <- sum(!is.na(xg))
      eligible <- n_valid >= min_years && n_valid >= MIN_N_ABSOLUTE
      yy <- yrs[!is.na(xg)]
      n_gaps <- if (length(yy) < 2L) 0L else sum(diff(yy) != 1L)
      
      base <- data.table(
        period = pn, period_role = PERIOD_ROLE[[pn]], index = ix, grid_id = GRIDS[g],
        n_period_years = n_period, min_years_required = min_years,
        n_valid = n_valid, completeness_pct = round(100 * n_valid / n_period, 2),
        n_gaps = n_gaps, has_gaps = n_gaps > 0L, eligible = eligible,
        units = UNITS[[ix]])
      
      if (!eligible) {
        ROWS[[kk]] <- cbind(base[rep(1L, length(LAG_VARIANTS))],
                            data.table(variant = names(LAG_VARIANTS),
                                       S = NA_real_, varS_uncorrected = NA_real_,
                                       lags_used = NA_integer_, n_lags_sig = NA_integer_,
                                       CF = NA_real_, cf_status = "insufficient_data",
                                       varS_corrected = NA_real_, Z = NA_real_, p = NA_real_,
                                       rho_lag1 = NA_real_, n_sig_rho_negative = NA_integer_,
                                       screen_rule = HR_SCREEN, screen_bound = NA_real_,
                                       sen_per_year = NA_real_, sen_per_decade = NA_real_,
                                       ci_lo_per_decade = NA_real_, ci_hi_per_decade = NA_real_))
        next
      }
      
      r <- mk_series(xg, tvec)
      if (!isTRUE(r$ok)) {
        ROWS[[kk]] <- cbind(base[rep(1L, length(LAG_VARIANTS))],
                            data.table(variant = names(LAG_VARIANTS),
                                       S = NA_real_, varS_uncorrected = NA_real_,
                                       lags_used = NA_integer_, n_lags_sig = NA_integer_,
                                       CF = NA_real_, cf_status = "insufficient_data",
                                       varS_corrected = NA_real_, Z = NA_real_, p = NA_real_,
                                       rho_lag1 = NA_real_, n_sig_rho_negative = NA_integer_,
                                       screen_rule = HR_SCREEN, screen_bound = NA_real_,
                                       sen_per_year = NA_real_, sen_per_decade = NA_real_,
                                       ci_lo_per_decade = NA_real_, ci_hi_per_decade = NA_real_))
        next
      }
      
      rr <- copy(r$rows)
      rr[, `:=`(S = r$S, varS_uncorrected = r$varS,
                sen_per_year = r$sen, sen_per_decade = r$sen * 10)]
      ci_lo <- ci_hi <- rep(NA_real_, nrow(rr))
      for (v in seq_len(nrow(rr))) {
        cc <- sen_ci(xg, tvec, rr$varS_corrected[v])
        ci_lo[v] <- cc[1] * 10; ci_hi[v] <- cc[2] * 10
      }
      rr[, `:=`(ci_lo_per_decade = ci_lo, ci_hi_per_decade = ci_hi)]
      ROWS[[kk]] <- cbind(base[rep(1L, nrow(rr))], rr)
    }
  }
  log_msg("Computed ", pn, " in ",
          sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

tr <- rbindlist(ROWS, use.names = TRUE, fill = TRUE)
rm(ROWS); invisible(gc(verbose = FALSE))
tr <- merge(tr, meta, by = "grid_id", all.x = TRUE, sort = FALSE)

EXP_TR_ROWS <- EXP_GRIDS * length(INDICES) * length(PERIODS) * length(LAG_VARIANTS)
record_check("trend_row_count", nrow(tr) == EXP_TR_ROWS,
             sprintf("observed %s vs expected %s (283 x 6 x 2 x 3)",
                     format(nrow(tr), big.mark = ","), format(EXP_TR_ROWS, big.mark = ",")))
record_check("no_duplicate_keys",
             !anyDuplicated(tr, by = c("period", "index", "grid_id", "variant")),
             sprintf("%d duplicated (period, index, grid_id, variant)",
                     sum(duplicated(tr, by = c("period", "index", "grid_id", "variant")))))
log_msg("Trend rows       : ", format(nrow(tr), big.mark = ","), " in ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("cf_status counts : ",
        paste(sprintf("%s=%d", names(table(tr$cf_status)), as.integer(table(tr$cf_status))),
              collapse = "  "))

## ---------------------------------------------------------------------------
## 6. BENJAMINI-HOCHBERG FDR
## ---------------------------------------------------------------------------
log_head("6. BENJAMINI-HOCHBERG FDR (family = index x period x variant)")

tr[, testable := is.finite(p)]
tr[, q := NA_real_]
tr[testable == TRUE, q := stats::p.adjust(p, method = "BH"),
   by = .(period, index, variant)]
tr[, n_tested_in_family := sum(testable), by = .(period, index, variant)]

record_check("q_defined_iff_testable",
             tr[testable == TRUE & is.na(q), .N] == 0L &&
               tr[testable == FALSE & !is.na(q), .N] == 0L,
             sprintf("%d testable rows without q; %d untestable rows with q",
                     tr[testable == TRUE & is.na(q), .N],
                     tr[testable == FALSE & !is.na(q), .N]))
record_check("q_never_below_p",
             tr[testable == TRUE & q < p - 1e-12, .N] == 0L,
             sprintf("%d rows with q < p (impossible under BH)",
                     tr[testable == TRUE & q < p - 1e-12, .N]))
record_check("fdr_families_are_per_period",
             tr[testable == TRUE, uniqueN(period), by = .(index, variant)][, all(V1 <= 2)] &&
               tr[testable == TRUE, .N, by = .(period, index, variant)][, all(N <= EXP_GRIDS)],
             "no family exceeds 283 grids and periods are never pooled")

## ---------------------------------------------------------------------------
## 7. TREND CLASSIFICATION
## ---------------------------------------------------------------------------
log_head("7. CLASSIFICATION (significance categorical; magnitude continuous)")

## Significance and magnitude are kept strictly separate. No strong/weak labels,
## no magnitude threshold. "significant_zero_slope" exists because WSDI can have
## S != 0 while the median pairwise slope is exactly 0 (68% zero years).
tr[, trend_class := fifelse(
  eligible == FALSE, "insufficient_data",
  fifelse(!testable, "test_undefined",
          fifelse(q <= Q_PRIMARY & sen_per_decade > 0, "significant_increase",
                  fifelse(q <= Q_PRIMARY & sen_per_decade < 0, "significant_decrease",
                          fifelse(q <= Q_PRIMARY & sen_per_decade == 0, "significant_zero_slope",
                                  "no_significant_trend")))))]
tr[, trend_class := factor(trend_class, levels = CLASS_LEVELS)]
tr[, significant_fdr := trend_class %in%
     c("significant_increase", "significant_decrease", "significant_zero_slope")]
tr[, significant_raw := testable & p < ALPHA_RAW]
tr[, q_threshold := Q_PRIMARY]
tr[, is_primary := period == PRIMARY_PERIOD & variant == PRIMARY_VARIANT]

record_check("classification_exhaustive", tr[is.na(trend_class), .N] == 0L,
             sprintf("%d unclassified rows", tr[is.na(trend_class), .N]))
record_check("classification_counts_sum",
             tr[, .N, by = .(period, index, variant)][, all(N == EXP_GRIDS)],
             "every (period, index, variant) covers exactly 283 grids")
record_check("significant_only_when_q_le_threshold",
             tr[significant_fdr == TRUE & (is.na(q) | q > Q_PRIMARY), .N] == 0L,
             sprintf("%d rows flagged significant without q <= %.2f",
                     tr[significant_fdr == TRUE & (is.na(q) | q > Q_PRIMARY), .N], Q_PRIMARY))
record_check("test_undefined_has_no_p_q_ci",
             tr[trend_class == "test_undefined" &
                  (!is.na(p) | !is.na(q) | !is.na(ci_lo_per_decade)), .N] == 0L,
             sprintf("%d test_undefined rows carry p, q or CI",
                     tr[trend_class == "test_undefined" &
                          (!is.na(p) | !is.na(q) | !is.na(ci_lo_per_decade)), .N]))
record_check("test_undefined_matches_CF_non_positive",
             tr[trend_class == "test_undefined" &
                  !(cf_status %in% c("non_positive", "non_finite")), .N] == 0L,
             sprintf("%d test_undefined rows not explained by CF <= 0",
                     tr[trend_class == "test_undefined" &
                          !(cf_status %in% c("non_positive", "non_finite")), .N]))
record_check("ci_brackets_sen_slope",
             tr[is.finite(ci_lo_per_decade) & is.finite(sen_per_decade) &
                  (sen_per_decade < ci_lo_per_decade - 1e-9 |
                     sen_per_decade > ci_hi_per_decade + 1e-9), .N] == 0L,
             sprintf("%d rows where the Sen slope falls outside its own CI",
                     tr[is.finite(ci_lo_per_decade) & is.finite(sen_per_decade) &
                          (sen_per_decade < ci_lo_per_decade - 1e-9 |
                             sen_per_decade > ci_hi_per_decade + 1e-9), .N]))
record_check("primary_rows_present",
             tr[is_primary == TRUE, .N] == EXP_GRIDS * length(INDICES),
             sprintf("%d primary rows vs expected %d (283 x 6)",
                     tr[is_primary == TRUE, .N], EXP_GRIDS * length(INDICES)))
record_check("units_populated", tr[is.na(units) | units == "", .N] == 0L,
             sprintf("%d rows without a units label", tr[is.na(units) | units == "", .N]))
record_check("completeness_rule_applied",
             tr[eligible == TRUE & n_valid < min_years_required, .N] == 0L,
             sprintf("all eligible rows meet the >= %.0f%% rule", COMPLETENESS_FRAC * 100))

## ---------------------------------------------------------------------------
## 8. SUMMARIES
## ---------------------------------------------------------------------------
log_head("8. SUMMARY TABLES")

summ <- tr[, .(
  n_grids = .N,
  n_eligible = sum(eligible),
  n_testable = sum(testable),
  n_test_undefined = sum(trend_class == "test_undefined"),
  n_insufficient = sum(trend_class == "insufficient_data"),
  n_expected_false_at_raw = round(ALPHA_RAW * sum(testable), 1),
  n_sig_raw = sum(significant_raw, na.rm = TRUE),
  n_sig_fdr = sum(significant_fdr),
  n_increase = sum(trend_class == "significant_increase"),
  n_decrease = sum(trend_class == "significant_decrease"),
  n_zero_slope_sig = sum(trend_class == "significant_zero_slope"),
  n_no_trend = sum(trend_class == "no_significant_trend"),
  median_sen_per_decade = stats::median(sen_per_decade, na.rm = TRUE),
  min_sen_per_decade = suppressWarnings(min(sen_per_decade, na.rm = TRUE)),
  max_sen_per_decade = suppressWarnings(max(sen_per_decade, na.rm = TRUE)),
  median_sen_sig_only = stats::median(sen_per_decade[significant_fdr], na.rm = TRUE),
  min_q = suppressWarnings(min(q, na.rm = TRUE)),
  units = units[1]),
  by = .(period, period_role, index, variant)]
for (cc in c("min_sen_per_decade", "max_sen_per_decade", "min_q"))
  set(summ, which(!is.finite(summ[[cc]])), cc, NA_real_)
setorder(summ, period, index, variant)

log_msg("PRIMARY RESULT | ", PRIMARY_PERIOD, " | ", PRIMARY_VARIANT,
        " | BH q <= ", Q_PRIMARY)
log_msg(sprintf("  %-5s %-12s %8s %8s %8s %6s %6s %10s",
                "index", "units", "testable", "raw<.05", "exp.false", "FDR+", "FDR-", "med slope"))
for (ix in INDICES) {
  z <- summ[period == PRIMARY_PERIOD & variant == PRIMARY_VARIANT & index == ix]
  if (!nrow(z)) next
  log_msg(sprintf("  %-5s %-12s %8d %8d %8.1f %6d %6d %10.4f",
                  z$index, z$units, z$n_testable, z$n_sig_raw,
                  z$n_expected_false_at_raw, z$n_increase, z$n_decrease,
                  z$median_sen_per_decade))
  if (z$n_test_undefined > 0L || z$n_insufficient > 0L)
    log_msg(sprintf("        (test_undefined %d ; insufficient_data %d)",
                    z$n_test_undefined, z$n_insufficient))
}

## sensitivity: agreement of the three variants on the primary period
sens <- dcast(tr[period == PRIMARY_PERIOD,
                 .(index, grid_id, variant, trend_class = as.character(trend_class), q)],
              index + grid_id ~ variant, value.var = c("trend_class", "q"))
sens_summ <- rbindlist(lapply(setdiff(names(LAG_VARIANTS), PRIMARY_VARIANT), function(v) {
  ca <- paste0("trend_class_", PRIMARY_VARIANT); cb <- paste0("trend_class_", v)
  sens[, .(comparison = paste0(PRIMARY_VARIANT, " vs ", v),
           n_grids = .N,
           n_same_class = sum(get(ca) == get(cb)),
           pct_same_class = round(100 * mean(get(ca) == get(cb)), 2),
           n_sig_primary = sum(grepl("^significant", get(ca))),
           n_sig_alt = sum(grepl("^significant", get(cb)))),
       by = index]
}))
setorder(sens_summ, index, comparison)
log_msg("Sensitivity agreement on the primary period (class identical, of 283):")
for (i in seq_len(nrow(sens_summ))) {
  z <- sens_summ[i]
  log_msg(sprintf("  %-5s %-26s : %3d same (%.1f%%) ; significant %3d vs %3d",
                  z$index, z$comparison, z$n_same_class, z$pct_same_class,
                  z$n_sig_primary, z$n_sig_alt))
}

## WSDI zero-inflation, reported explicitly as locked
wsdi_zi <- rbindlist(lapply(names(PERIODS), function(pn) {
  yr <- PERIODS[[pn]]
  s <- d6[year >= yr[1] & year <= yr[2],
          .(n_valid = sum(!is.na(WSDI)),
            n_zero_years = sum(WSDI == 0L, na.rm = TRUE),
            pct_zero_years = round(100 * mean(WSDI == 0L, na.rm = TRUE), 2),
            n_distinct_values = uniqueN(WSDI[!is.na(WSDI)]),
            max_WSDI = suppressWarnings(max(WSDI, na.rm = TRUE))), by = grid_id]
  s[!is.finite(max_WSDI), max_WSDI := NA_integer_]
  s[, period := pn][]
}))
wsdi_zi <- merge(wsdi_zi,
                 tr[index == "WSDI" & variant == PRIMARY_VARIANT,
                    .(period, grid_id, sen_per_decade, p, q,
                      trend_class = as.character(trend_class))],
                 by = c("period", "grid_id"), all.x = TRUE)
wsdi_zi[, sen_exactly_zero := is.finite(sen_per_decade) & sen_per_decade == 0]
setorder(wsdi_zi, period, grid_id)
for (pn in names(PERIODS)) {
  z <- wsdi_zi[period == pn]
  log_msg(sprintf("WSDI %-16s : mean %% zero years %.1f%% ; Sen slope exactly 0 on %d of %d grids",
                  pn, mean(z$pct_zero_years, na.rm = TRUE),
                  sum(z$sen_exactly_zero, na.rm = TRUE), nrow(z)))
}
record_check("wsdi_zero_inflation_reported", TRUE,
             sprintf("WSDI Sen slope is exactly zero on %d of %d grids at %s; report as a property of the index, not a null test",
                     wsdi_zi[period == PRIMARY_PERIOD & sen_exactly_zero == TRUE, .N],
                     EXP_GRIDS, PRIMARY_PERIOD),
             severity = "WARN")

## ---------------------------------------------------------------------------
## 9. MAP-READY PRODUCTS
## ---------------------------------------------------------------------------
log_head("9. MAP-READY PRODUCTS")

map_cols <- c("period", "period_role", "index", "grid_id", "lat", "lon",
              "sen_per_decade", "ci_lo_per_decade", "ci_hi_per_decade", "units",
              "S", "Z", "p", "q", "q_threshold", "trend_class", "significant_fdr",
              "significant_raw", "cf_status", "CF", "n_valid", "completeness_pct",
              "n_tested_in_family", "variant")
maps_long <- tr[variant == PRIMARY_VARIANT, ..map_cols]
setorder(maps_long, period, index, grid_id)

prim <- tr[is_primary == TRUE]
setorder(prim, index, grid_id)

wide <- dcast(prim[, .(grid_id, lat, lon, index,
                       slope = sen_per_decade, q = q,
                       class = as.character(trend_class))],
              grid_id + lat + lon ~ index, value.var = c("slope", "q", "class"))
setorder(wide, grid_id)

record_check("maps_long_row_count",
             nrow(maps_long) == EXP_GRIDS * length(INDICES) * length(PERIODS),
             sprintf("%s rows = 283 x 6 x 2", format(nrow(maps_long), big.mark = ",")))
record_check("primary_table_row_count", nrow(prim) == EXP_GRIDS * length(INDICES),
             sprintf("%s rows = 283 x 6", format(nrow(prim), big.mark = ",")))
record_check("maps_wide_row_count", nrow(wide) == EXP_GRIDS,
             sprintf("%d rows = one per grid", nrow(wide)))
record_check("maps_have_coordinates",
             !anyNA(maps_long$lat) && !anyNA(maps_long$lon) &&
               !anyNA(wide$lat) && !anyNA(wide$lon),
             "all map products carry complete lat/lon")

## ---------------------------------------------------------------------------
## 10. WRITE OUTPUTS + READ-BACK
## ---------------------------------------------------------------------------
log_head("10. WRITE OUTPUTS")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
out_tr <- copy(tr)[, trend_class := as.character(trend_class)]
arrow::write_parquet(out_tr, F_ALL, compression = "snappy")
fwrite(copy(prim)[, trend_class := as.character(trend_class)], F_PRIM)
fwrite(copy(maps_long)[, trend_class := as.character(trend_class)], F_MAPL)
fwrite(wide,      F_MAPW)
fwrite(summ,      F_SUMM)
fwrite(sens_summ, F_SENS)
fwrite(wsdi_zi,   F_WSDI)

for (f in c(F_ALL, F_PRIM, F_MAPL, F_MAPW, F_SUMM, F_SENS, F_WSDI))
  log_msg("Saved            : ", basename(f), sprintf("  (%.2f MB)", file.size(f) / 1024^2))

rb <- as.data.table(arrow::read_parquet(F_ALL))
record_check("readback_all_variants",
             nrow(rb) == EXP_TR_ROWS && identical(names(rb), names(out_tr)) &&
               isTRUE(all.equal(rb$sen_per_decade, out_tr$sen_per_decade)) &&
               isTRUE(all.equal(rb$q, out_tr$q)),
             sprintf("%s rows, %d cols, slopes and q identical after round-trip",
                     format(nrow(rb), big.mark = ","), ncol(rb)))
record_check("readback_primary_preserved",
             rb[period == PRIMARY_PERIOD & variant == PRIMARY_VARIANT, .N] ==
               EXP_GRIDS * length(INDICES),
             sprintf("%d primary rows survived the round-trip",
                     rb[period == PRIMARY_PERIOD & variant == PRIMARY_VARIANT, .N]))
rm(rb, out_tr)

.mtime_after <- file.mtime(IN_FILE)
record_check("input_unchanged", identical(.mtime_before, .mtime_after),
             "STEP 06 parquet opened read-only; modification time unchanged")

qc <- rbindlist(CHECKS)
qc[, `:=`(run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
          script = "07_final_trend_analysis.R")]
qc <- rbind(qc, data.table(check = paste0("engine_", tt$test),
                           status = fifelse(tt$ok, "PASS", "FATAL"),
                           detail = sprintf("got %s, want %s | %s", tt$got, tt$want, tt$note),
                           run_time = format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
                           script = "07_final_trend_analysis.R"))
fwrite(qc, F_QC)
log_msg("QC summary       : ", basename(F_QC), "  (", nrow(qc), " entries)")

## ---------------------------------------------------------------------------
## 11. REPORTING NOTES
## ---------------------------------------------------------------------------
log_head("11. REPORTING NOTES (carry these into the manuscript)")
log_msg("PRIMARY result   : ", PRIMARY_PERIOD, ", ", PRIMARY_VARIANT,
        ", BH-FDR q <= ", Q_PRIMARY, ". 1976-2025 is supplementary only.")
log_msg("Maps             : plot sen_per_decade as a continuous diverging fill and")
log_msg("                   mark significance (q <= ", Q_PRIMARY,
        ") by stippling. Do NOT bin magnitude.")
log_msg("Language         : only q <= ", Q_PRIMARY,
        " grids may be called an increase or decrease;")
log_msg("                   unadjudicated signs are a 'non-significant tendency'.")
log_msg("Categories       : 'test_undefined' (Hamed-Rao CF <= 0) and")
log_msg("                   'insufficient_data' are distinct from 'no significant")
log_msg("                   trend' and must be drawn as separate map classes.")
log_msg("WSDI             : report the zero-inflation statistics alongside the map;")
log_msg("                   a near-zero Sen slope field is a property of the index")
log_msg("                   over this season, not a failure of the test.")
log_msg("Method sentence  : Hamed & Rao (1998) variance correction with tie-corrected")
log_msg("                   Var(S), autocorrelation of Sen-detrended ranks screened at")
log_msg("                   +/- z_0.975/sqrt(n) and truncated at lag 3; Theil-Sen slope")
log_msg("                   with Gilbert (1987) CI from the same corrected variance;")
log_msg("                   Benjamini-Hochberg FDR at q <= 0.10 per Wilks (2016),")
log_msg("                   applied separately to each index and period.")

log_head("RUN COMPLETE")
log_msg("Engine tests     : ", sum(tt$ok), " / ", nrow(tt), " passed")
log_msg("Checks passed    : ", qc[status == "PASS", .N], " / ", nrow(qc))
log_msg("Warnings         : ", qc[status == "WARN", .N])
log_msg("Fatal            : ", qc[status == "FATAL", .N])
log_msg("Rows written     : ", format(nrow(tr), big.mark = ","), " (all variants) ; ",
        format(nrow(prim), big.mark = ","), " primary")
log_msg("Elapsed          : ",
        sprintf("%.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(utils::sessionInfo()))
writeLines(.LOG, F_LOG)
cat("\nConsole log written to:", F_LOG, "\n")

invisible(list(trends = tr, primary = prim, summary = summ,
               sensitivity = sens_summ, wsdi = wsdi_zi))
###############################################################################
## END OF SCRIPT
###############################################################################