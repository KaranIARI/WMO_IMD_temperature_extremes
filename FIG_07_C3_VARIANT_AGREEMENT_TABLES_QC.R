###############################################################################
##  REVISION 2 | REVIEWER COMMENT 3 | SCRIPT 01 of 02
##  01_C3_VARIANT_AGREEMENT_TABLES_QC.R
##
##  PURPOSE
##    Re-express the variant-sensitivity results using ONE consistent measure:
##    the FINAL FDR-based trend classification (BH q <= 0.10) for all eight
##    indices, and evaluate agreement among cells significant under at least
##    one formulation (Reviewer Comment 3).
##
##  WHY THIS IS NEEDED
##    The original TABLE S7 / FIGURE 6 script read two different summaries:
##      * TXx, TXn, TNx, TNn, WSDI, DTR  <- AMJ_STEP07A_variant_comparison.csv
##        (pre-lock diagnostic: agreement in RAW p < 0.05 significance)
##      * TX90p, TN10p                   <- AMJ_STEP07B_sensitivity_variants.csv
##        (agreement in the FINAL FDR trend class)
##    This script uses the final FDR trend class for all eight indices.
##
##  WHAT IT DOES NOT DO
##    No index, trend, p-value, q-value or class is recalculated. It only
##    COUNTS the existing classes in the locked STEP 07 / 07B outputs.
##
##  INPUTS (read-only, never modified)
##    data/AMJ_STEP07_trends_all_variants.parquet   (TXx TXn TNx TNn WSDI DTR)
##    data/AMJ_STEP07B_trends_all_variants.parquet  (TX90p TN10p)
##    optional: draft TABLE_S13_variant_agreement_significant_cells_1996_2025.csv
##
##  OUTPUTS -> STEP08/REVISION2_C3_VARIANT_AGREEMENT/
##    C3_TABLE_S7_REVISED_FDR_CLASS_AGREEMENT_1996_2025.csv   (+ .docx/.xlsx if packages available)
##    C3_TABLE_S13_SIGNIFICANT_CELL_AGREEMENT_1996_2025.csv    (+ .docx/.xlsx if packages available)
##    C3_TABLE_S13_RETAINED_ALL_THREE_BY_INDEX_1996_2025.csv
##    C3_DISCORDANCE_BREAKDOWN_1996_2025.csv
##    C3_GRID_CLASS_AUDIT_1996_2025.csv   (one row per index x grid, 3 classes; used by script 02)
##    C3_QC_REPORT_01.csv  +  C3_LOG_01.txt
##
##  RUN: open in RStudio -> Session > Restart R -> Source (Ctrl+Shift+S).
###############################################################################

rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE, scipen = 999)
SCRIPT_NAME <- "01_C3_VARIANT_AGREEMENT_TABLES_QC.R"
SCRIPT_VERSION <- "v1 (2026-10-15)"
cat("\n>>> Running ", SCRIPT_NAME, "  ", SCRIPT_VERSION, "\n", sep = "")

## ---------------------------------------------------------------------------
## 0. PACKAGES
##    Required : data.table, arrow
##    Optional : officer + flextable (Word tables), openxlsx (Excel tables)
##    install.packages(c("data.table", "arrow", "officer", "flextable", "openxlsx"))
## ---------------------------------------------------------------------------
for (p in c("data.table", "arrow")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Required package '", p, "' is missing. Run: install.packages('", p, "')", call. = FALSE)
}
suppressPackageStartupMessages({ library(data.table); library(arrow) })
HAS_DOCX <- requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)
HAS_XLSX <- requireNamespace("openxlsx", quietly = TRUE)

## ---------------------------------------------------------------------------
## 1. PATHS  (edit BASE_DIR only if your project folder moves)
## ---------------------------------------------------------------------------
BASE_DIR <- "H:/WMO_IMD_R/WMO_IMD"
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "STEP08", "REVISION2_C3_VARIANT_AGREEMENT")

F_S07  <- file.path(DATA_DIR, "AMJ_STEP07_trends_all_variants.parquet")
F_S07B <- file.path(DATA_DIR, "AMJ_STEP07B_trends_all_variants.parquet")
## Draft S13 (optional). Put the draft CSV in OUT_DIR or DATA_DIR; first match is used.
DRAFT_S13_CANDIDATES <- c(
  file.path(OUT_DIR,  "TABLE_S13_variant_agreement_significant_cells_1996_2025.csv"),
  file.path(DATA_DIR, "TABLE_S13_variant_agreement_significant_cells_1996_2025.csv"))

O_S7     <- file.path(OUT_DIR, "C3_TABLE_S7_REVISED_FDR_CLASS_AGREEMENT_1996_2025")
O_S13    <- file.path(OUT_DIR, "C3_TABLE_S13_SIGNIFICANT_CELL_AGREEMENT_1996_2025")
O_S13R   <- file.path(OUT_DIR, "C3_TABLE_S13_RETAINED_ALL_THREE_BY_INDEX_1996_2025.csv")
O_DISC   <- file.path(OUT_DIR, "C3_DISCORDANCE_BREAKDOWN_1996_2025.csv")
O_AUDIT  <- file.path(OUT_DIR, "C3_GRID_CLASS_AUDIT_1996_2025.csv")
O_QC     <- file.path(OUT_DIR, "C3_QC_REPORT_01.csv")
O_LOG    <- file.path(OUT_DIR, "C3_LOG_01.txt")

## ---------------------------------------------------------------------------
## 2. CONSTANTS (verified against the STEP 07 / 07B scripts and the files)
## ---------------------------------------------------------------------------
PERIOD       <- "TREND_1996_2025"          # primary analysis period
PRIMARY      <- "HR_lag3"
ALTERNATIVES <- c(HR_all = "HR_all", MK_original = "MK_original")
VARIANTS     <- c(PRIMARY, unname(ALTERNATIVES))
COMP_LABEL   <- c(HR_all = "HR_lag3_vs_HR_all", MK_original = "HR_lag3_vs_MK_original")
INDEX_ORDER  <- c("TXx", "TXn", "TNx", "TNn", "TX90p", "TN10p", "WSDI", "DTR")
IDX_S07      <- c("TXx", "TXn", "TNx", "TNn", "WSDI", "DTR")
IDX_S07B     <- c("TX90p", "TN10p")
EXP_GRIDS    <- 283L
Q_FDR        <- 0.10
CLASS_LEVELS <- c("significant_increase", "significant_decrease", "significant_zero_slope",
                  "no_significant_trend", "test_undefined", "insufficient_data")
SIG_CLASSES  <- c("significant_increase", "significant_decrease", "significant_zero_slope")
REQ_COLS     <- c("grid_id", "period", "index", "variant", "trend_class",
                  "significant_fdr", "q", "sen_per_decade", "eligible", "testable")

## Reference values from the earlier draft: CHECKS ONLY, never used in calculations.
REF_OVERALL <- data.table(
  index = rep(INDEX_ORDER, 2),
  comparison = rep(unname(COMP_LABEL), each = 8),
  ref_agreement_pct_1dp = c(98.6, 98.6, 97.9, 91.5, 99.6, 95.8, 98.2, 99.6,
                            97.9, 97.2, 97.9, 99.6, 98.9, 98.6, 100.0, 100.0))
REF_POOLED <- data.table(
  comparison = unname(COMP_LABEL),
  ref_n_either = c(90L, 56L), ref_n_same = c(46L, 31L),
  ref_gained = c(35L, 1L), ref_lost = c(9L, 24L), ref_reversals = c(0L, 0L))
REF_RETAINED <- data.table(index = c("TNx", "TN10p"), ref_primary = c(22L, 17L), ref_all3 = c(16L, 14L))

## ---------------------------------------------------------------------------
## 3. LOGGING + QC HELPERS
## ---------------------------------------------------------------------------
.LOG <- character(0)
log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), format = "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(line, "\n", sep = ""); .LOG <<- c(.LOG, line); invisible(line)
}
log_head <- function(txt) { bar <- strrep("-", 78); cat("\n", bar, "\n", txt, "\n", bar, "\n", sep = "")
.LOG <<- c(.LOG, "", bar, txt, bar) }
QC <- list()
qc <- function(check, passed, detail, severity = c("FATAL", "WARN")) {
  severity <- match.arg(severity)
  ok <- length(passed) == 1L && isTRUE(passed)
  status <- if (ok) "PASS" else severity
  QC[[length(QC) + 1L]] <<- data.table(check = check, status = status, detail = detail)
  log_msg(sprintf("  %-5s | %-48s | %s", status, check, detail))
  if (!ok && severity == "FATAL") {
    try({ dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
      fwrite(rbindlist(QC), O_QC); writeLines(.LOG, O_LOG) }, silent = TRUE)
    stop("FATAL QC failure: ", check, " -- ", detail, call. = FALSE)
  }
  invisible(ok)
}
pct <- function(num, den) ifelse(den > 0, 100 * num / den, NA_real_)

## ---------------------------------------------------------------------------
## 4. INPUT CHECKS + CHECKSUMS (locked files are never written)
## ---------------------------------------------------------------------------
log_head("4. INPUT FILES")
for (f in c(F_S07, F_S07B))
  if (!file.exists(f)) stop("Input file not found:\n  ", f,
                            "\nCheck BASE_DIR ('", BASE_DIR, "') and the file name.", call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
qc("output_folder_separate_from_data",
   normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE) !=
     normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE),
   paste0("outputs go to ", OUT_DIR))
md5_before <- tools::md5sum(c(F_S07, F_S07B))
log_msg("Input checksums recorded (re-checked at the end).")

read_trend <- function(f, label) {
  x <- as.data.table(arrow::read_parquet(f))
  miss <- setdiff(REQ_COLS, names(x))
  if (length(miss)) stop("File ", basename(f), " is missing column(s): ", paste(miss, collapse = ", "),
                         "\nPresent: ", paste(names(x), collapse = ", "), call. = FALSE)
  x <- x[, ..REQ_COLS]
  x[, `:=`(grid_id = as.character(grid_id), period = as.character(period),
           index = as.character(index), variant = as.character(variant),
           trend_class = as.character(trend_class), significant_fdr = as.logical(significant_fdr),
           source_file = label)]
  x
}
s07  <- read_trend(F_S07,  "STEP07")
s07b <- read_trend(F_S07B, "STEP07B")
log_msg("STEP07  rows: ", nrow(s07),  " | indices: ", paste(sort(unique(s07$index)), collapse = ", "))
log_msg("STEP07B rows: ", nrow(s07b), " | indices: ", paste(sort(unique(s07b$index)), collapse = ", "))
qc("STEP07_has_expected_indices",  setequal(unique(s07$index),  IDX_S07),
   paste(sort(unique(s07$index)), collapse = ", "))
qc("STEP07B_has_expected_indices", setequal(unique(s07b$index), IDX_S07B),
   paste(sort(unique(s07b$index)), collapse = ", "))
qc("no_index_in_both_files", length(intersect(unique(s07$index), unique(s07b$index))) == 0L,
   "each index comes from exactly one file")

d <- rbind(s07, s07b, use.names = TRUE)
qc("period_present", PERIOD %in% d$period, paste("periods found:", paste(unique(d$period), collapse = ", ")))
d <- d[period == PERIOD]
log_msg("Rows in primary period ", PERIOD, ": ", nrow(d))

## ---------------------------------------------------------------------------
## 5. STRUCTURE AND CLASSIFICATION-RULE CHECKS
## ---------------------------------------------------------------------------
log_head("5. STRUCTURE AND CLASSIFICATION CHECKS")
qc("all_eight_indices_present", setequal(unique(d$index), INDEX_ORDER),
   sprintf("%d indices: %s", uniqueN(d$index), paste(INDEX_ORDER[INDEX_ORDER %in% d$index], collapse = ", ")))
qc("all_three_variants_present", setequal(unique(d$variant), VARIANTS),
   paste(sort(unique(d$variant)), collapse = ", "))
dups <- d[, .N, by = .(index, variant, grid_id)][N > 1L]
qc("no_duplicate_index_variant_grid", nrow(dups) == 0L,
   sprintf("%d duplicated (index, variant, grid_id) records", nrow(dups)))
cnt <- d[, .(n = uniqueN(grid_id)), by = .(index, variant)]
qc("283_grids_every_index_variant", all(cnt$n == EXP_GRIDS) && nrow(cnt) == 8L * 3L,
   sprintf("grid counts per index x variant: min %d, max %d (%d blocks)", min(cnt$n), max(cnt$n), nrow(cnt)))
gsets <- d[, .(g = paste(sort(unique(grid_id)), collapse = ",")), by = .(index, variant)]
qc("same_grid_set_everywhere", uniqueN(gsets$g) == 1L, "identical 283-cell roster in all 24 blocks")
qc("no_missing_trend_class", !anyNA(d$trend_class) && all(nzchar(d$trend_class)),
   sprintf("%d missing classes", sum(is.na(d$trend_class) | !nzchar(d$trend_class))))
bad_lab <- setdiff(unique(d$trend_class), CLASS_LEVELS)
qc("only_known_class_labels", length(bad_lab) == 0L,
   if (length(bad_lab)) paste("unexpected:", paste(bad_lab, collapse = ", ")) else
     paste("labels:", paste(sort(unique(d$trend_class)), collapse = ", ")))
## The stored class must follow the locked STEP 07 rule (q <= 0.10 & sign of slope).
qc("significant_fdr_flag_matches_class",
   d[significant_fdr != (trend_class %in% SIG_CLASSES), .N] == 0L,
   sprintf("%d rows where significant_fdr disagrees with trend_class",
           d[significant_fdr != (trend_class %in% SIG_CLASSES), .N]))
qc("significant_classes_have_q_le_0.10",
   d[trend_class %in% SIG_CLASSES & !(is.finite(q) & q <= Q_FDR), .N] == 0L,
   sprintf("%d significant rows with q missing or > %.2f",
           d[trend_class %in% SIG_CLASSES & !(is.finite(q) & q <= Q_FDR), .N], Q_FDR))
qc("non_significant_classes_have_q_gt_0.10",
   d[trend_class == "no_significant_trend" & is.finite(q) & q <= Q_FDR, .N] == 0L,
   "no 'no_significant_trend' row has q <= 0.10")
qc("direction_matches_slope_sign",
   d[(trend_class == "significant_increase" & !(sen_per_decade > 0)) |
       (trend_class == "significant_decrease" & !(sen_per_decade < 0)) |
       (trend_class == "significant_zero_slope" & sen_per_decade != 0), .N] == 0L,
   "increase -> slope > 0 ; decrease -> slope < 0 ; zero_slope -> slope = 0")
qc("test_undefined_is_untestable",
   d[trend_class == "test_undefined" & (testable == TRUE | is.finite(q)), .N] == 0L,
   sprintf("test_undefined rows: %d (CF <= 0; no p/q, as in STEP 07/07B)", d[trend_class == "test_undefined", .N]))
qc("no_insufficient_data_in_primary_period", d[trend_class == "insufficient_data", .N] == 0L,
   sprintf("%d insufficient_data rows", d[trend_class == "insufficient_data", .N]), severity = "WARN")
undef_tab <- d[trend_class == "test_undefined", .N, by = .(index, variant)][order(index, variant)]
if (nrow(undef_tab))
  log_msg("test_undefined by index x variant: ",
          paste(sprintf("%s/%s=%d", undef_tab$index, undef_tab$variant, undef_tab$N), collapse = "; "))

## ---------------------------------------------------------------------------
## 6. ONE ROW PER INDEX x GRID WITH THE THREE CLASSES (audit table)
##    Classification rules used below:
##      * "same classification" = identical trend_class label
##      * "significant"         = trend_class in SIG_CLASSES (significant_fdr == TRUE)
##      * test_undefined counts as NON-significant and as its own class; it is
##        never dropped, so every comparison uses all 283 cells.
## ---------------------------------------------------------------------------
log_head("6. GRID-LEVEL CLASS AUDIT TABLE")
w <- dcast(d, index + grid_id ~ variant, value.var = "trend_class")
setnames(w, VARIANTS, paste0("class_", VARIANTS))
qc("audit_table_rows", nrow(w) == 8L * EXP_GRIDS && !anyNA(w[, paste0("class_", VARIANTS), with = FALSE]),
   sprintf("%d rows (expected %d), no missing classes", nrow(w), 8L * EXP_GRIDS))
for (v in VARIANTS) w[, (paste0("sig_", v)) := get(paste0("class_", v)) %in% SIG_CLASSES]
w[, index := factor(index, levels = INDEX_ORDER)]
setorder(w, index, grid_id)
w[, index := as.character(index)]

## ---------------------------------------------------------------------------
## 7. TABLE S7 (REVISED): overall agreement over all 283 cells
##    agreement_pct = n_same_class / n_compared * 100, n_compared = 283
## ---------------------------------------------------------------------------
log_head("7. TABLE S7 (REVISED) - FDR-class agreement over all cells")
S7 <- rbindlist(lapply(names(ALTERNATIVES), function(a) {
  pa <- paste0("class_", PRIMARY); pb <- paste0("class_", a)
  w[, .(comparison = COMP_LABEL[[a]],
        n_compared = .N,
        n_same_class = sum(get(pa) == get(pb)),
        n_discordant = sum(get(pa) != get(pb)),
        n_test_undefined_primary = sum(get(pa) == "test_undefined"),
        n_test_undefined_alternative = sum(get(pb) == "test_undefined")),
    by = index]
}))
S7[, `:=`(agreement_pct = round(pct(n_same_class, n_compared), 2),
          discordance_pct = round(pct(n_discordant, n_compared), 2),
          agreement_pct_1dp = round(pct(n_same_class, n_compared), 1),
          classification_basis = "final FDR trend class (BH q <= 0.10), STEP 07/07B")]
S7[, index := factor(index, levels = INDEX_ORDER)][, comparison := factor(comparison, levels = unname(COMP_LABEL))]
setorder(S7, comparison, index)
S7[, `:=`(index = as.character(index), comparison = as.character(comparison))]
qc("S7_has_16_rows", nrow(S7) == 16L, sprintf("%d rows = 8 indices x 2 comparisons", nrow(S7)))
qc("S7_counts_add_up", S7[, all(n_same_class + n_discordant == n_compared & n_compared == EXP_GRIDS)],
   "same + discordant = 283 in every row")
qc("S7_percentages_reproduce_from_counts",
   S7[, all(abs(agreement_pct - 100 * n_same_class / n_compared) < 0.006 &
              abs(agreement_pct + discordance_pct - 100) < 0.011)],
   "agreement% = 100*same/283 and agreement% + discordance% = 100")
for (cc in unname(COMP_LABEL)) {
  z <- S7[comparison == cc]
  log_msg(sprintf("%-24s : %s", cc, paste(sprintf("%s %.1f", z$index, z$agreement_pct_1dp), collapse = " | ")))
  log_msg(sprintf("%-24s   range %.1f - %.1f %%", "", min(z$agreement_pct_1dp), max(z$agreement_pct_1dp)))
}

## discordance breakdown (what kind of disagreement)
DISC <- rbindlist(lapply(names(ALTERNATIVES), function(a) {
  pa <- paste0("class_", PRIMARY); pb <- paste0("class_", a)
  sa <- paste0("sig_", PRIMARY);   sb <- paste0("sig_", a)
  w[, .(comparison = COMP_LABEL[[a]],
        n_discordant = sum(get(pa) != get(pb)),
        lost_sig_to_nonsig = sum(get(sa) & !get(sb)),
        gained_nonsig_to_sig = sum(!get(sa) & get(sb)),
        reversal_increase_decrease = sum((get(pa) == "significant_increase" & get(pb) == "significant_decrease") |
                                           (get(pa) == "significant_decrease" & get(pb) == "significant_increase")),
        other_sig_class_change = sum(get(sa) & get(sb) & get(pa) != get(pb) &
                                       !((get(pa) == "significant_increase" & get(pb) == "significant_decrease") |
                                           (get(pa) == "significant_decrease" & get(pb) == "significant_increase"))),
        nonsig_vs_test_undefined = sum(!get(sa) & !get(sb) & get(pa) != get(pb))),
    by = index]
}))
qc("discordance_breakdown_is_exhaustive",
   DISC[, all(n_discordant == lost_sig_to_nonsig + gained_nonsig_to_sig + reversal_increase_decrease +
                other_sig_class_change + nonsig_vs_test_undefined)],
   "every discordant cell falls in exactly one category")

## ---------------------------------------------------------------------------
## 8. TABLE S13: agreement among cells significant under >= 1 formulation
##    Unit = index-cell combination; each counted once per pairwise comparison.
##    Denominator (n_either) = cells significant under the primary OR the
##    alternative formulation. Same classification = identical trend_class.
## ---------------------------------------------------------------------------
log_head("8. TABLE S13 - significant-cell agreement")
pair_counts <- function(dd, a) {
  pa <- paste0("class_", PRIMARY); pb <- paste0("class_", a)
  sa <- paste0("sig_", PRIMARY);   sb <- paste0("sig_", a)
  dd[, .(n_sig_primary = sum(get(sa)), n_sig_alternative = sum(get(sb)),
         n_either = sum(get(sa) | get(sb)), n_both = sum(get(sa) & get(sb)),
         n_same_class = sum((get(sa) | get(sb)) & get(pa) == get(pb)),
         lost = sum(get(sa) & !get(sb)), gained = sum(!get(sa) & get(sb)),
         reversals = sum((get(pa) == "significant_increase" & get(pb) == "significant_decrease") |
                           (get(pa) == "significant_decrease" & get(pb) == "significant_increase")))]
}
S13 <- rbindlist(lapply(names(ALTERNATIVES), function(a) {
  per <- rbindlist(lapply(INDEX_ORDER, function(ix) cbind(index = ix, pair_counts(w[index == ix], a))))
  all <- cbind(index = "All indices", pair_counts(w, a))
  cbind(comparison = COMP_LABEL[[a]], rbind(per, all))
}))
S13[, pct_same_among_either := round(pct(n_same_class, n_either), 1)]
qc("S13_either_equals_both_plus_lost_plus_gained",
   S13[, all(n_either == n_both + lost + gained)], "n_either = n_both + lost + gained")
qc("S13_pooled_row_equals_sum_of_indices",
   S13[index != "All indices", .(e = sum(n_either), s = sum(n_same_class), l = sum(lost), g = sum(gained)), by = comparison][
     S13[index == "All indices", .(comparison, e2 = n_either, s2 = n_same_class, l2 = lost, g2 = gained)], on = "comparison"][
       , all(e == e2 & s == s2 & l == l2 & g == g2)],
   "pooled counts = sum over the 8 indices (no double counting)")
for (cc in unname(COMP_LABEL)) {
  z <- S13[comparison == cc & index == "All indices"]
  log_msg(sprintf("%-24s : either %d | same %d (%.1f%%) | lost %d | gained %d | reversals %d",
                  cc, z$n_either, z$n_same_class, z$pct_same_among_either, z$lost, z$gained, z$reversals))
}

## retained under all three formulations (per index)
RET <- w[, .(n_sig_primary = sum(sig_HR_lag3),
             n_sig_any_of_three = sum(sig_HR_lag3 | sig_HR_all | sig_MK_original),
             n_primary_retained_all_three = sum(sig_HR_lag3 & sig_HR_all & sig_MK_original),
             n_primary_retained_all_three_same_class = sum(sig_HR_lag3 & sig_HR_all & sig_MK_original &
                                                             class_HR_lag3 == class_HR_all &
                                                             class_HR_lag3 == class_MK_original)), by = index]
RET <- rbind(RET, RET[, .(index = "All indices", n_sig_primary = sum(n_sig_primary),
                          n_sig_any_of_three = sum(n_sig_any_of_three),
                          n_primary_retained_all_three = sum(n_primary_retained_all_three),
                          n_primary_retained_all_three_same_class = sum(n_primary_retained_all_three_same_class))])
log_msg("Primary detections retained under all three: ",
        paste(sprintf("%s %d/%d", RET$index, RET$n_primary_retained_all_three, RET$n_sig_primary), collapse = " | "))

## ---------------------------------------------------------------------------
## 9. VERIFY AGAINST REFERENCE VALUES AND THE DRAFT S13 (report, never force)
## ---------------------------------------------------------------------------
log_head("9. VERIFICATION AGAINST REFERENCE VALUES / DRAFT S13")
ro <- merge(S7[, .(index, comparison, agreement_pct_1dp)], REF_OVERALL, by = c("index", "comparison"))
ro_bad <- ro[abs(agreement_pct_1dp - ref_agreement_pct_1dp) > 0.05]
qc("overall_agreement_matches_reference", nrow(ro_bad) == 0L,
   if (nrow(ro_bad)) paste("DIFFERS:", paste(sprintf("%s %s computed %.1f vs ref %.1f", ro_bad$index,
                                                     ro_bad$comparison, ro_bad$agreement_pct_1dp,
                                                     ro_bad$ref_agreement_pct_1dp), collapse = "; "),
                           "-> check that period, class labels and test_undefined handling match") else
                             "all 16 values reproduce the reference (1 dp)", severity = "WARN")
rp <- merge(S13[index == "All indices", .(comparison, n_either, n_same_class, gained, lost, reversals)],
            REF_POOLED, by = "comparison")
rp_bad <- rp[n_either != ref_n_either | n_same_class != ref_n_same | gained != ref_gained |
               lost != ref_lost | reversals != ref_reversals]
qc("pooled_S13_matches_reference", nrow(rp_bad) == 0L,
   if (nrow(rp_bad)) paste("DIFFERS:", paste(sprintf("%s either %d/%d same %d/%d gained %d/%d lost %d/%d rev %d/%d",
                                                     rp_bad$comparison, rp_bad$n_either, rp_bad$ref_n_either,
                                                     rp_bad$n_same_class, rp_bad$ref_n_same, rp_bad$gained,
                                                     rp_bad$ref_gained, rp_bad$lost, rp_bad$ref_lost,
                                                     rp_bad$reversals, rp_bad$ref_reversals), collapse = "; "),
                           "(computed/reference)") else
                             "90/46/35/9 and 56/31/1/24, no reversals: reproduced", severity = "WARN")
rr <- merge(RET[, .(index, n_sig_primary, n_primary_retained_all_three)], REF_RETAINED, by = "index")
qc("retained_all_three_matches_reference",
   rr[, all(n_sig_primary == ref_primary & n_primary_retained_all_three == ref_all3)],
   paste(sprintf("%s %d/%d (ref %d/%d)", rr$index, rr$n_primary_retained_all_three, rr$n_sig_primary,
                 rr$ref_all3, rr$ref_primary), collapse = "; "), severity = "WARN")
qc("primary_counts_match_manuscript_Table3",
   identical(RET[match(INDEX_ORDER, index), n_sig_primary], c(6L, 5L, 22L, 1L, 3L, 17L, 0L, 1L)),
   paste("primary FDR cells:", paste(sprintf("%s %d", INDEX_ORDER, RET[match(INDEX_ORDER, index), n_sig_primary]),
                                     collapse = ", ")), severity = "WARN")

draft_file <- DRAFT_S13_CANDIDATES[file.exists(DRAFT_S13_CANDIDATES)][1]
if (!is.na(draft_file)) {
  dr <- fread(draft_file, check.names = FALSE)
  need <- c("Index", "HR_all: significant in either", "HR_all: same class among these (n)",
            "HR_all: lost (primary only)", "HR_all: gained (HR_all only)", "HR_all: opposite direction",
            "MK: significant in either", "MK: same class among these (n)", "MK: lost (primary only)",
            "MK: gained (MK only)", "MK: opposite direction", "Significant in all 3")
  if (all(need %in% names(dr))) {
    mine <- merge(S13[comparison == COMP_LABEL[["HR_all"]], .(index, e1 = n_either, s1 = n_same_class, l1 = lost, g1 = gained, r1 = reversals)],
                  S13[comparison == COMP_LABEL[["MK_original"]], .(index, e2 = n_either, s2 = n_same_class, l2 = lost, g2 = gained, r2 = reversals)],
                  by = "index")
    mine <- merge(mine, RET[, .(index, a3 = n_sig_any_of_three, r3 = n_primary_retained_all_three)], by = "index")
    cmpd <- merge(mine, dr[, .(index = Index,
                               de1 = get("HR_all: significant in either"), ds1 = get("HR_all: same class among these (n)"),
                               dl1 = get("HR_all: lost (primary only)"), dg1 = get("HR_all: gained (HR_all only)"),
                               dr1 = get("HR_all: opposite direction"),
                               de2 = get("MK: significant in either"), ds2 = get("MK: same class among these (n)"),
                               dl2 = get("MK: lost (primary only)"), dg2 = get("MK: gained (MK only)"),
                               dr2 = get("MK: opposite direction"), dall3 = get("Significant in all 3"))],
                  by = "index")
    cmpd[, mismatch := !(e1 == de1 & s1 == ds1 & l1 == dl1 & g1 == dg1 & r1 == dr1 &
                           e2 == de2 & s2 == ds2 & l2 == dl2 & g2 == dg2 & r2 == dr2 & r3 == dall3)]
    qc("draft_S13_counts_reproduced", nrow(cmpd) == 9L && !any(cmpd$mismatch),
       if (any(cmpd$mismatch)) paste("MISMATCH for:", paste(cmpd[mismatch == TRUE, index], collapse = ", "),
                                     "-> see C3_QC_REPORT_01.csv; check class labels / period") else
                                       sprintf("all %d rows of %s reproduced exactly", nrow(cmpd), basename(draft_file)), severity = "WARN")
  } else {
    qc("draft_S13_columns_recognised", FALSE,
       paste("draft found but columns differ; missing:", paste(setdiff(need, names(dr)), collapse = ", ")),
       severity = "WARN")
  }
} else {
  qc("draft_S13_available", FALSE, "draft S13 CSV not found in OUT_DIR or DATA_DIR; comparison skipped",
     severity = "WARN")
}

## ---------------------------------------------------------------------------
## 10. WRITE OUTPUTS (CSV always; Word / Excel when packages are installed)
## ---------------------------------------------------------------------------
log_head("10. WRITE OUTPUTS")
S7_out <- S7[, .(Index = index, Comparison = comparison, `Cells compared` = n_compared,
                 `Same classification` = n_same_class, Discordant = n_discordant,
                 `Agreement (%)` = agreement_pct_1dp, `Discordance (%)` = round(discordance_pct, 1),
                 `test_undefined (primary)` = n_test_undefined_primary,
                 `test_undefined (alternative)` = n_test_undefined_alternative,
                 agreement_pct_2dp = agreement_pct, classification_basis)]
S13_out <- S13[, .(Comparison = comparison, Index = index,
                   `Significant, primary` = n_sig_primary, `Significant, alternative` = n_sig_alternative,
                   `Significant in either (denominator)` = n_either, `Significant in both` = n_both,
                   `Same classification` = n_same_class, `Same classification (%)` = pct_same_among_either,
                   `Lost (primary only)` = lost, `Gained (alternative only)` = gained,
                   `Direction reversals` = reversals)]
RET_out <- RET[, .(Index = index, `Significant, primary (HR lag-3)` = n_sig_primary,
                   `Significant under >= 1 of 3` = n_sig_any_of_three,
                   `Primary detections retained under all 3` = n_primary_retained_all_three,
                   `Retained with identical class` = n_primary_retained_all_three_same_class)]
fwrite(S7_out, paste0(O_S7, ".csv"))
fwrite(S13_out, paste0(O_S13, ".csv"))
fwrite(RET_out, O_S13R)
fwrite(DISC, O_DISC)
fwrite(w, O_AUDIT)

if (HAS_DOCX) {
  make_docx <- function(tabs, titles, notes, file) {
    doc <- officer::read_docx()
    for (k in seq_along(tabs)) {
      doc <- officer::body_add_par(doc, titles[k], style = "heading 2")
      ft <- flextable::autofit(flextable::theme_booktabs(flextable::flextable(as.data.frame(tabs[[k]]))))
      doc <- flextable::body_add_flextable(doc, ft)
      doc <- officer::body_add_par(doc, notes[k], style = "Normal")
    }
    print(doc, target = file)
  }
  make_docx(list(S7_out[, 1:9]), "Table S7 (revised)",
            "Agreement in final FDR-based trend classification (q <= 0.10) over all 283 grid cells, AMJ 1996-2025.",
            paste0(O_S7, ".docx"))
  make_docx(list(S13_out, RET_out), c("Table S13a", "Table S13b"),
            c("Agreement among index-cell combinations significant under at least one formulation.",
              "Primary detections retained under all three formulations."),
            paste0(O_S13, ".docx"))
  log_msg("Word tables written (officer + flextable).")
} else log_msg("NOTE: officer/flextable not installed -> Word tables skipped (CSV written).")
if (HAS_XLSX) {
  openxlsx::write.xlsx(list(Table_S7_revised = S7_out), paste0(O_S7, ".xlsx"), overwrite = TRUE)
  openxlsx::write.xlsx(list(Table_S13a = S13_out, Table_S13b = RET_out, Discordance = DISC),
                       paste0(O_S13, ".xlsx"), overwrite = TRUE)
  log_msg("Excel tables written (openxlsx).")
} else log_msg("NOTE: openxlsx not installed -> Excel tables skipped (CSV written).")

md5_after <- tools::md5sum(c(F_S07, F_S07B))
qc("locked_inputs_unchanged", identical(unname(md5_before), unname(md5_after)),
   "MD5 of both parquet files identical before and after")

qct <- rbindlist(QC)
qct[, `:=`(script = SCRIPT_NAME, version = SCRIPT_VERSION, run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))]
fwrite(qct, O_QC)

log_head("SUMMARY - SCRIPT 01")
log_msg("QC: ", qct[status == "PASS", .N], " PASS | ", qct[status == "WARN", .N], " WARN | ",
        qct[status == "FATAL", .N], " FATAL")
if (qct[status == "WARN", .N]) log_msg("WARN checks: ", paste(qct[status == "WARN", check], collapse = ", "))
for (f in list.files(OUT_DIR, pattern = "^C3_", full.names = TRUE)) log_msg("  wrote ", basename(f))
log_msg("NEXT: run 02_C3_FIGURE6_REVISED.R")
.LOG <- c(.LOG, "", "--- sessionInfo() ---", capture.output(sessionInfo()))
writeLines(.LOG, O_LOG)
cat("\nLog written to: ", O_LOG, "\n", sep = "")






###############################################################################
##  REVISION 2 | REVIEWER COMMENT 3 | SCRIPT 02 of 02
##  02_C3_FIGURE6_REVISED.R
##
##  PURPOSE
##    Redraw Figure 6 (variant agreement) from the SAME FDR-class data used for
##    the revised Table S7, keeping the original FIG6 layout (grey bars, two
##    panels, 0-105 % axis, value labels).
##
##  INPUTS  (written by script 01; run 01 first)
##    STEP08/REVISION2_C3_VARIANT_AGREEMENT/C3_TABLE_S7_REVISED_FDR_CLASS_AGREEMENT_1996_2025.csv
##    STEP08/REVISION2_C3_VARIANT_AGREEMENT/C3_GRID_CLASS_AUDIT_1996_2025.csv
##
##  OUTPUTS -> same folder
##    C3_FIGURE_6_REVISED_VARIANT_AGREEMENT_AMJ_1996_2025_600dpi.png
##    C3_FIGURE_6_REVISED_VARIANT_AGREEMENT_AMJ_1996_2025.pdf
##    C3_FIGURE_6_REVISED_DATA.csv
##    C3_QC_REPORT_02.csv
##
##  The original STEP08 Figure 6 files are NOT overwritten (new file names).
##  RUN: Session > Restart R -> Source (Ctrl+Shift+S).
##  Packages: install.packages(c("data.table", "ggplot2"))
###############################################################################

rm(list = ls())
options(warn = 1, stringsAsFactors = FALSE)
SCRIPT_NAME <- "02_C3_FIGURE6_REVISED.R"; SCRIPT_VERSION <- "v1 (2026-10-15)"
cat("\n>>> Running ", SCRIPT_NAME, "  ", SCRIPT_VERSION, "\n", sep = "")
for (p in c("data.table", "ggplot2"))
  if (!requireNamespace(p, quietly = TRUE))
    stop("Required package '", p, "' is missing. Run: install.packages('", p, "')", call. = FALSE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

## ---------------------------------------------------------------------------
## 1. PATHS
## ---------------------------------------------------------------------------
BASE_DIR <- "H:/WMO_IMD_R/WMO_IMD"
OUT_DIR  <- file.path(BASE_DIR, "STEP08", "REVISION2_C3_VARIANT_AGREEMENT")
F_S7     <- file.path(OUT_DIR, "C3_TABLE_S7_REVISED_FDR_CLASS_AGREEMENT_1996_2025.csv")
F_AUDIT  <- file.path(OUT_DIR, "C3_GRID_CLASS_AUDIT_1996_2025.csv")
O_PNG    <- file.path(OUT_DIR, "C3_FIGURE_6_REVISED_VARIANT_AGREEMENT_AMJ_1996_2025_600dpi.png")
O_PDF    <- file.path(OUT_DIR, "C3_FIGURE_6_REVISED_VARIANT_AGREEMENT_AMJ_1996_2025.pdf")
O_DATA   <- file.path(OUT_DIR, "C3_FIGURE_6_REVISED_DATA.csv")
O_QC     <- file.path(OUT_DIR, "C3_QC_REPORT_02.csv")

INDEX_ORDER <- c("TXx", "TXn", "TNx", "TNn", "TX90p", "TN10p", "WSDI", "DTR")
COMP_ORDER  <- c("HR_lag3_vs_HR_all", "HR_lag3_vs_MK_original")
COMP_NICE   <- c(HR_lag3_vs_HR_all = "HR lag-3 vs HR all-lag",
                 HR_lag3_vs_MK_original = "HR lag-3 vs original MK")

QC <- list()
qc <- function(check, passed, detail) {
  ok <- length(passed) == 1L && isTRUE(passed)
  QC[[length(QC) + 1L]] <<- data.table(check = check, status = if (ok) "PASS" else "FATAL", detail = detail)
  cat(sprintf("  %-5s | %-44s | %s\n", if (ok) "PASS" else "FATAL", check, detail))
  if (!ok) { try(fwrite(rbindlist(QC), O_QC), silent = TRUE)
    stop("FATAL QC failure: ", check, " -- ", detail, call. = FALSE) }
}

## ---------------------------------------------------------------------------
## 2. READ SCRIPT 01 OUTPUTS
## ---------------------------------------------------------------------------
for (f in c(F_S7, F_AUDIT))
  if (!file.exists(f)) stop("File not found (run script 01 first):\n  ", f, call. = FALSE)
s7 <- fread(F_S7, check.names = FALSE)
au <- fread(F_AUDIT, colClasses = list(character = "grid_id"))
need_s7 <- c("Index", "Comparison", "Cells compared", "Same classification", "Agreement (%)", "agreement_pct_2dp")
need_au <- c("index", "grid_id", "class_HR_lag3", "class_HR_all", "class_MK_original")
if (length(setdiff(need_s7, names(s7)))) stop("Table S7 file lacks: ", paste(setdiff(need_s7, names(s7)), collapse = ", "))
if (length(setdiff(need_au, names(au)))) stop("Audit file lacks: ", paste(setdiff(need_au, names(au)), collapse = ", "))

## ---------------------------------------------------------------------------
## 3. QC: FIGURE DATA = TABLE S7 = RECOUNT FROM THE GRID-LEVEL CLASSES
## ---------------------------------------------------------------------------
cat("\nQC checks\n")
qc("S7_has_16_rows_8_indices_2_comparisons",
   nrow(s7) == 16L && setequal(s7$Index, INDEX_ORDER) && setequal(s7$Comparison, COMP_ORDER),
   sprintf("%d rows", nrow(s7)))
qc("audit_has_8x283_rows", nrow(au) == 8L * 283L && all(au[, .N, by = index]$N == 283L),
   sprintf("%d rows", nrow(au)))
recount <- rbind(
  au[, .(Comparison = "HR_lag3_vs_HR_all", n = .N, same = sum(class_HR_lag3 == class_HR_all)), by = .(Index = index)],
  au[, .(Comparison = "HR_lag3_vs_MK_original", n = .N, same = sum(class_HR_lag3 == class_MK_original)), by = .(Index = index)])
chk <- merge(s7[, .(Index, Comparison, n_s7 = `Cells compared`, same_s7 = `Same classification`,
                    pct2 = agreement_pct_2dp, pct1 = `Agreement (%)`)], recount, by = c("Index", "Comparison"))
qc("S7_counts_equal_recount_from_grid_classes",
   nrow(chk) == 16L && chk[, all(n_s7 == n & same_s7 == same)],
   "Table S7 and Figure 6 derive from the same FDR classes")
qc("S7_percent_reproduces_from_counts",
   chk[, all(abs(pct2 - round(100 * same / n, 2)) < 1e-9 & abs(pct1 - round(100 * same / n, 1)) < 1e-9)],
   "agreement % = 100 x same / 283")

## ---------------------------------------------------------------------------
## 4. FIGURE 6 (same layout as the original FIG6 script)
## ---------------------------------------------------------------------------
plot_data <- s7[, .(index = factor(Index, levels = INDEX_ORDER),
                    comparison = Comparison,
                    comparison_label = factor(COMP_NICE[Comparison], levels = COMP_NICE),
                    n_compared = `Cells compared`, n_same_class = `Same classification`,
                    agreement_pct = `Agreement (%)`)]
p6 <- ggplot(plot_data, aes(x = index, y = agreement_pct)) +
  geom_col(width = 0.65, fill = "grey35", colour = "black") +
  geom_text(aes(label = sprintf("%.1f", agreement_pct)), vjust = -0.35, size = 3.1) +
  facet_wrap(~ comparison_label, nrow = 1) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  labs(title = "Robustness of grid-level trend classification",
       subtitle = "AMJ 1996\u20132025 |  FDR-based classes (q \u2264 0.10), all 283 grid cells",
       x = NULL, y = "Classification agreement (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(colour = "black", size = 9.5),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold", size = 11),
        panel.grid = element_blank(), legend.position = "none")

ggsave(O_PNG, p6, width = 10, height = 6, units = "in", dpi = 600, bg = "white")
pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
ggsave(O_PDF, p6, width = 10, height = 6, units = "in", device = pdf_dev, bg = "white")
fwrite(plot_data, O_DATA)
qc("figure_files_written", file.exists(O_PNG) && file.exists(O_PDF), "PNG (600 dpi) and PDF saved")

qct <- rbindlist(QC)[, `:=`(script = SCRIPT_NAME, version = SCRIPT_VERSION)]
fwrite(qct, O_QC)
cat("\n================ SUMMARY - SCRIPT 02 ================\n")
cat("QC: ", qct[status == "PASS", .N], " PASS | ", qct[status == "FATAL", .N], " FATAL\n", sep = "")
for (cc in COMP_ORDER) {
  z <- plot_data[comparison == cc]
  cat(sprintf("%-26s range %.1f - %.1f %%\n", COMP_NICE[cc], min(z$agreement_pct), max(z$agreement_pct)))
}
cat("Files:\n  ", O_PNG, "\n  ", O_PDF, "\n  ", O_DATA, "\n  ", O_QC, "\n", sep = "")
if (interactive()) tryCatch(print(p6), error = function(e) message("Preview skipped: ", conditionMessage(e)))



R.version.string
