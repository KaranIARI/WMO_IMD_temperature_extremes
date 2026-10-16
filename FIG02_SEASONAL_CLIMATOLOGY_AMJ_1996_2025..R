###############################################################################
# FIGURE 2  — v2 PUBLICATION LAYOUT
# Seasonal climatology of temperature extremes during the hot weather season
# (April–June), 1996–2025, 283 IMD 1° × 1° grid cells
#
# WHAT CHANGED FROM THE ORIGINAL FIG01 SCRIPT (design only)
#   * Seven separate panels assembled with patchwork (4 + 3) so each panel can
#     carry its own y-axis title ("Temperature (°C)" / "Percentage of days (%)").
#   * Panel titles "(a) TXx (°C)" ... "(g) DTR (°C)".
#   * Navy mean line, white-filled circular markers, subtle light-blue IQR band.
#   * One shared legend in the empty 8th slot (no per-panel legends).
#   * Subtle dashed grey grid, white background, consistent journal font sizes.
#   * Index-specific y-limits computed FROM THE DATA (never truncating any
#     Mean/Q1/Q3 value) with headroom so the value labels never clip.
#     Percentage panels start at 0.
#   * Figure sized for a full-page Scientific Reports width (183 mm).
#   * Exports: 600-dpi PNG, vector PDF, 600-dpi LZW TIFF, plotted-data CSV, QC CSV.
#   * Title "pre-monsoon" replaced by "hot weather season" (manuscript wording).
#
# WHAT DID NOT CHANGE (scientific content)
#   * Same input files, same columns, same index-name cleaning, same
#     7 indices x 3 months, same Mean / Q1 / Q3 columns plotted as-is.
#   * No recalculation, smoothing or interpolation. Lines join the three monthly
#     means exactly as before; the IQR band joins the monthly Q1/Q3 as before.
#   * TX90p / TN10p remain the direct monthly exceedance climatology
#     (locked 1981–2010 thresholds, no bootstrap), exactly as before.
#
# The only path change: F:/ -> H:/ (your current project drive). Edit BASE_DIR
# if needed. Original FIG01 outputs are NOT overwritten (new file names, _v2).
###############################################################################

rm(list = ls())
SCRIPT_VERSION <- "FIG02 v2 (2026-10-15)"
cat("\n>>> Running ", SCRIPT_VERSION, "\n", sep = "")

# =============================================================================
# 1. PACKAGES  — install.packages(c("data.table", "ggplot2", "scales", "patchwork"))
# =============================================================================
for (p in c("data.table", "ggplot2", "scales", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package not installed: ", p, call. = FALSE)
}
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(patchwork)
})

# =============================================================================
# 2. PATHS (same files as the original script; only the drive letter changed)
# =============================================================================
BASE_DIR   <- "H:/WMO_IMD_R/WMO_IMD"
OUTPUT_DIR <- file.path(BASE_DIR, "STEP08")
CLIM_FILE  <- file.path(OUTPUT_DIR, "01_descriptive",
                        "STEP08_01_DESCRIPTIVE_MONTHLY_CLIMATOLOGY_1996_2025.csv")
PERCENTILE_FILE <- file.path(OUTPUT_DIR, "01_descriptive",
                             "STEP08_01_DESCRIPTIVE_MONTHLY_TX90p_TN10p_DIRECT_1996_2025.csv")
FIG_DIR <- file.path(OUTPUT_DIR, "12_figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Original plotted-data CSV (if it exists) is used ONLY to prove values are unchanged.
ORIGINAL_PLOTTED_DATA <- file.path(FIG_DIR, "FIG01_SEASONAL_CLIMATOLOGY_AMJ_1996_2025_7INDEX_PLOTTED_DATA.csv")

OUT_BASE <- file.path(FIG_DIR, "FIG02_SEASONAL_CLIMATOLOGY_AMJ_1996_2025_v2")

# =============================================================================
# 3. LOCKED SETTINGS (unchanged) + DESIGN SETTINGS (new)
# =============================================================================
MONTH_ORDER <- c("April", "May", "June")
INDEX_ORDER <- c("TXx", "TXn", "TNx", "TNn", "TX90p", "TN10p", "DTR")
DEG_C <- "\u00b0C"   # written as a Unicode escape so it prints correctly on every system
INDEX_UNITS <- c(TXx = DEG_C, TXn = DEG_C, TNx = DEG_C, TNn = DEG_C,
                 TX90p = "%", TN10p = "%", DTR = DEG_C)
PANEL_LETTER <- setNames(paste0("(", letters[1:7], ")"), INDEX_ORDER)

# ---- design (change freely; none of these affect the data) ------------------
SHOW_TITLE  <- TRUE      # Scientific Reports normally puts the title in the caption:
# set FALSE for the submission version.
FONT_FAMILY <- "sans"    # "sans" = Arial/Helvetica (recommended by Nature journals); "serif" = Times
COL_LINE    <- "#1B2F4B" # dark navy
COL_BAND    <- "#BFD3EA" # light blue
BAND_ALPHA  <- 0.55
COL_GRID    <- "grey88"
BASE_PT     <- 8         # axis text size (pt) at final print size
WIDTH_MM    <- 183       # full-page width
HEIGHT_MM   <- 118
LABEL_PT    <- 7         # value-label size (pt)

# =============================================================================
# 4. READ INPUTS (unchanged logic)
# =============================================================================
for (f in c(CLIM_FILE, PERCENTILE_FILE))
  if (!file.exists(f)) stop("\nInput file not found:\n", f, "\nCheck BASE_DIR.", call. = FALSE)

clim <- fread(CLIM_FILE)
pct  <- fread(PERCENTILE_FILE)

req_clim <- c("index", "month", "unit", "N", "Mean", "Median", "Q1", "Q3")
req_pct  <- c("index", "month", "scope", "unit", "period", "N", "Mean", "Median", "Q1", "Q3")
if (length(setdiff(req_clim, names(clim)))) stop("Climatology file missing: ", paste(setdiff(req_clim, names(clim)), collapse = ", "))
if (length(setdiff(req_pct,  names(pct))))  stop("TX90p/TN10p file missing: ", paste(setdiff(req_pct, names(pct)), collapse = ", "))

clim[, index := sub("_monthly_climatology$", "", as.character(index))]
pct[,  index := sub("_direct_monthly$",      "", as.character(index))]
clim[, month := as.character(month)]
pct[,  month := as.character(month)]

base <- clim[index %in% c("TXx", "TXn", "TNx", "TNn", "DTR") & month %in% MONTH_ORDER,
             .(index, month, unit, N, Mean, Median, Q1, Q3)]
percentile <- pct[index %in% c("TX90p", "TN10p") & month %in% MONTH_ORDER,
                  .(index, month, unit, N, Mean, Median, Q1, Q3)]
fig <- rbind(base, percentile, fill = TRUE)

# =============================================================================
# 5. DATA QC (same checks as the original, now recorded with real outcomes)
# =============================================================================
QC <- list()
qc <- function(check, ok, detail) {
  QC[[length(QC) + 1L]] <<- data.table(check = check, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail)
  cat(sprintf("  %-4s | %-38s | %s\n", if (isTRUE(ok)) "PASS" else "FAIL", check, detail))
  if (!isTRUE(ok)) { fwrite(rbindlist(QC), paste0(OUT_BASE, "_QC.csv")); stop("QC FAILED: ", check, call. = FALSE) }
}
cat("\nQC\n")
qc("seven_indices", setequal(unique(fig$index), INDEX_ORDER), paste(sort(unique(fig$index)), collapse = ", "))
qc("21_rows", nrow(fig) == 21L, sprintf("%d rows = 7 indices x 3 months", nrow(fig)))
qc("three_months_each", all(fig[, uniqueN(month), by = index]$V1 == 3L), "April, May, June for every index")
qc("no_duplicates", !anyDuplicated(fig, by = c("index", "month")), "one row per index x month")
qc("finite_plotted_values", fig[!is.finite(Mean) | !is.finite(Q1) | !is.finite(Q3), .N] == 0L, "Mean, Q1, Q3 finite")
qc("quartile_order", fig[Median < Q1 | Median > Q3 | Q1 > Q3, .N] == 0L, "Q1 <= Median <= Q3")
fig[, expected_unit := INDEX_UNITS[index]]
## compare units ignoring the degree sign's text encoding ("°C" vs "C"; "%" kept)
## "°C", "degC", "deg C" and "C" all mean degrees Celsius; "%" stays "%"
unit_key <- function(u) {
  k <- tolower(gsub("[^A-Za-z%]", "", enc2utf8(as.character(u))))
  sub("^deg", "", k)          # "degc" -> "c"
}
qc("units_match", fig[unit_key(unit) != unit_key(expected_unit), .N] == 0L,
   paste(unique(fig[, paste0(index, "=", unit_key(unit))]), collapse = ", "))

# Keep an untouched copy of exactly what is plotted
fig[, `:=`(index = factor(index, levels = INDEX_ORDER), month = factor(month, levels = MONTH_ORDER))]
setorder(fig, index, month)
fig[, month_num := as.integer(month)]
plotted <- copy(fig)

# Values must equal the original figure's plotted data, if that file exists
if (file.exists(ORIGINAL_PLOTTED_DATA)) {
  orig <- fread(ORIGINAL_PLOTTED_DATA)
  cmp <- merge(plotted[, .(index = as.character(index), month = as.character(month), Mean, Q1, Q3)],
               orig[, .(index = as.character(index), month = as.character(month),
                        Mean_o = Mean, Q1_o = Q1, Q3_o = Q3)], by = c("index", "month"))
  qc("values_identical_to_original_figure",
     nrow(cmp) == 21L && cmp[, max(abs(c(Mean - Mean_o, Q1 - Q1_o, Q3 - Q3_o)))] < 1e-9,
     sprintf("21 Mean/Q1/Q3 triplets compared with %s", basename(ORIGINAL_PLOTTED_DATA)))
} else {
  QC[[length(QC) + 1L]] <- data.table(check = "values_identical_to_original_figure", result = "NOT RUN",
                                      detail = "original plotted-data CSV not found; values come from the same input files")
  cat("  NOTE | original plotted-data CSV not found; comparison skipped\n")
}

# =============================================================================
# 6. Y-AXIS LIMITS PER PANEL (display only; computed from the data, never clip)
#    lower = min(Q1, Mean) - 10% of range ; upper = max(Q3, Mean) + 22% (label room)
#    Percentage panels start at 0. Breaks via scales::pretty_breaks.
# =============================================================================
lims <- fig[, {
  lo <- min(Q1, Mean); hi <- max(Q3, Mean); r <- hi - lo; if (r <= 0) r <- 1
  l <- lo - 0.14 * r; u <- hi + 0.22 * r
  if (INDEX_UNITS[[as.character(index[1])]] == "%") l <- 0
  .(ymin = l, ymax = u)
}, by = index]
chk <- merge(fig[, .(lo = min(Q1, Mean), hi = max(Q3, Mean)), by = index], lims, by = "index")
qc("y_limits_contain_all_data", chk[, all(ymin <= lo & ymax >= hi)], "no Mean/Q1/Q3 value outside its axis range")

# =============================================================================
# 6b. VALUE-LABEL POSITION (display only): each label goes where the mean line
#     is NOT, so labels never sit on the line or get clipped at panel edges.
#       end point (Apr/Jun), line leaves upwards      -> label centred BELOW
#       middle point, line rises to the right           -> label above-LEFT
#       middle point, line falls in from the left       -> label above-RIGHT
#       local minimum                                   -> label centred BELOW
#       peak / flat / end point with line going down    -> label centred ABOVE
# =============================================================================
setorder(fig, index, month_num)
fig[, `:=`(up_right = shift(Mean, -1L) > Mean, up_left = shift(Mean, 1L) > Mean), by = index]
fig[is.na(up_right), up_right := FALSE][is.na(up_left), up_left := FALSE]
fig[, is_end := month_num %in% c(1L, 3L)]
fig[, lab_pos := fifelse((up_right & up_left) | (is_end & (up_right | up_left)), "below",
                         fifelse(up_right & !up_left, "above_left",
                                 fifelse(up_left & !up_right, "above_right", "above")))]
fig[, `:=`(lab_hjust = fifelse(lab_pos == "above_left", 1, fifelse(lab_pos == "above_right", 0, 0.5)),
           lab_vjust = fifelse(lab_pos == "below", 2.1, -1.0),
           lab_nudge = fifelse(lab_pos == "above_left", -0.07, fifelse(lab_pos == "above_right", 0.07, 0)))]

# =============================================================================
# 7. PANEL BUILDER (one ggplot per index)
# =============================================================================
pt_to_mm <- 1 / .pt
theme_panel <- function() {
  theme_classic(base_size = BASE_PT, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", size = BASE_PT + 1.5, hjust = 0,
                                margin = margin(b = 3)),
      plot.title.position = "panel",
      axis.title.y = element_text(size = BASE_PT, colour = "grey15", margin = margin(r = 3)),
      axis.text = element_text(size = BASE_PT - 0.5, colour = "grey15"),
      axis.line = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "pt"),
      panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.25, linetype = "dashed"),
      panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.3),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 6, 4, 4)
    )
}
make_panel <- function(ix) {
  d  <- fig[index == ix]
  lm <- lims[index == ix]
  unit_lab <- INDEX_UNITS[[ix]]
  ylab <- if (unit_lab == "%") "Percentage of days (%)" else paste0("Temperature (", DEG_C, ")")
  ggplot(d, aes(x = month_num)) +
    geom_ribbon(aes(ymin = Q1, ymax = Q3), fill = COL_BAND, alpha = BAND_ALPHA) +   # IQR band
    geom_line(aes(y = Mean), colour = COL_LINE, linewidth = 0.6) +                  # domain mean
    geom_point(aes(y = Mean), shape = 21, fill = "white", colour = COL_LINE,
               size = 1.9, stroke = 0.6) +
    geom_text(aes(x = month_num + lab_nudge, y = Mean, label = sprintf("%.1f", Mean),   # value labels
                  hjust = lab_hjust, vjust = lab_vjust),
              colour = COL_LINE, size = LABEL_PT * pt_to_mm, fontface = "bold", family = FONT_FAMILY) +
    scale_x_continuous(breaks = 1:3, labels = c("Apr", "May", "Jun"),
                       expand = expansion(add = 0.28)) +
    scale_y_continuous(limits = c(lm$ymin, lm$ymax), breaks = pretty_breaks(n = 5),
                       expand = expansion(mult = 0)) +
    labs(title = paste0(PANEL_LETTER[[ix]], " ", ix, " (", unit_lab, ")"), x = NULL, y = ylab) +
    theme_panel()
}
panels <- lapply(INDEX_ORDER, make_panel)

# =============================================================================
# 8. SHARED LEGEND (drawn manually in the empty 8th slot; one legend only)
# =============================================================================
legend_panel <- ggplot() +
  annotate("rect", xmin = 0.08, xmax = 0.26, ymin = 0.33, ymax = 0.45, fill = COL_BAND, alpha = BAND_ALPHA) +
  annotate("segment", x = 0.08, xend = 0.26, y = 0.64, yend = 0.64, colour = COL_LINE, linewidth = 0.6) +
  annotate("point", x = 0.17, y = 0.64, shape = 21, fill = "white", colour = COL_LINE, size = 1.9, stroke = 0.6) +
  annotate("text", x = 0.31, y = 0.64, label = "Domain mean", hjust = 0,
           size = BASE_PT * pt_to_mm, family = FONT_FAMILY) +
  annotate("text", x = 0.31, y = 0.39, label = "Interquartile range\n(across grid-years)",
           hjust = 0, lineheight = 0.95, size = BASE_PT * pt_to_mm, family = FONT_FAMILY) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(4, 4, 4, 4))

# =============================================================================
# 9. ASSEMBLE 4 + 3 (+ legend)
# =============================================================================
FIG <- wrap_plots(c(panels, list(legend_panel)), ncol = 4, nrow = 2, byrow = TRUE)
if (SHOW_TITLE) {
  FIG <- FIG + plot_annotation(
    title = "Seasonal climatology of temperature extremes during the hot weather season",
    subtitle = "Domain mean across 283 IMD 1\u00b0 \u00d7 1\u00b0 grid cells, April\u2013June 1996\u20132025",
    theme = theme(plot.title = element_text(face = "bold", size = BASE_PT + 4, hjust = 0.5, family = FONT_FAMILY),
                  plot.subtitle = element_text(size = BASE_PT + 1, hjust = 0.5, colour = "grey25",
                                               family = FONT_FAMILY, margin = margin(b = 4)),
                  plot.background = element_rect(fill = "white", colour = NA)))
} else {
  FIG <- FIG + plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

# =============================================================================
# 10. EXPORT (PNG 600 dpi, vector PDF, TIFF 600 dpi LZW)
# =============================================================================
ggsave(paste0(OUT_BASE, ".png"), FIG, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600, bg = "white")
pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
ggsave(paste0(OUT_BASE, ".pdf"), FIG, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", device = pdf_dev, bg = "white")
ggsave(paste0(OUT_BASE, ".tiff"), FIG, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
qc("files_written", all(file.exists(paste0(OUT_BASE, c(".png", ".pdf", ".tiff")))), "PNG, PDF and TIFF saved")

# exact plotted values (identical columns to the original plotted-data CSV + y-limits)
pd <- merge(copy(plotted)[, `:=`(index = as.character(index), month = as.character(month), expected_unit = NULL)],
            lims[, .(index = as.character(index), ymin_display = ymin, ymax_display = ymax)], by = "index")
fwrite(pd[order(match(index, INDEX_ORDER), month_num)], paste0(OUT_BASE, "_PLOTTED_DATA.csv"))

QC_TAB <- rbind(rbindlist(QC),
                data.table(check = c("WSDI handling", "TX90p/TN10p handling", "Raw inputs modified", "Shading meaning"),
                           result = c("NOT IN FIGURE", "DIRECT MONTHLY CLIMATOLOGY", "NO", "IQR (not CI/SE)"),
                           detail = c("continuous 91-day AMJ WSDI; no monthly reset (unchanged)",
                                      "locked 1981-2010 calendar-day thresholds; no Zhang bootstrap (unchanged)",
                                      "input CSVs are read only", "25th-75th percentile across grid-years")))
fwrite(QC_TAB, paste0(OUT_BASE, "_QC.csv"))

cat("\n================ FIGURE 2 (v2) DONE ================\n")
cat("QC: ", sum(QC_TAB$result == "PASS"), " PASS, ", sum(QC_TAB$result == "FAIL"), " FAIL\n", sep = "")
cat("Files:\n ", paste0(OUT_BASE, c(".png", ".pdf", ".tiff", "_PLOTTED_DATA.csv", "_QC.csv"), collapse = "\n  "), "\n")
if (interactive()) tryCatch(print(FIG), error = function(e) message("Preview skipped: ", conditionMessage(e)))