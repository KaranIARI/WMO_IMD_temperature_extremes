# =====================================================================
# SUPPLEMENTARY FIGURE S1 — SPATIAL TREND MAPS, AMJ 1976–2025
# Longer-record sensitivity analysis (Reviewer comment: report 1976–2025)
#
# 283 IMD 1° x 1° land grid cells | 8 indices
#   TXx, TXn, TNx, TNn, TX90p, TN10p, WSDI, DTR
#
# Same style and class rule as the main spatial-trend figure (old Fig. 3),
# only the PERIOD differs: FULL_1976_2025 instead of TREND_1996_2025.
#
# Trend        : locked Theil–Sen slope per decade (STEP 07 / 07B)
# Significance : locked BH-FDR q <= 0.10, FDR family = index x period
# Variant      : HR_lag3 (primary test)
#
# NO trend, p, q or FDR value is calculated, refitted or adjusted here.
# The map class is DERIVED only from two locked fields:
#     sign(sen_per_decade)  and  significant_fdr
# and cross-checked against the locked trend_class in the QC block.
#
# INPUTS (read-only)
#   AMJ_STEP07_maps_long.csv   : TXx, TXn, TNx, TNn, WSDI, DTR (both periods)
#   AMJ_STEP07B_maps_long.csv  : TX90p, TN10p (both periods)
# =====================================================================

rm(list = ls())
SCRIPT_VERSION <- "v3 (single shared legend strip; no guides = 'collect')"
cat("\n>>> FIG S1 SCRIPT VERSION:", SCRIPT_VERSION, "<<<\n")

# =====================================================================
# 1. PACKAGES
# =====================================================================
for (p in c("data.table", "sf", "ggplot2", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package not installed: ", p, call. = FALSE)
}
suppressPackageStartupMessages({
  library(data.table); library(sf); library(ggplot2); library(patchwork)
})

# =====================================================================
# 2. PATHS (edit here only; use forward slashes on Windows)
# =====================================================================
BASE_DIR <- "H:/WMO_IMD_R/WMO_IMD"

FILE_07  <- file.path(BASE_DIR, "data", "AMJ_STEP07_maps_long.csv")
FILE_07B <- file.path(BASE_DIR, "data", "AMJ_STEP07B_maps_long.csv")

# Same SOI boundary as the new Figure 1 (verified to reproduce the 283 cells)
SHAPEFILE <- file.path(BASE_DIR, "shapefile/SOI/State_District_Subdistrict_PAN INDIA",
                       "State_District_Subdistrict_PAN INDIA/State Boundary/State Boundary.shp")

OUT_DIR <- file.path(BASE_DIR, "STEP08/00_Supp_final_maps")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(FILE_07, FILE_07B, SHAPEFILE))
  if (!file.exists(f)) stop("Input not found: ", f, "\nEdit the path in Section 2.", call. = FALSE)

# =====================================================================
# 3. LOCKED SETTINGS
# =====================================================================
PERIOD_CODE     <- "FULL_1976_2025"
START_YEAR      <- 1976L
END_YEAR        <- 2025L
PRIMARY_VARIANT <- "HR_lag3"
FDR_Q           <- 0.10
EXPECTED_GRIDS  <- 283L

INDEX_ORDER <- c("TXx", "TXn", "TNx", "TNn", "TX90p", "TN10p", "WSDI", "DTR")
EXPECTED_UNITS <- c(TXx = "degC/decade", TXn = "degC/decade", TNx = "degC/decade",
                    TNn = "degC/decade", TX90p = "%/decade", TN10p = "%/decade",
                    WSDI = "days/decade", DTR = "degC/decade")

CLASS_LEVELS <- c("Significant decrease",
                  "Decreasing tendency (not significant)",
                  "No trend (zero slope)",
                  "Increasing tendency (not significant)",
                  "Significant increase")
CLASS_COLOURS <- c("Significant decrease"                  = "#0571B0",
                   "Decreasing tendency (not significant)" = "#92C5DE",
                   "No trend (zero slope)"                 = "#FFFFFF",
                   "Increasing tendency (not significant)" = "#F4A582",
                   "Significant increase"                  = "#B2182B")

COL_BOUNDARY  <- "grey10"; LW_BOUNDARY  <- 0.42
COL_CELL_EDGE <- "grey55"; LW_CELL_EDGE <- 0.08

MAP_XLIM <- c(66.5, 98.5)
MAP_YLIM <- c(6.0,  39.6)
TITLE_Y  <- 38.55
NORTH_X  <- 96.30
NORTH_Y  <- 35.40

# =====================================================================
# 4. READ LOCKED STEP 07 / 07B RESULTS
# =====================================================================
cat("\n=====================================================\n")
cat("READING LOCKED TREND RESULTS (", PERIOD_CODE, ")\n", sep = "")
cat("=====================================================\n")

need <- c("period", "index", "grid_id", "lat", "lon", "variant", "units",
          "sen_per_decade", "p", "q", "trend_class", "significant_fdr")

read_part <- function(f) {
  x <- fread(f)
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("Missing columns in ", basename(f), ": ",
                         paste(miss, collapse = ", "), call. = FALSE)
  x[, ..need]
}
d <- rbind(read_part(FILE_07), read_part(FILE_07B))
d <- d[period == PERIOD_CODE & variant == PRIMARY_VARIANT & index %in% INDEX_ORDER]
d[, significant_fdr := as.logical(significant_fdr)]
d[, testable := is.finite(p)]
cat("Rows after filtering: ", nrow(d), " (expected ", 8 * EXPECTED_GRIDS, ")\n", sep = "")

# =====================================================================
# 5. DERIVE MAP CLASS (no statistic recomputed) — identical rule to Fig. 3
# =====================================================================
d[, sen_sign := fifelse(!is.finite(sen_per_decade), NA_integer_,
                        fifelse(sen_per_decade > 0, 1L,
                                fifelse(sen_per_decade < 0, -1L, 0L)))]
d[, map_class := fcase(
  trend_class == "insufficient_data", NA_character_,
  trend_class == "test_undefined",    NA_character_,
  is.na(sen_sign),                    NA_character_,
  sen_sign ==  0L,                    "No trend (zero slope)",
  sen_sign == -1L &  significant_fdr, "Significant decrease",
  sen_sign == -1L & !significant_fdr, "Decreasing tendency (not significant)",
  sen_sign ==  1L &  significant_fdr, "Significant increase",
  sen_sign ==  1L & !significant_fdr, "Increasing tendency (not significant)")]
d[, map_class := factor(map_class, levels = CLASS_LEVELS)]

# =====================================================================
# 6. DATA QC (stops on any failure)
# =====================================================================
qc_stop <- function(ok, msg) if (!isTRUE(ok)) stop("QC FAILED: ", msg, call. = FALSE)

cnt <- d[, .(n = .N, n_grid = uniqueN(grid_id)), by = index]
qc_stop(setequal(cnt$index, INDEX_ORDER), "the 8 indices are not all present.")
qc_stop(all(cnt$n == EXPECTED_GRIDS & cnt$n_grid == EXPECTED_GRIDS),
        "every index must contain exactly 283 unique grid cells.")
qc_stop(d[, .N, by = .(index, grid_id)][N > 1, .N] == 0L, "duplicate index x grid rows.")
qc_stop(all(d[, unique(units), by = index][, V1 == EXPECTED_UNITS[index]]),
        "units do not match the expected unit for each index.")
qc_stop(d[testable == TRUE, all(significant_fdr == (q <= FDR_Q))],
        "significant_fdr disagrees with q <= 0.10.")
qc_stop(d[trend_class == "significant_increase" & map_class != "Significant increase", .N] +
          d[trend_class == "significant_decrease" & map_class != "Significant decrease", .N] == 0L,
        "derived map_class contradicts the locked trend_class.")
qc_stop(d[map_class %in% c("Significant increase", "Significant decrease") &
            significant_fdr != TRUE, .N] == 0L,
        "a cell would be coloured significant without significant_fdr.")
qc_stop(d[is.na(map_class) & !trend_class %in% c("test_undefined", "insufficient_data"), .N] == 0L,
        "unclassified cells without a declared reason.")

class_tab <- dcast(d[, .N, by = .(index, map_class)], index ~ map_class,
                   value.var = "N", fill = 0L)
cat("\nClass counts per index (1976–2025):\n"); print(class_tab)

cat("\nsignificant_zero_slope cells (q <= 0.10, slope exactly 0; shown as zero-slope class):\n")
print(d[trend_class == "significant_zero_slope", .N, by = index])
cat("\nNon-plottable cells (test_undefined / insufficient_data):\n")
print(d[trend_class %in% c("test_undefined", "insufficient_data"), .N, by = .(index, trend_class)])
cat("\nData QC: PASS\n")

# =====================================================================
# 7. STATE BOUNDARY
# =====================================================================
bnd <- st_read(SHAPEFILE, quiet = TRUE)
if (is.na(st_crs(bnd))) stop("QC FAILED: boundary has no CRS.", call. = FALSE)
if (!all(st_is_valid(bnd))) bnd <- st_make_valid(bnd)
bnd <- st_transform(bnd, 4326)

# =====================================================================
# 8. THEME, NORTH ARROW, MAP BUILDER (same as Fig. 3)
# =====================================================================
theme_map <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.65),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key       = element_rect(fill = "white", colour = NA),
          legend.text      = element_text(size = base_size - 1),
          plot.margin      = margin(4, 4, 4, 4))
}

north_arrow <- function(p, x = NORTH_X, y = NORTH_Y, h = 2.10, w = 0.62, lab_size = 4.4) {
  p +
    annotate("text", x = x, y = y + h + 0.72, label = "N",
             fontface = "bold", size = lab_size, colour = "black") +
    annotate("polygon", x = c(x, x - w, x, x + w), y = c(y + h, y, y + 0.30 * h, y),
             fill = "black", colour = "black", linewidth = 0.30)
}

make_map <- function(dd, index_name, panel_label) {
  undefined    <- dd[trend_class == "test_undefined"]
  insufficient <- dd[trend_class == "insufficient_data"]
  plottable    <- dd[!is.na(map_class)]
  
  p <- ggplot() +
    geom_tile(data = plottable, aes(x = lon, y = lat, fill = map_class),
              width = 1, height = 1, colour = COL_CELL_EDGE, linewidth = LW_CELL_EDGE) +
    scale_fill_manual(values = CLASS_COLOURS, limits = CLASS_LEVELS, drop = FALSE,
                      na.value = "white", guide = "none")
  
  if (nrow(insufficient) > 0L)
    p <- p + geom_tile(data = insufficient, aes(x = lon, y = lat), width = 1, height = 1,
                       fill = "white", colour = "grey25", linewidth = 0.38, linetype = "dotted")
  if (nrow(undefined) > 0L)
    p <- p + geom_tile(data = undefined, aes(x = lon, y = lat), width = 1, height = 1,
                       fill = NA, colour = "black", linewidth = 0.42, linetype = "solid")
  
  p <- p +
    geom_sf(data = bnd, inherit.aes = FALSE, fill = NA,
            colour = COL_BOUNDARY, linewidth = LW_BOUNDARY) +
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE, datum = NA) +
    theme_map(11) +
    annotate("text", x = mean(MAP_XLIM), y = TITLE_Y,
             label = paste0(index_name, " Trend, AMJ ", START_YEAR, "\u2013", END_YEAR),
             fontface = "bold", size = 4.0, colour = "black") +
    annotate("text", x = MAP_XLIM[1] + 0.55, y = MAP_YLIM[2] - 0.55,
             label = panel_label, hjust = 0, vjust = 1,
             fontface = "bold", size = 4.4, colour = "black")
  north_arrow(p)
}

# =====================================================================
# 9. BUILD 2 x 4 COMPOSITE
#    Row 1 : TXx | TXn | TNx | TNn      Row 2 : TX90p | TN10p | WSDI | DTR
# =====================================================================
panels <- lapply(seq_along(INDEX_ORDER), function(j)
  make_map(d[index == INDEX_ORDER[j]], INDEX_ORDER[j], LETTERS[j]))

n_undef <- d[trend_class == "test_undefined", .N]
FOOT <- paste0(
  "Sensitivity analysis for the longer record (n = 50 years). Identical method to the primary ",
  "1996\u20132025 analysis: Hamed\u2013Rao modified Mann\u2013Kendall test (lag 3), Theil\u2013Sen slope, ",
  "\u226580% completeness, BH-FDR q \u2264 0.10 applied separately for each index and period.",
  if (n_undef > 0L) paste0("\nSolid black outline: test undefined (non-positive Hamed\u2013Rao variance; ",
                           n_undef, " cells). ") else "\n",
  "Cells significant with an exactly zero slope are shown in the zero-slope class. ",
  "The 1976\u20132025 period contains the 1996\u20132025 years and is therefore not independent of the primary analysis.")

# One shared legend, drawn directly (independent of ggplot2/patchwork version,
# so it can never be duplicated or clipped).
leg <- data.table(cls = factor(CLASS_LEVELS, levels = CLASS_LEVELS),
                  x = c(0, 1.05, 2.45, 3.45, 4.85))
LEGEND_STRIP <- ggplot(leg) +
  geom_tile(aes(x = x, y = 0, fill = cls), width = 0.16, height = 0.55,
            colour = "grey30", linewidth = 0.3) +
  geom_text(aes(x = x + 0.12, y = 0, label = cls), hjust = 0, size = 4.3) +
  scale_fill_manual(values = CLASS_COLOURS, guide = "none") +
  coord_cartesian(xlim = c(-0.15, 6.0), ylim = c(-0.5, 0.5), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

FIG_S1 <- (wrap_plots(panels, nrow = 2, ncol = 4, byrow = TRUE) / LEGEND_STRIP) +
  plot_layout(heights = c(1, 0.045)) +
  plot_annotation(
    title    = "Spatial distribution of temperature-extreme trends during AMJ, 1976\u20132025",
    subtitle = paste0("283 IMD 1\u00B0 \u00D7 1\u00B0 land grid cells | Theil\u2013Sen slope per decade, ",
                      "significance from BH-FDR q \u2264 0.10 | native 1\u00B0 cells, no interpolation"),
    caption  = FOOT,
    theme = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 19),
                  plot.subtitle = element_text(hjust = 0.5, size = 11, colour = "grey25"),
                  plot.caption  = element_text(hjust = 0, size = 10, colour = "grey20",
                                               lineheight = 1.15),
                  plot.background = element_rect(fill = "white", colour = NA)))

# =====================================================================
# 10. SAVE
# =====================================================================
BASE <- file.path(OUT_DIR, "FIG_S1_SPATIAL_TRENDS_ALL8_CLASS_AMJ_1976_2025_2x4_v3")
ggsave(paste0(BASE, "_600dpi.png"),  FIG_S1, width = 22, height = 14.2, units = "in",
       dpi = 600, bg = "white", limitsize = FALSE)
ggsave(paste0(BASE, "_600dpi.tiff"), FIG_S1, width = 22, height = 14.2, units = "in",
       dpi = 600, bg = "white", compression = "lzw", limitsize = FALSE)
pdf_dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"
ggsave(paste0(BASE, ".pdf"),         FIG_S1, width = 22, height = 14.2, units = "in",
       device = pdf_dev, bg = "white", limitsize = FALSE)

fwrite(class_tab, file.path(OUT_DIR, "FIG_S1_CLASS_COUNTS_AMJ_1976_2025.csv"))
fwrite(d[, .(index, grid_id, lat, lon, sen_per_decade, p, q, significant_fdr,
             trend_class, map_class)],
       file.path(OUT_DIR, "FIG_S1_PLOTTED_DATA_AMJ_1976_2025.csv"))

cat("\n=====================================================\n")
cat("FIGURE S1 COMPLETED | period ", PERIOD_CODE, " | variant ", PRIMARY_VARIANT,
    " | q <= ", FDR_Q, "\n", sep = "")
cat("Statistics recomputed: NONE | Interpolation: NONE\n")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("Saved figure : ", BASE, "_600dpi.png\n", sep = "")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("=====================================================\n")

if (interactive()) tryCatch({ grDevices::graphics.off(); print(FIG_S1) },
                            error = function(e) message("Preview skipped; open the saved PNG."))