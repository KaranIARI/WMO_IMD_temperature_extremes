# =====================================================================
# FIGURE 4 (manuscript numbering) - PRESENTATION REVISION v2
# SPATIAL TREND MAPS | AMJ 1996-2025 | 283 IMD 1 deg x 1 deg LAND GRIDS
#
# TXx, TXn, TNx, TNn, TX90p, TN10p, WSDI, DTR
#
# Trend        : locked Theil-Sen slope per decade (STEP 07)
# Significance : locked BH-FDR q <= 0.10          (STEP 07)
#
# WHAT CHANGED IN v2 (presentation only - see section 11 onward)
#   * Sections 1-10 (packages, paths, locked settings, reading, filtering,
#     map-class derivation, data QC, boundary handling, continuous limits)
#     are copied unchanged from the original script.
#   * Composite redesigned at the Scientific Reports double-column width
#     (183 mm), so fonts and lines are sized for the printed figure.
#   * New title, no subtitle; panel letters A-H; index name above each map.
#   * One shared two-row legend that also explains the outlined
#     test-undefined cells (TXn), which the old legend did not show.
#   * New output file names (suffix _v2) - original files are NOT overwritten.
#   * Extra presentation QC: plotted cells re-read from the built ggplot
#     objects and compared with the locked classes, and (if present) with
#     the previously exported FIGURE_3_PLOTTED_DATA file.
#
# NO trend, p, q or FDR value is calculated, refitted or adjusted here.
# The map class is DERIVED from two locked fields only:
#     sign(sen_per_decade)  and  significant_fdr
# and is cross-checked against the locked trend_class in the QC block.
#
# MAP_STYLE = "class"      -> 5-class categorical map (main Figure 4)
# MAP_STYLE = "continuous" -> continuous Theil-Sen colourbar (Supplementary)
# =====================================================================

rm(list = ls())

# =====================================================================
# 1. PACKAGES
# =====================================================================

library(data.table)
library(sf)
library(ggplot2)
library(patchwork)
library(grid)

# =====================================================================
# 2. PATHS
# =====================================================================

TREND_FILE <- paste0(
  "H:/WMO_IMD_R/WMO_IMD/data/",
  "AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv"
)

SHAPEFILE <- paste0(
  "H:/WMO_IMD_R/WMO_IMD/shapefile/",
  "STATE_BOUNDARY.shp"
)

OUT_DIR <- "H:/WMO_IMD_R/WMO_IMD/STEP08/00_Supp_final_maps"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 3. LOCKED SETTINGS
# =====================================================================

START_YEAR      <- 1996L
END_YEAR        <- 2025L
PRIMARY_VARIANT <- "HR_lag3"
FDR_Q           <- 0.10
EXPECTED_GRIDS  <- 283L

# "class" for the main figure, "continuous" for the supplementary version
MAP_STYLE <- "class"

INDEX_ORDER <- c("TXx", "TXn", "TNx", "TNn",
                 "TX90p", "TN10p", "WSDI", "DTR")

INDEX_UNITS <- c(
  TXx   = "\u00B0C decade\u207B\u00B9",
  TXn   = "\u00B0C decade\u207B\u00B9",
  TNx   = "\u00B0C decade\u207B\u00B9",
  TNn   = "\u00B0C decade\u207B\u00B9",
  TX90p = "% decade\u207B\u00B9",
  TN10p = "% decade\u207B\u00B9",
  WSDI  = "days decade\u207B\u00B9",
  DTR   = "\u00B0C decade\u207B\u00B9"
)

# ---- five-class scheme -----------------------------------------------------
CLASS_LEVELS <- c(
  "Significant decrease",
  "Decreasing tendency (not significant)",
  "No trend (zero slope)",
  "Increasing tendency (not significant)",
  "Significant increase"
)

CLASS_COLOURS <- c(
  "Significant decrease"                  = "#0571B0",  # dark blue
  "Decreasing tendency (not significant)" = "#92C5DE",  # light blue
  "No trend (zero slope)"                 = "#FFFFFF",  # white
  "Increasing tendency (not significant)" = "#F4A582",  # light red
  "Significant increase"                  = "#B2182B"   # dark red
)

COL_BOUNDARY <- "grey10"   # dark, clearly visible, not heavy
LW_BOUNDARY  <- 0.42
COL_CELL_EDGE <- "grey55"  # thin cell outline keeps 1 deg support visible
LW_CELL_EDGE  <- 0.08

# =====================================================================
# 4. MAP EXTENT
# =====================================================================
# v2: the extent is now set in section 11 from the boundary and the grid
# cells themselves (identical for every panel). The original fixed limits
# with title headroom are no longer needed because titles sit outside
# the map frame. CRS, cell positions and boundaries are unchanged.

# =====================================================================
# 5. READ LOCKED STEP07 RESULTS
# =====================================================================

cat("\n=====================================================\n")
cat("READING LOCKED STEP07 TREND RESULTS\n")
cat("=====================================================\n")

d <- fread(TREND_FILE)
cat("Rows read: ", nrow(d), "\n", sep = "")

required_cols <- c("grid_id", "period", "index", "variant", "lat", "lon",
                   "sen_per_decade", "p", "q", "trend_class",
                   "testable", "significant_fdr")

missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0L)
  stop(paste0("Missing required columns:\n",
              paste(missing_cols, collapse = ", ")))

# =====================================================================
# 6. FILTER TO THE LOCKED PRIMARY RESULT
# =====================================================================

d <- d[period == "TREND_1996_2025" &
         variant == PRIMARY_VARIANT &
         index %in% INDEX_ORDER]

if ("is_primary" %in% names(d)) d <- d[is_primary == TRUE]

cat("Rows after primary filtering: ", nrow(d), "\n", sep = "")

# =====================================================================
# 7. DERIVE THE MAP CLASS  (no statistic is recomputed)
# =====================================================================
# sign of the locked slope + locked significant_fdr flag, nothing else.

d[, sen_sign := fifelse(!is.finite(sen_per_decade), NA_integer_,
                        fifelse(sen_per_decade > 0,  1L,
                                fifelse(sen_per_decade < 0, -1L, 0L)))]

d[, map_class := fcase(
  trend_class == "insufficient_data",                       NA_character_,
  trend_class == "test_undefined",                          NA_character_,
  is.na(sen_sign),                                          NA_character_,
  sen_sign ==  0L,                        "No trend (zero slope)",
  sen_sign == -1L &  significant_fdr,     "Significant decrease",
  sen_sign == -1L & !significant_fdr,     "Decreasing tendency (not significant)",
  sen_sign ==  1L &  significant_fdr,     "Significant increase",
  sen_sign ==  1L & !significant_fdr,     "Increasing tendency (not significant)"
)]

d[, map_class := factor(map_class, levels = CLASS_LEVELS)]

# =====================================================================
# 8. DATA QC
# =====================================================================

index_counts <- d[, .(n = .N, n_grid = uniqueN(grid_id)), by = index]
print(index_counts)

if (nrow(index_counts) != length(INDEX_ORDER))
  stop("QC FAILED: expected exactly 8 indices.")

if (any(index_counts$n_grid != EXPECTED_GRIDS))
  stop("QC FAILED: every index must contain exactly 283 grid cells.")

if (nrow(d[, .N, by = .(index, grid_id)][N > 1]) > 0L)
  stop("QC FAILED: duplicate index x grid_id combinations detected.")

# significant_fdr must agree with the locked q
if (!isTRUE(d[testable == TRUE & is.finite(q), all(significant_fdr == (q <= FDR_Q))]))
  stop("QC FAILED: significant_fdr disagrees with q <= 0.10.")

# the derived class must never contradict the locked trend_class
bad_up <- d[trend_class == "significant_increase" &
              map_class != "Significant increase", .N]
bad_dn <- d[trend_class == "significant_decrease" &
              map_class != "Significant decrease", .N]
if (bad_up + bad_dn > 0L)
  stop("QC FAILED: derived map_class contradicts the locked trend_class.")

# no cell may be shown as significant unless the locked flag says so
if (d[map_class %in% c("Significant increase", "Significant decrease") &
      significant_fdr != TRUE, .N] > 0L)
  stop("QC FAILED: a cell is coloured significant without significant_fdr.")

# every cell is either classified or is one of the two non-plottable states
unclass_bad <- d[is.na(map_class) &
                   !trend_class %in% c("test_undefined", "insufficient_data"), .N]
if (unclass_bad > 0L)
  stop("QC FAILED: ", unclass_bad, " cells are unclassified for no declared reason.")

cat("\nClass composition per index (derived from locked fields only):\n")
print(dcast(d[, .N, by = .(index, map_class)],
            index ~ map_class, value.var = "N", fill = 0L))

cat("\nsignificant_zero_slope cells (q <= 0.10 with an exactly zero slope), ",
    "shown in the zero-slope class:\n", sep = "")
print(d[trend_class == "significant_zero_slope", .N, by = index])

cat("\nNon-plottable cells:\n")
print(d[trend_class %in% c("test_undefined", "insufficient_data"),
        .N, by = .(index, trend_class)])

cat("\nData QC: PASS\n")

# =====================================================================
# 9. STATE BOUNDARY
# =====================================================================

cat("\n=====================================================\n")
cat("READING STATE BOUNDARY\n")
cat("=====================================================\n")

bnd <- st_read(SHAPEFILE, quiet = TRUE)
cat("Boundary features: ", nrow(bnd), "\n", sep = "")

if (is.na(st_crs(bnd)))
  stop("QC FAILED: boundary has no CRS.")

cat("Original boundary CRS: ", st_crs(bnd)$input, "\n", sep = "")

if (!all(st_is_valid(bnd))) {
  cat("Invalid geometries detected - repairing with st_make_valid()...\n")
  bnd <- st_make_valid(bnd)
}

bnd <- st_transform(bnd, crs = 4326)
cat("Boundary transformed to EPSG:4326\n")

# Multipart geometries are preserved; geom_sf handles them directly.
# The boundary is drawn as an overlay only - it clips nothing.

# =====================================================================
# 10. CONTINUOUS COLOUR LIMITS  (used only when MAP_STYLE == "continuous")
# =====================================================================

INDEX_LIMITS <- list()

for (ix in INDEX_ORDER) {
  vals <- d[index == ix & is.finite(sen_per_decade), sen_per_decade]
  lim  <- suppressWarnings(max(abs(vals), na.rm = TRUE))
  if (!is.finite(lim) || lim == 0) lim <- 1
  base <- 10^(floor(log10(lim)) - 1)
  INDEX_LIMITS[[ix]] <- ceiling(lim / base) * base   # round UP: never squished
}

# #####################################################################
# EVERYTHING BELOW THIS LINE IS PRESENTATION ONLY (v2).
# It reads the objects d, bnd and INDEX_LIMITS created above and never
# changes them.
# #####################################################################

# =====================================================================
# 11. PRESENTATION SETTINGS (v2)
# =====================================================================

FIG_VERSION   <- "v2"
FIG_TITLE     <- paste0("Spatial Distribution of Temperature-Extreme Trends ",
                        "over India during AMJ (", START_YEAR, "–",
                        END_YEAR, ")")

# Final printed size: Scientific Reports double-column width = 183 mm.
# All font sizes below are in points AT THIS PRINTED SIZE.
FIG_WIDTH_MM  <- 183
SINGLE_COL_MM <- 89            # width used for the individual maps
BASE_FAMILY   <- "sans"        # Arial / Helvetica

FS_TITLE      <- 10            # figure title
FS_PANEL      <- 8.5           # index name above each map
FS_TAG        <- 9.5           # panel letters A-H
FS_LEGEND     <- 7.5           # legend labels
FS_NORTH      <- 6.5           # "N" of the north arrow

# Same five class colours as the original. The only change: the zero-slope
# class is drawn in a very light grey instead of pure white, so that cells
# with an exactly zero slope (e.g. 282 of 283 WSDI cells) are visible as
# analysed cells rather than looking like empty map. The class meaning is
# unchanged. Set to FALSE to restore pure white.
ZERO_FILL_LIGHT_GREY <- TRUE
CLASS_COLOURS_PLOT   <- CLASS_COLOURS
if (ZERO_FILL_LIGHT_GREY)
  CLASS_COLOURS_PLOT["No trend (zero slope)"] <- "#EBEBEB"

# Line weights re-scaled for the 183 mm print (original values were set
# for a 22-inch-wide canvas and become hairlines or blobs when reduced).
COL_BOUNDARY  <- "grey15";  LW_BOUNDARY  <- 0.22
COL_CELL_EDGE <- "grey60";  LW_CELL_EDGE <- 0.06
LW_FRAME      <- 0.35
LW_UNDEFINED  <- 0.38          # outline of test-undefined cells
LW_INSUFF     <- 0.30          # dotted outline of insufficient-data cells

# North arrow: "all" panels (as in the original) or "first" panel only.
NORTH_ARROW_PANELS <- "all"

# ---- map extent: identical for every panel ----------------------------
# Derived from the boundary and the 1 deg cells so that the whole of India
# (including Ladakh, the northeast, Lakshadweep and the Andaman & Nicobar
# Islands) and every grid cell fit, with a small even margin.
bb_b <- st_bbox(bnd)
PAD  <- 0.35
MAP_XLIM <- c(min(bb_b[["xmin"]], min(d$lon) - 0.5) - PAD,
              max(bb_b[["xmax"]], max(d$lon) + 0.5) + PAD)
MAP_YLIM <- c(min(bb_b[["ymin"]], min(d$lat) - 0.5) - PAD,
              max(bb_b[["ymax"]], max(d$lat) + 0.5) + PAD)

# geographic aspect ratio applied by coord_sf for EPSG:4326
MAP_ASPECT <- diff(MAP_YLIM) /
  (diff(MAP_XLIM) * cos(mean(MAP_YLIM) * pi / 180))

cat("\nMap extent (lon): ", sprintf("%.2f to %.2f", MAP_XLIM[1], MAP_XLIM[2]),
    "\nMap extent (lat): ", sprintf("%.2f to %.2f", MAP_YLIM[1], MAP_YLIM[2]),
    "\n", sep = "")

# north arrow in the empty Arabian Sea corner (bottom-left), clear of
# Gujarat and Lakshadweep
NORTH_X <- MAP_XLIM[1] + 1.25
NORTH_Y <- MAP_YLIM[1] + 0.80

# cells that are shown only as outlines
UNDEF_INDICES  <- d[trend_class == "test_undefined",   sort(unique(as.character(index)))]
INSUFF_INDICES <- d[trend_class == "insufficient_data", sort(unique(as.character(index)))]

# =====================================================================
# 12. SHARED THEME AND NORTH ARROW (v2)
# =====================================================================

theme_map_v2 <- function() {
  theme_void(base_size = FS_PANEL, base_family = BASE_FAMILY) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "black",
                                      linewidth = LW_FRAME),
      panel.grid       = element_blank(),
      plot.title       = element_text(size = FS_PANEL, face = "bold",
                                      hjust = 0.5,
                                      margin = margin(0, 0, 1.2, 0, "mm")),
      plot.tag         = element_text(size = FS_TAG, face = "bold"),
      plot.tag.position = c(0.035, 0.985),
      legend.position  = "none",
      plot.margin      = margin(0.8, 1.2, 0.8, 1.2, "mm")
    )
}

north_arrow_v2 <- function(p, x = NORTH_X, y = NORTH_Y, h = 2.2, w = 0.55) {
  p +
    annotate("polygon",
             x = c(x, x - w, x, x + w),
             y = c(y + h, y, y + 0.32 * h, y),
             fill = "black", colour = "black", linewidth = 0.15) +
    annotate("text", x = x, y = y + h + 0.75, label = "N",
             fontface = "bold", family = BASE_FAMILY,
             size = FS_NORTH / .pt, colour = "black")
}

# =====================================================================
# 13. MAP BUILDER (v2) - same layers and data as the original
# =====================================================================

make_map_v2 <- function(dd, index_name, title = index_name, tag = NULL,
                        style = MAP_STYLE, show_north = TRUE) {
  
  undefined    <- dd[trend_class == "test_undefined"]
  insufficient <- dd[trend_class == "insufficient_data"]
  plottable    <- dd[!is.na(map_class)]
  
  p <- ggplot()
  
  # ---- 1 deg x 1 deg cells: no interpolation, no regridding, no smoothing
  if (identical(style, "class")) {
    p <- p +
      geom_tile(data = plottable,
                aes(x = lon, y = lat, fill = map_class),
                width = 1, height = 1,
                colour = COL_CELL_EDGE, linewidth = LW_CELL_EDGE) +
      scale_fill_manual(values = CLASS_COLOURS_PLOT, limits = CLASS_LEVELS,
                        drop = FALSE, na.value = "white", guide = "none")
  } else {
    lim <- INDEX_LIMITS[[index_name]]
    p <- p +
      geom_tile(data = dd, aes(x = lon, y = lat, fill = sen_per_decade),
                width = 1, height = 1,
                colour = COL_CELL_EDGE, linewidth = LW_CELL_EDGE) +
      scale_fill_gradient2(
        low = "#0571B0", mid = "#FFFFFF", high = "#CA0020",
        midpoint = 0, limits = c(-lim, lim),
        name = INDEX_UNITS[[index_name]], na.value = "white",
        guide = guide_colourbar(
          title.position = "top", title.hjust = 0.5,
          barheight = unit(18, "mm"), barwidth = unit(1.8, "mm"),
          frame.colour = "grey20", ticks.colour = "grey20"))
  }
  
  # ---- insufficient data: white fill, dotted outline (as original)
  if (nrow(insufficient) > 0L)
    p <- p + geom_tile(data = insufficient, aes(x = lon, y = lat),
                       width = 1, height = 1, fill = "white",
                       colour = "grey25", linewidth = LW_INSUFF,
                       linetype = "dotted")
  
  # ---- test undefined: solid black outline, no fill (as original)
  if (nrow(undefined) > 0L)
    p <- p + geom_tile(data = undefined, aes(x = lon, y = lat),
                       width = 1, height = 1, fill = NA,
                       colour = "black", linewidth = LW_UNDEFINED,
                       linetype = "solid")
  
  # ---- state boundary overlay (clips nothing)
  p <- p +
    geom_sf(data = bnd, inherit.aes = FALSE, fill = NA,
            colour = COL_BOUNDARY, linewidth = LW_BOUNDARY) +
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE, datum = NA) +
    labs(title = title, tag = tag) +
    theme_map_v2()
  
  if (!identical(style, "class"))
    p <- p + theme(legend.position = "right",
                   legend.title = element_text(size = FS_LEGEND - 0.5),
                   legend.text  = element_text(size = FS_LEGEND - 1))
  
  if (isTRUE(show_north)) p <- north_arrow_v2(p)
  p
}

# =====================================================================
# 14. SHARED LEGEND (v2) - drawn explicitly so every map symbol is listed
# =====================================================================
# Two rows: decreases on the top row, increases on the bottom row, so the
# significant / not-significant pairs line up in columns. The third column
# holds the zero-slope class and the outlined cells.

text_width_mm <- function(lbl, fs, face = "plain") {
  vapply(lbl, function(s)
    convertWidth(grobWidth(textGrob(s, gp = gpar(fontsize = fs,
                                                 fontfamily = BASE_FAMILY,
                                                 fontface = face))),
                 "mm", valueOnly = TRUE), numeric(1))
}

legend_items <- function() {
  it <- data.table(
    label = c("Significant decrease",
              "Decreasing tendency (not significant)",
              "No trend (zero slope)",
              "Significant increase",
              "Increasing tendency (not significant)"),
    row   = c(1, 1, 1, 2, 2),
    col   = c(1, 2, 3, 1, 2))
  it[, fill := CLASS_COLOURS_PLOT[label]]
  it[, `:=`(edge = "grey30", lty = "solid", lw = 0.25)]
  if (length(UNDEF_INDICES))
    it <- rbind(it, data.table(
      label = paste0("Test undefined (", paste(UNDEF_INDICES, collapse = ", "), ")"),
      row = 2, col = 3, fill = "white", edge = "black", lty = "solid",
      lw = LW_UNDEFINED))
  if (length(INSUFF_INDICES))
    it <- rbind(it, data.table(
      label = paste0("Insufficient data (", paste(INSUFF_INDICES, collapse = ", "), ")"),
      row = if (length(UNDEF_INDICES)) 3 else 2, col = 3, fill = "white",
      edge = "grey25", lty = "dotted", lw = LW_INSUFF))
  it
}

build_legend_v2 <- function(width_mm, fs = FS_LEGEND, one_col = FALSE) {
  it <- copy(legend_items())
  if (one_col) {                       # narrow layout for single-column maps
    ord <- c(match(CLASS_LEVELS, it$label),
             setdiff(seq_len(nrow(it)), match(CLASS_LEVELS, it$label)))
    it <- it[ord][, `:=`(row = seq_len(.N), col = 1)]
  }
  key_w <- 4.2; key_h <- 2.6; gap_txt <- 1.3; gap_col <- 4.5; row_h <- 4.2
  it[, tw := text_width_mm(label, fs)]
  colw <- it[, .(w = max(tw)), by = col][order(col)]
  colw[, x0 := cumsum(c(0, head(key_w + gap_txt + w + gap_col, -1)))]
  total_w <- colw[, max(x0 + key_w + gap_txt + w)]
  off <- (width_mm - total_w) / 2
  it <- merge(it, colw[, .(col, x0)], by = "col")
  it[, `:=`(x0 = x0 + off)]
  n_rows <- max(it$row)
  head_txt <- "Trend class (significance: BH-FDR q ≤ 0.10)"
  height <- 3.4 + n_rows * row_h + 1
  it[, yc := height - 3.4 - (row - 0.5) * row_h]
  
  g <- ggplot() +
    annotate("text", x = width_mm / 2, y = height - 1.3, label = head_txt,
             size = fs / .pt, family = BASE_FAMILY, fontface = "bold",
             vjust = 0.5)
  for (i in seq_len(nrow(it))) {
    r <- it[i]
    g <- g +
      annotate("rect", xmin = r$x0, xmax = r$x0 + key_w,
               ymin = r$yc - key_h / 2, ymax = r$yc + key_h / 2,
               fill = r$fill, colour = r$edge, linewidth = r$lw,
               linetype = r$lty) +
      annotate("text", x = r$x0 + key_w + gap_txt, y = r$yc, label = r$label,
               hjust = 0, vjust = 0.5, size = fs / .pt,
               family = BASE_FAMILY)
  }
  g <- g +
    coord_cartesian(xlim = c(0, width_mm), ylim = c(0, height),
                    expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", colour = NA))
  list(plot = g, height_mm = height, labels = it$label,
       total_w = total_w)
}

# =====================================================================
# 15. PANELS AND INDIVIDUAL MAPS
# =====================================================================

cat("\n=====================================================\n")
cat("GENERATING PANELS | style = ", MAP_STYLE, " | ", FIG_VERSION, "\n", sep = "")
cat("=====================================================\n")

sfx <- if (identical(MAP_STYLE, "class")) "CLASS" else "SENSLOPE"

panel_list <- list()
for (j in seq_along(INDEX_ORDER)) {
  ix <- INDEX_ORDER[j]
  show_n <- identical(NORTH_ARROW_PANELS, "all") || j == 1L
  panel_list[[ix]] <- make_map_v2(d[index == ix], ix, title = ix,
                                  tag = LETTERS[j], show_north = show_n)
}

# ---- device helpers (ragg if installed, otherwise base devices) -------
save_png <- function(file, plot, w_mm, h_mm) {
  if (requireNamespace("ragg", quietly = TRUE))
    ggsave(file, plot, width = w_mm, height = h_mm, units = "mm", dpi = 600,
           bg = "white", device = ragg::agg_png)
  else
    ggsave(file, plot, width = w_mm, height = h_mm, units = "mm", dpi = 600,
           bg = "white")
}
save_tiff <- function(file, plot, w_mm, h_mm) {
  if (requireNamespace("ragg", quietly = TRUE))
    ggsave(file, plot, width = w_mm, height = h_mm, units = "mm", dpi = 600,
           bg = "white", device = ragg::agg_tiff, compression = "lzw")
  else
    ggsave(file, plot, width = w_mm, height = h_mm, units = "mm", dpi = 600,
           bg = "white", device = "tiff", compression = "lzw")
}
save_pdf <- function(file, plot, w_mm, h_mm) {
  dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file, plot, width = w_mm, height = h_mm, units = "mm",
         device = dev, bg = "white")
}

# ---- individual maps (single-column width, one-column legend below) ---
TITLE_BAND_MM <- 5.0            # height reserved for a panel title line
ind_map_w  <- SINGLE_COL_MM - 2.4
ind_leg    <- build_legend_v2(SINGLE_COL_MM, one_col = TRUE)
ind_h      <- ind_map_w * MAP_ASPECT + TITLE_BAND_MM + 1.6 + ind_leg$height_mm

for (j in seq_along(INDEX_ORDER)) {
  ix <- INDEX_ORDER[j]
  m  <- make_map_v2(d[index == ix], ix,
                    title = paste0(ix, " trend, AMJ ", START_YEAR, "–",
                                   END_YEAR),
                    tag = NULL, show_north = TRUE)
  if (identical(MAP_STYLE, "class"))
    m <- m / ind_leg$plot +
    plot_layout(heights = c(ind_h - ind_leg$height_mm, ind_leg$height_mm))
  fbase <- file.path(OUT_DIR, paste0("MAP_", ix, "_", sfx,
                                     "_AMJ_1996_2025_", FIG_VERSION))
  save_png (paste0(fbase, "_600dpi.png"),  m, SINGLE_COL_MM, ind_h)
  save_tiff(paste0(fbase, "_600dpi.tiff"), m, SINGLE_COL_MM, ind_h)
  save_pdf (paste0(fbase, ".pdf"),         m, SINGLE_COL_MM, ind_h)
  cat("Saved individual map: ", basename(fbase), "\n", sep = "")
}

# =====================================================================
# 16. COMPOSITE: 2 ROWS x 4 COLUMNS (183 mm wide)
# =====================================================================
#  Row 1 : TXx | TXn | TNx | TNn
#  Row 2 : TX90p | TN10p | WSDI | DTR

cat("\n=====================================================\n")
cat("CREATING 2 x 4 COMPOSITE FIGURE\n")
cat("=====================================================\n")

panel_w_mm  <- FIG_WIDTH_MM / 4
map_w_mm    <- panel_w_mm - 2.4                    # minus side margins
panel_h_mm  <- map_w_mm * MAP_ASPECT + TITLE_BAND_MM + 1.6
FIG_TITLE_MM <- 7

comp_leg <- if (identical(MAP_STYLE, "class"))
  build_legend_v2(FIG_WIDTH_MM) else NULL
leg_h_mm <- if (is.null(comp_leg)) 0 else comp_leg$height_mm

maps_block <- wrap_plots(panel_list[INDEX_ORDER], nrow = 2, ncol = 4,
                         byrow = TRUE)

panel_fig <- if (is.null(comp_leg)) maps_block else
  (maps_block / comp_leg$plot) +
  plot_layout(heights = c(2 * panel_h_mm, leg_h_mm))

panel_fig <- panel_fig +
  plot_annotation(
    title = FIG_TITLE,
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = FS_TITLE,
                                family = BASE_FAMILY,
                                margin = margin(1.5, 0, 1.5, 0, "mm")),
      plot.margin = margin(1, 1.5, 1.5, 1.5, "mm"),
      plot.background = element_rect(fill = "white", colour = NA)))

FIG_HEIGHT_MM <- FIG_TITLE_MM + 2 * panel_h_mm + leg_h_mm + 3

# =====================================================================
# 17. SAVE THE COMPOSITE (new file names - originals are not touched)
# =====================================================================

PANEL_BASE <- file.path(
  OUT_DIR,
  paste0("FIGURE_4_SPATIAL_TRENDS_ALL8_", sfx, "_AMJ_1996_2025_2x4_",
         FIG_VERSION))

PANEL_PNG  <- paste0(PANEL_BASE, "_600dpi.png")
PANEL_TIFF <- paste0(PANEL_BASE, "_600dpi.tiff")
PANEL_PDF  <- paste0(PANEL_BASE, ".pdf")

save_png (PANEL_PNG,  panel_fig, FIG_WIDTH_MM, FIG_HEIGHT_MM)
save_tiff(PANEL_TIFF, panel_fig, FIG_WIDTH_MM, FIG_HEIGHT_MM)
save_pdf (PANEL_PDF,  panel_fig, FIG_WIDTH_MM, FIG_HEIGHT_MM)

# =====================================================================
# 18. PRESENTATION QC (v2)
# =====================================================================

cat("\n=====================================================\n")
cat("PRESENTATION QC\n")
cat("=====================================================\n")

QC <- list()
qc <- function(check, passed, detail) {
  ok <- isTRUE(passed)
  QC[[length(QC) + 1L]] <<- data.table(check = check,
                                       status = if (ok) "PASS" else "FAIL",
                                       detail = detail)
  cat(sprintf("  %-4s | %-46s | %s\n", if (ok) "PASS" else "FAIL",
              check, detail))
}

# (1) plotted cells = locked classes, read back from the built plots
mism <- 0L; n_checked <- 0L
for (ix in INDEX_ORDER) {
  b   <- ggplot_build(panel_list[[ix]])
  ld  <- as.data.table(b$data[[1]])[, .(x, y, fill)]
  exp <- d[index == ix & !is.na(map_class)]
  if (identical(MAP_STYLE, "class")) {
    exp <- exp[, .(x = lon, y = lat,
                   fill = unname(CLASS_COLOURS_PLOT[as.character(map_class)]))]
    cmp <- merge(ld, exp, by = c("x", "y"), all = TRUE,
                 suffixes = c("_plot", "_locked"))
    mism <- mism + cmp[is.na(fill_plot) | is.na(fill_locked) |
                         toupper(fill_plot) != toupper(fill_locked), .N]
    n_checked <- n_checked + nrow(cmp)
  } else {
    n_checked <- n_checked + nrow(ld)
    mism <- mism + abs(nrow(ld) - nrow(d[index == ix]))
  }
}
qc("plotted_cells_equal_locked_classes", mism == 0L,
   sprintf("%d cells checked across 8 panels, %d mismatches", n_checked, mism))

# (2) outlined cells drawn exactly where the locked file says
n_undef_plot <- sum(vapply(INDEX_ORDER, function(ix) {
  dd <- d[index == ix]
  if (!nrow(dd[trend_class == "test_undefined"])) return(0L)
  li <- 2L + as.integer(nrow(dd[trend_class == "insufficient_data"]) > 0L)
  ld <- as.data.table(ggplot_build(panel_list[[ix]])$data[[li]])
  exp <- dd[trend_class == "test_undefined", .(x = lon, y = lat)]
  if (all(ld$colour == "black") && nrow(ld) == nrow(exp) &&
      nrow(merge(ld[, .(x, y)], exp, by = c("x", "y"))) == nrow(exp))
    nrow(ld) else -1000L }, integer(1)))
qc("test_undefined_cells_outlined",
   n_undef_plot == d[trend_class == "test_undefined", .N],
   sprintf("%d test-undefined cell(s) outlined: %s", n_undef_plot,
           if (length(UNDEF_INDICES)) paste(UNDEF_INDICES, collapse = ", ") else "none"))

# (3) comparison with the previously exported plotted data, if present
OLD_PLOTTED <- file.path(OUT_DIR, paste0("FIGURE_3_PLOTTED_DATA_", sfx,
                                         "_AMJ_1996_2025.csv"))
if (file.exists(OLD_PLOTTED)) {
  old <- fread(OLD_PLOTTED, colClasses = list(character = "grid_id"))
  new <- d[, .(index = as.character(index), grid_id = as.character(grid_id),
               sen_new = sen_per_decade, cls_new = as.character(map_class))]
  cmp <- merge(old[, .(index = as.character(index),
                       grid_id = as.character(grid_id),
                       sen_old = sen_per_decade,
                       cls_old = as.character(map_class))],
               new, by = c("index", "grid_id"), all = TRUE)
  cmp[is.na(cls_old) | cls_old == "", cls_old := NA_character_]
  cmp[is.na(cls_new) | cls_new == "", cls_new := NA_character_]
  cmp[, same_cls := (is.na(cls_old) & is.na(cls_new)) |
        (!is.na(cls_old) & !is.na(cls_new) & cls_old == cls_new)]
  n_bad <- cmp[is.na(sen_old) != is.na(sen_new) |
                 (!is.na(sen_old) & !is.na(sen_new) &
                    abs(sen_old - sen_new) > 1e-9) |
                 !same_cls, .N]
  qc("identical_to_previous_figure_data", nrow(cmp) == nrow(d) && n_bad == 0L,
     sprintf("%d rows compared with %s, %d differences", nrow(cmp),
             basename(OLD_PLOTTED), n_bad))
} else {
  qc("identical_to_previous_figure_data", TRUE,
     "previous plotted-data CSV not found - comparison skipped")
}

# (4) eight panels, correct titles and letters
titles <- vapply(panel_list[INDEX_ORDER], function(p) p$labels$title, "")
tags   <- vapply(panel_list[INDEX_ORDER], function(p) p$labels$tag, "")
qc("eight_panels_titled_and_lettered",
   length(panel_list) == 8L && identical(unname(titles), INDEX_ORDER) &&
     identical(unname(tags), LETTERS[1:8]),
   paste(paste0(tags, "=", titles), collapse = " "))

# (5) legend lists every class used on the maps
if (!is.null(comp_leg)) {
  leg_ok <- all(CLASS_LEVELS %in% comp_leg$labels) &&
    (length(UNDEF_INDICES) == 0L || any(grepl("^Test undefined", comp_leg$labels))) &&
    (length(INSUFF_INDICES) == 0L || any(grepl("^Insufficient data", comp_leg$labels)))
  qc("legend_matches_map_symbols", leg_ok,
     paste0(length(comp_leg$labels), " entries: ",
            paste(comp_leg$labels, collapse = " | ")))
  qc("legend_fits_figure_width", comp_leg$total_w <= FIG_WIDTH_MM - 3,
     sprintf("legend %.1f mm wide within %d mm", comp_leg$total_w, FIG_WIDTH_MM))
}

# (6) print size and smallest text
qc("print_size_and_min_font",
   min(FS_TITLE, FS_PANEL, FS_TAG, FS_LEGEND, FS_NORTH) >= 6.5,
   sprintf("%d x %.0f mm; fonts %.1f-%.1f pt at print size; each map %.1f mm wide",
           FIG_WIDTH_MM, FIG_HEIGHT_MM,
           min(FS_TITLE, FS_PANEL, FS_TAG, FS_LEGEND, FS_NORTH),
           max(FS_TITLE, FS_PANEL, FS_TAG, FS_LEGEND), map_w_mm))

# (7) files written
qc("composite_files_written",
   all(file.exists(c(PANEL_PNG, PANEL_TIFF, PANEL_PDF))),
   "PNG (600 dpi), TIFF (600 dpi, LZW), PDF (vector)")

qct <- rbindlist(QC)[, version := FIG_VERSION]
fwrite(qct, paste0(PANEL_BASE, "_QC.csv"))

# =====================================================================
# 19. FINAL QC AND PROVENANCE (as original, new file names)
# =====================================================================

cat("\n=====================================================\n")
cat("FINAL MAP QC\n")
cat("=====================================================\n")
cat("Indices                : 8\n")
cat("Grid cells per index   : ", EXPECTED_GRIDS, "\n", sep = "")
cat("Period                 : AMJ ", START_YEAR, "-", END_YEAR, "\n", sep = "")
cat("Variant                : ", PRIMARY_VARIANT, "\n", sep = "")
cat("Trend                  : locked Theil-Sen slope per decade\n")
cat("Significance           : locked BH-FDR q <= ", FDR_Q, "\n", sep = "")
cat("Map style              : ", MAP_STYLE, "\n", sep = "")
cat("Class rule             : sign(sen_per_decade) + significant_fdr only\n")
cat("Statistics recomputed  : NONE\n")
cat("Interpolation          : None\n")
cat("Regridding             : None\n")
cat("Smoothing              : None\n")
cat("Clipping               : None (boundary is an overlay)\n")
cat("Boundary CRS           : EPSG:4326\n")
cat("Composite layout       : 2 rows x 4 columns, ", FIG_WIDTH_MM, " mm wide\n", sep = "")
cat("Presentation QC        : ", qct[status == "PASS", .N], " PASS | ",
    qct[status == "FAIL", .N], " FAIL\n", sep = "")

class_tab <- dcast(d[, .N, by = .(index, map_class)],
                   index ~ map_class, value.var = "N", fill = 0L)
fwrite(class_tab, file.path(OUT_DIR,
                            paste0("FIGURE_4_CLASS_COUNTS_", sfx, "_AMJ_1996_2025_", FIG_VERSION, ".csv")))

fwrite(d[, .(index, grid_id, lat, lon, sen_per_decade, p, q,
             significant_fdr, trend_class, map_class)],
       file.path(OUT_DIR, paste0("FIGURE_4_PLOTTED_DATA_", sfx,
                                 "_AMJ_1996_2025_", FIG_VERSION, ".csv")))

for (f in c(PANEL_PNG, PANEL_TIFF, PANEL_PDF))
  cat(sprintf("%-6s : %s (%.1f MB)\n", toupper(tools::file_ext(f)),
              basename(f), file.size(f) / 1024^2))

cat("\nOutput directory:\n", OUT_DIR, "\n", sep = "")
if (qct[status == "FAIL", .N] > 0L)
  warning("One or more presentation QC checks FAILED - see the QC CSV.")
cat("\n=====================================================\n")
cat("FINAL MAP GENERATION COMPLETED (", FIG_VERSION, ")\n", sep = "")
cat("=====================================================\n")

if (interactive()) print(panel_fig)