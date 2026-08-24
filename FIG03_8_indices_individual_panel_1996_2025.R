# =====================================================================
# FINAL JOURNAL-READY SPATIAL TREND MAPS
# AMJ 1996-2025 | 283 IMD 1 deg x 1 deg LAND GRIDS
#
# TXx, TXn, TNx, TNn, TX90p, TN10p, WSDI, DTR
#
# Trend        : locked Theil-Sen slope per decade (STEP 07)
# Significance : locked BH-FDR q <= 0.10          (STEP 07)
#
# NO trend, p, q or FDR value is calculated, refitted or adjusted here.
# The map class is DERIVED from two locked fields only:
#     sign(sen_per_decade)  and  significant_fdr
# and is cross-checked against the locked trend_class in the QC block.
#
# MAP_STYLE = "class"      -> 5-class categorical map (main Figure 3)
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
  "F:/WMO_IMD_R/WMO_IMD/data/",
  "AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv"
)

SHAPEFILE <- paste0(
  "F:/WMO_IMD_R/WMO_IMD/shapefile/",
  "STATE_BOUNDARY.shp"
)

OUT_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08/00_Supp_final_maps"

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
# Headroom at the top is reserved for the title strip and north arrow,
# matching the reference layout.

MAP_XLIM <- c(66.5, 98.5)
MAP_YLIM <- c(6.0,  39.6)

TITLE_Y  <- 38.55          # centre of the title strip
NORTH_X  <- 96.30          # north arrow centre
NORTH_Y  <- 35.40

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

# =====================================================================
# 11. SHARED THEME
# =====================================================================
# White everywhere: panel, plot and NA fill. No axes, no gridlines,
# no graticule - the reference layout.

theme_map <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "black",
                                      linewidth = 0.65),
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key       = element_rect(fill = "white", colour = NA),
      legend.title     = element_text(face = "bold", size = base_size - 1.5),
      legend.text      = element_text(size = base_size - 2),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# =====================================================================
# 12. NORTH ARROW  (clean dart, reference style)
# =====================================================================

north_arrow <- function(p, x = NORTH_X, y = NORTH_Y,
                        h = 2.10, w = 0.62, lab_size = 4.4) {
  p +
    annotate("text", x = x, y = y + h + 0.72, label = "N",
             fontface = "bold", size = lab_size, colour = "black") +
    annotate("polygon",
             x = c(x, x - w, x,          x + w),
             y = c(y + h, y, y + 0.30 * h, y),
             fill = "black", colour = "black", linewidth = 0.30)
}


# =====================================================================
# 14. MAP BUILDER
# =====================================================================

make_map <- function(dd, index_name, panel_label = NULL,
                     style = MAP_STYLE, inset_legend = TRUE,
                     base_size = 10) {
  
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
      scale_fill_manual(values = CLASS_COLOURS,
                        limits = CLASS_LEVELS,
                        drop   = FALSE,
                        na.value = "white",
                        name   = NULL) +
      guides(fill = guide_legend(
        ncol = 1, byrow = TRUE, reverse = FALSE,
        keywidth  = unit(0.62, "cm"),
        keyheight = unit(0.40, "cm"),
        override.aes = list(colour = "grey30", linewidth = 0.25)))
  } else {
    lim <- INDEX_LIMITS[[index_name]]
    p <- p +
      geom_tile(data = dd,
                aes(x = lon, y = lat, fill = sen_per_decade),
                width = 1, height = 1) +
      scale_fill_gradient2(
        low = "#0571B0", mid = "#FFFFFF", high = "#CA0020",
        midpoint = 0, limits = c(-lim, lim),
        name = INDEX_UNITS[[index_name]],
        na.value = "white",
        guide = guide_colourbar(
          title.position = "top", title.hjust = 0.5,
          barheight = unit(3.4, "cm"), barwidth = unit(0.40, "cm"),
          frame.colour = "grey20", ticks.colour = "grey20"))
  }
  
  # ---- insufficient data: white fill, dotted outline (no grey background)
  if (nrow(insufficient) > 0L)
    p <- p + geom_tile(data = insufficient, aes(x = lon, y = lat),
                       width = 1, height = 1, fill = "white",
                       colour = "grey25", linewidth = 0.38, linetype = "dotted")
  
  # ---- test undefined: slope exists, no valid test -> solid outline
  if (nrow(undefined) > 0L)
    p <- p + geom_tile(data = undefined, aes(x = lon, y = lat),
                       width = 1, height = 1, fill = NA,
                       colour = "black", linewidth = 0.42, linetype = "solid")
  
  # NOTE: no stippling layer. In the class style, BH-FDR significance is
  # carried by the fill colour itself, so a stipple would encode the same
  # information twice and is the source of the stray dots. Significance is
  # never shown by two channels at once.
  
  # ---- state boundary overlay (clips nothing)
  p <- p +
    geom_sf(data = bnd, inherit.aes = FALSE, fill = NA,
            colour = COL_BOUNDARY, linewidth = LW_BOUNDARY)
  
  p <- p +
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE, datum = NA) +
    theme_map(base_size = base_size)
  
  # ================================================================
  # TITLE — NO EXTRA BOX
  # ================================================================
  
  p <- p +
    annotate(
      "text",
      x = mean(MAP_XLIM),
      y = TITLE_Y,
      label = paste0(
        index_name,
        " Trend, AMJ ",
        START_YEAR,
        "\u2013",
        END_YEAR
      ),
      fontface = "bold",
      size = 4.0,
      colour = "black"
    )
  p <- north_arrow(p)
  
  if (!is.null(panel_label))
    p <- p + annotate("text",
                      x = MAP_XLIM[1] + 0.55, y = MAP_YLIM[2] - 0.55,
                      label = panel_label, hjust = 0, vjust = 1,
                      fontface = "bold", size = 4.4, colour = "black")
  
  if (isTRUE(inset_legend)) {
    p <- p + theme(
      legend.position    = c(0.685, 0.235),
      legend.justification = c(0, 0.5),
      legend.background  = element_rect(fill = "white", colour = NA),
      legend.margin      = margin(2, 4, 2, 2),
      legend.key.spacing.y = unit(1.2, "pt")
    )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  p
}

# =====================================================================
# 15. INDIVIDUAL MAPS
# =====================================================================

cat("\n=====================================================\n")
cat("GENERATING INDIVIDUAL MAPS | style = ", MAP_STYLE, "\n", sep = "")
cat("=====================================================\n")

sfx <- if (identical(MAP_STYLE, "class")) "CLASS" else "SENSLOPE"

map_list  <- list()   # with inset legend, for standalone use
panel_list <- list()  # legend suppressed, for the composite

for (j in seq_along(INDEX_ORDER)) {
  
  ix <- INDEX_ORDER[j]
  cat("Creating ", ix, "...\n", sep = "")
  
  dd <- d[index == ix]
  
  map_list[[ix]]   <- make_map(dd, ix, panel_label = LETTERS[j],
                               inset_legend = TRUE,  base_size = 10)
  panel_list[[ix]] <- make_map(dd, ix, panel_label = LETTERS[j],
                               inset_legend = FALSE, base_size = 11)
  
  fbase <- file.path(OUT_DIR, paste0("MAP_", ix, "_", sfx,
                                     "_AMJ_1996_2025"))
  
  ggsave(paste0(fbase, "_600dpi.png"),  map_list[[ix]],
         width = 7.2, height = 9.4, units = "in", dpi = 600, bg = "white")
  
  ggsave(paste0(fbase, "_600dpi.tiff"), map_list[[ix]],
         width = 7.2, height = 9.4, units = "in", dpi = 600, bg = "white",
         compression = "lzw")
  
  ggsave(paste0(fbase, ".pdf"),         map_list[[ix]],
         width = 7.2, height = 9.4, units = "in", device = cairo_pdf,
         bg = "white")
}

# =====================================================================
# 16. COMPOSITE: 2 ROWS x 4 COLUMNS
# =====================================================================
#  Row 1 : TXx | TXn | TNx | TNn
#  Row 2 : TX90p | TN10p | WSDI | DTR

cat("\n=====================================================\n")
cat("CREATING 2 x 4 COMPOSITE FIGURE\n")
cat("=====================================================\n")

# one shared legend for the whole composite, collected at the bottom
legend_donor <- make_map(d[index == INDEX_ORDER[1]], INDEX_ORDER[1],
                         panel_label = NULL, inset_legend = FALSE,
                         base_size = 12) +
  theme(legend.position = "bottom")

if (identical(MAP_STYLE, "class")) {
  legend_donor <- legend_donor +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, title = NULL,
                               keywidth  = unit(0.85, "cm"),
                               keyheight = unit(0.42, "cm"),
                               override.aes = list(colour = "grey30",
                                                   linewidth = 0.3)))
}

panel_fig <- wrap_plots(panel_list[INDEX_ORDER], nrow = 2, ncol = 4,
                        byrow = TRUE) /
  wrap_elements(full = cowplot_free_legend <- {
    g <- ggplotGrob(legend_donor)
    lg <- g$grobs[[which(vapply(g$grobs, function(z) z$name, "") == "guide-box")]]
    lg
  }) +
  plot_layout(heights = c(1, 0.075))

panel_fig <- panel_fig +
  plot_annotation(
    title    = "Spatial distribution of temperature-extreme trends during AMJ",
    subtitle = paste0("283 IMD 1\u00B0 \u00D7 1\u00B0 land grid cells, ",
                      START_YEAR, "\u2013", END_YEAR,
                      " | Theil\u2013Sen slope per decade, ",
                      "significance from BH-FDR q \u2264 0.10 | ",
                      "native 1\u00B0 cells, no interpolation"),
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 19),
      plot.subtitle = element_text(hjust = 0.5, size = 11, colour = "grey25"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

# =====================================================================
# 17. SAVE THE COMPOSITE
# =====================================================================

PANEL_BASE <- file.path(
  OUT_DIR,
  paste0("FIGURE_3_SPATIAL_TRENDS_ALL8_", sfx, "_AMJ_1996_2025_2x4")
)

PANEL_PNG  <- paste0(PANEL_BASE, "_600dpi.png")
PANEL_TIFF <- paste0(PANEL_BASE, "_600dpi.tiff")
PANEL_PDF  <- paste0(PANEL_BASE, ".pdf")

ggsave(PANEL_PNG,  panel_fig, width = 22, height = 13.6, units = "in",
       dpi = 600, bg = "white", limitsize = FALSE)

ggsave(PANEL_TIFF, panel_fig, width = 22, height = 13.6, units = "in",
       dpi = 600, bg = "white", compression = "lzw", limitsize = FALSE)

ggsave(PANEL_PDF,  panel_fig, width = 22, height = 13.6, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)

# =====================================================================
# 18. FINAL QC AND PROVENANCE
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
cat("Background             : white (panel, plot and NA fill)\n")
cat("Boundary CRS           : EPSG:4326\n")
cat("Composite layout       : 2 rows x 4 columns\n")

class_tab <- dcast(d[, .N, by = .(index, map_class)],
                   index ~ map_class, value.var = "N", fill = 0L)
fwrite(class_tab, file.path(OUT_DIR,
                            paste0("FIGURE_3_CLASS_COUNTS_", sfx, "_AMJ_1996_2025.csv")))

fwrite(d[, .(index, grid_id, lat, lon, sen_per_decade, p, q,
             significant_fdr, trend_class, map_class)],
       file.path(OUT_DIR,
                 paste0("FIGURE_3_PLOTTED_DATA_", sfx, "_AMJ_1996_2025.csv")))

for (f in c(PANEL_PNG, PANEL_TIFF, PANEL_PDF))
  cat(sprintf("%-6s : %s (%.1f MB)\n", toupper(tools::file_ext(f)),
              basename(f), file.size(f) / 1024^2))

cat("\nOutput directory:\n", OUT_DIR, "\n", sep = "")
cat("\n=====================================================\n")
cat("FINAL MAP GENERATION COMPLETED\n")
cat("=====================================================\n")

print(panel_fig)