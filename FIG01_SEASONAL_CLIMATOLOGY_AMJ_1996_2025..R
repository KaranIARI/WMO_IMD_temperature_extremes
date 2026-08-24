###############################################################################
# FIGURE 1
# Seasonal climatology of temperature extremes during AMJ, 1996–2025
#
# FINAL PANEL SET = 7 INDICES
#
#   TXx | TXn | TNx | TNn
#   TX90p | TN10p | DTR
#
# WSDI is NOT included here because the locked WSDI definition treats
# April–June as one continuous 91-day sequence. A month-reset WSDI would
# represent a different index.
#
# INPUT FILES ARE READ ONLY.
# NO raw/index values are recalculated here.
#
# Outputs:
#   600-dpi PNG
#   vector PDF
#   exact plotted-data CSV
#   QC CSV
###############################################################################


rm(list = ls())


# =============================================================================
# 1. PACKAGES
# =============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "scales"
)

for (p in required_packages) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    
    stop(
      "Package not installed: ",
      p
    )
  }
}

suppressPackageStartupMessages({
  
  library(data.table)
  library(ggplot2)
  library(scales)
  
})


# =============================================================================
# 2. PATHS
# =============================================================================

OUTPUT_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08"

# Existing monthly climatology:
CLIM_FILE <- file.path(
  OUTPUT_DIR,
  "01_descriptive",
  "STEP08_01_DESCRIPTIVE_MONTHLY_CLIMATOLOGY_1996_2025.csv"
)

# Existing monthly TX90p/TN10p direct-labelled climatology:
PERCENTILE_FILE <- file.path(
  OUTPUT_DIR,
  "01_descriptive",
  "STEP08_01_DESCRIPTIVE_MONTHLY_TX90p_TN10p_DIRECT_1996_2025.csv"
)

FIG_DIR <- file.path(
  OUTPUT_DIR,
  "12_figures"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. LOCKED SETTINGS
# =============================================================================

START_YEAR <- 1996L
END_YEAR   <- 2025L

N_GRIDS <- 283L
N_YEARS <- 30L

MONTH_ORDER <- c(
  "April",
  "May",
  "June"
)

# -------------------------------------------------------------------------
# FINAL 7 FIGURE-1 INDICES
# -------------------------------------------------------------------------

INDEX_ORDER <- c(
  "TXx",
  "TXn",
  "TNx",
  "TNn",
  "TX90p",
  "TN10p",
  "DTR"
)

# -------------------------------------------------------------------------
# Units
# -------------------------------------------------------------------------

INDEX_UNITS <- c(
  
  TXx   = "°C",
  TXn   = "°C",
  TNx   = "°C",
  TNn   = "°C",
  
  TX90p = "%",
  TN10p = "%",
  
  DTR   = "°C"
  
)


# =============================================================================
# 4. CHECK INPUT FILES
# =============================================================================

if (!file.exists(CLIM_FILE)) {
  
  stop(
    "\nMonthly climatology file not found:\n",
    CLIM_FILE,
    "\n\nRun STEP08.1 first."
  )
  
}

if (!file.exists(PERCENTILE_FILE)) {
  
  stop(
    "\nTX90p/TN10p monthly climatology file not found:\n",
    PERCENTILE_FILE,
    "\n\nRun STEP08.1 with MONTHLY_TX90P_MODE = 'direct_labelled'."
  )
  
}


cat("\n============================================================\n")
cat("FIGURE 1 — FINAL 7-INDEX SEASONAL CLIMATOLOGY\n")
cat("============================================================\n\n")


# =============================================================================
# 5. READ EXISTING LOCKED OUTPUTS
# =============================================================================

cat("Reading existing STEP08 monthly climatology...\n")

clim <- fread(
  CLIM_FILE
)

cat(
  "Rows in general monthly climatology:",
  nrow(clim),
  "\n"
)


cat("\nReading existing TX90p/TN10p direct monthly climatology...\n")

pct <- fread(
  PERCENTILE_FILE
)

cat(
  "Rows in percentile monthly climatology:",
  nrow(pct),
  "\n"
)


# =============================================================================
# 6. CHECK REQUIRED COLUMNS
# =============================================================================

required_clim <- c(
  "index",
  "month",
  "unit",
  "N",
  "Mean",
  "Median",
  "Q1",
  "Q3"
)

missing_clim <- setdiff(
  required_clim,
  names(clim)
)

if (length(missing_clim) > 0L) {
  
  stop(
    "General climatology file missing columns:\n",
    paste(missing_clim, collapse = ", ")
  )
  
}


required_pct <- c(
  "index",
  "month",
  "scope",
  "unit",
  "period",
  "N",
  "Mean",
  "Median",
  "Q1",
  "Q3"
)

missing_pct <- setdiff(
  required_pct,
  names(pct)
)

if (length(missing_pct) > 0L) {
  
  stop(
    "TX90p/TN10p file missing columns:\n",
    paste(missing_pct, collapse = ", ")
  )
  
}


# =============================================================================
# 7. CLEAN INDEX NAMES
# =============================================================================

clim[
  ,
  index := sub(
    "_monthly_climatology$",
    "",
    as.character(index)
  )
]

pct[
  ,
  index := sub(
    "_direct_monthly$",
    "",
    as.character(index)
  )
]


clim[
  ,
  month := as.character(month)
]

pct[
  ,
  month := as.character(month)
]


# =============================================================================
# 8. EXTRACT FIVE ORIGINAL MONTHLY-DERIVABLE INDICES
# =============================================================================

BASE_INDICES <- c(
  "TXx",
  "TXn",
  "TNx",
  "TNn",
  "DTR"
)

base <- clim[
  index %in% BASE_INDICES &
    month %in% MONTH_ORDER,
  .(
    index,
    month,
    unit,
    N,
    Mean,
    Median,
    Q1,
    Q3
  )
]


# =============================================================================
# 9. EXTRACT TX90p + TN10p
# =============================================================================

percentile <- pct[
  index %in% c(
    "TX90p",
    "TN10p"
  ) &
    month %in% MONTH_ORDER,
  .(
    index,
    month,
    unit,
    N,
    Mean,
    Median,
    Q1,
    Q3
  )
]


# =============================================================================
# 10. COMBINE
# =============================================================================

fig <- rbind(
  base,
  percentile,
  fill = TRUE
)


# =============================================================================
# 11. FINAL QC
# =============================================================================

cat("\n============================================================\n")
cat("QUALITY CONTROL\n")
cat("============================================================\n")


# ---- 11.1 Exactly seven indices ---------------------------------------------

found_indices <- unique(
  fig$index
)

cat(
  "\nIndices found:\n",
  paste(
    found_indices,
    collapse = ", "
  ),
  "\n"
)


if (
  !setequal(
    found_indices,
    INDEX_ORDER
  )
) {
  
  stop(
    "\nQC FAILED: Figure 1 does not contain exactly the required 7 indices.\n",
    "Expected: ",
    paste(INDEX_ORDER, collapse = ", "),
    "\nFound: ",
    paste(found_indices, collapse = ", ")
  )
  
}


# ---- 11.2 Exactly 21 rows ---------------------------------------------------

expected_rows <- 7L * 3L

if (nrow(fig) != expected_rows) {
  
  stop(
    "\nQC FAILED: Expected ",
    expected_rows,
    " rows = 7 indices × 3 months.\n",
    "Found: ",
    nrow(fig)
  )
  
}


# ---- 11.3 Every index has April, May, June -------------------------------

coverage <- fig[
  ,
  .(
    n_months = uniqueN(month)
  ),
  by = index
]

print(coverage)

if (
  any(
    coverage$n_months != 3L
  )
) {
  
  stop(
    "\nQC FAILED: One or more indices do not contain all three months."
  )
  
}


# ---- 11.4 No missing plotted statistics ------------------------------------

bad_values <- fig[
  !is.finite(Mean) |
    !is.finite(Q1) |
    !is.finite(Q3),
  .N
]

if (bad_values > 0L) {
  
  stop(
    "\nQC FAILED: ",
    bad_values,
    " rows contain missing/non-finite Mean, Q1 or Q3."
  )
  
}


# ---- 11.5 Q1 <= Median <= Q3 -----------------------------------------------

bad_order <- fig[
  Median < Q1 |
    Median > Q3 |
    Q1 > Q3,
  .N
]

if (bad_order > 0L) {
  
  stop(
    "\nQC FAILED: Quartile/median ordering problem in ",
    bad_order,
    " rows."
  )
  
}


# ---- 11.6 Units -------------------------------------------------------------

fig[
  ,
  expected_unit := INDEX_UNITS[index]
]

bad_units <- fig[
  unit != expected_unit,
  .N
]

if (bad_units > 0L) {
  
  cat(
    "\nWARNING: ",
    bad_units,
    " unit labels differ from expected labels.\n"
  )
  
} else {
  
  cat(
    "\nUnit QC: PASS\n"
  )
  
}


# =============================================================================
# 12. LOCK ORDER
# =============================================================================

fig[
  ,
  index := factor(
    index,
    levels = INDEX_ORDER
  )
]

fig[
  ,
  month := factor(
    month,
    levels = MONTH_ORDER
  )
]

setorder(
  fig,
  index,
  month
)


# =============================================================================
# 13. PANEL LABELS
# =============================================================================

PANEL_LABELS <- c(
  
  TXx   = "(a) TXx",
  TXn   = "(b) TXn",
  TNx   = "(c) TNx",
  TNn   = "(d) TNn",
  
  TX90p = "(e) TX90p",
  TN10p = "(f) TN10p",
  
  DTR   = "(g) DTR"
  
)
# =============================================================================
# 14. CREATE FIGURE
# =============================================================================

# Numeric month axis for ribbon/line geometry
fig[
  ,
  month_num := match(
    as.character(month),
    MONTH_ORDER
  )
]

# Make sure final index order is locked
fig[
  ,
  index := factor(
    as.character(index),
    levels = INDEX_ORDER
  )
]

# Panel labels
PANEL_LABELS <- c(
  TXx   = "(a) TXx",
  TXn   = "(b) TXn",
  TNx   = "(c) TNx",
  TNn   = "(d) TNn",
  TX90p = "(e) TX90p",
  TN10p = "(f) TN10p",
  DTR   = "(g) DTR"
)

fig[
  ,
  panel_label := PANEL_LABELS[
    as.character(index)
  ]
]


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

p <- ggplot(
  fig,
  aes(
    x = month_num,
    y = Mean,
    group = 1
  )
) +
  
  # ---------------------------------------------------------------------------
# IQR ribbon
# IMPORTANT:
# x is explicitly supplied here because inherit.aes = FALSE
# ---------------------------------------------------------------------------

geom_ribbon(
  data = fig,
  aes(
    x = month_num,
    ymin = Q1,
    ymax = Q3,
    group = 1
  ),
  fill = "#9FB3C8",
  alpha = 0.35,
  inherit.aes = FALSE
) +
  
  # ---------------------------------------------------------------------------
# Q1 boundary
# ---------------------------------------------------------------------------

geom_line(
  aes(
    x = month_num,
    y = Q1,
    group = 1
  ),
  colour = "#9FB3C8",
  linewidth = 0.35
) +
  
  # ---------------------------------------------------------------------------
# Q3 boundary
# ---------------------------------------------------------------------------

geom_line(
  aes(
    x = month_num,
    y = Q3,
    group = 1
  ),
  colour = "#9FB3C8",
  linewidth = 0.35
) +
  
  # ---------------------------------------------------------------------------
# Mean
# ---------------------------------------------------------------------------

geom_line(
  aes(
    x = month_num,
    y = Mean,
    group = 1
  ),
  colour = "#1F3B57",
  linewidth = 0.9
) +
  
  geom_point(
    aes(
      x = month_num,
      y = Mean
    ),
    shape = 21,
    fill = "white",
    colour = "#1F3B57",
    size = 2.8,
    stroke = 0.9
  ) +
  
  # ---------------------------------------------------------------------------
# Mean labels
# ---------------------------------------------------------------------------

geom_text(
  aes(
    x = month_num,
    y = Mean,
    label = sprintf("%.1f", Mean)
  ),
  vjust = -1.2,
  colour = "#1F3B57",
  size = 3.1
) +
  
  # ---------------------------------------------------------------------------
# Seven panels
# ---------------------------------------------------------------------------

facet_wrap(
  ~ panel_label,
  ncol = 4,
  nrow = 2,
  scales = "free_y"
) +
  
  # ---------------------------------------------------------------------------
# Month axis
# ---------------------------------------------------------------------------

scale_x_continuous(
  breaks = c(1, 2, 3),
  labels = c(
    "Apr",
    "May",
    "Jun"
  ),
  expand = expansion(
    mult = c(
      0.12,
      0.12
    )
  )
) +
  
  # ---------------------------------------------------------------------------
# Y axis
# ---------------------------------------------------------------------------

scale_y_continuous(
  breaks = scales::pretty_breaks(
    n = 5
  ),
  expand = expansion(
    mult = c(
      0.05,
      0.12
    )
  )
) +
  
  # ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------

labs(
  title =
    "Seasonal climatology of temperature extremes during the pre-monsoon season",
  
  subtitle =
    "Domain mean over 283 IMD 1° land grid cells, April–June 1996–2025; shading spans the interquartile range across grid-years",
  
  x = NULL,
  y = NULL
) +
  
  # ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------

theme_minimal(
  base_size = 11
) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 17,
        colour = "grey10",
        hjust = 0.5
      ),
    
    plot.subtitle =
      element_text(
        size = 10.5,
        colour = "grey35",
        hjust = 0.5,
        margin = margin(
          b = 12
        )
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11.5,
        colour = "grey15"
      ),
    
    panel.grid =
      element_blank(),
    
    panel.border =
      element_rect(
        colour = "grey80",
        fill = NA,
        linewidth = 0.4
      ),
    
    axis.text =
      element_text(
        colour = "grey25",
        size = 9.5
      ),
    
    axis.title =
      element_text(
        colour = "grey25"
      ),
    
    axis.ticks =
      element_line(
        colour = "grey70",
        linewidth = 0.3
      ),
    
    plot.margin =
      margin(
        10,
        12,
        10,
        12
      )
  )


# =============================================================================
# 15. DISPLAY
# =============================================================================

print(p)

fig[
  ,
  panel_label := PANEL_LABELS[
    as.character(index)
  ]
]



# =============================================================================
# 15. DISPLAY
# =============================================================================

print(p)


# =============================================================================
# 16. SAVE PNG
# =============================================================================

PNG_FILE <- file.path(
  FIG_DIR,
  "FIG01_SEASONAL_CLIMATOLOGY_AMJ_1996_2025_7INDEX.png"
)

ggsave(
  filename = PNG_FILE,
  plot = p,
  width = 12,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# =============================================================================
# 17. SAVE PDF
# =============================================================================

PDF_FILE <- file.path(
  FIG_DIR,
  "FIG01_SEASONAL_CLIMATOLOGY_AMJ_1996_2025_7INDEX.pdf"
)

ggsave(
  filename = PDF_FILE,
  plot = p,
  width = 12,
  height = 7.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 18. SAVE EXACT PLOTTED DATA
# =============================================================================

PLOT_DATA_FILE <- file.path(
  FIG_DIR,
  "FIG01_SEASONAL_CLIMATOLOGY_AMJ_1996_2025_7INDEX_PLOTTED_DATA.csv"
)

plot_data <- copy(fig)

plot_data[
  ,
  index := as.character(index)
]

plot_data[
  ,
  month := as.character(month)
]

plot_data[
  ,
  expected_unit := NULL
]

fwrite(
  plot_data,
  PLOT_DATA_FILE
)


# =============================================================================
# 19. SAVE QC
# =============================================================================

QC_FILE <- file.path(
  FIG_DIR,
  "FIG01_SEASONAL_CLIMATOLOGY_AMJ_1996_2025_7INDEX_QC.csv"
)

QC <- data.table(
  
  check = c(
    "Number of indices",
    "Number of rows",
    "Months per index",
    "Missing plotted values",
    "Quartile ordering",
    "WSDI handling",
    "TX90p/TN10p handling",
    "Raw inputs modified"
  ),
  
  result = c(
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "BLOCKED FROM FIGURE 1",
    "DIRECT MONTHLY CLIMATOLOGY",
    "NO"
  ),
  
  detail = c(
    
    "7 indices",
    
    "21 rows = 7 indices × 3 months",
    
    "April, May and June for every index",
    
    "Mean, Q1 and Q3 are finite for all plotted rows",
    
    "Q1 <= Median <= Q3 for all rows",
    
    "Continuous 91-day AMJ WSDI; no monthly reset",
    
    "Locked 1981–2010 calendar-day thresholds; direct monthly exceedance; no Zhang bootstrap",
    
    "All input STEP08 files are read only"
    
  )
  
)

fwrite(
  QC,
  QC_FILE
)


# =============================================================================
# 20. FINAL OUTPUT
# =============================================================================

cat("\n============================================================\n")
cat("FIGURE 1 SUCCESSFULLY CREATED\n")
cat("============================================================\n")

cat("\nIndices included:\n")

cat(
  paste(
    INDEX_ORDER,
    collapse = " | "
  ),
  "\n"
)

cat(
  "\nTotal panels: 7\n"
)

cat(
  "Total plotted observations: 21 (7 × 3 months)\n"
)

cat(
  "\nWSDI: NOT included — remains in Figure 5\n"
)

cat(
  "TX90p/TN10p: included using existing direct monthly climatology\n"
)

cat(
  "\nPNG:\n",
  PNG_FILE,
  "\n"
)

cat(
  "\nPDF:\n",
  PDF_FILE,
  "\n"
)

cat(
  "\nPlotted data:\n",
  PLOT_DATA_FILE,
  "\n"
)

cat(
  "\nQC:\n",
  QC_FILE,
  "\n"
)

cat("\n============================================================\n")