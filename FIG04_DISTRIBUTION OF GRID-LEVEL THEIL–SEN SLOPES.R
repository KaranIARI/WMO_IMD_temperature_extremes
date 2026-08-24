# =====================================================================
# FIGURE 4 — DISTRIBUTION OF GRID-LEVEL THEIL–SEN SLOPES
# AMJ 1996–2025 | 283 IMD 1° × 1° grids
#
# FINAL REVIEWED VERSION
#
# Layout:
#   TXx      TXn      TNx      TNn
#   TX90p    TN10p    WSDI     DTR
#
# Data:
#   TXx, TXn, TNx, TNn, WSDI, DTR
#       -> AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv
#
#   TX90p, TN10p
#       -> AMJ_STEP07B_trends_PRIMARY_1996_2025.csv
#
# IMPORTANT:
#   - No trend is recalculated.
#   - No p/q/FDR is recalculated.
#   - Only locked primary Theil–Sen slopes are visualized.
#   - TXn has 3 test-undefined grids.
#   - Therefore TXn has 280 testable slopes.
#   - All other indices have 283 testable slopes.
# =====================================================================


# =====================================================================
# 0. CLEAN SESSION
# =====================================================================

rm(list = ls())


# =====================================================================
# 1. PACKAGES
# =====================================================================

library(data.table)
library(ggplot2)


# =====================================================================
# 2. PATHS
# =====================================================================

DATA_DIR <- "F:/WMO_IMD_R/WMO_IMD/data"

F1 <- file.path(
  DATA_DIR,
  "AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv"
)

F2 <- file.path(
  DATA_DIR,
  "AMJ_STEP07B_trends_PRIMARY_1996_2025.csv"
)

OUT_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =====================================================================
# 3. LOCKED SETTINGS
# =====================================================================

START_YEAR <- 1996L
END_YEAR   <- 2025L

PRIMARY_VARIANT <- "HR_lag3"

INDEX_ORDER <- c(
  "TXx",
  "TXn",
  "TNx",
  "TNn",
  "TX90p",
  "TN10p",
  "WSDI",
  "DTR"
)

INDEX_UNITS <- c(
  TXx   = "°C decade⁻¹",
  TXn   = "°C decade⁻¹",
  TNx   = "°C decade⁻¹",
  TNn   = "°C decade⁻¹",
  TX90p = "% decade⁻¹",
  TN10p = "% decade⁻¹",
  WSDI  = "days decade⁻¹",
  DTR   = "°C decade⁻¹"
)


# =====================================================================
# 4. EXPECTED TESTABLE GRID COUNTS
#
# TXn:
#   283 total grids
#   3 test-undefined
#   280 testable
#
# All others:
#   283 testable
# =====================================================================

EXPECTED_TESTABLE <- c(
  TXx   = 283L,
  TXn   = 280L,
  TNx   = 283L,
  TNn   = 283L,
  TX90p = 283L,
  TN10p = 283L,
  WSDI  = 283L,
  DTR   = 283L
)


# =====================================================================
# 5. INPUT FILE CHECK
# =====================================================================

if (!file.exists(F1)) {
  
  stop(
    "\nSTEP07 input file does not exist:\n",
    F1
  )
}

if (!file.exists(F2)) {
  
  stop(
    "\nSTEP07B input file does not exist:\n",
    F2
  )
}


# =====================================================================
# 6. READ DATA
# =====================================================================

cat("\n")
cat("============================================================\n")
cat("READING LOCKED PRIMARY TREND DATA\n")
cat("============================================================\n")

d1 <- fread(F1)
d2 <- fread(F2)

cat(
  "STEP07 rows  : ",
  nrow(d1),
  "\n",
  sep = ""
)

cat(
  "STEP07B rows : ",
  nrow(d2),
  "\n",
  sep = ""
)


# =====================================================================
# 7. REQUIRED COLUMN CHECK
# =====================================================================

REQUIRED_COLUMNS <- c(
  "grid_id",
  "period",
  "period_role",
  "index",
  "variant",
  "sen_per_decade",
  "testable",
  "is_primary"
)

missing_d1 <- setdiff(
  REQUIRED_COLUMNS,
  names(d1)
)

missing_d2 <- setdiff(
  REQUIRED_COLUMNS,
  names(d2)
)

if (length(missing_d1) > 0L) {
  
  stop(
    "\nSTEP07 is missing required column(s):\n",
    paste(
      missing_d1,
      collapse = ", "
    )
  )
}

if (length(missing_d2) > 0L) {
  
  stop(
    "\nSTEP07B is missing required column(s):\n",
    paste(
      missing_d2,
      collapse = ", "
    )
  )
}

cat(
  "Column QC: PASSED\n"
)


# =====================================================================
# 8. FILTER STEP07
#
# Six indices:
# TXx, TXn, TNx, TNn, WSDI, DTR
# =====================================================================

d1 <- d1[
  period == "TREND_1996_2025" &
    period_role == "PRIMARY" &
    variant == PRIMARY_VARIANT &
    is_primary == TRUE &
    index %in% c(
      "TXx",
      "TXn",
      "TNx",
      "TNn",
      "WSDI",
      "DTR"
    ) &
    testable == TRUE &
    is.finite(sen_per_decade)
]


# =====================================================================
# 9. FILTER STEP07B
#
# Two indices:
# TX90p, TN10p
# =====================================================================

d2 <- d2[
  period == "TREND_1996_2025" &
    period_role == "PRIMARY" &
    variant == PRIMARY_VARIANT &
    is_primary == TRUE &
    index %in% c(
      "TX90p",
      "TN10p"
    ) &
    testable == TRUE &
    is.finite(sen_per_decade)
]


# =====================================================================
# 10. COMBINE
# =====================================================================

d <- rbind(
  d1,
  d2,
  fill = TRUE
)


# =====================================================================
# 11. FORCE INDEX TO CHARACTER
#
# Keep character during QC.
# Convert to factor only immediately before plotting.
# =====================================================================

d[
  ,
  index := as.character(index)
]


# =====================================================================
# 12. CHECK EXACTLY EIGHT INDICES
# =====================================================================

PRESENT_INDICES <- sort(
  unique(d$index)
)

MISSING_INDICES <- setdiff(
  INDEX_ORDER,
  PRESENT_INDICES
)

EXTRA_INDICES <- setdiff(
  PRESENT_INDICES,
  INDEX_ORDER
)

if (length(MISSING_INDICES) > 0L) {
  
  stop(
    "\nQC FAILED — Missing index/indices:\n",
    paste(
      MISSING_INDICES,
      collapse = ", "
    )
  )
}

if (length(EXTRA_INDICES) > 0L) {
  
  stop(
    "\nQC FAILED — Unexpected index/indices:\n",
    paste(
      EXTRA_INDICES,
      collapse = ", "
    )
  )
}

if (length(PRESENT_INDICES) != 8L) {
  
  stop(
    "\nQC FAILED — Expected 8 indices, found ",
    length(PRESENT_INDICES)
  )
}

cat(
  "Index QC: PASSED — all 8 indices present\n"
)


# =====================================================================
# 13. GRID COUNT QC
# =====================================================================

cat("\n")
cat("============================================================\n")
cat("TESTABLE GRID COUNT QC\n")
cat("============================================================\n")

COUNT_TABLE <- d[
  ,
  .(
    n_slope = .N,
    n_grid = uniqueN(grid_id)
  ),
  by = index
]

COUNT_TABLE[
  ,
  expected_testable := EXPECTED_TESTABLE[index]
]

COUNT_TABLE[
  ,
  difference := n_grid - expected_testable
]

COUNT_TABLE <- COUNT_TABLE[
  order(
    match(
      index,
      INDEX_ORDER
    )
  )
]

print(COUNT_TABLE)


# Check each index individually
for (idx in INDEX_ORDER) {
  
  actual_n <- COUNT_TABLE[
    index == idx,
    n_grid
  ]
  
  expected_n <- EXPECTED_TESTABLE[idx]
  
  if (length(actual_n) != 1L) {
    
    stop(
      "\nQC FAILED — count missing for index: ",
      idx
    )
  }
  
  if (actual_n != expected_n) {
    
    stop(
      "\nQC FAILED — ",
      idx,
      ": expected ",
      expected_n,
      " testable grids but found ",
      actual_n
    )
  }
}

cat(
  "\nGrid-count QC: PASSED\n"
)

cat(
  "TXn correctly contains 280 testable grids.\n"
)

cat(
  "The other seven indices correctly contain 283 testable grids.\n"
)


# =====================================================================
# 14. DUPLICATE GRID-INDEX QC
# =====================================================================

DUPLICATES <- d[
  ,
  .N,
  by = .(
    index,
    grid_id
  )
][
  N > 1
]

if (nrow(DUPLICATES) > 0L) {
  
  cat("\nDuplicate records detected:\n")
  print(DUPLICATES)
  
  stop(
    "\nQC FAILED — duplicate index-grid combinations detected."
  )
}

cat(
  "Duplicate QC: PASSED\n"
)


# =====================================================================
# 15. SLOPE FINITENESS QC
# =====================================================================

if (
  any(
    !is.finite(
      d$sen_per_decade
    )
  )
) {
  
  stop(
    "\nQC FAILED — non-finite Sen slope detected."
  )
}

cat(
  "Slope finiteness QC: PASSED\n"
)


# =====================================================================
# 16. CREATE INDEX FACTOR
#
# This is done ONLY after all character-based QC.
# =====================================================================

d[
  ,
  index := factor(
    index,
    levels = INDEX_ORDER
  )
]


# =====================================================================
# 17. CALCULATE MEDIAN SEN SLOPE
#
# No new trend calculation.
# This is only a descriptive summary of existing slopes.
# =====================================================================

MEDIAN_TABLE <- d[
  ,
  .(
    n_grid = uniqueN(grid_id),
    median_sen_slope = median(
      sen_per_decade,
      na.rm = TRUE
    ),
    min_sen_slope = min(
      sen_per_decade,
      na.rm = TRUE
    ),
    max_sen_slope = max(
      sen_per_decade,
      na.rm = TRUE
    )
  ),
  by = index
]

MEDIAN_TABLE[
  ,
  index := factor(
    index,
    levels = INDEX_ORDER
  )
]

MEDIAN_TABLE[
  ,
  unit := INDEX_UNITS[
    as.character(index)
  ]
]

MEDIAN_TABLE <- MEDIAN_TABLE[
  order(index)
]

cat("\n")
cat("============================================================\n")
cat("SEN SLOPE SUMMARY\n")
cat("============================================================\n")

print(
  MEDIAN_TABLE
)


# =====================================================================
# 18. PREPARE MEDIAN LINES FOR FIGURE
# =====================================================================

MEDIAN_LINES <- MEDIAN_TABLE[
  ,
  .(
    index = index,
    median_slope = median_sen_slope
  )
]


# =====================================================================
# 19. FIGURE 4
#
# 2 × 4:
#
# TXx      TXn      TNx      TNn
# TX90p    TN10p    WSDI     DTR
#
# Each panel has its own x-axis scale because the units/magnitudes
# differ between temperature, percentile and duration indices.
# =====================================================================

p <- ggplot(
  d,
  aes(
    x = sen_per_decade
  )
) +
  
  # ---------------------------------------------------------------
# Density distribution
# ---------------------------------------------------------------
geom_density(
  fill = "grey85",
  colour = "black",
  linewidth = 0.45,
  adjust = 1.05,
  na.rm = TRUE
) +
  
  # ---------------------------------------------------------------
# Zero-slope reference
# ---------------------------------------------------------------
geom_vline(
  xintercept = 0,
  linetype = "dashed",
  linewidth = 0.45
) +
  
  # ---------------------------------------------------------------
# Median Sen slope
# ---------------------------------------------------------------
geom_vline(
  data = MEDIAN_LINES,
  aes(
    xintercept = median_slope
  ),
  linewidth = 0.70,
  inherit.aes = FALSE
) +
  
  # ---------------------------------------------------------------
# 2 × 4 layout
# ---------------------------------------------------------------
facet_wrap(
  ~ index,
  ncol = 4,
  nrow = 2,
  scales = "free_x"
) +
  
  # ---------------------------------------------------------------
# Labels
# ---------------------------------------------------------------
labs(
  title = paste0(
    "Distribution of grid-level Theil–Sen slopes during AMJ, ",
    START_YEAR,
    "–",
    END_YEAR
  ),
  x = "Theil–Sen slope per decade",
  y = "Density"
) +
  
  # ---------------------------------------------------------------
# Theme
# ---------------------------------------------------------------
theme_classic(
  base_size = 11
) +
  
  theme(
    
    # Main title
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    # Panel titles
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    # Axis titles
    axis.title.x = element_text(
      face = "bold",
      size = 10.5
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 10.5
    ),
    
    # Axis text
    axis.text = element_text(
      colour = "black",
      size = 9
    ),
    
    # Panel borders
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    
    # No grid
    panel.grid = element_blank(),
    
    # Panel spacing
    panel.spacing = grid::unit(
      0.9,
      "lines"
    ),
    
    # Margins
    plot.margin = margin(
      10,
      12,
      10,
      12
    )
  )


# =====================================================================
# 20. SAVE 600-DPI PNG
# =====================================================================

PNG_FILE <- file.path(
  OUT_DIR,
  "FIGURE_4_SEN_SLOPE_DISTRIBUTIONS_AMJ_1996_2025.png"
)

ggsave(
  filename = PNG_FILE,
  plot = p,
  width = 10,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


# =====================================================================
# 21. SAVE 600-DPI TIFF
# =====================================================================

TIFF_FILE <- file.path(
  OUT_DIR,
  "FIGURE_4_SEN_SLOPE_DISTRIBUTIONS_AMJ_1996_2025.tiff"
)

ggsave(
  filename = TIFF_FILE,
  plot = p,
  width = 10,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =====================================================================
# 22. SAVE VECTOR PDF
# =====================================================================

PDF_FILE <- file.path(
  OUT_DIR,
  "FIGURE_4_SEN_SLOPE_DISTRIBUTIONS_AMJ_1996_2025.pdf"
)

ggsave(
  filename = PDF_FILE,
  plot = p,
  width = 10,
  height = 7.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =====================================================================
# 23. SAVE SUMMARY CSV
# =====================================================================

SUMMARY_FILE <- file.path(
  OUT_DIR,
  "FIGURE_4_SEN_SLOPE_SUMMARY_1996_2025.csv"
)

SUMMARY_OUT <- copy(
  MEDIAN_TABLE
)

SUMMARY_OUT[
  ,
  index := as.character(index)
]

setcolorder(
  SUMMARY_OUT,
  c(
    "index",
    "unit",
    "n_grid",
    "median_sen_slope",
    "min_sen_slope",
    "max_sen_slope"
  )
)

fwrite(
  SUMMARY_OUT,
  SUMMARY_FILE
)


# =====================================================================
# 24. DISPLAY FIGURE
# =====================================================================

print(p)


# =====================================================================
# 25. FINAL QC / OUTPUT MESSAGE
# =====================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 4 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat(
  "Period       : AMJ ",
  START_YEAR,
  "–",
  END_YEAR,
  "\n",
  sep = ""
)

cat(
  "Layout       : 2 × 4\n"
)

cat(
  "Indices      : 8\n"
)

cat(
  "TXn          : 280 testable / 3 test-undefined\n"
)

cat(
  "Other 7      : 283 testable each\n"
)

cat("------------------------------------------------------------\n")

cat(
  "PNG  : ",
  PNG_FILE,
  "\n",
  sep = ""
)

cat(
  "TIFF : ",
  TIFF_FILE,
  "\n",
  sep = ""
)

cat(
  "PDF  : ",
  PDF_FILE,
  "\n",
  sep = ""
)

cat(
  "DATA : ",
  SUMMARY_FILE,
  "\n",
  sep = ""
)

cat("============================================================\n")