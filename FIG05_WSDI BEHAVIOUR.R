# =====================================================================
# FIGURE 5 — WSDI BEHAVIOUR
# =====================================================================
# AMJ 1996–2025 | 283 IMD 1° × 1° land grids
#
# PANEL A
#   WSDI zero-inflation / frequency distribution
#
# PANEL B
#   Distribution of grid-level Theil–Sen slopes
#   expressed as percentage of grid cells
#
# DATA
#   WSDI observations:
#     AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet
#
#   Locked WSDI trends:
#     AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv
#
# IMPORTANT
#   - No WSDI values are recalculated.
#   - No Sen slopes are recalculated.
#   - No p-values are recalculated.
#   - No q-values are recalculated.
#   - No FDR procedure is rerun.
#   - Figure uses only locked STEP06 / STEP07 results.
#
# OUTPUT
#   PDF  -> separate PDF folder
#   TIFF -> separate TIFF folder
# =====================================================================


# =====================================================================
# 0. CLEAN SESSION
# =====================================================================

rm(list = ls())

gc()


# =====================================================================
# 1. PACKAGES
# =====================================================================

library(data.table)
library(ggplot2)
library(arrow)
library(patchwork)


# =====================================================================
# 2. PATHS
# =====================================================================

DATA_DIR <- "F:/WMO_IMD_R/WMO_IMD/data"

# -------------------------------------------------------------
# STEP06 WSDI observation file
# -------------------------------------------------------------

WSDI_FILE <- file.path(
  DATA_DIR,
  "AMJ_STEP06_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.parquet"
)

# -------------------------------------------------------------
# STEP07 locked trend file
# -------------------------------------------------------------

TREND_FILE <- file.path(
  DATA_DIR,
  "AMJ_STEP07_ALL8_trends_PRIMARY_1996_2025.csv"
)

# -------------------------------------------------------------
# Output folders
# -------------------------------------------------------------

PDF_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08/Figure5/PDF"

TIF_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08/Figure5/TIFF"


dir.create(
  PDF_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TIF_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =====================================================================
# 3. LOCKED SETTINGS
# =====================================================================

START_YEAR <- 1996L
END_YEAR   <- 2025L

EXPECTED_GRIDS <- 283L

PRIMARY_VARIANT <- "HR_lag3"


# =====================================================================
# 4. CHECK INPUT FILES
# =====================================================================

if (!file.exists(WSDI_FILE)) {
  
  stop(
    "\nWSDI observation file not found:\n",
    WSDI_FILE
  )
}


if (!file.exists(TREND_FILE)) {
  
  stop(
    "\nSTEP07 trend file not found:\n",
    TREND_FILE
  )
}


# =====================================================================
# 5. READ STEP06 WSDI OBSERVATIONS
# =====================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 5 — WSDI BEHAVIOUR\n")
cat("============================================================\n")

cat("\nReading STEP06 WSDI observations...\n")


wsdi_raw <- as.data.table(
  read_parquet(WSDI_FILE)
)


cat(
  "Rows read:",
  nrow(wsdi_raw),
  "\n"
)


cat(
  "Available columns:\n"
)

print(
  names(wsdi_raw)
)


# =====================================================================
# 6. IDENTIFY WSDI COLUMN
# =====================================================================

wsdi_candidates <- names(wsdi_raw)[
  toupper(names(wsdi_raw)) == "WSDI"
]


if (length(wsdi_candidates) != 1L) {
  
  wsdi_candidates <- names(wsdi_raw)[
    grepl(
      "^WSDI$|WSDI",
      names(wsdi_raw),
      ignore.case = TRUE
    )
  ]
}


if (length(wsdi_candidates) != 1L) {
  
  stop(
    "\nCould not uniquely identify WSDI column.\n",
    "Candidates found:\n",
    paste(
      wsdi_candidates,
      collapse = ", "
    )
  )
}


WSDI_COL <- wsdi_candidates[1]


cat(
  "\nWSDI column:",
  WSDI_COL,
  "\n"
)


# =====================================================================
# 7. IDENTIFY YEAR COLUMN
# =====================================================================

year_candidates <- names(wsdi_raw)[
  tolower(names(wsdi_raw)) %in%
    c("year", "yr")
]


if (length(year_candidates) < 1L) {
  
  year_candidates <- names(wsdi_raw)[
    grepl(
      "year|yr",
      names(wsdi_raw),
      ignore.case = TRUE
    )
  ]
}


if (length(year_candidates) < 1L) {
  
  stop(
    "\nCould not identify year column."
  )
}


YEAR_COL <- year_candidates[1]


cat(
  "Year column:",
  YEAR_COL,
  "\n"
)


# =====================================================================
# 8. STANDARDISE WSDI DATA
# =====================================================================

wsdi <- wsdi_raw[
  ,
  .(
    year = as.integer(
      get(YEAR_COL)
    ),
    
    WSDI = as.numeric(
      get(WSDI_COL)
    )
  )
]


# =====================================================================
# 9. FILTER AMJ 1996–2025
# =====================================================================

wsdi <- wsdi[
  year >= START_YEAR &
    year <= END_YEAR
]


# Keep valid WSDI observations only

wsdi <- wsdi[
  is.finite(WSDI)
]


cat(
  "\nValid WSDI grid-year observations:",
  nrow(wsdi),
  "\n"
)


# =====================================================================
# 10. WSDI BASIC QC
# =====================================================================

if (nrow(wsdi) == 0L) {
  
  stop(
    "\nQC FAILED: No valid WSDI observations."
  )
}


if (any(wsdi$WSDI < 0)) {
  
  stop(
    "\nQC FAILED: Negative WSDI values detected."
  )
}


# WSDI is a number of days

non_integer <- wsdi[
  WSDI != round(WSDI)
]


if (nrow(non_integer) > 0L) {
  
  warning(
    "\n",
    nrow(non_integer),
    " non-integer WSDI observations detected."
  )
}


# =====================================================================
# 11. ZERO-INFLATION SUMMARY
# =====================================================================

N_VALID <- nrow(wsdi)


N_ZERO <- sum(
  wsdi$WSDI == 0,
  na.rm = TRUE
)


PCT_ZERO <- (
  N_ZERO / N_VALID
) * 100


WSDI_MEAN <- mean(
  wsdi$WSDI,
  na.rm = TRUE
)


WSDI_MEDIAN <- median(
  wsdi$WSDI,
  na.rm = TRUE
)


WSDI_IQR <- IQR(
  wsdi$WSDI,
  na.rm = TRUE
)


WSDI_MIN <- min(
  wsdi$WSDI,
  na.rm = TRUE
)


WSDI_MAX <- max(
  wsdi$WSDI,
  na.rm = TRUE
)


cat("\n")
cat("------------------------------------------------------------\n")
cat("WSDI ZERO-INFLATION SUMMARY\n")
cat("------------------------------------------------------------\n")

cat(
  "Valid observations :",
  N_VALID,
  "\n"
)

cat(
  "Zero observations  :",
  N_ZERO,
  "\n"
)

cat(
  "Zero percentage    :",
  round(PCT_ZERO, 2),
  "%\n"
)

cat(
  "Mean               :",
  round(WSDI_MEAN, 2),
  "days\n"
)

cat(
  "Median             :",
  WSDI_MEDIAN,
  "days\n"
)

cat(
  "IQR                :",
  WSDI_IQR,
  "days\n"
)

cat(
  "Range              :",
  WSDI_MIN,
  "–",
  WSDI_MAX,
  "days\n"
)


# =====================================================================
# 12. PANEL A DATA
# =====================================================================

wsdi_freq <- wsdi[
  ,
  .N,
  by = WSDI
][
  order(WSDI)
]


wsdi_freq[
  ,
  percentage := (
    N / sum(N)
  ) * 100
]


# Create complete integer sequence

all_days <- data.table(
  WSDI = seq(
    floor(min(wsdi_freq$WSDI)),
    ceiling(max(wsdi_freq$WSDI)),
    by = 1
  )
)


wsdi_freq <- merge(
  all_days,
  wsdi_freq,
  by = "WSDI",
  all.x = TRUE
)


wsdi_freq[
  is.na(N),
  `:=`(
    N = 0L,
    percentage = 0
  )
]


# =====================================================================
# 13. PANEL A — WSDI ZERO-INFLATION
# =====================================================================

zero_label_y <- wsdi_freq[
  WSDI == 0,
  percentage
]


pA <- ggplot(
  wsdi_freq,
  aes(
    x = WSDI,
    y = percentage
  )
) +
  
  # All WSDI values
  geom_col(
    fill = "grey70",
    colour = "black",
    linewidth = 0.25,
    width = 0.82
  ) +
  
  # Explicitly highlight WSDI = 0
  geom_col(
    data = wsdi_freq[
      WSDI == 0
    ],
    fill = "grey35",
    colour = "black",
    linewidth = 0.3,
    width = 0.82
  ) +
  
  # Zero percentage annotation
  annotate(
    "text",
    x = 0.7,
    y = zero_label_y +
      max(wsdi_freq$percentage) * 0.04,
    
    label = paste0(
      round(PCT_ZERO, 1),
      "% of observations = 0 days"
    ),
    
    fontface = "bold",
    size = 4.0,
    hjust = 0
  ) +
  
  labs(
    title = "A. WSDI zero-inflation",
    x = "WSDI (days)",
    y = "Grid-year observations (%)"
  ) +
  
  scale_x_continuous(
    breaks = pretty(
      wsdi_freq$WSDI,
      n = 8
    ),
    
    expand = expansion(
      mult = c(
        0.01,
        0.03
      )
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.10
      )
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      colour = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )


# =====================================================================
# 14. READ LOCKED STEP07 TREND RESULTS
# =====================================================================

cat("\n")
cat("Reading locked STEP07 WSDI trend results...\n")


trend <- fread(
  TREND_FILE
)


required_trend_cols <- c(
  "grid_id",
  "period",
  "variant",
  "index",
  "sen_per_decade",
  "testable",
  "is_primary"
)


missing_trend <- setdiff(
  required_trend_cols,
  names(trend)
)


if (length(missing_trend) > 0L) {
  
  stop(
    "\nMissing STEP07 columns:\n",
    paste(
      missing_trend,
      collapse = ", "
    )
  )
}


# =====================================================================
# 15. FILTER LOCKED PRIMARY WSDI RESULTS
# =====================================================================

wsdi_trend <- trend[
  period == "TREND_1996_2025" &
    variant == PRIMARY_VARIANT &
    index == "WSDI" &
    is_primary == TRUE &
    testable == TRUE &
    is.finite(sen_per_decade)
]


cat(
  "Testable WSDI trend rows:",
  nrow(wsdi_trend),
  "\n"
)


N_TREND_GRIDS <- uniqueN(
  wsdi_trend$grid_id
)


cat(
  "Unique WSDI trend grids:",
  N_TREND_GRIDS,
  "\n"
)


# =====================================================================
# 16. TREND QC
# =====================================================================

if (N_TREND_GRIDS != EXPECTED_GRIDS) {
  
  warning(
    "\nExpected ",
    EXPECTED_GRIDS,
    " testable grids but found ",
    N_TREND_GRIDS,
    "."
  )
}


# Check duplicate grid IDs

dup_grid <- wsdi_trend[
  ,
  .N,
  by = grid_id
][
  N > 1
]


if (nrow(dup_grid) > 0L) {
  
  stop(
    "\nQC FAILED: Duplicate WSDI grid-level trend records detected."
  )
}


# =====================================================================
# 17. SEN-SLOPE SUMMARY
# =====================================================================

N_ZERO_SLOPE <- sum(
  wsdi_trend$sen_per_decade == 0,
  na.rm = TRUE
)


N_POSITIVE_SLOPE <- sum(
  wsdi_trend$sen_per_decade > 0,
  na.rm = TRUE
)


N_NEGATIVE_SLOPE <- sum(
  wsdi_trend$sen_per_decade < 0,
  na.rm = TRUE
)


MEDIAN_SLOPE <- median(
  wsdi_trend$sen_per_decade,
  na.rm = TRUE
)


MIN_SLOPE <- min(
  wsdi_trend$sen_per_decade,
  na.rm = TRUE
)


MAX_SLOPE <- max(
  wsdi_trend$sen_per_decade,
  na.rm = TRUE
)


cat("\n")
cat("------------------------------------------------------------\n")
cat("WSDI SEN-SLOPE SUMMARY\n")
cat("------------------------------------------------------------\n")

cat(
  "Zero slopes     :",
  N_ZERO_SLOPE,
  "\n"
)

cat(
  "Positive slopes :",
  N_POSITIVE_SLOPE,
  "\n"
)

cat(
  "Negative slopes :",
  N_NEGATIVE_SLOPE,
  "\n"
)

cat(
  "Median slope    :",
  MEDIAN_SLOPE,
  "days decade^-1\n"
)

cat(
  "Range           :",
  MIN_SLOPE,
  "to",
  MAX_SLOPE,
  "days decade^-1\n"
)


# =====================================================================
# 18. PANEL B DATA
# =====================================================================
#
# IMPORTANT:
#
# We use percentage of grid cells rather than density.
#
# This makes the WSDI zero-slope dominance directly interpretable.
#
# Expected:
#
#   282 / 283 = 99.65% zero slope
#     1 / 283 =  0.35% positive slope
#
# =====================================================================

slope_freq <- wsdi_trend[
  ,
  .N,
  by = sen_per_decade
][
  order(sen_per_decade)
]


slope_freq[
  ,
  percentage := (
    N / sum(N)
  ) * 100
]


# =====================================================================
# 19. PANEL B — GRID-LEVEL SEN-SLOPE DISTRIBUTION
# =====================================================================

# Separate zero slopes from non-zero slopes.
#
# This is important because 282 grids have exactly zero slope.
# A normal histogram can visually hide this point mass.


zero_slope_pct <- slope_freq[
  sen_per_decade == 0,
  percentage
]


pB <- ggplot(
  slope_freq,
  aes(
    x = sen_per_decade,
    y = percentage
  )
) +
  
  # Zero slope point mass
  geom_col(
    data = slope_freq[
      sen_per_decade == 0
    ],
    fill = "grey35",
    colour = "black",
    width = 0.04
  ) +
  
  # Non-zero slope values
  geom_point(
    data = slope_freq[
      sen_per_decade != 0
    ],
    aes(
      size = percentage
    ),
    shape = 21,
    fill = "grey70",
    colour = "black",
    stroke = 0.5
  ) +
  
  # Zero reference
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # Median line only if median differs from zero
  {
    if (
      is.finite(MEDIAN_SLOPE) &&
      abs(MEDIAN_SLOPE) > .Machine$double.eps
    ) {
      geom_vline(
        xintercept = MEDIAN_SLOPE,
        linewidth = 0.9
      )
    } else {
      NULL
    }
  } +
  
  # Zero-slope annotation
  annotate(
    "text",
    x = 0.04,
    y = zero_slope_pct * 0.94,
    
    label = paste0(
      N_ZERO_SLOPE,
      " grids (",
      round(zero_slope_pct, 2),
      "%)"
    ),
    
    fontface = "bold",
    size = 4.0,
    hjust = 0
  ) +
  
  labs(
    title = "B. Grid-level Sen-slope distribution",
    x = "Theil–Sen slope (days decade⁻¹)",
    y = "Grid cells (%)"
  ) +
  
  scale_size_continuous(
    guide = "none"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      colour = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )


# =====================================================================
# 20. COMBINE PANELS
# =====================================================================

FIG5 <- (
  pA +
    pB
) +
  
  plot_layout(
    ncol = 2
  ) +
  
  plot_annotation(
    
    title =
      "WSDI behaviour during AMJ, 1996–2025",
    
    subtitle =
      "283 IMD 1° × 1° land grids | zero-inflation and grid-level Theil–Sen slopes",
    
    theme =
      theme(
        
        plot.title = element_text(
          face = "bold",
          size = 18,
          hjust = 0.5
        ),
        
        plot.subtitle = element_text(
          size = 11,
          hjust = 0.5
        )
      )
  )


# =====================================================================
# 21. DISPLAY
# =====================================================================

print(
  FIG5
)


# =====================================================================
# 22. OUTPUT FILE NAMES
# =====================================================================

PDF_FILE <- file.path(
  PDF_DIR,
  "Figure_5_WSDI_behaviour_AMJ_1996_2025.pdf"
)


TIF_FILE <- file.path(
  TIF_DIR,
  "Figure_5_WSDI_behaviour_AMJ_1996_2025.tif"
)


# =====================================================================
# 23. SAVE PDF
# =====================================================================

ggsave(
  
  filename = PDF_FILE,
  
  plot = FIG5,
  
  width = 11,
  
  height = 5.8,
  
  units = "in",
  
  device = cairo_pdf,
  
  bg = "white"
)


# =====================================================================
# 24. SAVE TIFF
# =====================================================================

ggsave(
  
  filename = TIF_FILE,
  
  plot = FIG5,
  
  width = 11,
  
  height = 5.8,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white"
)


# =====================================================================
# 25. FINAL REPORT
# =====================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 5 COMPLETE\n")
cat("============================================================\n")

cat(
  "\nPanel A — WSDI zero inflation\n"
)

cat(
  "Valid grid-year observations :",
  N_VALID,
  "\n"
)

cat(
  "Zero WSDI observations       :",
  N_ZERO,
  "\n"
)

cat(
  "Zero WSDI percentage         :",
  round(PCT_ZERO, 2),
  "%\n"
)

cat(
  "Mean WSDI                    :",
  round(WSDI_MEAN, 2),
  "days\n"
)

cat(
  "Median WSDI                  :",
  WSDI_MEDIAN,
  "days\n"
)

cat(
  "IQR                          :",
  WSDI_IQR,
  "days\n"
)

cat(
  "Range                        :",
  WSDI_MIN,
  "–",
  WSDI_MAX,
  "days\n"
)


cat(
  "\nPanel B — Sen-slope distribution\n"
)

cat(
  "Testable grids               :",
  N_TREND_GRIDS,
  "\n"
)

cat(
  "Zero Sen slopes              :",
  N_ZERO_SLOPE,
  " (",
  round(
    100 * N_ZERO_SLOPE /
      N_TREND_GRIDS,
    2
  ),
  "%)\n",
  sep = ""
)

cat(
  "Positive Sen slopes          :",
  N_POSITIVE_SLOPE,
  "\n"
)

cat(
  "Negative Sen slopes          :",
  N_NEGATIVE_SLOPE,
  "\n"
)

cat(
  "Median Sen slope             :",
  MEDIAN_SLOPE,
  " days decade^-1\n"
)

cat(
  "Slope range                  :",
  MIN_SLOPE,
  " to ",
  MAX_SLOPE,
  " days decade^-1\n",
  sep = ""
)


cat(
  "\nPDF output:\n",
  PDF_FILE,
  "\n"
)

cat(
  "\nTIFF output:\n",
  TIF_FILE,
  "\n"
)

cat("\n============================================================\n")