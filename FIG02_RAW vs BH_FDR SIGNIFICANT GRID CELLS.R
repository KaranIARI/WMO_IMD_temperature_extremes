# ================================================================
# FIGURE 2 — RAW vs BH-FDR SIGNIFICANT GRID CELLS
# Section 3.2
#
# PRIMARY ANALYSIS:
#   Period              : 1996–2025
#   Variant             : Hamed–Rao MMK, lag-3
#   Completeness rule   : 80% (>=24 of 30 years)
#   FDR threshold       : q <= 0.10
#
# Expected primary FDR percentages:
#   TXx    = 2.1%
#   TXn    = 1.8%
#   TNx    = 7.8%
#   TNn    = 0.4%
#   TX90p  = 1.1%
#   TN10p  = 6.0%
#   WSDI   = 0.0%
#   DTR    = 0.4%
# ================================================================

library(data.table)
library(ggplot2)

# ================================================================
# 1. PATHS
# ================================================================

FILE_07 <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP07A_fdr_comparison.csv"

FILE_07B <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP07B_trends_PRIMARY_1996_2025.csv"

OUT_DIR <- "F:/WMO_IMD_R/WMO_IMD"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ================================================================
# 2. LOCK PRIMARY ANALYSIS SETTINGS
# ================================================================

PRIMARY_PERIOD <- "TREND_1996_2025"

PRIMARY_VARIANT <- "HR_lag3"

PRIMARY_COMPLETENESS <- 0.8

PRIMARY_Q <- 0.10

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

# ================================================================
# 3. READ STEP 07A FDR SUMMARY
# ================================================================

d07_all <- fread(FILE_07)

# ------------------------------------------------
# Required columns
# ------------------------------------------------

required_07 <- c(
  "period",
  "index",
  "variant",
  "completeness_rule",
  "q_level",
  "n_tested",
  "n_sig_raw",
  "n_sig_bh"
)

missing_07 <- setdiff(
  required_07,
  names(d07_all)
)

if (length(missing_07) > 0L) {
  
  stop(
    "STEP07A missing required columns: ",
    paste(missing_07, collapse = ", ")
  )
  
}

# ------------------------------------------------
# Select ONLY primary configuration
# ------------------------------------------------

d07 <- d07_all[
  period == PRIMARY_PERIOD &
    variant == PRIMARY_VARIANT &
    completeness_rule == PRIMARY_COMPLETENESS &
    q_level == PRIMARY_Q
]

# ------------------------------------------------
# STEP07A should contain exactly six non-percentile
# indices for the primary configuration
# ------------------------------------------------

STEP07A_INDEXES <- c(
  "TXx",
  "TXn",
  "TNx",
  "TNn",
  "WSDI",
  "DTR"
)

d07 <- d07[
  index %in% STEP07A_INDEXES
]

# ------------------------------------------------
# Strict duplicate check
# ------------------------------------------------

dup_07 <- d07[
  ,
  .N,
  by = index
][N != 1L]

if (nrow(dup_07) > 0L) {
  
  stop(
    "QC FAILED: duplicate STEP07A rows remain for primary configuration: ",
    paste(dup_07$index, collapse = ", ")
  )
  
}

# ------------------------------------------------
# Check expected number of indices
# ------------------------------------------------

if (nrow(d07) != length(STEP07A_INDEXES)) {
  
  stop(
    "QC FAILED: expected ",
    length(STEP07A_INDEXES),
    " STEP07A indices, found ",
    nrow(d07)
  )
  
}

# ------------------------------------------------
# Standardize names
# ------------------------------------------------

d07 <- d07[
  ,
  .(
    index = index,
    n_tested = n_tested,
    n_raw_sig = n_sig_raw,
    n_fdr_sig = n_sig_bh
  )
]

# ================================================================
# 4. READ STEP 07B GRID-LEVEL PRIMARY RESULTS
# ================================================================

d07b_all <- fread(FILE_07B)

# ------------------------------------------------
# Required columns
# ------------------------------------------------

required_07b <- c(
  "grid_id",
  "index",
  "p",
  "q",
  "testable",
  "significant_raw",
  "significant_fdr",
  "n_tested_in_family",
  "variant",
  "is_primary"
)

missing_07b <- setdiff(
  required_07b,
  names(d07b_all)
)

if (length(missing_07b) > 0L) {
  
  stop(
    "STEP07B missing required columns: ",
    paste(missing_07b, collapse = ", ")
  )
  
}

# ------------------------------------------------
# Primary HR lag-3 only
# ------------------------------------------------

d07b <- d07b_all[
  is_primary == TRUE &
    variant == PRIMARY_VARIANT
]

# ------------------------------------------------
# Keep only percentile indices
# ------------------------------------------------

d07b <- d07b[
  index %in% c(
    "TX90p",
    "TN10p"
  )
]

# ================================================================
# 5. VALIDATE STEP07B PRIMARY RESULTS
# ================================================================

if (nrow(d07b) == 0L) {
  
  stop(
    "QC FAILED: no primary TX90p/TN10p records found in STEP07B."
  )
  
}

# ------------------------------------------------
# Calculate summary
# ------------------------------------------------

d07b_summary <- d07b[
  ,
  .(
    n_tested = sum(
      testable,
      na.rm = TRUE
    ),
    
    n_raw_sig = sum(
      significant_raw & testable,
      na.rm = TRUE
    ),
    
    n_fdr_sig = sum(
      significant_fdr & testable,
      na.rm = TRUE
    )
  ),
  by = index
]

# ------------------------------------------------
# Check exactly two percentile indices
# ------------------------------------------------

EXPECTED_07B <- c(
  "TX90p",
  "TN10p"
)

if (
  !setequal(
    d07b_summary$index,
    EXPECTED_07B
  )
) {
  
  stop(
    "QC FAILED: STEP07B did not return exactly TX90p and TN10p."
  )
  
}

# ------------------------------------------------
# Check duplicate summary rows
# ------------------------------------------------

if (
  anyDuplicated(
    d07b_summary$index
  ) > 0L
) {
  
  stop(
    "QC FAILED: duplicate percentile-index summaries."
  )
  
}

# ================================================================
# 6. COMBINE STEP07A + STEP07B
# ================================================================

d <- rbindlist(
  list(
    d07,
    d07b_summary
  ),
  use.names = TRUE,
  fill = TRUE
)

# ================================================================
# 7. CHECK EXACTLY 8 INDICES
# ================================================================

d <- d[
  index %in% INDEX_ORDER
]

if (nrow(d) != 8L) {
  
  stop(
    "QC FAILED: expected 8 indices, found ",
    nrow(d)
  )
  
}

# ------------------------------------------------
# Check every index occurs exactly once
# ------------------------------------------------

dup_final <- d[
  ,
  .N,
  by = index
][N != 1L]

if (nrow(dup_final) > 0L) {
  
  stop(
    "QC FAILED: duplicate final index rows: ",
    paste(
      dup_final$index,
      collapse = ", "
    )
  )
  
}

# ------------------------------------------------
# Order indices
# ------------------------------------------------

d[
  ,
  index := factor(
    index,
    levels = INDEX_ORDER
  )
]

# ================================================================
# 8. FINAL DATA QC
# ================================================================

if (
  any(
    d$n_fdr_sig > d$n_raw_sig,
    na.rm = TRUE
  )
) {
  
  stop(
    "QC FAILED: FDR significant count exceeds raw significant count."
  )
  
}

if (
  any(
    d$n_fdr_sig > d$n_tested,
    na.rm = TRUE
  )
) {
  
  stop(
    "QC FAILED: FDR significant count exceeds tested-grid count."
  )
  
}

if (
  any(
    d$n_raw_sig > d$n_tested,
    na.rm = TRUE
  )
) {
  
  stop(
    "QC FAILED: raw significant count exceeds tested-grid count."
  )
  
}

if (
  any(
    d$n_tested <= 0,
    na.rm = TRUE
  )
) {
  
  stop(
    "QC FAILED: invalid tested-grid count."
  )
  
}

# ================================================================
# 9. CALCULATE PERCENTAGES
# ================================================================

d[
  ,
  raw_pct := 100 * n_raw_sig / n_tested
]

d[
  ,
  fdr_pct := 100 * n_fdr_sig / n_tested
]

# ================================================================
# 10. PRINT FINAL NUMBERS
# ================================================================

cat("\n")
cat("====================================================\n")
cat("FIGURE 2 — FINAL PRIMARY SIGNIFICANCE SUMMARY\n")
cat("====================================================\n")

cat(
  "Period            : ",
  PRIMARY_PERIOD,
  "\n",
  sep = ""
)

cat(
  "Variant           : ",
  PRIMARY_VARIANT,
  "\n",
  sep = ""
)

cat(
  "Completeness      : ",
  PRIMARY_COMPLETENESS * 100,
  "%\n",
  sep = ""
)

cat(
  "FDR threshold     : q <= ",
  PRIMARY_Q,
  "\n",
  sep = ""
)

cat("====================================================\n")

print(
  d[
    ,
    .(
      index,
      n_tested,
      n_raw_sig,
      raw_pct,
      n_fdr_sig,
      fdr_pct
    )
  ][
    order(index)
  ]
)

cat("====================================================\n")

# ================================================================
# 11. EXPECTED PRIMARY FDR VALUES
# ================================================================

EXPECTED_FDR <- data.table(
  index = INDEX_ORDER,
  
  expected_n_fdr = c(
    6,
    5,
    22,
    1,
    3,
    17,
    0,
    1
  ),
  
  expected_fdr_pct = c(
    2.1,
    1.8,
    7.8,
    0.4,
    1.1,
    6.0,
    0.0,
    0.4
  )
)

# Convert factor to character for comparison

CHECK_DATA <- copy(d)

CHECK_DATA[
  ,
  index := as.character(index)
]

CHECK_DATA <- merge(
  CHECK_DATA,
  EXPECTED_FDR,
  by = "index",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------
# Check expected FDR counts
# ------------------------------------------------

if (
  any(
    CHECK_DATA$n_fdr_sig !=
    CHECK_DATA$expected_n_fdr,
    na.rm = TRUE
  )
) {
  
  print(
    CHECK_DATA[
      n_fdr_sig != expected_n_fdr,
      .(
        index,
        actual = n_fdr_sig,
        expected = expected_n_fdr
      )
    ]
  )
  
  stop(
    "QC FAILED: primary FDR counts do not match expected locked results."
  )
  
}

# ------------------------------------------------
# Check expected percentages
# ------------------------------------------------

if (
  any(
    abs(
      CHECK_DATA$fdr_pct -
      CHECK_DATA$expected_fdr_pct
    ) > 0.05,
    na.rm = TRUE
  )
) {
  
  print(
    CHECK_DATA[
      abs(
        fdr_pct -
          expected_fdr_pct
      ) > 0.05,
      .(
        index,
        actual = fdr_pct,
        expected = expected_fdr_pct
      )
    ]
  )
  
  stop(
    "QC FAILED: primary FDR percentages do not match expected locked results."
  )
  
}

cat(
  "PRIMARY FDR QC: PASSED\n"
)

# ================================================================
# 12. SAVE MASTER FIGURE DATA
# ================================================================

MASTER_DATA <- file.path(
  OUT_DIR,
  "FIGURE_2_RAW_vs_FDR_SIGNIFICANCE_1996_2025_DATA.csv"
)

fwrite(
  d,
  MASTER_DATA
)

# ================================================================
# 13. LONG FORMAT FOR PLOTTING
# ================================================================

plot_data <- melt(
  d,
  id.vars = c(
    "index",
    "n_tested"
  ),
  measure.vars = c(
    "raw_pct",
    "fdr_pct"
  ),
  variable.name = "criterion",
  value.name = "percentage"
)

# ------------------------------------------------
# Rename criteria
# ------------------------------------------------

plot_data[
  criterion == "raw_pct",
  criterion := "Raw p < 0.05"
]

plot_data[
  criterion == "fdr_pct",
  criterion := "BH-FDR q ≤ 0.10"
]

plot_data[
  ,
  criterion := factor(
    criterion,
    levels = c(
      "Raw p < 0.05",
      "BH-FDR q ≤ 0.10"
    )
  )
]

# ------------------------------------------------
# Ensure index order
# ------------------------------------------------

plot_data[
  ,
  index := factor(
    index,
    levels = INDEX_ORDER
  )
]

# ================================================================
# 14. FIGURE
# ================================================================

p <- ggplot(
  plot_data,
  aes(
    x = index,
    y = percentage,
    fill = criterion
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.72
    ),
    width = 0.62
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        percentage
      )
    ),
    position = position_dodge(
      width = 0.72
    ),
    vjust = -0.35,
    size = 3.1
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      max(
        plot_data$percentage,
        na.rm = TRUE
      ) * 1.18
    ),
    expand = c(
      0,
      0
    ),
    breaks = seq(
      0,
      100,
      by = 10
    ),
    name = "Significant grid cells (%)"
  ) +
  
  scale_fill_manual(
    values = c(
      "Raw p < 0.05" = "grey65",
      "BH-FDR q ≤ 0.10" = "grey20"
    ),
    name = NULL
  ) +
  
  labs(
    title =
      "Statistical significance of temperature-extreme trends",
    
    subtitle =
      "AMJ 1996–2025 | Primary Hamed–Rao MMK specification",
    
    x = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10.5,
      colour = "grey30"
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text.x = element_text(
      face = "bold",
      colour = "black",
      size = 10
    ),
    
    axis.text.y = element_text(
      colour = "black"
    ),
    
    panel.grid.major = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      12,
      15,
      12,
      15
    )
  )

# ================================================================
# 15. SAVE 600-DPI PNG
# ================================================================

PNG_FILE <- file.path(
  OUT_DIR,
  "FIGURE_2_RAW_vs_FDR_SIGNIFICANCE_1996_2025.png"
)

ggsave(
  filename = PNG_FILE,
  plot = p,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ================================================================
# 16. SAVE VECTOR PDF
# ================================================================

PDF_FILE <- file.path(
  OUT_DIR,
  "FIGURE_2_RAW_vs_FDR_SIGNIFICANCE_1996_2025.pdf"
)

ggsave(
  filename = PDF_FILE,
  plot = p,
  width = 10,
  height = 6.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ================================================================
# 17. DISPLAY
# ================================================================

print(p)

# ================================================================
# 18. FINAL MESSAGE
# ================================================================

cat("\n")
cat("====================================================\n")
cat("FIGURE 2 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

cat(
  "PNG : ",
  PNG_FILE,
  "\n",
  sep = ""
)

cat(
  "PDF : ",
  PDF_FILE,
  "\n",
  sep = ""
)

cat(
  "DATA: ",
  MASTER_DATA,
  "\n",
  sep = ""
)

cat("====================================================\n")
cat("PRIMARY CONFIGURATION LOCKED:\n")
cat("  Completeness = 80%\n")
cat("  FDR q        = 0.10\n")
cat("  Variant      = HR_lag3\n")
cat("  Period       = 1996–2025\n")
cat("====================================================\n")