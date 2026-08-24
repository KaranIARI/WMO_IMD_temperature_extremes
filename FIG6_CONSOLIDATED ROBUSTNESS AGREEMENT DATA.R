###############################################################################
# TABLE S7 + FIGURE 6
# Variant sensitivity of grid-level trend classification
#
# Primary method : Hamed-Rao modified Mann-Kendall, lag-3
# Period         : TREND_1996_2025
#
# Comparisons:
#   1. HR_lag3 vs HR_all
#   2. HR_lag3 vs MK_original
#
# Indices:
#   TXx, TXn, TNx, TNn, TX90p, TN10p, WSDI, DTR
#
# STEP07A:
#   TXx, TXn, TNx, TNn, WSDI, DTR
#
# STEP07B:
#   TX90p, TN10p
#
# FINAL:
#   8 indices x 2 comparisons = 16 rows
###############################################################################

rm(list = ls())

# ============================================================================
# 1. PACKAGES
# ============================================================================

library(data.table)
library(ggplot2)

# ============================================================================
# 2. FILE PATHS
# ============================================================================

STEP07A_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP07A_variant_comparison.csv"

STEP07B_FILE <- "F:/WMO_IMD_R/WMO_IMD/data/AMJ_STEP07B_sensitivity_variants.csv"

OUT_DIR <- "F:/WMO_IMD_R/WMO_IMD/STEP08"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================================
# 3. CONSTANTS
# ============================================================================

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

INDEX_STEP07A <- c(
  "TXx",
  "TXn",
  "TNx",
  "TNn",
  "WSDI",
  "DTR"
)

INDEX_STEP07B <- c(
  "TX90p",
  "TN10p"
)

COMPARISON_ORDER <- c(
  "HR_lag3_vs_HR_all",
  "HR_lag3_vs_MK_original"
)

# ============================================================================
# 4. CHECK INPUT FILES
# ============================================================================

if (!file.exists(STEP07A_FILE)) {
  stop(
    "ERROR: STEP07A file not found:\n",
    STEP07A_FILE
  )
}

if (!file.exists(STEP07B_FILE)) {
  stop(
    "ERROR: STEP07B file not found:\n",
    STEP07B_FILE
  )
}

# ============================================================================
# 5. READ STEP07A
# ============================================================================

d07a <- fread(STEP07A_FILE)

cat("\n============================================================\n")
cat("STEP07A LOADED\n")
cat("============================================================\n")

cat("Rows:", nrow(d07a), "\n")

# ============================================================================
# 6. CHECK STEP07A COLUMNS
# ============================================================================

required_07a <- c(
  "period",
  "index",
  "variant_a",
  "variant_b",
  "n_grids",
  "n_discordant"
)

missing_07a <- setdiff(
  required_07a,
  names(d07a)
)

if (length(missing_07a) > 0L) {
  stop(
    "ERROR: STEP07A missing columns: ",
    paste(missing_07a, collapse = ", ")
  )
}

# ============================================================================
# 7. FILTER STEP07A TO PRIMARY PERIOD
# ============================================================================

d07a <- d07a[
  period == "TREND_1996_2025"
]

# Keep only the six indices available in STEP07A

d07a <- d07a[
  index %in% INDEX_STEP07A
]

# ============================================================================
# 8. IDENTIFY REQUIRED COMPARISONS IN STEP07A
# ============================================================================

d07a[
  ,
  comparison := fifelse(
    
    (
      variant_a == "HR_lag3" &
        variant_b == "HR_all"
    ) |
      (
        variant_a == "HR_all" &
          variant_b == "HR_lag3"
      ),
    
    "HR_lag3_vs_HR_all",
    
    fifelse(
      
      (
        variant_a == "HR_lag3" &
          variant_b == "MK_original"
      ) |
        (
          variant_a == "MK_original" &
            variant_b == "HR_lag3"
        ),
      
      "HR_lag3_vs_MK_original",
      
      NA_character_
    )
  )
]

d07a <- d07a[
  !is.na(comparison)
]

# ============================================================================
# 9. CHECK STEP07A REQUIRED COMPARISONS
# ============================================================================

if (
  !all(
    d07a$comparison %in%
    COMPARISON_ORDER
  )
) {
  stop(
    "ERROR: Unexpected comparison found in STEP07A."
  )
}

# ============================================================================
# 10. CREATE STEP07A SUMMARY
# ============================================================================

s07a <- d07a[
  ,
  .(
    n_compared = first(
      as.integer(n_grids)
    ),
    
    n_discordant = first(
      as.integer(n_discordant)
    )
  ),
  
  by = .(
    index,
    comparison
  )
]

# ============================================================================
# 11. CHECK STEP07A DUPLICATES
# ============================================================================

dup07a <- s07a[
  ,
  .N,
  by = .(
    index,
    comparison
  )
][N > 1L]

if (nrow(dup07a) > 0L) {
  
  cat("\nDuplicate STEP07A combinations:\n")
  print(dup07a)
  
  stop(
    "ERROR: Duplicate index × comparison combinations in STEP07A."
  )
}

# ============================================================================
# 12. CALCULATE STEP07A SAME / DISCORDANT
# ============================================================================

s07a[
  ,
  `:=`(
    n_same_class =
      n_compared - n_discordant
  )
]

s07a[
  ,
  `:=`(
    
    agreement_pct =
      round(
        100 *
          n_same_class /
          n_compared,
        2
      ),
    
    discordance_pct =
      round(
        100 *
          n_discordant /
          n_compared,
        2
      )
  )
]

# ============================================================================
# 13. STEP07A INTERNAL QC
# ============================================================================

if (
  any(
    s07a$n_compared <= 0L
  )
) {
  stop(
    "ERROR: Invalid n_compared in STEP07A."
  )
}

if (
  any(
    s07a$n_discordant < 0L |
    s07a$n_discordant >
    s07a$n_compared
  )
) {
  stop(
    "ERROR: Invalid n_discordant in STEP07A."
  )
}

if (
  any(
    s07a$n_same_class +
    s07a$n_discordant !=
    s07a$n_compared
  )
) {
  stop(
    "ERROR: STEP07A count identity failed."
  )
}

if (
  nrow(s07a) != 12L
) {
  stop(
    "ERROR: STEP07A should contain 12 rows ",
    "(6 indices × 2 comparisons), found ",
    nrow(s07a)
  )
}

# ============================================================================
# 14. PRINT STEP07A
# ============================================================================

cat("\n============================================================\n")
cat("STEP07A SUMMARY\n")
cat("============================================================\n")

print(
  s07a[
    ,
    .(
      index,
      comparison,
      n_compared,
      n_same_class,
      n_discordant,
      agreement_pct,
      discordance_pct
    )
  ]
)

# ============================================================================
# 15. READ STEP07B
# ============================================================================

d07b <- fread(STEP07B_FILE)

cat("\n============================================================\n")
cat("STEP07B LOADED\n")
cat("============================================================\n")

cat("Rows:", nrow(d07b), "\n")

cat("Columns:\n")
print(names(d07b))

# ============================================================================
# 16. CHECK STEP07B ACTUAL COLUMNS
# ============================================================================

required_07b <- c(
  "index",
  "comparison",
  "n_grids",
  "n_same_class"
)

missing_07b <- setdiff(
  required_07b,
  names(d07b)
)

if (length(missing_07b) > 0L) {
  stop(
    "ERROR: STEP07B missing columns: ",
    paste(missing_07b, collapse = ", ")
  )
}

# IMPORTANT:
# STEP07B DOES NOT HAVE A 'period' COLUMN.
# Therefore no period filter is applied here.

# ============================================================================
# 17. KEEP ONLY TX90p AND TN10p
# ============================================================================

d07b <- d07b[
  index %in% INDEX_STEP07B
]

# ============================================================================
# 18. STANDARDIZE STEP07B COMPARISON NAMES
# ============================================================================

d07b[
  ,
  comparison := fifelse(
    
    comparison ==
      "HR_lag3 vs HR_all",
    
    "HR_lag3_vs_HR_all",
    
    fifelse(
      
      comparison ==
        "HR_lag3 vs MK_original",
      
      "HR_lag3_vs_MK_original",
      
      NA_character_
    )
  )
]

d07b <- d07b[
  !is.na(comparison)
]

# ============================================================================
# 19. CHECK STEP07B COMPARISONS
# ============================================================================

if (
  !all(
    d07b$comparison %in%
    COMPARISON_ORDER
  )
) {
  stop(
    "ERROR: Unexpected comparison found in STEP07B."
  )
}

# ============================================================================
# 20. CHECK STEP07B DUPLICATES
# ============================================================================

dup07b <- d07b[
  ,
  .N,
  by = .(
    index,
    comparison
  )
][N > 1L]

if (nrow(dup07b) > 0L) {
  
  cat("\nDuplicate STEP07B combinations:\n")
  print(dup07b)
  
  stop(
    "ERROR: Duplicate index × comparison combinations in STEP07B."
  )
}

# ============================================================================
# 21. CREATE STEP07B SUMMARY
# ============================================================================

s07b <- d07b[
  ,
  .(
    n_compared =
      first(
        as.integer(n_grids)
      ),
    
    n_same_class =
      first(
        as.integer(n_same_class)
      )
  ),
  
  by = .(
    index,
    comparison
  )
]

# ============================================================================
# 22. CALCULATE STEP07B DISCORDANCE
# ============================================================================

s07b[
  ,
  `:=`(
    n_discordant =
      n_compared - n_same_class
  )
]

s07b[
  ,
  `:=`(
    
    agreement_pct =
      round(
        100 *
          n_same_class /
          n_compared,
        2
      ),
    
    discordance_pct =
      round(
        100 *
          n_discordant /
          n_compared,
        2
      )
  )
]

# ============================================================================
# 23. STEP07B QC
# ============================================================================

if (
  nrow(s07b) != 4L
) {
  stop(
    "ERROR: STEP07B should contain 4 rows ",
    "(2 indices × 2 comparisons), found ",
    nrow(s07b)
  )
}

if (
  any(
    s07b$n_same_class +
    s07b$n_discordant !=
    s07b$n_compared
  )
) {
  stop(
    "ERROR: STEP07B count identity failed."
  )
}

# ============================================================================
# 24. PRINT STEP07B
# ============================================================================

cat("\n============================================================\n")
cat("STEP07B SUMMARY\n")
cat("============================================================\n")

print(
  s07b[
    ,
    .(
      index,
      comparison,
      n_compared,
      n_same_class,
      n_discordant,
      agreement_pct,
      discordance_pct
    )
  ]
)

# ============================================================================
# 25. COMBINE STEP07A + STEP07B
# ============================================================================

s7 <- rbindlist(
  list(
    s07a,
    s07b
  ),
  use.names = TRUE,
  fill = TRUE
)

# ============================================================================
# 26. KEEP ONLY THE FINAL 8 INDICES
# ============================================================================

s7 <- s7[
  index %in% INDEX_ORDER
]

s7 <- s7[
  comparison %in%
    COMPARISON_ORDER
]

# ============================================================================
# 27. CHECK FINAL DUPLICATES
# ============================================================================

dup_final <- s7[
  ,
  .N,
  by = .(
    index,
    comparison
  )
][N > 1L]

if (nrow(dup_final) > 0L) {
  
  cat("\nFinal duplicate combinations:\n")
  print(dup_final)
  
  stop(
    "ERROR: Final table contains duplicate index × comparison combinations."
  )
}

# ============================================================================
# 28. CHECK EXACT EXPECTED COMBINATIONS
# ============================================================================

expected <- CJ(
  index = INDEX_ORDER,
  comparison = COMPARISON_ORDER
)

missing_final <- expected[
  !s7,
  on = .(
    index,
    comparison
  )
]

if (nrow(missing_final) > 0L) {
  
  cat("\nMissing final combinations:\n")
  print(missing_final)
  
  stop(
    "ERROR: Required final combinations are missing."
  )
}

# ============================================================================
# 29. FINAL 16-ROW CHECK
# ============================================================================

if (
  nrow(s7) != 16L
) {
  stop(
    "ERROR: Expected exactly 16 rows ",
    "(8 indices × 2 comparisons), found ",
    nrow(s7)
  )
}

# ============================================================================
# 30. FINAL COUNT QC
# ============================================================================

if (
  any(
    s7$n_same_class +
    s7$n_discordant !=
    s7$n_compared
  )
) {
  stop(
    "ERROR: Final count identity failed."
  )
}

if (
  any(
    s7$n_compared <= 0L
  )
) {
  stop(
    "ERROR: Invalid n_compared value."
  )
}

if (
  any(
    s7$n_same_class < 0L |
    s7$n_discordant < 0L
  )
) {
  stop(
    "ERROR: Negative classification count detected."
  )
}

# ============================================================================
# 31. RECALCULATE FINAL PERCENTAGES
# ============================================================================

s7[
  ,
  `:=`(
    
    agreement_pct =
      round(
        100 *
          n_same_class /
          n_compared,
        2
      ),
    
    discordance_pct =
      round(
        100 *
          n_discordant /
          n_compared,
        2
      )
  )
]

# ============================================================================
# 32. PERCENTAGE QC
# ============================================================================

if (
  any(
    abs(
      s7$agreement_pct +
      s7$discordance_pct -
      100
    ) > 0.01
  )
) {
  stop(
    "ERROR: Agreement + discordance does not equal 100%."
  )
}

# ============================================================================
# 33. FINAL ORDER
# ============================================================================

s7[
  ,
  index := factor(
    index,
    levels = INDEX_ORDER
  )
]

s7[
  ,
  comparison := factor(
    comparison,
    levels = COMPARISON_ORDER
  )
]

setorder(
  s7,
  comparison,
  index
)

# ============================================================================
# 34. FINAL TABLE S7
# ============================================================================

TABLE_S7 <- s7[
  ,
  .(
    index,
    comparison,
    n_compared,
    n_same_class,
    n_discordant,
    agreement_pct,
    discordance_pct
  )
]

cat("\n\n")
cat("################################################################\n")
cat("# FINAL TABLE S7\n")
cat("################################################################\n")

print(TABLE_S7)

# ============================================================================
# 35. SAVE TABLE S7
# ============================================================================

S7_FILE <- file.path(
  OUT_DIR,
  "TABLE_S7_VARIANT_SENSITIVITY_AMJ_1996_2025.csv"
)

fwrite(
  TABLE_S7,
  S7_FILE
)

# ============================================================================
# 36. PREPARE FIGURE 6 DATA
# ============================================================================

plot_data <- copy(
  TABLE_S7
)

plot_data[
  ,
  index := factor(
    as.character(index),
    levels = INDEX_ORDER
  )
]

plot_data[
  ,
  comparison_label := fifelse(
    
    as.character(comparison) ==
      "HR_lag3_vs_HR_all",
    
    "HR lag-3 vs HR all-lag",
    
    "HR lag-3 vs original MK"
  )
]

plot_data[
  ,
  comparison_label := factor(
    comparison_label,
    levels = c(
      "HR lag-3 vs HR all-lag",
      "HR lag-3 vs original MK"
    )
  )
]

# ============================================================================
# 37. FIGURE 6
# ============================================================================

p6 <- ggplot(
  plot_data,
  aes(
    x = index,
    y = agreement_pct
  )
) +
  
  geom_col(
    width = 0.65,
    fill = "grey35",
    colour = "black"
  ) +
  
  geom_text(
    aes(
      label =
        sprintf(
          "%.1f",
          agreement_pct
        )
    ),
    vjust = -0.35,
    size = 3.1
  ) +
  
  facet_wrap(
    ~ comparison_label,
    nrow = 1
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      105
    ),
    
    breaks = seq(
      0,
      100,
      by = 10
    ),
    
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    
    title =
      "Robustness of grid-level trend classification",
    
    subtitle =
      "AMJ 1996–2025",
    
    x = NULL,
    
    y =
      "Classification agreement (%)"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        hjust = 0.5,
        size = 10,
        colour = "grey30"
      ),
    
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11
      ),
    
    axis.text.x =
      element_text(
        colour = "black",
        size = 9.5
      ),
    
    axis.text.y =
      element_text(
        colour = "black"
      ),
    
    axis.title.y =
      element_text(
        face = "bold",
        size = 11
      ),
    
    panel.grid.major =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "none"
  )

# ============================================================================
# 38. SAVE FIGURE 6 — PNG
# ============================================================================

FIG6_PNG <- file.path(
  OUT_DIR,
  "FIGURE_6_VARIANT_AGREEMENT_AMJ_1996_2025_600dpi.png"
)

ggsave(
  filename = FIG6_PNG,
  plot = p6,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================================
# 39. SAVE FIGURE 6 — PDF
# ============================================================================

FIG6_PDF <- file.path(
  OUT_DIR,
  "FIGURE_6_VARIANT_AGREEMENT_AMJ_1996_2025.pdf"
)

ggsave(
  filename = FIG6_PDF,
  plot = p6,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ============================================================================
# 40. SAVE FIGURE 6 DATA
# ============================================================================

FIG6_DATA <- file.path(
  OUT_DIR,
  "FIGURE_6_VARIANT_AGREEMENT_DATA.csv"
)

fwrite(
  plot_data[
    ,
    .(
      index,
      comparison,
      n_compared,
      n_same_class,
      n_discordant,
      agreement_pct,
      discordance_pct
    )
  ],
  FIG6_DATA
)

# ============================================================================
# 41. FINAL QC REPORT
# ============================================================================

QC_REPORT <- data.table(
  
  check = c(
    "Final row count = 16",
    "Eight indices present",
    "Two comparisons present",
    "No duplicate combinations",
    "No missing combinations",
    "Counts internally consistent",
    "Percentages internally consistent"
  ),
  
  result = c(
    
    nrow(TABLE_S7) == 16L,
    
    uniqueN(
      TABLE_S7$index
    ) == 8L,
    
    uniqueN(
      TABLE_S7$comparison
    ) == 2L,
    
    nrow(
      dup_final
    ) == 0L,
    
    nrow(
      missing_final
    ) == 0L,
    
    all(
      TABLE_S7$n_same_class +
        TABLE_S7$n_discordant ==
        TABLE_S7$n_compared
    ),
    
    all(
      abs(
        TABLE_S7$agreement_pct +
          TABLE_S7$discordance_pct -
          100
      ) <= 0.01
    )
  )
)

QC_FILE <- file.path(
  OUT_DIR,
  "TABLE_S7_FIGURE_6_QC.csv"
)

fwrite(
  QC_REPORT,
  QC_FILE
)

# ============================================================================
# 42. DISPLAY FIGURE
# ============================================================================

print(p6)

# ============================================================================
# 43. FINAL OUTPUT
# ============================================================================

cat("\n\n")
cat("################################################################\n")
cat("# TABLE S7 + FIGURE 6 COMPLETED\n")
cat("################################################################\n")

cat(
  "Final rows       : ",
  nrow(TABLE_S7),
  " / 16\n",
  sep = ""
)

cat(
  "Indices           : ",
  uniqueN(TABLE_S7$index),
  " / 8\n",
  sep = ""
)

cat(
  "Comparisons       : ",
  uniqueN(TABLE_S7$comparison),
  " / 2\n",
  sep = ""
)

cat(
  "Table S7 CSV      : ",
  S7_FILE,
  "\n",
  sep = ""
)

cat(
  "Figure 6 PNG      : ",
  FIG6_PNG,
  "\n",
  sep = ""
)

cat(
  "Figure 6 PDF      : ",
  FIG6_PDF,
  "\n",
  sep = ""
)

cat(
  "Figure 6 data     : ",
  FIG6_DATA,
  "\n",
  sep = ""
)

cat(
  "QC report         : ",
  QC_FILE,
  "\n",
  sep = ""
)

cat("################################################################\n")