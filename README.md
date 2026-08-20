# WMO_IMD_temperature_extremes

## Project overview

This repository contains the R scripts used for the analysis of
pre-monsoon temperature extremes over India using daily gridded
temperature data from the India Meteorological Department (IMD).

The analysis covers the April–June (AMJ) season and uses daily maximum
temperature (Tmax) and minimum temperature (Tmin) data on the original
1° × 1° IMD grid. The analysis was conducted for 283 land grid cells.

Eight temperature-extreme indices were analysed:

- TXx – annual maximum of daily maximum temperature
- TXn – annual minimum of daily maximum temperature
- TNx – annual maximum of daily minimum temperature
- TNn – annual minimum of daily minimum temperature
- TX90p – percentage of warm days
- TN10p – percentage of cold nights
- WSDI – Warm Spell Duration Index
- DTR – Diurnal Temperature Range

## Analysis workflow

The analysis was performed in the following sequence:

01. IMD data extraction, AMJ selection and quality control
02. Calendar-day percentile threshold calculation for TX90 and TN10
03. Calculation of annual TX90p and TN10p
04. Bootstrap procedure for TX90p and TN10p
05. Validation, flagging and preparation of the locked analysis-ready dataset
06. Calculation of WSDI, TXx, TXn, TNx, TNn and DTR
07. Pre-lock diagnostics for the trend methodology
07A. Additional pre-lock diagnostics for the trend methodology
07B. Final trend analysis of TX90p and TN10p
08. Final results synthesis, tables, figures and quality control

## Repository scripts

1. `01_IMD_Extraction_to_AMJ_TX90p_TN10p_valid_observation_QC_1981_2010.R`
2. `02_AMJ_calendar_day_percentile_thresholds_TX90_TN10.R`
3. `03_AMJ_TX90p_TN10p_1976_2025.R`
4. `04_AMJ_TX90p_TN10p_bootstrap_1976_2025.R`
5. `05_Validation_flagging_and_LOCKED_analysis_ready_dataset.R`
6. `06_AMJ_WSDI_TXx_TXn_TNx_TNn_DTR_1976_2025.R`
7. `07_PRE-LOCK_DIAGNOSTIC_for_trend_methodology.R`
8. `07A_PRE-LOCK_DIAGNOSTIC_for_trend_methodology.R`
9. `07B_FINAL_TREND_ANALYSIS_for_the_percentile_indices_TX90p_TN10p.R`
10. `08_FINAL_RESULTS_SYNTHESIS_TABLES_FIGURES_AND_QC.R`
