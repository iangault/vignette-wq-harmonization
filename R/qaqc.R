# QAQC functions for multi-lab water quality datasets.
#
# Each function returns a data frame (either an annotated version of the input
# or a summary table) so results can be displayed directly in the vignette.

library(dplyr)
library(tidyr)

# Flag results reported below the detection limit and substitute half-DL.
# Labs prefix below-DL results with "<"; this parses that notation into a logical
# flag (dl.f) and replaces the value with 0.5 × DL for downstream calculations.
# The flag must be preserved so censored values can be excluded from parametric
# statistics where appropriate.
flag_below_dl <- function(df) {
  df %>%
    mutate(
      dl        = as.numeric(dl),
      dl.impute = dl,
      unit.impute = unit,
      dl.f      = grepl("<", result),
      # Strip "<" before coercion so below-DL strings convert cleanly.
      result.impute = if_else(dl.f, dl * 0.5, as.numeric(gsub("<", "", result)))
    )
}

# Parse lab holding-time symbols to a consistent "T" (exceeded) / "F" (ok) flag.
# Historical lab exports use Unicode checkmark/cross glyphs (ü/û), the letters
# T/F, or logical TRUE/FALSE (introduced when readr guesses column type).
# Converting to character first ensures case_when always operates on one type.
parse_holding_time <- function(df) {
  df %>%
    mutate(
      hold.time = as.character(hold.time),
      hold.time = case_when(
        hold.time %in% c("\xfb", "T", "exceedance", "TRUE")            ~ "T",
        hold.time %in% c("\xfc", "F", "FALSE") | is.na(hold.time)      ~ "F",
        TRUE ~ hold.time
      )
    )
}

# Screen blank samples for detections exceeding 5× the detection limit.
# A result > 5× DL in a field or travel blank indicates probable contamination
# during sampling or transport; these analytes need to be flagged in the report.
screen_blanks <- function(df) {
  df %>%
    filter(grepl("Blank", sample.type), !is.na(dl.impute)) %>%
    mutate(exceeds_5x_dl = (result.impute / dl.impute) > 5) %>%
    select(program, sample.type, site.impute, sample.date, analysis,
           analyte.impute, result.impute, dl.impute, unit.impute, exceeds_5x_dl)
}

# Pivot blank results to wide format for tabular reporting.
# One column per blank event (site × date × sample ID), rows are analytes.
build_blank_table <- function(blank_df, blank_type = "Field Blank") {
  blank_df %>%
    filter(sample.type == blank_type) %>%
    pivot_wider(
      id_cols    = c("analysis", "analyte.impute", "dl.impute", "unit.impute"),
      names_from = c("site.impute", "sample.date"),
      values_from = "result.impute"
    ) %>%
    arrange(analysis, analyte.impute)
}

# Compute RPD for field duplicate pairs.
# Each "blind" duplicate row carries the real site name in the `duplicate` column;
# this maps it back to the original sample for pairing.
# RPD = 2|A - B| / (A + B) × 100
# Pairs where both results are below DL → "<DL"
# Pairs where either result is below 5× DL → "NC" (not calculable; precision
# cannot be assessed when concentrations approach the detection threshold)
compute_rpd <- function(df) {
  df_dup <- df %>%
    mutate(
      dup.hash = if_else(
        !is.na(duplicate),
        paste(duplicate,   sample.date, sep = "|"),  # blind dup → real site hash
        paste(site.impute, sample.date, sep = "|")   # original → its own hash
      )
    )

  flagged_hashes <- df_dup %>%
    filter(!is.na(duplicate)) %>%
    pull(dup.hash) %>%
    unique()

  df_dup %>%
    filter(dup.hash %in% flagged_hashes) %>%
    group_by(dup.hash, analysis, analyte.impute, unit.impute, dl.impute) %>%
    filter(n() == 2) %>%
    summarise(
      result_a = min(result.impute),
      result_b = max(result.impute),
      rpd = case_when(
        all(dl.f)                              ~ "<DL",
        any(result.impute < 5 * dl.impute)    ~ "NC",
        TRUE ~ paste0(
          round(2 * abs(result_b - result_a) / (result_a + result_b) * 100, 1), "%"
        )
      ),
      .groups = "drop"
    ) %>%
    separate(dup.hash, into = c("site.impute", "sample.date"), sep = "\\|") %>%
    select(site.impute, sample.date, analysis, analyte.impute,
           unit.impute, dl.impute, result_a, result_b, rpd)
}

# Summarize holding-time exceedances by program and analyte.
holding_time_summary <- function(df) {
  df %>%
    group_by(program, analysis, analyte.impute) %>%
    mutate(total_obs = n()) %>%
    ungroup() %>%
    filter(hold.time == "T") %>%
    group_by(program, analysis, analyte.impute, total_obs) %>%
    summarise(n_exceed = n(), .groups = "drop") %>%
    mutate(pct_exceed = round(n_exceed / total_obs * 100, 1)) %>%
    arrange(program, desc(pct_exceed)) %>%
    select(Program = program, Analysis = analysis, Analyte = analyte.impute,
           `N Exceeded` = n_exceed, `Total Obs` = total_obs, `% Exceeded` = pct_exceed)
}

# Flag dissolved:total metals ratios > 1.5.
# A ratio substantially above 1.0 is physically implausible because dissolved
# concentration cannot greatly exceed total. Values > 1.5 most often reflect
# a labelling error, a sample swap, or a matrix effect in one of the digests.
check_dt_ratio <- function(df) {
  df %>%
    filter(analysis %in% c("Total Metals", "Dissolved Metals")) %>%
    select(program, sample.date, site.impute, sample.id, water.type,
           analysis, analyte.impute, result.impute, unit.impute, qualifier, source) %>%
    pivot_wider(names_from = "analysis", values_from = "result.impute") %>%
    mutate(
      dtratio    = `Dissolved Metals` / `Total Metals`,
      ratio_flag = dtratio > 1.5
    ) %>%
    filter(ratio_flag) %>%
    arrange(desc(dtratio)) %>%
    select(Program = program, Site = site.impute, Date = sample.date,
           `Sample ID` = sample.id, `Water Type` = water.type,
           Analyte = analyte.impute, Unit = unit.impute,
           `Dissolved` = `Dissolved Metals`, `Total` = `Total Metals`, Ratio = dtratio)
}
