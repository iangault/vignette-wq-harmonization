# Database assembly functions: reshape, deduplicate, label, and prepare for
# guideline comparison.

library(dplyr)
library(tidyr)

# Reshape wide in situ field columns to long format so they can be merged with
# lab results. Column names encode both the analyte abbreviation and unit
# (e.g. "sp.cond.uS.cm"); these are split on "." to recover each.
reshape_field_insitu <- function(df) {
  insitu_cols <- c("pH", "sp.cond.uS.cm", "sal.ppt", "temp.C",
                   "diss.O2.mg.L", "diss.O2.percent", "redox.ORP", "turbidity.NTU")

  lab_cols <- c("analysis", "analyte.impute", "test.code", "result.impute",
                "unit.impute", "hold.time", "dl.impute", "dl.f", "source", "lab.id")

  df_long <- df %>%
    select(-any_of(lab_cols)) %>%
    pivot_longer(
      cols      = all_of(insitu_cols),
      names_to  = "analyte.impute",
      values_to = "result.impute"
    ) %>%
    mutate(
      # Rename before split to avoid ambiguous dot splits
      analyte.impute = case_when(
        analyte.impute == "sp.cond.uS.cm"   ~ "spcond.uScm",
        analyte.impute == "diss.O2.mg.L"    ~ "disso2.mgl",
        analyte.impute == "diss.O2.percent" ~ "dissO2.percent",
        TRUE                               ~ analyte.impute
      ),
      analysis = "Field Measured"
    ) %>%
    distinct() %>%
    # pH has no unit suffix; fill = "right" leaves unit.impute NA, resolved below.
    separate(analyte.impute, into = c("analyte.impute", "unit.impute"),
             sep = "\\.", fill = "right") %>%
    mutate(
      unit.impute = case_when(
        unit.impute == "uScm"                  ~ "µS/cm",
        unit.impute == "C"                     ~ "°C",
        unit.impute == "mgl"                   ~ "mg/L",
        unit.impute == "percent"               ~ "%",
        grepl("pH", analyte.impute)            ~ "pH units",
        unit.impute == "ORP"                   ~ "mV",
        unit.impute == "NTU"                   ~ "NTU",
        TRUE                                  ~ unit.impute
      ),
      analyte.impute = case_when(
        analyte.impute == "spcond"    ~ "Specific Conductivity",
        analyte.impute %in% c("dissO2", "disso2") ~ "Dissolved Oxygen",
        analyte.impute == "redox"     ~ "Oxidation-Reduction Potential",
        analyte.impute == "turbidity" ~ "Turbidity",
        analyte.impute == "sal"       ~ "Salinity",
        analyte.impute == "temp"      ~ "Temperature",
        TRUE                         ~ analyte.impute
      )
    )

  df_long
}

# Remove duplicate analyte observations by keeping the row with the lowest
# detection limit. When the same analyte appears twice under the same sample
# (e.g. submitted to two labs, or under two test codes), the result at the
# finer detection threshold is more informative.
deduplicate_wq <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    arrange(dl.impute, .by_group = TRUE) %>%
    distinct(across(all_of(group_vars)), .keep_all = TRUE) %>%
    ungroup()
}

# Assign sequential week labels numbered from the first sampling date.
# Week boundaries run Sunday-to-Saturday to match typical field schedules.
label_weeks <- function(df, n_weeks = 5) {
  start_date <- min(df$sample.date, na.rm = TRUE)
  breaks     <- seq(start_date, by = "week", length.out = n_weeks + 1)

  df %>%
    mutate(
      week = cut(sample.date, breaks = breaks,
                 labels = paste("Week", seq_len(n_weeks)),
                 include.lowest = TRUE)
    )
}

# Join water-chemistry cofactors (pH, hardness, DOC, etc.) as wide columns onto
# each row. Required for hardness-dependent guideline calculations, where the
# criterion value is a function of in-sample hardness rather than a fixed threshold.
attach_cofactors <- function(df, cofactors = c("pH", "Temperature", "Chloride",
                                               "Dissolved Organic Carbon",
                                               "Hardness, as CaCO3", "Salinity")) {
  cofactor_wide <- df %>%
    filter(analyte.impute %in% cofactors) %>%
    select(program, site.impute, sample.id, sample.date, actual.test,
           analyte.impute, result.impute) %>%
    pivot_wider(names_from = analyte.impute, values_from = result.impute)

  df %>%
    left_join(cofactor_wide,
              by = c("program", "site.impute", "sample.id", "actual.test", "sample.date"))
}
