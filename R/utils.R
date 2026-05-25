# Shared utilities for the WQ harmonization pipeline.

# Magnitude-based rounding: finer precision at sub-ppm concentrations, coarser at
# high concentrations. Avoids over-reporting spurious digits in lab results.
make_rspec <- function() {
  table.glue::round_spec() %>%
    table.glue::round_using_magnitude(
      digits = c(6, 5, 4, 3, 2, 1, 0),
      breaks = c(0.00001, 0.0001, 0.001, 0.01, 1, 10, Inf)
    )
}

# Apply magnitude-based rounding to every numeric column in a data frame.
round_table_cols <- function(df, rspec = make_rspec()) {
  for (i in seq_along(df)) {
    if (is.numeric(df[[i]])) {
      df[[i]] <- table.glue::table_glue("{df[[i]]}", rspec = rspec)
    }
  }
  df
}

# Ordered factor levels for analysis groups; ensures consistent table/plot ordering.
analysis_levels <- c(
  "Field Measured",
  "Conventional Parameters",
  "Major Ions",
  "Nutrients and Biological Indicators",
  "Total Metals",
  "Dissolved Metals"
)
