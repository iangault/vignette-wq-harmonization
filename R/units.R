# Unit standardization for lab water quality results.
#
# Labs report in whatever unit the method dictates. For cross-analyte comparability:
#   - Mercury species: µg/L → ng/L  (low DL requires finer precision reporting)
#   - All other metals: µg/L → mg/L (aligns with nutrient and conventional units)
#   - Total Phosphorus: µg/L → mg/L (submitted as metals but reported as nutrient)
#   - Major Ions: µg/L → mg/L
#   - Conductivity: µmhos/cm → µS/cm (numerically identical; label only)
#
# Detection limits are scaled by the same factor as results so that DL comparisons
# remain valid after conversion.

library(dplyr)

standardize_units <- function(df) {
  df %>%
    mutate(
      result.impute = case_when(
        unit == "µg/L" & grepl("mercury", analyte.impute, ignore.case = TRUE) ~ result.impute * 1000,
        unit == "µg/L" & grepl("metal",   analysis,       ignore.case = TRUE) ~ result.impute / 1000,
        unit == "µg/L" & analyte.impute == "Total Phosphorus"                 ~ result.impute / 1000,
        unit == "µg/L" & analysis == "Major Ions"                             ~ result.impute / 1000,
        TRUE                                                                  ~ result.impute
      ),
      dl.impute = case_when(
        unit == "µg/L" & grepl("mercury", analyte.impute, ignore.case = TRUE) ~ dl.impute * 1000,
        unit == "µg/L" & grepl("metal",   analysis,       ignore.case = TRUE) ~ dl.impute / 1000,
        unit == "µg/L" & analyte.impute == "Total Phosphorus"                 ~ dl.impute / 1000,
        unit == "µg/L" & analysis == "Major Ions"                             ~ dl.impute / 1000,
        TRUE                                                                  ~ dl.impute
      ),
      unit.impute = case_when(
        unit == "µg/L" & grepl("mercury", analyte.impute, ignore.case = TRUE) ~ "ng/L",
        unit == "µg/L" & grepl("metal",   analysis,       ignore.case = TRUE) ~ "mg/L",
        unit == "µg/L" & analyte.impute == "Total Phosphorus"                 ~ "mg/L",
        unit == "µg/L" & analysis == "Major Ions"                             ~ "mg/L",
        unit == "µmhos/cm"                                                    ~ "µS/cm",
        TRUE                                                                  ~ unit.impute
      )
    )
}
