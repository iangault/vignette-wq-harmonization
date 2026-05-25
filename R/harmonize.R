# Analyte name harmonization and analysis group classification.
#
# Different lab submission codes produce different names for the same compound
# (e.g. E100 reports "Ammonia, total (as N)" while E235 reports
# "Nitrogen - Ammonia as N"). The functions here resolve those variants to a
# single canonical name before any joining or comparison.

library(dplyr)

# Canonical name lookup. Keys are raw lab names; values are the harmonized names
# used throughout the pipeline. Names not in this map pass through unchanged.
analyte_name_map <- c(
  "Ammonia, total (as N)"                     = "Total Ammonia",
  "Nitrogen - Ammonia as N"                   = "Total Ammonia",
  "Ammonia, un-ionized (as N), field"         = "Un-ionized Ammonia Nitrogen",
  "Carbon, dissolved organic [DOC]"           = "Dissolved Organic Carbon",
  "Carbon, total organic [TOC]"               = "Total Organic Carbon",
  "Hardness (as CaCO3), from total Ca/Mg"     = "Hardness, as CaCO3",
  "Hardness as CaCO3"                         = "Hardness, as CaCO3",
  "Kjeldahl nitrogen, total [TKN]"            = "Total Kjeldahl Nitrogen",
  "Nitrogen - Total Kjeldahl"                 = "Total Kjeldahl Nitrogen",
  "Nitrate (as N)"                            = "Nitrate",
  "Nitrogen - Nitrate as N"                   = "Nitrate",
  "Nitrite (as N)"                            = "Nitrite",
  "Nitrogen - Nitrite as N"                   = "Nitrite",
  "Nitrogen, total organic"                   = "Total Nitrogen",
  "Solids, total dissolved [TDS]"             = "Total Dissolved Solids",
  "Solids, total suspended [TSS]"             = "Total Suspended Solids",
  "Methylmercury (as MeHg), total"            = "Methylmercury",
  "Hardness (as CaCO3), dissolved"            = "Dissolved Hardness, as CaCO3",
  "Silicate (as SiO2)"                        = "Reactive Silica, as SiO2",
  "Volatile Suspended Solids"                 = "Volatile TSS",
  "Biochemical oxygen demand [BOD]"           = "Biochemical Oxygen Demand",
  "pH"                                        = "pH, Laboratory",
  "Conductivity"                              = "Conductivity, Laboratory",
  "Turbidity"                                 = "Turbidity, Laboratory",
  "Sulfur"                                    = "Sulphur"
)

# Replace raw analyte names with canonical names. Unknown names pass through.
harmonize_analyte_names <- function(df, name_map = analyte_name_map) {
  df %>%
    mutate(analyte.impute = dplyr::coalesce(name_map[analyte], analyte))
}

# Parse ", total" / ", dissolved" suffixes into a separate analysis column, then
# strip the suffix from the analyte name. This is necessary because dissolved and
# total fractions of the same element (e.g. Copper) require different handling
# downstream, and the suffix encodes the fraction, not the analyte identity.
strip_fraction_suffix <- function(df) {
  df %>%
    mutate(
      analysis = case_when(
        grepl(", total$| total$", analyte.impute, ignore.case = TRUE)    ~ "Total Metals",
        grepl(", dissolved$| dissolved$", analyte.impute, ignore.case = TRUE) ~ "Dissolved Metals",
        TRUE ~ NA_character_
      ),
      analyte.impute = gsub(
        ", total$|, dissolved$| total$| dissolved$", "",
        analyte.impute, ignore.case = TRUE
      )
    )
}

# Assign each analyte to a reporting group. Groups drive table structure and
# determine which guidelines apply.
#
# phosphorus_nutrient_codes: test codes that identify the low-DL nutrient-panel
# phosphorus submission, which maps to "Total Phosphorus" in the Nutrients group
# rather than to the higher-DL metals-panel "Phosphorus".
classify_analysis_group <- function(df, phosphorus_nutrient_codes = NULL) {

  conventional <- c(
    "Alkalinity, Total", "Alkalinity, Bicarbonate as CaCO3",
    "Alkalinity, Carbonate as CaCO3", "pH, Laboratory", "Sulfate",
    "Total Organic Carbon", "Conductivity, Laboratory", "Dissolved Oxygen",
    "Oxidation-Reduction Potential", "Temperature", "Turbidity, Laboratory",
    "Dissolved Organic Carbon", "Total Dissolved Solids", "Total Suspended Solids",
    "Volatile TSS", "Hardness, as CaCO3", "Dissolved Hardness, as CaCO3"
  )

  nutrients <- c(
    "Nitrate", "Nitrite", "Reactive Silica, as SiO2", "Total Ammonia",
    "Total Phosphorous", "Total Phosphorus", "Chloride",
    "Total Kjeldahl Nitrogen", "Total Nitrogen", "Un-ionized Ammonia Nitrogen",
    "Chlorophyll a", "Biochemical Oxygen Demand"
  )

  df <- df %>%
    mutate(
      analysis = case_when(
        analyte.impute %in% conventional                                          ~ "Conventional Parameters",
        analyte.impute == "Methylmercury"                                         ~ "Total Metals",
        analyte.impute %in% nutrients                                             ~ "Nutrients and Biological Indicators",
        grepl("Reactive silica", analyte.impute, ignore.case = TRUE)              ~ "Major Ions",
        grepl("Chloride", analyte.impute)                                         ~ "Major Ions",
        analysis == "Dissolved Metals" &
          analyte.impute %in% c("Magnesium", "Sodium", "Calcium", "Potassium")   ~ "Major Ions",
        TRUE ~ analysis
      )
    )

  if (!is.null(phosphorus_nutrient_codes)) {
    df <- df %>%
      mutate(
        analysis      = if_else(test.code %in% phosphorus_nutrient_codes,
                                "Nutrients and Biological Indicators", analysis),
        analyte.impute = if_else(test.code %in% phosphorus_nutrient_codes,
                                 "Total Phosphorus", analyte.impute)
      )
  }

  df
}

# Classify each sample as Freshwater, Seawater, or Effluent.
#
# water.type — derived from field specific conductance (sp.cond.uS.cm) joined
#   from wq_insitu.csv before this function is called. The 4000 µS/cm threshold
#   separates fresh from saline water in this estuary. Blank sites lack a
#   conductance reading, so they are protected by the -FB/-TB check above.
#
# actual.test — reflects what the lab actually analysed, derived from the
#   submission code (E469S = seawater metals package, E420 = freshwater metals).
#   A site could have seawater conductance but have been submitted incorrectly
#   under a freshwater code; actual.test captures that distinction.
assign_water_type <- function(df) {
  df %>%
    mutate(
      water.type = case_when(
        sub.matrix == "Effluent"         ~ "Effluent",
        grepl("-FB|-TB", site.impute)    ~ "Freshwater",   # blanks lack sp.cond
        sp.cond.uS.cm < 4000            ~ "Freshwater",
        sp.cond.uS.cm >= 4000           ~ "Seawater",
        TRUE                            ~ NA_character_
      ),
      sample.type = case_when(
        grepl("-FB", site.impute)        ~ "Field Blank",
        grepl("-TB", site.impute)        ~ "Travel Blank",
        !is.na(duplicate)               ~ "Duplicate",
        TRUE                            ~ NA_character_
      )
    ) %>%
    group_by(site.impute, sample.date, sample.id) %>%
    mutate(
      actual.test = case_when(
        any(test.code == "E469S")        ~ "Seawater analysis",
        any(test.code == "E420")         ~ "Freshwater analysis",
        sub.matrix == "Effluent"         ~ "Effluent",
        TRUE                            ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    mutate(
      actual.test = case_when(
        is.na(actual.test) & water.type == "Freshwater" ~ "Freshwater analysis",
        is.na(actual.test) & water.type == "Seawater"   ~ "Seawater analysis",
        TRUE                                            ~ actual.test
      )
    )
}
