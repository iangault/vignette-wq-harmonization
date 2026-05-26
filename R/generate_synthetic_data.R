# Generates data/wq_raw.csv and data/wq_insitu.csv for the harmonization vignette.
#
# wq_raw.csv  — lab bottle results only; no field-instrument columns.
# wq_insitu.csv — one row per site × date with in-situ instrument readings;
#                 left-joined to the harmonized lab data downstream.
#
# Design choices that surface every pipeline challenge:
#
#   Naming inconsistency: River1 submits N-species under E100 ("Ammonia, total (as N)");
#     Exposure submits under E235 ("Nitrogen - Ammonia as N"). Same compound, two names.
#
#   Fraction suffix: metals arrive as "Copper, total" / "Copper, dissolved";
#     the suffix encodes the fraction and must be split from the analyte name.
#
#   Mixed units: most metals µg/L; Ca/Mg/Na/K in mg/L; mercury in µg/L but
#     reported as ng/L given ultra-low detection thresholds.
#
#   Lab split: mercury (all fractions) → Research Lab (CVAFS method);
#     all other metals and nutrients → Commercial Lab.
#
#   Source split: Effluent sub.matrix → "Client"; everything else → "Contractor".
#
#   Below-DL notation: results near DL are prefixed with "<" by the lab.
#
#   Sample types: field blanks (-FB), travel blanks (-TB), one field duplicate (-DUP).
#
#   Water type: EXP-02 is seawater (high conductance, uses E469S metals package).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(42)

dates    <- as.Date(c("2023-02-07", "2023-02-14", "2023-02-21"))
dup_date <- dates[1]

# ---- Metals analyte table ---------------------------------------------------
# One row per element; expanded into total and dissolved pairs below.
# lo / hi are typical freshwater ambient ranges used to simulate results.
# Effluent block uses elevated multipliers applied later.

metals_base <- tibble::tribble(
  ~element,      ~unit,    ~dl_t,    ~dl_d,    ~lo_t,    ~hi_t,    ~lo_d,    ~hi_d,
  "Aluminum",    "µg/L",   5.000,    2.000,    20,       800,       5,       150,
  "Antimony",    "µg/L",   0.200,    0.100,     0.2,      4.0,      0.1,      2.0,
  "Arsenic",     "µg/L",   0.200,    0.100,     0.5,      8.0,      0.3,      5.0,
  "Barium",      "µg/L",   1.000,    0.500,    20,       200,      15,       150,
  "Beryllium",   "µg/L",   0.050,    0.020,     0.05,     0.5,      0.02,     0.30,
  "Bismuth",     "µg/L",   0.100,    0.050,     0.1,      2.0,      0.05,     1.0,
  "Boron",       "µg/L",  10.000,    5.000,    50,       500,      40,       450,
  "Cadmium",     "µg/L",   0.050,    0.020,     0.05,     1.0,      0.02,     0.60,
  "Calcium",     "mg/L",   0.500,    0.200,    10,        80,       9,        75,
  "Chromium",    "µg/L",   0.500,    0.200,     1.0,     15.0,      0.5,     10.0,
  "Cobalt",      "µg/L",   0.100,    0.050,     0.2,      5.0,      0.1,      3.0,
  "Copper",      "µg/L",   0.500,    0.200,     1.0,     20.0,      0.5,     12.0,
  "Iron",        "µg/L",   5.000,    2.000,    50,      3000,      10,       200,
  "Lead",        "µg/L",   0.200,    0.100,     0.2,      5.0,      0.1,      2.0,
  "Magnesium",   "mg/L",   0.200,    0.100,     3.0,     40.0,      2.5,     38.0,
  "Manganese",   "µg/L",   1.000,    0.500,     5,       300,       2,        80,
  "Molybdenum",  "µg/L",   0.200,    0.100,     0.5,     10.0,      0.4,      8.0,
  "Nickel",      "µg/L",   0.500,    0.200,     1.0,     20.0,      0.5,     12.0,
  "Potassium",   "mg/L",   0.500,    0.200,     2.0,     20.0,      1.8,     18.0,
  "Selenium",    "µg/L",   0.200,    0.100,     0.2,      5.0,      0.1,      3.0,
  "Silver",      "µg/L",   0.050,    0.020,     0.05,     1.0,      0.02,     0.40,
  "Sodium",      "mg/L",   0.500,    0.200,     5.0,     80.0,      4.5,     75.0,
  "Strontium",   "µg/L",   1.000,    0.500,    50,       500,      45,       480,
  "Thallium",    "µg/L",   0.050,    0.020,     0.05,     0.5,      0.02,     0.30,
  "Tin",         "µg/L",   0.500,    0.200,     0.5,      5.0,      0.2,      3.0,
  "Titanium",    "µg/L",   1.000,    0.500,     2.0,     50.0,      0.5,     10.0,
  "Uranium",     "µg/L",   0.050,    0.020,     0.1,      5.0,      0.08,     4.0,
  "Vanadium",    "µg/L",   0.500,    0.200,     1.0,     15.0,      0.5,     10.0,
  "Zinc",        "µg/L",   1.000,    0.500,     5.0,     80.0,      3.0,     50.0
)

# Expand to total + dissolved analyte rows (ALS lab, ICPMS method)
metals_total <- metals_base %>%
  transmute(
    analyte    = paste0(element, ", total"),
    unit       = unit, dl = dl_t,
    lo = lo_t, hi = hi_t,
    test.code  = "ICPMS", lab.id = "Commercial Lab"
  )

metals_dissolved <- metals_base %>%
  transmute(
    analyte    = paste0(element, ", dissolved"),
    unit       = unit, dl = dl_d,
    lo = lo_d, hi = hi_d,
    test.code  = "ICPMS", lab.id = "Commercial Lab"
  )

# Mercury: separate Research Lab (ultra-low DL via CVAFS)
mercury_analytes <- tibble::tribble(
  ~analyte,                          ~unit,   ~dl,         ~lo,          ~hi,          ~test.code,  ~lab.id,
  "Mercury, total",                  "µg/L",  0.000050,    0.00010,      0.00500,      "CVAFS",     "Research Lab",
  "Mercury, dissolved",              "µg/L",  0.000020,    0.00005,      0.00300,      "CVAFS",     "Research Lab",
  "Methylmercury (as MeHg), total",  "µg/L",  0.000010,    0.00002,      0.00100,      "METHYL",    "Research Lab"
)

# Nutrients and conventional parameters (ALS lab)
river1_nutrients <- tibble::tribble(
  ~analyte,                                ~test.code, ~lab.id, ~unit,   ~dl,    ~lo,    ~hi,
  "Ammonia, total (as N)",                 "E100",     "Commercial Lab",   "mg/L",  0.010,  0.05,   0.80,
  "Nitrate (as N)",                        "E100",     "Commercial Lab",   "mg/L",  0.050,  0.20,   2.50,
  "Solids, total suspended [TSS]",         "E100",     "Commercial Lab",   "mg/L",  2.000,  4.0,   60.0,
  "Carbon, dissolved organic [DOC]",       "E100",     "Commercial Lab",   "mg/L",  0.100,  2.0,   12.0,
  "Phosphorus, total",                     "E372-U",   "Commercial Lab",   "mg/L",  0.003,  0.01,   0.15,
  "Biochemical oxygen demand [BOD]",       "E100",     "Commercial Lab",   "mg/L",  1.000,  2.0,   15.0
)

exposure_nutrients <- tibble::tribble(
  ~analyte,                                ~test.code, ~lab.id, ~unit,   ~dl,    ~lo,    ~hi,
  "Nitrogen - Ammonia as N",               "E235",     "Commercial Lab",   "mg/L",  0.005,  0.05,   0.80,
  "Nitrogen - Nitrate as N",               "E235",     "Commercial Lab",   "mg/L",  0.010,  0.20,   2.50,
  "Solids, total suspended [TSS]",         "E100",     "Commercial Lab",   "mg/L",  2.000,  4.0,   60.0,
  "Carbon, dissolved organic [DOC]",       "E100",     "Commercial Lab",   "mg/L",  0.100,  2.0,   12.0,
  "Phosphorus, total",                     "E372-U",   "Commercial Lab",   "mg/L",  0.003,  0.01,   0.15,
  "Biochemical oxygen demand [BOD]",       "E100",     "Commercial Lab",   "mg/L",  1.000,  2.0,   15.0
)

effluent_nutrients <- tibble::tribble(
  ~analyte,                                ~test.code, ~lab.id, ~unit,   ~dl,    ~lo,    ~hi,
  "Ammonia, total (as N)",                 "E100",     "Commercial Lab",   "mg/L",  0.010,  0.10,   5.00,
  "Solids, total suspended [TSS]",         "E100",     "Commercial Lab",   "mg/L",  2.000, 10.0,  120.0,
  "Phosphorus, total",                     "E372-U",   "Commercial Lab",   "mg/L",  0.010,  0.10,   1.50,
  "Biochemical oxygen demand [BOD]",       "E100",     "Commercial Lab",   "mg/L",  1.000,  5.0,   50.0
)

# Effluent metals: subset of elements, elevated concentrations
effluent_metals_base <- metals_base %>%
  filter(element %in% c("Aluminum", "Arsenic", "Cadmium", "Chromium", "Copper",
                         "Iron", "Lead", "Manganese", "Nickel", "Selenium", "Zinc"))

effluent_metals_total <- effluent_metals_base %>%
  transmute(
    analyte = paste0(element, ", total"),
    unit    = unit, dl = dl_t,
    lo = lo_t * 2, hi = hi_t * 4,   # elevated vs ambient
    test.code = "ICPMS", lab.id = "Commercial Lab"
  )

effluent_metals_dissolved <- effluent_metals_base %>%
  transmute(
    analyte = paste0(element, ", dissolved"),
    unit    = unit, dl = dl_d,
    lo = lo_d * 2, hi = hi_d * 4,
    test.code = "ICPMS", lab.id = "Commercial Lab"
  )

effluent_metals <- bind_rows(effluent_metals_total, effluent_metals_dissolved)

effluent_mercury <- mercury_analytes %>%
  mutate(lo = lo * 3, hi = hi * 5)  # elevated mercury in effluent

# Seawater metals use E469S submission code
sw_metals_total <- metals_total %>%
  filter(!grepl("^(Calcium|Magnesium|Potassium|Sodium)", analyte)) %>%
  mutate(test.code = "E469S")

sw_metals_dissolved <- metals_dissolved %>%
  filter(!grepl("^(Calcium|Magnesium|Potassium|Sodium)", analyte)) %>%
  mutate(test.code = "E469S")

# ---- Result simulator -------------------------------------------------------

sim_result <- function(lo, hi, dl, n) {
  raw <- runif(n, lo, hi)
  ifelse(raw < dl, paste0("<", dl), as.character(signif(raw, 4)))
}

# ---- Helper: build a sampling block -----------------------------------------

make_block <- function(program, sites, sub_matrix, analytes_df, dates_vec,
                       source_label, duplicate_val = NA_character_,
                       hold_exceedance_prob = 0.07) {
  crossing(
    tibble(site = sites, sub.matrix = sub_matrix),
    analytes_df %>% select(analyte, test.code, lab.id, unit, dl, lo, hi),
    sample.date = dates_vec
  ) %>%
    mutate(
      program     = program,
      site.impute = site,
      sample.id   = paste0(site, "-", format(sample.date, "%m%d")),
      sample.time = "09:00",
      matrix      = "Water",
      source      = source_label,
      duplicate   = duplicate_val,
      qualifier   = NA_character_,
      hold.time   = sample(
        c("F", "T"), n(), replace = TRUE,
        prob = c(1 - hold_exceedance_prob, hold_exceedance_prob)
      ),
      result = mapply(sim_result, lo, hi, dl, 1)
    ) %>%
    select(-lo, -hi)
}

# ---- Build each sample group ------------------------------------------------

river1_analytes <- bind_rows(
  metals_total, metals_dissolved, mercury_analytes %>% select(-lo, -hi),
  river1_nutrients %>% select(-lo, -hi)
)
# add lo/hi back for make_block to use (merge from source tables)
river1_analytes_full <- bind_rows(
  metals_total,
  metals_dissolved,
  mercury_analytes,
  river1_nutrients
)

exposure_analytes_full <- bind_rows(
  metals_total,
  metals_dissolved,
  mercury_analytes,
  exposure_nutrients
)

exposure_sw_analytes_full <- bind_rows(
  sw_metals_total,
  sw_metals_dissolved,
  mercury_analytes,
  exposure_nutrients %>%
    mutate(test.code = if_else(
      analyte %in% c("Nitrogen - Ammonia as N", "Nitrogen - Nitrate as N",
                     "Solids, total suspended [TSS]", "Carbon, dissolved organic [DOC]",
                     "Phosphorus, total", "Biochemical oxygen demand [BOD]"),
      test.code, "E469S"
    ))
)

effluent_analytes_full <- bind_rows(
  effluent_metals,
  effluent_mercury,
  effluent_nutrients
)

river1_sites <- make_block(
  program        = "River1",
  sites          = c("SITE-01", "SITE-02", "SITE-REF"),
  sub_matrix     = "Freshwater",
  analytes_df    = river1_analytes_full,
  dates_vec      = dates,
  source_label   = "Contractor"
)

exposure_fw <- make_block(
  program        = "Exposure",
  sites          = c("EXP-01", "REF-01"),
  sub_matrix     = "Freshwater",
  analytes_df    = exposure_analytes_full,
  dates_vec      = dates,
  source_label   = "Contractor"
)

exposure_sw <- make_block(
  program        = "Exposure",
  sites          = "EXP-02",
  sub_matrix     = "Seawater",
  analytes_df    = exposure_sw_analytes_full,
  dates_vec      = dates,
  source_label   = "Contractor"
)

effluent_rows <- make_block(
  program        = "Exposure",
  sites          = "Effluent-T1",
  sub_matrix     = "Effluent",
  analytes_df    = effluent_analytes_full,
  dates_vec      = dates,
  source_label   = "Client",           # effluent collected by Client
  hold_exceedance_prob = 0
)

# Blanks: almost all <DL; deliberate lead exceedance in the FB for QAQC demo
blank_analytes_full <- river1_analytes_full

blank_rows <- make_block(
  program        = "River1",
  sites          = c("SITE-01-FB", "SITE-01-TB"),
  sub_matrix     = "Freshwater",
  analytes_df    = blank_analytes_full,
  dates_vec      = dates[1],
  source_label   = "Contractor",
  hold_exceedance_prob = 0
) %>%
  mutate(
    result = case_when(
      analyte == "Lead, total" & site == "SITE-01-FB" ~
        as.character(signif(dl * 7.5, 4)),
      TRUE ~ paste0("<", dl)
    )
  )

# Duplicate: same analytes and date as SITE-01; duplicate column holds parent site
dup_rows <- make_block(
  program        = "River1",
  sites          = "SITE-01-DUP",
  sub_matrix     = "Freshwater",
  analytes_df    = river1_analytes_full,
  dates_vec      = dup_date,
  source_label   = "Contractor",
  duplicate_val  = "SITE-01"
) %>%
  mutate(sample.time = "09:30")

# ---- Assemble lab dataset (no in-situ columns) ------------------------------

wq_raw <- bind_rows(
  river1_sites, exposure_fw, exposure_sw,
  effluent_rows, blank_rows, dup_rows
) %>%
  select(
    program, site, site.impute, sample.id, sample.date, sample.time,
    matrix, sub.matrix,
    analyte, test.code, lab.id, result, dl, unit,
    qualifier, hold.time, source, duplicate
  ) %>%
  arrange(program, site, sample.date, analyte)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
write_csv(wq_raw, "data/raw/wq_raw.csv")
message(sprintf("Wrote %d rows to data/raw/wq_raw.csv", nrow(wq_raw)))

# ---- In-situ field measurements (one row per site × date) -------------------
# Saved separately; left-joined to harmonized lab data downstream.
# Blank and duplicate sites share readings with their parent or are NA.

site_base <- tibble::tribble(
  ~site,          ~sp.cond.uS.cm, ~pH,   ~sal.ppt, ~temp.C, ~diss.O2.mg.L, ~diss.O2.percent, ~redox.ORP, ~turbidity.NTU, ~e.coordinate, ~n.coordinate,
  "SITE-01",            352,       7.4,    0.2,     9.8,      10.2,          88,               240,          3.1,          502300,        5452100,
  "SITE-02",            418,       7.3,    0.3,     9.5,       9.8,          84,               235,          5.4,          502800,        5451600,
  "SITE-REF",           284,       7.6,    0.2,    10.1,      10.8,          91,               255,          1.8,          501800,        5452800,
  "EXP-01",             601,       7.2,    0.4,     9.2,       9.5,          82,               228,          8.2,          497200,        5448400,
  "EXP-02",           19200,       7.6,   14.8,     8.9,       8.3,          86,               220,         12.1,          496700,        5447900,
  "REF-01",             308,       7.5,    0.2,    10.3,      10.9,          92,               258,          2.2,          497700,        5449000,
  "Effluent-T1",       2150,       7.1,    1.6,    12.4,       7.8,          67,               195,         28.6,          496200,        5447400
)

# Expand to all sampling dates with slight random variation per date
set.seed(123)
wq_insitu <- crossing(site_base, sample.date = dates) %>%
  mutate(
    # small temporal variation to make each date distinct
    sp.cond.uS.cm  = round(sp.cond.uS.cm  * runif(n(), 0.96, 1.04)),
    pH             = round(pH              + runif(n(), -0.1,  0.1), 1),
    temp.C         = round(temp.C          + runif(n(), -0.5,  0.5), 1),
    diss.O2.mg.L   = round(diss.O2.mg.L   + runif(n(), -0.3,  0.3), 1),
    diss.O2.percent= round(diss.O2.percent + runif(n(), -2,    2)),
    redox.ORP      = round(redox.ORP       + runif(n(), -5,    5)),
    turbidity.NTU  = round(turbidity.NTU   * runif(n(), 0.90, 1.10), 1),
    depth.m        = 0.5
  ) %>%
  arrange(site, sample.date)

write_csv(wq_insitu, "data/raw/wq_insitu.csv")
message(sprintf("Wrote %d rows to data/raw/wq_insitu.csv", nrow(wq_insitu)))
