# Water Quality Data Harmonization

Multi-lab water quality data arrive with inconsistent analyte names, mixed units, and detection limit notations. This vignette demonstrates a production-grade pipeline for standardizing raw lab exports into a clean, analysis-ready dataset.

The pipeline reflects work across multiple environmental monitoring programs. All program names, site identifiers, and laboratory references have been replaced with generic labels; synthetic data is used so the vignette runs standalone.

## What this covers

- Analyte name harmonization across lab submission codes (e.g. E100 vs E235 naming for N-species)
- Dissolved/total fraction parsing from analyte name suffixes
- Analysis group classification (Conventional Parameters, Nutrients, Major Ions, Total/Dissolved Metals)
- Freshwater vs. seawater classification by specific conductance
- Detection limit flagging and ½-DL substitution
- Unit standardization (µg/L → mg/L; mercury µg/L → ng/L) with DL scaling
- QAQC: field/travel blank screening at 5× DL
- QAQC: field duplicate RPD with censoring rules
- QAQC: holding-time exceedance flags
- QAQC: dissolved:total metals ratio check

## Structure

```
analysis/
  wq_harmonization.Rmd        # vignette — the reader-facing document
R/
  harmonize.R                 # analyte name map, group classification, water type
  units.R                     # unit conversion with DL propagation
  qaqc.R                      # blank screening, RPD, holding time, D:T ratio
  database.R                  # field data reshape, deduplication, week labels, cofactors
  utils.R                     # rounding spec, shared constants
data/
  generate_synthetic_data.R   # builds wq_raw.csv
  wq_raw.csv                  # 183-row synthetic dataset
outputs/
  figures/
```

## Running the vignette

Open `vignette-wq-harmonization.Rproj` in RStudio, then knit `analysis/wq_harmonization.Rmd`. The setup chunk sets the working directory to the project root so all relative paths resolve correctly.

The synthetic data file (`data/wq_raw.csv`) is included. If it is missing, the Rmd will regenerate it automatically by sourcing `data/generate_synthetic_data.R`.

## Key packages

`dplyr`, `tidyr`, `readr`, `lubridate`, `knitr`, `kableExtra`, `table.glue`
