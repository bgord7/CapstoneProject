# nass_pull.R -- USDA-NASS honey bee colony data ingestion for the capstone
# Loops by year to stay under the 50,000-record API cap, filters to bee-colony
# records, cleans NASS withheld codes, and caches to disk so you don't re-hit
# the API every time you knit.

library(rnassqs)
library(dplyr)
library(purrr)
library(readr)
library(stringr)

API_KEY <- Sys.getenv("NASSQS_TOKEN")
nassqs_auth(key = API_KEY)                      # reads NASSQS_TOKEN from .Renviron
dir.create("data", showWarnings = FALSE)

years <- 2015:2026

# One request per year keeps each call well under the 50k cap.
pull_year <- function(yr) {
  params <- list(
    commodity_desc = "HONEY",      # <- the key filter; without it you pull everything
    agg_level_desc = "STATE",
    source_desc    = "SURVEY",
    year           = yr
  )
  message("Pulling ", yr)
  nassqs(params)
}

raw <- map_dfr(years, pull_year)

raw <- as_tibble(raw)

# --- Inspect before filtering -------------------------------------------------
# Run this once to see the exact data-item strings NASS uses, so you can pick
# your response, offset, and stressor predictors deliberately:

inspection <- raw %>% 
  distinct(short_desc) %>% 
  arrange(short_desc)

print(inspection, n = Inf)

# Keep only bee-colony records (inventory / loss / renovated / added / stressors).
# The annual honey production, price, stocks, and yield rows drop out here.
# Note: CCD-symptom loss is published at the U.S. level only, so it naturally
# does not appear in this STATE-level pull.
colonies <- raw %>%
  filter(str_detect(short_desc, "BEE COLONIES"))

# --- Clean the Value column ---------------------------------------------------
# NASS encodes withheld / rounded cells as (D), (Z), (X), (NA). Strip thousands
# separators and coerce; the parenthetical codes become NA.
clean_value <- function(x) {
  x <- str_remove_all(x, ",")
  suppressWarnings(as.numeric(x))
}

colonies <- colonies %>%
  mutate(
    value_num = clean_value(Value),
    withheld  = is.na(value_num) & str_detect(Value, "\\(")   # flag suppressed cells
  )

# --- Cache --------------------------------------------------------------------
saveRDS(colonies, "data/nass_bee_colonies_2015_2026.rds")
write_csv(colonies, "data/nass_bee_colonies_2015_2026.csv")


message(nrow(colonies), " colony records | ",
        n_distinct(colonies$state_alpha), " states | ",
        min(colonies$year), "-", max(colonies$year))
