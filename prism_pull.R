# prism_pull.R -- PRISM daily weather ingestion + state-quarterly moment
# aggregation for the capstone honey bee colony loss model.
#
# Produces a tidy state x year x quarter panel of the 8 weather moments
# from the proposal's Wijk vector:
#   tmin, tmin_sd, tmin_skew, tmin_kurt,
#   prec_mean, prec_sd, prec_skew, prec_kurt
#
# Requires: prism, terra, exactextractr, tigris, sf, dplyr, tidyr, purrr,
#           moments, lubridate, readr

library(prism)
library(terra)
library(exactextractr)
library(tigris)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(moments)   # skewness() / kurtosis()
library(lubridate)
library(readr)

options(tigris_use_cache = TRUE)
dir.create("data/prism", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Download daily PRISM rasters -----------------------------------
prism_set_dl_dir("data/prism")

# NOTE: this is a LOT of data -- ~11 years x 2 variables x 365 days of
# CONUS 4km daily rasters -- pulled one file per day from an external web
# service. A single multi-year call is fragile: any transient network blip
# (SSL hiccup, DNS timeout, a brief rate-limit from the PRISM service) kills
# the whole request, and you lose everything after the point of failure --
# which is exactly what happened partway through April 2017.
#
# This loops year-by-year (same pattern as your NASS pull) and retries each
# year a few times with a backoff pause before giving up, so one bad request
# only costs you that year's retries, not the whole run. It's also safe to
# just re-run this whole block if it dies again -- get_prism_dailys checks
# your archive folder and skips days it already has, so you won't burn time
# or bandwidth re-downloading anything through April 2017. (Worth a quick
# `prism_archive_ls()` afterward to confirm nothing's silently missing.)
years <- 2015:2026
last_date <- "2026-06-30"  # adjust to your last complete quarter

download_year_robust <- function(type, year, max_attempts = 4) {
  min_date <- sprintf("%d-01-01", year)
  max_date <- if (year == max(years)) last_date else sprintf("%d-12-31", year)
  
  for (attempt in seq_len(max_attempts)) {
    ok <- tryCatch({
      get_prism_dailys(type = type, minDate = min_date, maxDate = max_date, keepZip = FALSE)
      TRUE
    }, error = function(e) {
      message("  attempt ", attempt, "/", max_attempts, " failed (", type, " ", year, "): ",
              conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) return(invisible(TRUE))
    if (attempt < max_attempts) Sys.sleep(15 * attempt)  # back off before retrying
  }
  warning("Giving up on ", type, " ", year, " after ", max_attempts, " attempts -- ",
          "re-run this year manually once your connection is stable.")
  invisible(FALSE)
}

for (yr in years) {
  message("Downloading tmin ", yr, "...")
  download_year_robust("tmin", yr)
  message("Downloading ppt ", yr, "...")
  download_year_robust("ppt", yr)
}

# ---- 2. State boundaries (CONUS, matching your 46-state NASS panel) -----
# Drop AK/HI/territories/DC (not CONUS states) and the four states your
# proposal notes are absent from the bee survey (NV, RI, DE, NH), so the
# zonal extraction only runs where you actually have colony-loss data.
not_conus   <- c("AK", "HI", "PR", "VI", "GU", "AS", "MP", "DC")
no_bee_data <- c("NV", "RI", "DE", "NH")

states_sf <- tigris::states(cb = TRUE, resolution = "20m") %>%
  filter(!STUSPS %in% c(not_conus, no_bee_data))

# Pin to PRISM's real CRS once, rather than assuming 4326
sample_file <- pd_to_file(prism_archive_subset("tmin", "daily", resolution = "4km"))[1]
states_sf <- st_transform(states_sf, crs(rast(sample_file)))      

st_crs(states_sf)
st_crs(rast(sample_file))
# ---- 3. Zonal extraction: one value per state per day per variable ------
# NOTE: check `?prism_archive_subset` / `?pd_to_file` for your installed
# prism package version -- these helper names have shifted across releases.
extract_daily <- function(prism_type) {
  files <- prism_archive_subset(prism_type, "daily", resolution = "4km")
  paths <- pd_to_file(files)
  dates <- pd_get_date(files)
  
  map2_dfr(paths, dates, function(path, dt) {
    r <- rast(path)
    vals <- exact_extract(r, states_sf, "mean")  # area-weighted zonal mean
    tibble(state = states_sf$STUSPS, date = dt, value = vals)
  })
}

tmin_daily <- extract_daily("tmin") %>% rename(tmin = value)
ppt_daily   <- extract_daily("ppt")   %>% rename(ppt = value)

daily <- full_join(tmin_daily, ppt_daily, by = c("state", "date"))
saveRDS(daily, "data/prism_state_daily.rds")


# ---- 4. Map calendar dates to NASS's quarterly reporting periods --------
# These strings are my best guess at NASS's `reference_period_desc` values --
# verify against the real distinct values in your NASS pull before joining,
# and adjust here if they don't match exactly.
daily <- daily %>%
  mutate(
    year   = year(date),
    month  = month(date),
    period = case_when(
      month %in% 1:3   ~ "Q1",
      month %in% 4:6   ~ "Q2",
      month %in% 7:9   ~ "Q3",
      month %in% 10:12 ~ "Q4"
    )
  )

# ---- 5. Collapse daily -> quarterly, first four moments per variable -----
# NOTE: moments::kurtosis() returns Pearson's kurtosis (normal = 3), not
# excess kurtosis (normal = 0). Subtract 3 below if your model / Insolia et
# al.'s convention expects excess kurtosis -- worth checking their paper.
quarterly_moments <- daily %>%
  group_by(state, year, period) %>%
  summarise(
    tmin_mean = mean(tmin, na.rm = TRUE),
    tmin_sd   = sd(tmin, na.rm = TRUE),
    tmin_skew = skewness(tmin, na.rm = TRUE),
    tmin_kurt = kurtosis(tmin, na.rm = TRUE),
    prec_mean = mean(ppt, na.rm = TRUE),
    prec_sd   = sd(ppt, na.rm = TRUE),
    prec_skew = skewness(ppt, na.rm = TRUE),
    prec_kurt = kurtosis(ppt, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(quarterly_moments, "data/prism_quarterly_moments_2015_2026.csv")
saveRDS(quarterly_moments, "data/prism_quarterly_moments_2015_2026.rds")

message(nrow(quarterly_moments), " state-quarter weather records | ",
        n_distinct(quarterly_moments$state), " states")
raw_files <- list.files("data/prism", pattern = "tmin", recursive = TRUE, full.names = TRUE)
length(raw_files)          # how many tmin files actually exist on disk
length(tmin_files)         # what prism_archive_subset() thinks exists
range(tmin_dates)
length(tmin_dates)
length(unique(tmin_dates))
sum(expected_dates %in% tmin_dates)
expected_dates <- seq(as.Date("2015-01-01"), as.Date(last_date), by = "day")

tmin_files <- prism_archive_subset("tmin", "daily", resolution = "4km")
tmin_dates <- pd_get_date(tmin_files) |>
  as.Date()
class(expected_dates); head(expected_dates)
class(tmin_dates);     head(tmin_dates)
missing_tmin <- as.Date(setdiff(expected_dates, tmin_dates), origin = "1970-01-01")
missing_tmin   # 2015-06-28 should show up here — check if there are others too

missing_tmin_chr <- format(missing_tmin, "%Y-%m-%d")
for (d in missing_tmin_chr) {
  get_prism_dailys(type = "tmin", minDate = as.character(d), maxDate = as.character(d), keepZip = FALSE)
}