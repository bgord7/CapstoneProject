library(dplyr)
# install.packages("here")
library(here)
library(readr)
library(tidyverse)

usda_data  <- readRDS(here("data", "nass_model_panel.rds"))
prism_data <- readRDS(here("data", "prism_quarterly_moments_2015_2026.rds"))

dim(usda_data)
dim(prism_data)

# Need to change the values in the period column to Q1, Q2, Q3, Q4 for the USDA data
usda_data <- usda_data |>
  mutate(
    period = as.character(fct_recode(period,
      "Q1" = "JAN THRU MAR",
      "Q2" = "APR THRU JUN" ,
      "Q3" = "JUL THRU SEP",
      "Q4" = "OCT THRU DEC"    
      )
  )
)
  


# sanity-check the join keys before merging - type mismatches here
# (e.g. "Q1" vs 1L, full state name vs abbreviation) cause silent drops
str(usda_data[c("state", "year", "period")])
str(prism_data[c("state", "year", "period")])

class(usda_data$year)
class(prism_data$year) # There is a class misatch. Convert usda year to integer

usda_data <- usda_data |>
  mutate(year = as.integer(year))

anyNA(usda_data$year)

# check for unmatched combinations on both sides
unmatched_in_prism <- anti_join(usda_data, prism_data, by = c("state", "year", "period"))
unmatched_in_usda  <- anti_join(prism_data, usda_data, by = c("state", "year", "period"))

nrow(unmatched_in_prism)  # NASS state/qtr/yr with no PRISM weather match
nrow(unmatched_in_usda)   # PRISM state/qtr/yr with no NASS record

unmatched_in_usda |> count(state, sort = TRUE)
unmatched_in_usda |> count(year, period)

reduced_states <- unmatched_in_usda |> 
  filter(year >= 2024) |> 
  distinct(state) |>
  pull(state)


# full join keeps everything, filling NAs where one side is missing -
# use inner_join() instead if you'd rather only keep complete matches
period_levels <- c("Q1", "Q2", "Q3", "Q4")  # confirm exact NASS wording
bee_data <- full_join(usda_data, prism_data, by = c("state", "year", "period")) |>
  filter(!is.na(lost_colonies)) |>
  mutate(
  quarter = factor(period, 
                   levels = period_levels, 
                   ordered = FALSE),
  ) |>
  select(-period) |>
  filter(year < 2024 | (year == 2024 & quarter %in% c("Q1", "Q2"))) |>
  filter(!is.na(varroa))

bee_data <- bee_data |>
  mutate(year_z = year)

colnames(bee_data)

scale_cols <- c("varroa", "diseases", "pesticides", "year_z",
                grep("^(tmin|prec)_", names(bee_data), value = TRUE))

scaling <- data.frame(
  variable = scale_cols,
  center   = sapply(bee_data[scale_cols], mean, na.rm = TRUE),
  scale    = sapply(bee_data[scale_cols], sd,   na.rm = TRUE)
)
saveRDS(scaling, here("data", "scaling_constants.rds"))

scaling <- readRDS(here("data", "scaling_constants.rds"))
scaling[scaling$variable == "varroa",]

bee_data <- bee_data |>
  mutate(across(all_of(scale_cols), ~ as.numeric(scale(.x))))

sum(is.na(bee_data$quarter))

dim(bee_data)

saveRDS(bee_data, file = here("data", "bee_data.rds"))
write_csv(bee_data, here("data", "bee_data.csv"))
bee_data |> filter(state == "OK", year == 2019, quarter == "Q1")
bee_data |> filter(state == "VT", year == 2021, quarter == "Q3")

# Need to look at missing disease and pesticide reporting 
# is it concentrated in particular states, or broad?
disease_rate <- bee_data |> filter(year >= 2020) |> 
  group_by(state) |> 
  summarise(diseases_na_rate = mean(is.na(diseases)), .groups = "drop") |> 
  arrange(desc(diseases_na_rate))

print(disease_rate, n=44)

sum(is.na(bee_data$varroa))
sum(is.na(bee_data$diseases))
sum(is.na(bee_data$pesticides))
# does it correlate with the reduced-reporting states, or is it a separate phenomenon?
bee_data |> filter(year >= 2020) |> 
  group_by(series_truncated) |> 
  summarise(diseases_na_rate = mean(is.na(diseases)), pesticides_na_rate = mean(is.na(pesticides)))

bee_data |> 
  summarise(both_missing = sum(is.na(diseases) & is.na(pesticides)),
            only_diseases = sum(is.na(diseases) & !is.na(pesticides)),
            only_pesticides = sum(!is.na(diseases) & is.na(pesticides)))
