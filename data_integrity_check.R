# Check OK Q1 loss for data integrity 
ok_2019 <- bee_data |>
  filter(year == 2019, state == "OK", quarter == "Q1") |>
  select(c(year, state, quarter, y, t))

print(ok_2019)

raw_nass <- read_rds(here("data", "nass_model_panel.rds"))


raw_nass_ok_2019 <- raw_nass |>
  filter(year == 2019, state == "OK")

print(raw_nass_ok_2019, n = 20)

bee_data |> filter(state == "OK") |>
  select(year, quarter, max_colonies, lost_colonies, t) |> print(n = 50)

max_colony_swings <- bee_data |>
  group_by(state) |>
  arrange(year, quarter) |>
  mutate(max_ratio = max_colonies / lag(max_colonies)) |>
  ungroup() |>
  filter(y < -5 | max_ratio > 3 | max_ratio < 1/3)

print(max_colony_swings, n = 100)

# Verify count of daily moments and check for null values
prism_daily <- read_rds(here("data", "prism_state_daily.rds"))
nrow(prism_daily)
colSums(is.na(prism_daily))