# nass_parse.R -- trim + pivot the raw USDA-NASS bee colony pull down to the
# variables in the capstone proposal's model formula, and compute the
# response variable (log-odds of proportional colony loss).
#
# Input:  data/nass_bee_colonies_2015_2026.rds   (from nass_pull.R)
# Output: data/nass_model_panel.csv / .rds -- one row per state x year x
#         quarter: state, year, period, varroa, diseases, pesticides,
#         max_colonies, lost_colonies, t, y

library(dplyr)
library(tidyr)
library(readr)

colonies <- readRDS("data/nass_bee_colonies_2015_2026.rds")

# ---- 1. Drop non-CONUS / non-state groupings -----------------------------
# HI is present in your pull (46 distinct state_alpha values) but Hawaii is
# not part of the contiguous US -- your proposal's own scope -- so it comes
# out here even though it isn't one of the four states you flagged as absent.
# "OT" (location_desc "OTHER STATES") is a NASS confidentiality aggregate:
# small-production states rolled together to avoid disclosing individual
# figures. Which states get folded into it isn't fixed across years, so it
# can't be a hierarchical grouping level -- drop it, don't treat it as a
# 45th/46th state.
colonies <- colonies %>%
  filter(!state_alpha %in% c("HI", "OT"))

# Sanity check: your data also has zero rows for AK, DE, NH, NV, and RI --
# five states, not the four your proposal text currently lists (NV, RI, DE,
# NH). Worth a one-line fix to the proposal, or confirming AK really doesn't
# report to this survey.
n_states <- n_distinct(colonies$state_alpha)
message(n_states, " states remain after dropping HI/OT (expect 44).")

# ---- 2. Pull out exactly the short_desc items the model formula needs ----
item_map <- c(
  "HONEY, BEE COLONIES - INVENTORY, MAX, MEASURED IN COLONIES"                              = "max_colonies",
  "HONEY, BEE COLONIES - LOSS, DEADOUT, MEASURED IN COLONIES"                                = "lost_colonies",
  "HONEY, BEE COLONIES, AFFECTED BY VARROA MITES - INVENTORY, MEASURED IN PCT OF COLONIES"   = "varroa",
  "HONEY, BEE COLONIES, AFFECTED BY DISEASE - INVENTORY, MEASURED IN PCT OF COLONIES"        = "diseases",
  "HONEY, BEE COLONIES, AFFECTED BY PESTICIDES - INVENTORY, MEASURED IN PCT OF COLONIES"     = "pesticides"
)

panel <- colonies %>%
  filter(short_desc %in% names(item_map)) %>%
  mutate(variable = recode(short_desc, !!!item_map)) %>%
  select(state = state_alpha, year, period = reference_period_desc, variable, value_num) %>%
  pivot_wider(names_from = variable, values_from = value_num)

# ---- 3. Response variable: log-odds of proportional colony loss ---------
# t_ijk = lost_colonies / max_colonies ; y_ijk = log(t / (1 - t))
panel <- panel %>%
  mutate(
    t = lost_colonies / max_colonies,
    y = log(t / (1 - t))
  )

# In this pull there are no t == 0 / t == 1 rows (checked directly), so no
# continuity correction is needed for the response variable itself. Re-check
# this after any year you add, since a single 0% or 100%-loss state-quarter
# would send y to +/-Inf.
n_undefined <- sum(is.infinite(panel$y) | is.nan(panel$y))
if (n_undefined > 0) {
  warning(n_undefined, " rows have undefined log-odds (t = 0 or 1) -- ",
          "apply a continuity correction, e.g. (lost+0.5)/(max+1), before modeling.")
}

# ---- 4. Missing stressor coverage ----------------------------------------
# ~18% of state-quarter rows are missing at least one of varroa/diseases/
# pesticides (NASS non-response, not structural zero), spread across most
# states rather than concentrated in a few. This needs a deliberate call
# before fitting: Insolia et al.'s own pipeline treated missing stressor
# cells as 0, which you could replicate for comparability -- but given
# you're already in a Bayesian framework, imputing them (e.g. via brms's
# missing-data support) is arguably the more defensible choice methodologically.
missing_any <- with(panel, is.na(varroa) | is.na(diseases) | is.na(pesticides))
message(sum(missing_any), " / ", nrow(panel), " rows (",
        round(100 * mean(missing_any), 1), "%) missing at least one stressor, across ",
        n_distinct(panel$state[missing_any]), " states -- decide: impute, drop, or treat as 0.")

write_csv(panel, "data/nass_model_panel.csv")
saveRDS(panel, "data/nass_model_panel.rds")

message(nrow(panel), " state-quarter rows | ", n_distinct(panel$state), " states | ",
        min(panel$year), "-", max(panel$year))
