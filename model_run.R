library(tidybayes)
# install.packages("brms")
library(here)
library(readr)
library(tidyverse)
library(rstanarm)
library(broom.mixed)
library(brms)
library(bayesplot)
library(shinystan)
install.packages("performance")
library(performance)


list.files(here("data"))

bee_data <- read_rds(here("data", "bee_data.rds"))
dir.create(here("fits"), showWarnings = FALSE)
# and use file = here("fits", "test_loss_normal") so paths match your data loading

summary(bee_data$y)

class(bee_data$year)
class(bee_data$year_z)
range(bee_data$t)

# Gaussian Model definition
bform <- bf(y ~ varroa + mi(diseases) + 
              mi(pesticides) + tmin_mean + tmin_sd + tmin_skew + 
              tmin_kurt + prec_mean + prec_sd + prec_skew + prec_kurt +
              quarter + year_z + (varroa | state), family = gaussian()) + 
  bf(diseases | mi() ~ quarter + year_z + (1 | state)) + 
  bf(pesticides | mi() ~ quarter + year_z + (1 | state)) + 
  set_rescor(rescor = FALSE)

priors <- c(
  set_prior("normal(0, 2.5)", class = "b", resp = "y"),
  set_prior("normal(0, 1)",   class = "sd",  
            coef = "Intercept", group = "state", resp = "y"),
  set_prior("normal(0, 0.5)", class = "sd",  
            coef = "varroa", group = "state", resp = "y"),
  set_prior("lkj(2)", class = "cor", 
            group = "state"),
  set_prior("normal(0, 1)", class = "b", resp = "pesticides"),
  set_prior("normal(0, 1)", class = "b", resp = "diseases")
)

cols <- c("varroa", "diseases", "pesticides", "year_z",
          grep("^(tmin|prec)_", names(bee_data), value = TRUE))
sapply(bee_data[cols], function(x) c(mean = mean(x, na.rm = TRUE),
                                     sd = sd(x, na.rm = TRUE)))

validate_prior(formula = bform, prior = priors, data = bee_data)


test_fit_normal <- brm(bform, data = bee_data, prior = priors,
                   warmup = 500, iter = 1000, chains = 2, cores = 4,
                   control = list(adapt_delta = 0.95, max_treedepth = 12),
                   save_pars = save_pars(all = TRUE),
                   seed = 2026, file = here("fits", "test_fit_normal"),
                   file_refit = "on_change")

# summary(test_fit_normal)
# rstan::check_hmc_diagnostics(test_fit_normal$fit)
# rhat_vals_normal <- rhat(test_fit_normal)
# sort(rhat_vals_normal, decreasing = TRUE)[1:10]

# pp_check(test_fit_normal, resp = "y", type = "loo_pit_overlay")
# mcmc_trace(test_fit, size = 0.1)
# mcmc_dens_overlay(test_fit)
# pp_check(test_fit_normal, resp = "y")
# pp_check(test_fit_normal, resp = "y", type = "stat_grouped", group = "state")

#tidy(test_fit, effects = "fixed",
#   conf.int = TRUE, conf.level = 0.95)


#Student t model
bform1 <- bf(y ~ varroa + mi(diseases) + 
              mi(pesticides) + tmin_mean + tmin_sd + tmin_skew + 
              tmin_kurt + prec_mean + prec_sd + prec_skew + prec_kurt +
              quarter + year_z + (varroa | state), family = student()) + 
  bf(diseases | mi() ~ quarter + year_z + varroa +  (1 | state)) + 
  bf(pesticides | mi() ~ quarter + year_z + varroa + (1 | state)) + 
  set_rescor(rescor = FALSE)


validate_prior(formula = bform1, data = bee_data, prior = priors)

test_fit_student <- brm(bform1, data = bee_data, prior = priors,
                        warmup = 500, iter = 1000, chains = 2, cores = 4,
                        control = list(adapt_delta = 0.95, max_treedepth = 12),
                        save_pars = save_pars(all = TRUE),
                        seed = 2026, file = here("fits", "test_fit_student"),
                        file_refit = "on_change")

#summary(test_fit_student)
#rstan::check_hmc_diagnostics(test_fit_student$fit)
#rhat_vals_student <- rhat(test_fit_student)
#sort(rhat_vals_student, decreasing = TRUE)[1:10]
#pp_check(test_fit_student, resp = "y")
#pp_check(test_fit_student, resp = "y", type = "loo_pit_overlay")
#pp_check(test_fit_student, resp = "y", type = "stat_grouped", group = "state")

loo_normal  <- loo(test_fit_normal,  resp = "y")
loo_student <- loo(test_fit_student, resp = "y")
loo_compare(loo_normal, loo_student)

colony_loss <- brm(bform1, data = bee_data, prior = priors,
                   warmup = 2500, iter = 5000, chains = 4, cores = 4,
                   control = list(adapt_delta = 0.95, max_treedepth = 12),
                   save_pars = save_pars(all = TRUE),
                   seed = 2026, file = here("fits", "colony_loss_full"), 
                   file_refit = "on_change")

summary(colony_loss)
icc(colony_loss)
draws   <- as_draws_df(colony_loss)
tau0    <- draws$sd_state__y_Intercept
tau1    <- draws$sd_state__y_varroa
sigma   <- draws$sigma_y
nu      <- draws$nu_y
res_var <- sigma^2 * nu / (nu - 2)

icc_int <- tau0^2 / (tau0^2 + res_var)
icc_all <- (tau0^2 + tau1^2) / (tau0^2 + tau1^2 + res_var)
quantile(icc_int, c(.025, .5, .975))
quantile(icc_all, c(.025, .5, .975))
print(icc_int)
print(icc_all)
rstan::check_hmc_diagnostics(colony_loss$fit)
loo_full <- loo(colony_loss, resp = "y")
print(loo_full)
bayes_R2(colony_loss, resp = "y")
bayes_R2(colony_loss, resp = "y", re_formula = NA)
bayes_R2(colony_loss, resp = "y", re_formula = ~ (1 | state))
pp_check(colony_loss, resp = "y")
pp_check(colony_loss, resp = "y", type = "loo_pit_overlay")
pp_check(colony_loss, resp = "y", type = "stat_grouped", group = "state")
# mcmc_trace(colony_loss, regex_pars = "^b_y_|^sd_state__y|^sigma_y|^nu_y")
p <- mcmc_trace(colony_loss, regex_pars = "^b_y_|^sd_state__y|^sigma_y|^nu_y",
                facet_args = list(ncol = 3)) +
  scale_x_continuous(breaks = c(0, 1250, 2500))
ggsave("mcmc_trace.pdf", p, width = 10, height = 12, device = cairo_pdf)
# mcmc_dens_overlay
launch_shinystan(colony_loss)

priors_intercept <- c(
  set_prior("normal(0, 2.5)", class = "b", resp = "y"),
  set_prior("normal(0, 1)",   class = "sd",  
            coef = "Intercept", group = "state", resp = "y"),
  set_prior("normal(0, 1)", class = "b", resp = "pesticides"),
  set_prior("normal(0, 1)", class = "b", resp = "diseases")
)

bform2 <- bf(y ~ varroa + mi(diseases) + 
               mi(pesticides) + tmin_mean + tmin_sd + tmin_skew + 
               tmin_kurt + prec_mean + prec_sd + prec_skew + prec_kurt +
               quarter + year_z + (1 | state), family = student()) + 
  bf(diseases | mi() ~ quarter + year_z + (1 | state)) + 
  bf(pesticides | mi() ~ quarter + year_z + (1 | state)) + 
  set_rescor(rescor = FALSE)

intercept_model <- brm(bform2, data = bee_data, prior = priors_intercept,
                       warmup = 2500, iter = 5000, chains = 4, cores = 4,
                       control = list(adapt_delta = 0.95, max_treedepth = 12),
                       save_pars = save_pars(all = TRUE),
                       seed = 2026, file = here("fits", "colony_loss_intercept"), 
                       file_refit = "on_change")

bform3 <- bf(y ~ varroa + mi(diseases) + 
               mi(pesticides) + tmin_mean + tmin_sd + prec_mean + prec_sd + 
               quarter + year_z + (varroa | state), family = student()) + 
  bf(diseases | mi() ~ quarter + year_z + (1 | state)) + 
  bf(pesticides | mi() ~ quarter + year_z + (1 | state)) + 
  set_rescor(rescor = FALSE)

reduced_model <- brm(bform3, data = bee_data, prior = priors_intercept,
                       warmup = 2500, iter = 5000, chains = 4, cores = 4,
                       control = list(adapt_delta = 0.95, max_treedepth = 12),
                       save_pars = save_pars(all = TRUE),
                       seed = 2026, file = here("fits", "colony_loss_reduced"), 
                       file_refit = "on_change")

loo_full <- loo(colony_loss, resp = "y")
loo_intercept <- loo(intercept_model, resp = "y")
loo_reduced <- loo(reduced_model, resp = "y")
loo_compare(loo_full, loo_intercept, loo_reduced)
print(loo_full)
print(loo_reduced)
print(loo_intercept)
fixef(intercept_model)["y_varroa",]
fixef(colony_loss)["y_varroa",]

s <- readRDS(here("data", "scaling_constants.rds"))
sd_v <- s$scale[s$variable == "varroa"]
print(s)

b100 <- as_draws_df(colony_loss)$b_y_varroa * 100 / sd_v
c(estimate = mean(b100), quantile(b100, c(0.025, 0.975)))
for (v in c("diseases", "pesticides")) {
  sd_x <- s$scale[s$variable == v]
  d <- as_draws_df(colony_loss)[[paste0("bsp_y_mi", v)]] * 100 / sd_x
  print(c(v, round(c(mean(d), quantile(d, c(.025, .975))), 3)))
}

np <- brms::nuts_params(colony_loss)
sum(np$Value[np$Parameter == "divergent__"])
