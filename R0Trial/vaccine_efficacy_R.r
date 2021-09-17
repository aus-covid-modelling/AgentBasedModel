# compute the expected vaccine efficacy against breakthrough infections that is
# due to a reduction in the case attributable fraction

# do this so that we can compare it against the estimated overall vaccine
# efficacy against breakthrough infections

# then compute correction factors to the infectiousness for vaccinated cases in
# each of the 4 combinations of vaccine type and  number of doses

# p_s = proportion symptomatic among unvaccinated infections

# p_s_star = proportion symptomatic among vaccinated infections

# VE_s = vaccine efficacy (% reduction) against symptomatic disease

# VE_o = vaccine efficacy (% reduction) against onward infectiousness of
# breakthrough infections

# VE_o(s) = vaccine efficacy (% reduction) against onward infectiousness of
# breakthrough infections, that is mediated by reductions in the fraction
# symptomatic

# absolute infectiousness (z is a transmission probability) of all unvaccinated
# people (symptomatic or asymptomatic)
# z * (p_s + (1 - p_s) / 2)

# absolute infectiousness (z is a transmission probability) of all vaccinated
# people (symptomatic or asymptomatic)
# z * (p_s_star + (1 - p_s_star) / 2)

# relative infectiousness of vaccinated versus unvaccinated (breakthrough)
# infections
#  (z * (p_s_star + (1 - p_s_star) / 2)) / (z * (p_s + (1 - p_s) / 2))
#   = (p_s_star + (1 - p_s_star) / 2) / (p_s + (1 - p_s) / 2)
#   = (1 + p_s_star) / (1 + p_s)

# proportion symptomatic among vaccinated infections
# p_s_star = p_s * (1 - VE_s)

# correction is the multiplier for reduction in onward transmission we want,
# divided by the one expected via reduction in symptomatics
# correction = ((1 - VE_o) / (1 - VE_o_s))

# this ^ can be multiplied by the probability of infection (by age bin, vaccine
# type, vaccine doses) to adjust it ot the correct value

library(tidyverse)

efficacies <- tibble::tribble(
  ~vaccine, ~dose, ~ve_infection, ~ve_symptom, ~ve_onward,
  "pfizer",     1,            30,          33,         46,
  "pfizer",     2,            79,          83,         65,
      "AZ",     1,            18,          33,         48,
      "AZ",     2,            60,          61,         65
) %>%
  mutate(
    across(
      starts_with("ve"),
      ~./100
    )
  )

clinical_fractions <- structure(
  list(
    age_group_5y = c("0-4", "5-9", "10-14", "15-19", 
                     "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", 
                     "55-59", "60-64", "65-69", "70-74", "75-79", "80+"),
    clinical_fraction_mean = c(
      0.29, 
      0.29, 0.21, 0.21, 0.27, 0.27, 0.33, 0.33, 0.4, 0.4, 0.49, 0.49, 
      0.63, 0.63, 0.69, 0.69, 0.69
    )
  ),
  row.names = c(NA, -17L),
  class = c("tbl_df", "tbl", "data.frame")
)

clinical_fractions %>%
  left_join(
    efficacies,
    by = character()
  ) %>%
  rename(
    clinical_fraction_unvacc = clinical_fraction_mean
  ) %>%
  mutate(
    clinical_fraction_vacc = clinical_fraction_unvacc * (1 - ve_symptom),
    ve_onward_via_symptom = 1 - ((1 + clinical_fraction_vacc) / (1 + clinical_fraction_unvacc)),
    onward_correction = ((1 - ve_onward) / (1 - ve_onward_via_symptom))
  ) %>%
  relocate(
    clinical_fraction_unvacc, .after = clinical_fraction_vacc
  ) %>%
  relocate(
    ve_onward_via_symptom, .after = ve_onward
  )

# p_s <- seq(0.1, 0.9, by = 0.1)
# VE_s <- efficacies$ve_symptom[2]
# VE_o <- efficacies$ve_onward[2]
# p_s_star = p_s * (1 - VE_s)
# VE_o_s <- 1 - ((1 + p_s_star) / (1 + p_s))
# correction <- ((1 - VE_o) / (1 - VE_o_s))