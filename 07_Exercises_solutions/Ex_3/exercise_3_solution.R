# R Workshop, March 21-22 2024
# Ty Wagner, USGS

# Exercise 3: Analysis of CPE data
# Description: To gain experience fitting HLMs
# Dataset: Brown trout CPE


# Load R package
library(dplyr) # data management
library(tidyverse)
library(ggplot2) # plot
library(lubridate) # dates
library(stringr) # manipulate character stings
library(sf) # map creation
library(spData) # spatial data
library(car)  # logit function
library(qs) # save and read in large files
library(cmdstanr) # Bayesian estimation using stan
library(rstanarm) # Bayesian models using R functions and stan
library(mcmcplots) # as.mcmc function
library(ggrepel)
library(kableExtra) # tables
library(ggmap) # Google API data
library(bayesplot)


# Adult brown trout CPE ---------------------------------------------------

# Read in data
adult_dat <- read_csv("././02_Data/Delaware_tailwaters/NBE WB 2017 to 2021 TNP ver 121622 djp 122022.csv")

str(adult_dat)

adult_dat <- adult_dat %>% 
  mutate(SampleDate = mdy(SampleDate),
         Site = factor(Site),
         Wild = factor(Wild),
         Stage = factor(Stage)) %>%
  filter(Stage != "Yearling" & Stage != "Young of year") %>%
  filter(Wild != "Stocked") %>%
  filter(SpeciesCd == 328) %>% # select brown trout
  rename_with(tolower) # make all column names lower case
  
# Get summary statistics: total catch, effort, discharge, and 
# calculate cpe
catch_dat <-adult_dat %>% 
  group_by(site, year, month, riverside) %>% 
  summarize(total_catch = n(), effort_hrs = mean(timefish), 
            discharge = mean(discharge), cpe = total_catch/effort_hrs) %>% 
  ungroup()


# Make a table of summary statistics
catch_dat %>% 
  kbl(caption = "Adult summary catch by site, month, and year", digits=2,
      col.names = c("Site",
                    "Year",
                    "Month",
                    "River side",
                    "Total catch",
                    "Total effort (hrs)",
                    "Discharge",
                    "CPE")) %>%
  kable_classic(full_width = F, html_font = "Cambria")


# Plot data
ggplot(data = catch_dat, mapping = aes(x=year, y=cpe)) + 
  facet_wrap(~site) +
  geom_point(aes(color=factor(month))) + 
  theme_bw() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(size=13)) +
  theme(plot.margin = margin(0, .5, .5, .5, "cm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.text=element_text(size=11)) +
  labs(title="", y="CPE (fish/hr)", x="Year", color="Month") 



# Standardize discharge
catch_dat <- catch_dat %>% 
  mutate(z_discharge = as.numeric(scale(discharge)),
         year = factor(year))


m1 <- stan_glmer(formula = total_catch ~ 1 + z_discharge + year +
                   (1 | month) + (1|site) , family = neg_binomial_2(link='log'),
                 offset = log(effort_hrs),
                 data = catch_dat,
                 iter = 1500, chains = 3)
print(m1, digits=3)

sum1 <- summary(m1, probs = c(0.025, 0.975), digits=2)
# Summary for table
sum1[-c(18,21,22),c(1,3:7)] %>% 
  kbl(caption = "MCMC posterior summaries", booktabs = TRUE, digits=2,
      col.names = c("Parameter",
                    "Posterior mean",
                    "SD",
                    "Lower 95% CI",
                    "Upper 95% CI",
                    "Effective sample size",
                    "Rhat")) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  kable_styling(latex_options=c("striped","scale_down"))

posterior <- as.array(m1)

mcmc_trace(posterior, pars = c("(Intercept)", "z_discharge"))


# Grab mcmc draws for selected parameters
fits <- m1 %>%
  as_tibble() %>%
  rename(intercept = `(Intercept)`) %>%
  select(-starts_with("Sigma"),-reciprocal_dispersion ) # exclude the random effects for no and the dispersion parameter


############ Function to sum to get random intercept and slopes
shift_draws <- function(draws) {
  sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
}
###################################

# Extract intercept and add to beta year values
ints <- fits %>%
  select(intercept, starts_with("year"))
# Calculate  intercepts
year_effects <- shift_draws(as.matrix(ints))
year_effects <- data.frame(ints$intercept, year_effects)
# Rename columns based on year
colnames(year_effects) <- 2017:2021
head(year_effects)

# Put year effects on cpe scale
year_cpe <- exp(year_effects)
head(year_cpe)

# Posterior mean for each year
mean_cpe <- apply(year_cpe, 2, mean)
# Lower CI
lower <- apply(year_cpe, 2, quantile, 0.025)
# Upper CI
upper <- apply(year_cpe, 2, quantile, 0.975)
# Year indicator for plotting
year <- 2017:2021

# Mgt threshold/goal CPE (fish/hr)
mgt_goal <- 100

# Probability that a given year's CPE is > 100 fish/hr mgt goal
# Write simple function for passing through apply
prob_fn <- function(x){
  mean(x > mgt_goal)
}
# Probabilities of exceeding mgt goal
mgt_probs <- apply(year_cpe, 2, prob_fn)

plot_dat <- data.frame(year, mean_cpe, lower, upper, mgt_probs)

ggplot(plot_dat, aes(x=year, y=mean_cpe)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  geom_point(size=2.0) +
  geom_text(aes(label = round(mgt_probs, digits = 2)),
            position=position_dodge(width=0.9), vjust=-0.75, 
            hjust = 1.3,size=3.5) +
  ylim(0,400) +
  ylab('CPE (fish/hr)') +
  xlab('Year') +
  theme_bw() +
  geom_hline(yintercept = mgt_goal, linetype="dashed", color='red') +
  theme( axis.title = element_text(size=12), axis.text = element_text(size=12) )


#--------------------- Fig non-HM GLM

m_glm <- stan_glm(formula = total_catch ~ 1 + z_discharge + year, 
               family = neg_binomial_2(link='log'),
                 offset = log(effort_hrs),
                 data = catch_dat,
                 iter = 1500, chains = 3)
print(m_glm, digits=3)
summary(m_glm, probs = c(0.025, 0.975), digits=2)

# Grab mcmc draws for selected parameters
fits2 <- m_glm %>%
  as_tibble() %>%
  rename(intercept = `(Intercept)`) %>%
  select(-reciprocal_dispersion )


############ Function to sum to get random intercept and slopes
shift_draws <- function(draws) {
  sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
}
###################################

ints <- fits2 %>%
  select(intercept, starts_with("year"))
# Calculate intercepts
year_effects <- shift_draws(as.matrix(ints))
year_effects <- data.frame(ints$intercept, year_effects)
# Rename columns based on year
colnames(year_effects) <- 2017:2021
head(year_effects)

# Put year effects on cpe scale
year_cpe <- exp(year_effects)
head(year_cpe)

# Posterior mean for each year
mean_cpe <- apply(year_cpe, 2, mean)
# Lower CI
lower <- apply(year_cpe, 2, quantile, 0.025)
# Upper CI
upper <- apply(year_cpe, 2, quantile, 0.975)
# Year indicator for plotting
year <- 2017:2021

# Mgt goal CPE (fish/hr)
mgt_goal <- 100

# Probability that a given year's CPE is > 100 fish/hr mgt goal
# Write simple function for passing through apply
prob_fn <- function(x){
  mean(x > mgt_goal)
}
# Probabilities of exceeding mgt goal
mgt_probs <- apply(year_cpe, 2, prob_fn)

plot_dat <- data.frame(year, mean_cpe, lower, upper, mgt_probs)

ggplot(plot_dat, aes(x=year, y=mean_cpe)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  geom_point(size=2.0) +
  geom_text(aes(label = round(mgt_probs, digits = 2)),
            position=position_dodge(width=0.9), vjust=-0.75, 
            hjust = 1.3,size=3.5) +
  ylim(0,400) +
  ylab('CPE (fish/hr)') +
  xlab('Year') +
  theme_bw() +
  labs(title = "Standard GLM") +
  geom_hline(yintercept = mgt_goal, linetype="dashed", color='red') +
  theme( axis.title = element_text(size=12), axis.text = element_text(size=12) )


### OLD CODE

# Grab posterior samples as an array for use in 
# bayesplot::mcmc_trace
# posterior <- as.array(m2)
# 
# mcmc_trace(posterior, pars = c("(Intercept)"))
# 
# 
# # Grab mcmc draws for selected parameters
# fits <- m2 %>%
#   as_tibble() %>%
#   rename(intercept = `(Intercept)`) %>%
#   select(-starts_with("Sigma"),-reciprocal_dispersion )
# 
# 
# ############ Function to sum to get random intercept and slopes
# shift_draws <- function(draws) {
#   sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
# }
# ###################################
# 
# # Extract population average intercept and random intercepts for adding to get random intercept values
# ints <- fits %>%
#   select(intercept, starts_with("b[(Intercept) year"))
# # Calculate random intercepts
# year_effects <- shift_draws(as.matrix(ints))
# # Rename columns based on year
# colnames(year_effects) <- 2017:2021
# head(year_effects)
# 
# # Put year effects on cpe scale
# year_cpe <- exp(year_effects)
# head(year_cpe)
# 
# mean_cpe <- apply(year_cpe, 2, mean)
# 
# plot(2017:2021, mean_cpe, pch=16)

# Fit model
# m1 <- stan_glmer(formula = total_catch ~ 1 + 
#                    (1 | month) + (1|site) + (1|year), family = neg_binomial_2(link='log'), 
#                  offset = log(effort_hrs),
#                  data = catch_dat,
#                  iter = 1500, chains = 3)
# print(m1, digits=3)

# m2 <- stan_glmer(formula = total_catch ~ 1 + z_discharge +
#                    (1 | month) + (1|site) + (1|year), family = neg_binomial_2(link='log'), 
#                  offset = log(effort_hrs),
#                  data = catch_dat,
#                  iter = 2500, chains = 3)
# print(m2, digits=3)
# 
# summary(m2, probs = c(0.025, 0.975), digits=2)
# 
# posterior <- as.array(m2)
# 
# mcmc_trace(posterior, pars = c("(Intercept)"))
# 
# 
# # Grab mcmc draws for selected parameters
# fits <- m2 %>%
#   as_tibble() %>%
#   rename(intercept = `(Intercept)`) %>%
#   select(-starts_with("Sigma"),-reciprocal_dispersion )
# 
# 
# ############ Function to sum to get random intercept and slopes
# shift_draws <- function(draws) {
#   sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
# }
# ###################################
# 
# # Extract population average intercept and random intercepts for adding to get random intercept values
# ints <- fits %>%
#   select(intercept, starts_with("b[(Intercept) year"))
# # Calculate random intercepts
# year_effects <- shift_draws(as.matrix(ints))
# # Rename columns based on year
# colnames(year_effects) <- 2017:2021
# head(year_effects)
# 
# # Put year effects on cpe scale
# year_cpe <- exp(year_effects)
# head(year_cpe)
# 
# mean_cpe <- apply(year_cpe, 2, mean)
# 
# plot(2017:2021, mean_cpe, pch=16)
# # Posterior prediction
# pred_dat <- with(catch_dat, 
#                    expand.grid(year = levels(factor(year)),
#                                site=0,
#                                month=0,
#                                z_discharge=0))
# ppd_years <- posterior_predict(m2, newdata = pred_dat, offset = 0.15,
#                                draws = 500)
# str(ppd_years)
# 
# ppd_yr_means <-apply(ppd_years, 2, mean)
# plot(2017:2021, ppd_yr_means, pch=16)
