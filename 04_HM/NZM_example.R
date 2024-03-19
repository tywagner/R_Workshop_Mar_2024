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


# Read in NZM data
nzm_dat <- read_csv("././02_data/NZM_data/NZM Presence_Absence Dataset.csv")

str(nzm_dat)

# Clean data
nzm_dat <- nzm_dat %>%
  mutate(Present_Absent = factor(Present_Absent),
         observation_date = dmy(observation_date),
         year = year(observation_date)) %>%
  mutate(yvar = if_else(Present_Absent == "P", 1, 0)) %>%
  rename("long" = x,
         "lat" = y)


# Get summary statistics
t1nzm <- nzm_dat %>%
  group_by(year) %>%
  summarize(n = n(), Proportion_occupied = round(mean(yvar),2) )

# Make a table of summary statistics
t1nzm %>%
  kbl(caption = " ",
      col.names = c("Year",
                    "N", "Proportion occupied")) %>%
  kable_classic(full_width = F, html_font = "Cambria")


# nzm_dat %>%
#   group_by(year) %>%
#   summarize(n = n(), Proportion_occupied = round(mean(yvar),2) ) %>%
#   kbl(caption = " ",
#       col.names = c("Year",
#                     "N", "Proportion occupied")) %>%
#   kable_classic(full_width = F, html_font = "Cambria")


# Map sites
# Make NZM a sf object for plotting unique Lat/Long sites
nzm_map <- nzm_dat %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# Use Google API to grab map of PA
get_map(location = 'pennsylvania', zoom=7, color="bw", maptype = "terrain") %>%
  ggmap() +
  geom_sf(data = nzm_map, size = 3,aes(colour=Present_Absent), alpha=0.15, inherit.aes = FALSE) +
  geom_sf(data = nzm_map, size = 3,aes(colour=Present_Absent), shape=1, inherit.aes = FALSE) +
  ylab("Latitude") +
  xlab("Longitude") +
  labs(color='Present/absent')


# Stan fit - hierarchical model
# Simply looking at annual estimates of
# the probability of occurrence

m1 <- stan_glmer(formula = yvar ~  1 +
                   (1 | year),
                 family = binomial,
                 data = nzm_dat,
                 seed = 349, iter = 2000, chains = 3)

summary(m1)

# Grab mcmc draws for intercepts and slopes
fits <- m1 %>%
  as_tibble() %>%
  rename(intercept = `(Intercept)`) %>%
  select(-starts_with("Sigma") )

############ Function to sum to get random intercept and slopes
shift_draws <- function(draws) {
  sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
}
###################################

# Extract population average intercept and random intercepts for adding to get random intercept values
ints <- fits %>%
  select(intercept, starts_with("b[(Intercept)"))
# Calculate random intercepts
year_effects_logit <- shift_draws(as.matrix(ints))
head(year_effects_logit)

year_effects_probs <- apply(year_effects_logit, 2, plogis)

# Posterior summaries of year effects
post_mean <- as.vector(apply(year_effects_probs, 2, mean))
Lower <- as.vector(apply(year_effects_probs, 2, quantile, 0.025))
Upper <- as.vector(apply(year_effects_probs, 2, quantile, 0.975))
year <- t1nzm$year
n <- t1nzm$n
Proportion_occupied <- post_mean
Estimate <- "HLM"

hlm_ests <- data.frame(year, n, Proportion_occupied,
                       Estimate, Lower, Upper)

# Overall mean probability of occurrence
overall_mean <- mean(plogis(ints$intercept))

# Get table with same headings as post summaries for merging
t1nzm$Estimate <- "Raw"
t1nzm$Lower <- NA
t1nzm$Upper <- NA

# Combine raw estimates and hlm estimates for plotting
plot_dat <- rbind(t1nzm, hlm_ests)

# The points overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.18)

ggplot(plot_dat, aes(x=as.factor(year), y=Proportion_occupied,
                     color = Estimate)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1, position=pd) +
  geom_point(position=pd, size=3) +
  xlab('Year') +
  theme_bw() +
  geom_hline(yintercept=overall_mean ) +
  ylab('Probability of NZM occurrence') +
  theme(legend.position = "top", axis.title = element_text(size=12), axis.text = element_text(size=12),
        legend.title = element_text(size=20), legend.text = element_text(size=20))


