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


# Channel catfish length at age -------------------------------------------


# Read in Channel Catfish (CC) growth data
ccgrowth <- read_csv("././02_Data/Growth/LAKECCFlengthage.csv")

# str(ccgrowth)

# Clean data
ccgrowth <- ccgrowth %>% 
  mutate(site = factor(WATER)) %>% 
  rename("length" = LENGTH,
          "age" = AGE) 

# Get summary statistics
t1 <- ccgrowth %>% 
  group_by(site) %>% 
  summarize(n = n(), median_length = quantile(length, c(0.5)), min_length = min(length), max_length = max(length),
            median_age = quantile(age, c(0.5)), min_age = min(age), max_age = max(age))

# Make a table of summary statistics
t1 %>% 
  kbl(caption = "Channel catfish length-at-age summary statistics.",
      col.names = c("Site",
                    "N",
                    "Median length",
                    "Min length",
                    "Max length",
                    "Median age",
                    "Min age",
                    "Max age")) %>%
  kable_classic(full_width = F, html_font = "Cambria")

# Read in lake-level data
lake_dat <- read_csv("././02_Data/Growth/LakeCCF_locationinfo.csv")
str(lake_dat)

# Grab lat/long, lake area, and mean depth
lake_dat <- lake_dat %>% 
  select(site, Latitude, Longitude, AreaSurf, DepthMean) %>% 
  mutate(site = factor(site))

# Merge with growth data using site name
ccgrowth <- ccgrowth %>% 
  left_join(lake_dat, by=c("site")) %>%  # merge
  filter(site != "Lake Erie") %>% # Exclude Lake Erie
  droplevels()

# Map sites
# Make smbgrowth a sf object for plotting unique Lat/Long sites
ccmap.dat <- ccgrowth %>% 
  distinct(Longitude, Latitude, .keep_all = TRUE) %>% 
  select(site, AreaSurf, DepthMean, Latitude, Longitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 

# Use Google API to grab map of PA
get_map(location = 'pennsylvania', zoom=7, color="bw", maptype = "terrain") %>%
  ggmap() +
  geom_sf(data = ccmap.dat, size = 5,aes(colour=log(AreaSurf)), inherit.aes = FALSE) +
  ylab("Latitude") +
  xlab("Longitude") +
  geom_text_repel(data = ccmap.dat, aes(label = site, geometry = geometry),
                  stat = "sf_coordinates", size = 4, colour = "darkred",
                  inherit.aes = FALSE) +
  labs(color='log(Surface area)')

# Plot data
ggplot(ccgrowth, aes(x=age, y=length)) +
  geom_point() +
  facet_wrap(~site, ncol = 2) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  ylab("Length (mm)") +
  xlab("Age (yrs)") +
  theme(axis.text = element_text(size = 13), # increase axis numbers size
        axis.title = element_text(size = 13), # increase axis labels size
        strip.text.x = element_text(size=13)) # increase facet label size
# ggsave("figures/CC_length_age.pdf", height=8, width=8, units="in")


# Fit hierarchical growth model using stan and cmdstanr package
str(ccgrowth)

model_dat <- ccgrowth %>% 
  arrange(site, YEAR) %>% 
  mutate(siteID = as.numeric(site)) 

head(model_dat)
tail(model_dat)


# Load data
stan.data <- list(nFish = dim(model_dat)[1], # Number of fish
                  nSites = max(model_dat$siteID), # Number of sites
                  siteID = model_dat$siteID, # System ID for each observation
                  age = model_dat$age, # Fish age
                  length = model_dat$length,  # Length of fish
                  hp_tau = 1.5, # Scale
                  hp_sigma = 10, # Observation error
                  hp_omega = 2 # Correlation matrix

)

str(stan.data)

# identifying stan file path

stanfile <- "./09_stan_models/VBmodel_rawt0.stan"


# creating stan model
mod <- cmdstan_model(stanfile)

# sampling from posterior
out <- mod$sample(data = stan.data, 
                  seed = 1,
                  chains = 3,
                  step_size = 0.1,
                  parallel_chains = 3,
                  iter_warmup = 1000,
                  iter_sampling = 2000,
                  refresh = 10)

# qsave(out, './04_ModelOutput/stan_out_cc_growth.qs')
# out <- qread('./04_ModelOutput/stan_out_cc_growth.qs')

# Site-specific growth parameters
out$summary(
  variables = c("Linf", "K", "t0"), 
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)


post_table <- out$summary(
  variables = c("Linf", "K", "t0"),
  posterior::default_summary_measures()
)

post_table %>% 
  kbl(caption = "Posterior summaries.",digits=1,
      col.names = c("Parameter",
                    "Mean",
                    "Median",
                    "SD",
                    "MAD",
                    "Q5",
                    "Q95")) %>%
  kable_classic(full_width = F, html_font = "Cambria")

# draws_arr <- out$draws() # or format="array"
# str(draws_arr)
# mcmc_trace(draws_arr, regex_pars = c("Linf"))


# draws x variables data frame
draws_df <- out$draws(format = "df")
str(draws_df)


linfs <- draws_df %>% 
  select(starts_with("Linf")) %>% 
  data.frame()
head(linfs)

ks <- draws_df %>% 
  select(starts_with("K")) %>% 
  data.frame()
head(ks)

t0s <- draws_df %>% 
  select(starts_with("t0")) %>% 
  data.frame()
head(t0s)

## Begin code for plotting
# # age  value for each site 
age_range_J <- list()
J <- max(model_dat$siteID) # number of sites
n.age <- 50 # number of ages to predict across
for(i in 1:J){
  age_range_J[[i]] <- seq(min(model_dat$age[model_dat$siteID==i]), max(model_dat$age[model_dat$siteID==i]), length.out = n.age)
}

# Number of mcmc samples
nsim <- dim(draws_df)[1]

# Container for predicted values
est.lineB <- array(NA, c(nsim,n.age,J) )
dim(est.lineB)


for(j in 1:J){ # loop over groups (J)
    for(t in 1:n.age){
      est.lineB[,t,j] <-  linfs[,j] * 
        (1 -exp(-ks[,j] * 
                  (age_range_J[[j]][t] - t0s[,j]))) 
    }
  }


groupMean <- array(NA, c(n.age,J) )
upper.CIB <- array(NA, c(n.age,J) )
lower.CIB <- array(NA, c(n.age,J) )

for(i in 1:J){
  # Posterior means
  groupMean[,i] <- apply(est.lineB[,,i], 2, mean )
  # 95% CIs for fitted values
  upper.CIB[,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.95) )
  lower.CIB[,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.05) )
}



### GGPLOT
# Generate fake y-axis for creating plot
z <- seq(min(model_dat$length),max(model_dat$length),length=n.age) 

Ymin <- 0
Ymax <- max(model_dat$length)

#make a dataframe for lines
#xvalues in depth_std_range_J
for(i in 1:J){
  temp.data=data.frame("siteID"=i,"Zage"=age_range_J[[i]], length=groupMean[,i])
  if(i==1) line.plot=temp.data else line.plot=rbind(line.plot, temp.data)
  temp.CIs=data.frame("siteID"=i,"Zage"=age_range_J[[i]], "length.lower"=lower.CIB[,i],"length.upper"=upper.CIB[,i])
  if(i==1) ci.line.plot=temp.CIs else ci.line.plot=rbind(ci.line.plot,temp.CIs)
}

head(line.plot)
head(ci.line.plot)

site.labs <- as.character(unique(model_dat$site) )
names(site.labs) <- c("1", "2", "3", "4", "5", "6", "7", "8")

ggplot()+geom_point(data=model_dat, aes(age, length),colour="blue", size=.5)+
  facet_wrap(~siteID, ncol = 2, labeller = labeller(siteID = site.labs) )+
  theme_bw()+
  scale_y_continuous(limits = c(Ymin, Ymax)) +
  geom_line(data=line.plot, aes(Zage, length), lwd=1) +
  geom_ribbon(data=ci.line.plot, aes(x=Zage, ymax=length.upper, ymin=length.lower), fill="grey", alpha=.75)+
  theme(panel.grid = element_blank(),strip.background = element_blank())+
  xlab('Age (yrs)')+ylab("Length (mm)")

# Calculate omega (early growth in mm/year)
omegas <- ks * linfs
omegas <- as_tibble(omegas, .name_repair = ~paste0("omega", 1:dim(omegas)[2]))
head(omegas)

mean_omega <- apply(omegas, 2, mean)
lower_omega <- apply(omegas, 2, quantile, 0.025)
upper_omega <- apply(omegas, 2, quantile, 0.975)
site <- unique(model_dat$site)

omega_dat <- data.frame(site, mean_omega, lower_omega, upper_omega)
head(omega_dat)

ggplot(omega_dat, aes(site, mean_omega)) +
  geom_pointrange(aes(ymax=upper_omega, ymin=lower_omega),size=0.5) +
  labs(title = "", y = expression(paste(omega, ' (mm/yr)')), x = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Regress omega vs. lake area and depth

# covDat <- model_dat %>% 
#   group_by(site) %>% 
#   summarize(depth = mean(DepthMean), surfarea = mean(AreaSurf))
# 
# covDat <- data.frame(covDat, mean_omega)
# 
# plot(mean_omega ~ depth, data=covDat)
# plot(mean_omega ~ log(surfarea), data=covDat)

