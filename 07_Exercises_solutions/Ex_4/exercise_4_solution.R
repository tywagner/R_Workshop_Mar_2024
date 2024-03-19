# R Workshop, March 21-22 2024
# Ty Wagner, USGS

# Exercise 4: Temporal trends in fantail darter abundance
# Description: To gain experience fitting HLMs
# Dataset: Fantail darter abundance


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


# Fantail Darter data  ---------------------------------------------------

# Read in data
dart_dat <- read_csv("././02_Data/Fischer_data/fantaildarter.csv")

# Clean data
dart_dat <- dart_dat %>% 
  select(RecordID, Date, County, Basin, N_Detected, Lat, Long, HUC_8, HUC_12) %>% 
  mutate(Date = mdy(Date),
         year = year(Date),
         County = factor(County),
         Basin = factor(Basin),
         HUC_8 = factor(HUC_8),
         HUC_12 = factor(HUC_12),
         decade = floor(year/10)*10) %>% # Assign each date to a decade (not used in analysis)
  filter(!is.na(N_Detected)) %>% 
  rename_with(tolower) # make all column names lower case
  
str(dart_dat)

# Sample size per HUC 8
huc_subset <- dart_dat %>% 
  group_by(huc_8) %>% 
  summarize(n = n()) %>% 
  arrange(-n)

huc_subset$n = cell_spec(huc_subset$n, color = ifelse(huc_subset$n > 50, "red", "darkgray"))

# Create table
huc_subset %>%
  kbl(booktabs = TRUE, digits=1, escape = FALSE,
      col.names = c("HUC 8", "Number of observations")) %>%
  kable_paper(full_width = FALSE) 


# Summarize number of observations by HUC_8
# Select HUC_8's with at least 50 observations
dart_subset <- dart_dat %>% 
  group_by(huc_8) %>% 
  summarize(n = n()) %>% 
  filter(n > 50) %>%  # retain HUCs with > 50 observations
  arrange(-n)

dart_dat <- dart_dat %>% 
  subset(huc_8 %in% dart_subset$huc_8) %>% 
  mutate(log_n = log(n_detected+1)) %>% # log-transform N (add 1 for 0's)
  droplevels()

str(dart_dat)

# Map locations 

# Grab state boundaries from spData and
# transform the coordinate references system (crs)
# crs = 4326 = WGS84; WGS84 CRS is often used for lat and long positions
us_states2 <- st_transform(us_states, crs = 4326)
# Rename column
colnames(us_states2)[2] <- "State"
# Select state(s) of interest
selectStates <- c("Pennsylvania")
# Subset data for plotting
us_state_select <- us_states2[us_states2$State %in% selectStates, ]

map.dat <- st_as_sf(dart_dat, coords = c("long", "lat"), crs = 4326)

# Create simple map
ggplot() +
  geom_sf(data = us_state_select, color = "gray30", lwd=1, fill="grey80") +
  geom_sf(data=map.dat, shape=16, size = 1, aes(colour=log(n_detected+1))) +
  labs(title="Fantail Darter locations", y="Latitude", x="Longitude", color="log(N detected+1)") +
  scale_color_gradientn(colours = rainbow(10)) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

# Create data summary tables
# Summary of the mean number of detects by HUC 8
dart_tab1 <- dart_dat %>% 
  group_by(huc_8) %>% 
  summarize(mean_abundance = mean(n_detected))

# Summary of mean number of detects by HUC and year
dart_tab2 <- dart_dat %>% 
  group_by(huc_8, year) %>% 
  summarize(mean_logN = mean(log_n))

# Create a list that has the number detected by year for each HUC 8
dart_list <- split(dart_tab2$mean_logN, dart_tab2$huc_8)

# Add a new column to dart_tab1 that will contain our time trend of detections
# in the table
dart_tab1 <- data.frame(dart_tab1, "Trend" = " ")

# Create table
dart_tab1 %>%
  kbl(booktabs = TRUE, digits=1, 
      col.names = c("HUC 8", "Mean abundance", "Temporal trend (log[N])")) %>%
  kable_paper(full_width = FALSE) %>%
  column_spec(3, image = spec_plot(dart_list, same_lim = TRUE, minmax = list(),
                                   col = "black",width = 500,
                                   height = 80)) 

str(dart_dat)
# Plot data
ggplot(data = dart_dat, mapping = aes(x=year, y=log_n) )+ 
  facet_wrap(~huc_8) +
  geom_point(alpha=0.3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(size=13)) +
  theme(plot.margin = margin(0, .5, .5, .5, "cm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_x_continuous(breaks=seq(1934,2013,10)) +
  labs(title="", y="log(N)", x="Year") 


# Standardize year predictor
dart_dat <- dart_dat %>% 
  mutate(z_year = as.numeric(scale(year)))
  
######## ---------------------------------
# Fit varying intercept and slope model
m1 <- stan_glmer(formula = log_n ~ 1 + z_year + (1 + z_year | huc_8), 
                 family = gaussian,
                 data = dart_dat,
                 iter = 1500, chains = 3)
print(m1, digits=3)


# Prepare MCMC summary stats for a table
sum1 <- summary(m1, probs = c(0.025, 0.975), digits=2)
sum1 <- data.frame(sum1[-c(27, 28),c(1,3:5)] )
sum1$parameter <- rownames(sum1)
rownames(sum1) <- NULL 

# Summary for table
sum1[-c(27, 28),c(5,1:4)] %>% 
  kbl(caption = "MCMC posterior summaries", booktabs = TRUE, digits=2,
      col.names = c("Parameter",
                    "Posterior mean",
                    "SD",
                    "Lower 95% CI",
                    "Upper 95% CI")) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  kable_styling(latex_options=c("striped","scale_down"))

# Put m1 output as an array for use in mcmc_trace() function
posterior <- as.array(m1)

mcmc_trace(posterior, pars = c("(Intercept)", "z_year"))


# Grab mcmc draws for selected parameters
fits <- m1 %>%
  as_tibble() %>%
  rename(intercept = `(Intercept)`) %>%
  select(-starts_with("Sigma") ) # exclude the estimated variances


############ Function to sum to get random intercept and slopes
shift_draws <- function(draws) {
  sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
}
###################################

# Extract intercept and random huc_8 intercepts and add them together
ints <- fits %>%
  select(intercept, starts_with("b[(Intercept)"))
# Calculate  huc_8 specific intercepts
huc_ints <- shift_draws(as.matrix(ints))


# Extract slope and random huc_8 slopes and add them together
slopes <- fits %>%
  select(z_year, starts_with("b[z_year"))
# Calculate  huc_8 specific intercepts
huc_slopes <- shift_draws(as.matrix(slopes))


# Create predicted regression line and uncertainty for each HUC_8

##############################
# ### BEGIN plotting code
# Create numeric HUC indicator for plotting
dart_dat <- dart_dat %>% 
  mutate(huc_8_id = as.numeric(huc_8)) %>% 
  arrange(huc_8_id)

# Number of hucs
J <- max(dart_dat$huc_8_id)

# Get range of year predictor for each huc
z_year_range <- list()
for(i in 1:J){
  z_year_range[[i]] <- seq(min(dart_dat$z_year[dart_dat$huc_8_id==i]), max(dart_dat$z_year[dart_dat$huc_8_id==i]), length.out = 30)
}


# Container for predicted values for each group (i.e., each HUC)
linPredGroup  <- array(NA, c(dim(huc_ints)[1],length(z_year_range[[1]]),J)) 
dim(linPredGroup)

# Put each groups MCMC draws for all parameters in its own list
group.ints <- list()
group.slopes <- list()
for(m in 1:J){
  group.ints[[m]] <- huc_ints[,m]
  group.slopes[[m]] <- huc_slopes[,m]
}


for(p in 1:J){ # loop over groups (J)
    for(t in 1:length(z_year_range[[1]])){
      linPredGroup[ ,t,p] <- group.ints[[p]] + group.slopes[[p]] * z_year_range[[p]][t] 
    }	  
}

dim(linPredGroup)
# Create containers
# Store posterior means
meanProbGroup <- array(NA, c(length(z_year_range[[1]]),J) )
# Store CIs
upperCI.Group <- array(NA, c(length(z_year_range[[1]]),J) )
lowerCI.Group <- array(NA, c(length(z_year_range[[1]]),J) )

for(i in 1:J){
  # Means
  meanProbGroup[,i] <- apply(linPredGroup[,,i], 2, mean )
  # 95% CIs for fitted values
  upperCI.Group[,i] <- apply(linPredGroup[,,i], 2, quantile, probs=c(0.975) )
  lowerCI.Group[,i] <- apply(linPredGroup[,,i], 2, quantile, probs=c(0.025) )
}

#######################GGPLOT####################
# New data frame for plotting (not really necessary)
toplot <- dart_dat 
# If you use this scale it cuts out a few points from the graph
Ymin <- 0
Ymax <- max(dart_dat$log_n)

# Make a dataframe for lines
# x-values in z_year_range
for(i in 1:J){
  temp.data=data.frame("huc_8_id"=i,"z_year"=z_year_range[[i]], log_n=meanProbGroup[,i])
  if(i==1) line.plot=temp.data else line.plot=rbind(line.plot, temp.data)
  temp.CIs=data.frame("huc_8_id"=i,"z_year"=z_year_range[[i]], "lower"=lowerCI.Group[,i],"upper"=upperCI.Group[,i])
  if(i==1) ci.line.plot=temp.CIs else ci.line.plot=rbind(ci.line.plot,temp.CIs)
}

# New lables for facets
site.labs <- as.character(unique(dart_dat$huc_8) )
names(site.labs) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

# Function for passing through apply to get
# posterior probs of a decline
post_prob <- function(x){
  mean(x < 0)
}
# Posterior prob of negative trend
slope_probs <- apply(huc_slopes, 2, post_prob)
slope_probs <- as.vector(slope_probs)
# Create data frame required for adding probs in a facet for ggplot
slope_probs_plot <- data.frame(huc_8_id = 1:10, probs = slope_probs)

# Plot
ggplot() +
  geom_point(data=dart_dat, aes(z_year, log_n),colour="blue", size=.5) + 
  facet_wrap(~huc_8_id, labeller = labeller(huc_8_id = site.labs)) +
  theme_bw() +
  scale_y_continuous(limits = c(Ymin, Ymax)) +
  geom_line(data=line.plot, aes(z_year, log_n), lwd=1) +
  geom_ribbon(data=ci.line.plot, aes(x=z_year, ymax=upper, ymin=lower), fill="grey", alpha=.75) +
  theme(panel.grid = element_blank(),strip.background = element_blank()) +
  xlab("Year (standardized)")+
  ylab(expression(paste(log[e],'(N)' ))) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(size=13)) +
  geom_text(aes(x=-3, y=5, label = round(probs, digits = 2)), data=slope_probs_plot) # add post probs


