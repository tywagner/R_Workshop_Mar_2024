# rm(list=ls())
library(dplyr)
library(lme4)
library(rstanarm)
library(ggplot2)
library(tidyverse)


######## Data simulation function #####
#######################################
###    Function arguement definitions #
#######################################
# J = Number of groups (e.g., ecoregions)
# n.lower = lower number of sample units (e.g., lakes) per group
# n.upper = upper number of sample units (e.g., lakes) per group
# mu.a.true = Grand-mean intercept
# sigma.y.true = SD for level-1
# sigma.a.true = SD for intercepts
# sigma.b.true = SD for slopes
## CSI parameters
# gamma.b.0 <- Intercept term in relationship between level-2 predictor and level-1 slope parameters (i.e., slopes of x1-y regression)
# gamma.b.1 Slope-effect of z1 on slopes in relationship between x1 and y


################### FUNCTION FOR SIMULATING BIVARIATE NORMAL DATA
rbivariate <- function(mean.x = 70, sd.x=3, mean.y=162, sd.y=14, r=.50, iter=100) {
  z1 <- rnorm(iter)
  z2 <- rnorm(iter)
  x <- sqrt(1-r^2)*sd.x*z1 + r*sd.x*z2 + mean.x
  y <- sd.y*z2 + mean.y
  return(list(x,y))
}
################################################

DataSim <- function(J=9, n.lower=1, n.upper=20, mu.a=0.50, sigma.y=0.30, sigma.a=0.25, sigma.b=0.05, gamma.b.0=-0.2,
                    gamma.b.1=0, censorlevel=0.0){
  # Create a vector that repeats the number of lakes per region
  n <- round(runif(J, n.lower, n.upper))
  # Rep region-specific sample sizes
  site <- numeric()
  for(i in 1:J){
    a <- rep(i,n[i])
    site <- append(site,a)
  }
  # Level-2 covariate
  # z1 <- seq(-0.003, 0.40, length=J)
  # Level-1 covariate: varies randomly across regions and across range of predictor variable
  x1 <- numeric()
  for(i in 1:J){
    b <- runif(n[i],-0.83, 1.75)
    x1<-append(x1,b)
  }
  # Simulate varying coefficients from bivariate normal
  params <- rbivariate(mean.x = mu.a, sd.x = sigma.a,
                       mean.y = gamma.b.0, sd.y = sigma.b, r = -0.05, iter=J)
  # Simulate fake lognormally distributed data
  y.fake <- exp(rnorm(sum(n), params[[1]][site] + params[[2]][site] * x1 , sigma.y) )
  # z1.full <- z1[site]
  y.fake <- log(y.fake)
  dat <- data.frame(y.fake, site, x1)
  ## Create replicate y.fake for future analyses
  dat$y.fake2 <- dat$y.fake
  dat_list <- list(dat, params)
  # Left-censor data
  # censor <- quantile(dat$y.fake2, censorlevel)
  # censor2 <- rep(censor, sum(n))
  # dat <- data.frame(dat, censor2)
  # ### Censor data (data < 30th percentile is censored)
  # dat$detect <- ifelse(dat$y.fake < censor,0,1)
  # # Detection limit is same for all observations
  # dat$dl <- as.numeric(censor)
  # # If an observation was not detected then set to NA
  # dat$y.fake[dat$detect==0] <- NA
  # # Set dl to 1000 for detect == 1
  # dat$dl[dat$detect==1] <- 1000
  return( dat_list )

}


datSim <- DataSim()
dat <- datSim[[1]]
head(dat)
dim(dat)
### Fit varying intercept model using lmer function
# m1 <- lmer(y.fake2 ~ 1 + x1 + x1 + (1 + x1|site), data=dat)
# summary(m1)


# Raw data with lm fit
ggplot(data=dat, aes(x1,y.fake2))+
  geom_point( size=1)+
  theme_bw()+
  stat_smooth(method='lm') +
  facet_wrap(.~site, scales="free") +
  theme(panel.grid = element_blank(),strip.background = element_blank())+
  xlab('X') +
  ylab('Y') +
  theme(legend.position = "none", axis.title = element_text(size=12), axis.text = element_text(size=12))




# stan fit

m1 <- stan_lmer(y.fake2 ~ 1 + x1 + (1 + x1|site), data=dat,
                seed = 349, iter = 2000, chains = 3)
# print(m1)
# coef(m1)

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
rand.ints <- shift_draws(as.matrix(ints))
# Add new olumn names to random intercepts
rand.ints <- as_tibble(rand.ints, .name_repair = ~paste0("int", 1:dim(rand.ints)[2]))
head(rand.ints)

# Extract population average slope and random slopes for adding to get random intercept values
slopes <- fits %>%
  select(x1, starts_with("b[x1"))
# Calculate random slopes
rand.slopes <- shift_draws(as.matrix(slopes))
# Add new olumn names to random intercepts
rand.slopes <- as_tibble(rand.slopes, .name_repair = ~paste0("slope", 1:dim(rand.slopes)[2]))
head(rand.slopes)
apply(rand.slopes, 2, mean)
# Number of MCMC samples
n_samps <- dim(ints)[1]
# Create new predictor and design matrix
x_new <- data.frame(seq(min(dat$x1), max(dat$x1), length = 50) )
x_new$int <- 1
x_new <- cbind(x_new$int, x_new$seq.min.dat.x1...max.dat.x1...length...50.)
# Container to hold predictions
post_pred = array(NA, dim = c(n_samps, dim(x_new)[1], dim(rand.ints)[2]) )
for(j in 1:dim(rand.ints)[2]){
  BETA <- as.matrix(cbind(rand.ints[,j], rand.slopes[,j]))
  for(i in 1:n_samps){
    post_pred[i,, j] = x_new %*% BETA[i,]
  }
}

dim(post_pred)
str(post_pred)
# Posterior means
post_mean <- apply(post_pred, c(2,3), mean)
post_mean <- as_tibble(post_mean, .name_repair = ~paste0("site", 1:dim(post_mean)[2]))
# Add x for plotting
post_mean <- cbind(x_new[,2], post_mean)
post_mean <- post_mean %>%
  rename(X = 'x_new[, 2]')
dim(post_mean)
head(post_mean)

post_mean <- post_mean %>%
  pivot_longer(cols=c(starts_with("site")),
                    names_to='site',
                    values_to='Mean') %>%
  arrange(site)

post_mean$site <- rep(1:dim(rand.ints)[2], each = dim(x_new)[1])

# Lower 95% CI
lower <- apply(post_pred, c(2,3), quantile, c(0.025))
lower <- as_tibble(lower, .name_repair = ~paste0("site", 1:dim(lower)[2]))
dim(lower)
head(lower)

lower <- lower %>%
  pivot_longer(cols=c(starts_with("site")),
               names_to='site',
               values_to='lower') %>%
  arrange(site)

# Upper 95% CI
upper <- apply(post_pred, c(2,3), quantile, c(0.975))
upper <- as_tibble(upper, .name_repair = ~paste0("site", 1:dim(upper)[2]))
dim(upper)
head(upper)

upper <- upper %>%
  pivot_longer(cols=c(starts_with("site")),
               names_to='site',
               values_to='upper') %>%
  arrange(site)



plot.dat <- data.frame(post_mean, lower, upper)
head(plot.dat)

plot.dat$site <- as.factor(plot.dat$site)
dat$site <- as.factor(dat$site)

ggplot() +
  geom_point(data=dat, aes(x=x1, y = y.fake2)) +
  geom_ribbon(data=plot.dat, aes(x=X, ymax=upper, ymin=lower),  alpha=.15) +
  geom_line(data = plot.dat, aes(x = X, y = Mean), lwd=1) +
  facet_wrap(~site) +
  theme_bw() +
  theme(panel.grid = element_blank(),strip.background = element_blank())+
  xlab('X') +
  ylab('Y') +
  theme(legend.position = "none", axis.title = element_text(size=12), axis.text = element_text(size=12))



# Plot estimated and true regression parameters

# True parameters
beta0_true <- as.vector(unlist(datSim[[2]][1]))
beta1_true <- as.vector(unlist(datSim[[2]][2]))

# lm parameters
lm1 <- lm(y.fake2 ~ -1 + x1*site - x1, data=dat)
summary(lm1)
lm_ints <- as.vector(coef(lm1)[1:dim(rand.ints)[2]])
lm_slopes <- as.vector(coef(lm1)[(dim(rand.ints)[2]+1):(dim(rand.ints)[2]*2)])
# Parameter SEs
int_lm_ses <- as.vector(summary(lm1)$coefficients[, 2][1:dim(rand.ints)[2]])
slope_lm_ses <- as.vector(summary(lm1)$coefficients[, 2][(dim(rand.ints)[2]+1):(dim(rand.ints)[2]*2)])


# Bayes parameters
hm_ints <- as.vector(apply(rand.ints, 2, mean))
hm_slopes <- as.vector(apply(rand.slopes, 2, mean))
hm_ints_ses <- as.vector(apply(rand.ints, 2, sd))
hm_slopes_ses <- as.vector(apply(rand.slopes, 2, sd))

# Combine
intercepts <- c(beta0_true, lm_ints, hm_ints)
slopes <- c(beta1_true, lm_slopes, hm_slopes)
ses_ints <- c(rep(0, dim(rand.ints)[2]), int_lm_ses, hm_ints_ses)
ses_slopes <- c(rep(0, dim(rand.ints)[2]), slope_lm_ses, hm_slopes_ses)

estimates <- rep(c("Truth", "lm", "HM"), each = dim(rand.ints)[2])
parameter <- rep(c("Intercept", "Slope"), each = dim(rand.ints)[2] * 3)
site <- rep(1:dim(rand.ints)[2])

plot.dat2 <- data.frame(Model = c(estimates, estimates),
                        Parameter = parameter,
                        Value = c(intercepts, slopes),
                        SE = c(ses_ints, ses_slopes),
                        site = c(site, site, site, site, site, site))


ggplot(plot.dat2, aes(x = Value, y = site, xmin = Value-SE, xmax = Value+SE, color = Model)) +
  facet_wrap(~Parameter)+
  geom_point() +
  geom_errorbarh(height=.2) +
  theme_bw() +
  scale_y_continuous(breaks=seq(1,dim(rand.ints)[2],1)) +
  # theme(panel.grid = element_blank(),strip.background = element_blank())+
  xlab('Parameter values') +
  ylab('Site number') +
  theme(legend.position = "top", axis.title = element_text(size=12), axis.text = element_text(size=12),
        legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=4)))



pd <- position_dodge(0.8)
ggplot(plot.dat2, aes(x=site, y=Value,
                     color = Model)) +
  geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=.1, position=pd) +
  geom_point(position=pd, size=2.0) +
  facet_wrap(~Parameter) +
  scale_x_continuous(breaks=seq(1,dim(rand.ints)[2],1)) +
  ylab('Parameter values') +
  xlab('Site number') +
  theme_bw() +
  theme(legend.position = "top", axis.title = element_text(size=12), axis.text = element_text(size=12),
        legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=4)))
