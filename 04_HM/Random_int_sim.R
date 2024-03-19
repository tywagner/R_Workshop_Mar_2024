# rm(list=ls())
library(lme4)
library(ggplot2)
library(dplyr)
library(rstanarm)
library(tidyverse)

#####################################
### Simulate some data ###############
#####################################
#True parameters

# Mean for distribution of group means (these values come from actual CSI data for log10(TP))
g.0 <- 1.15
# Variance for level-1 sd
sigma.y <- 0.28
#Variance for level-2 sd (among-group sd)
sigma.a <- 0.20
# Intraclass correlation coefficient
sigma.a^2/(sigma.a^2 + sigma.y^2)
# Number of groups (e.g., lakes)
J = 50
#Number of samples per group (allow n to vary by group)
set.seed(12358)
n <- round(abs(rnorm(J, mean=2,sd=30))+1)

range(n)

# Create site indicator
site <- numeric()
for(i in 1:J){
	a<-rep(i,n[i])
	site<-append(site,a)
}
# site

################
#Simulate varying coefficients (alpha j's)
################
a.true <- rep(NA,J)
for(j in 1:J){
	a.true[j]<-rnorm(1, g.0, sigma.a)
}
a.true

####################
#Simulate fake data
####################
y.fake <- rep(NA, length(site))
for(i in 1:length(site)){
	y.fake[i] <- rnorm(1, a.true[site[i]], sigma.y)
}
# y.fake

##########################################
######## END Data simulation #############
##########################################

# Create dataframe of data
dat <- data.frame(y.fake,site)
dat

head(dat)
# Create site as a factor and call group
dat$group <- as.factor(dat$site)

# Fit varying intercept model using lmer function
m.lmer <- lmer(y.fake ~ 1 + (1|site), data=dat)
summary(m.lmer)

# Look at group 'shrinkage means'
coef(m.lmer)
range(coef(m.lmer)$site)

# Calculate ICC
sigma_y2 <- attr(VarCorr(m.lmer), "sc")^2 # Extract sigma y^2 from model (level-1 variance)
sigma_a2 <- as.numeric(VarCorr(m.lmer))   # Extract level-2 variance (sigam a^2)
ICC <- sigma_a2 / (sigma_a2 + sigma_y2)
ICC

samp.size <- dat %>%
  group_by(site) %>%
  summarise(n = n())

dat <- dat %>%
  left_join(samp.size, by = 'site') %>%
  mutate(n = as.numeric(n))

# Plot raw group means
ggplot(dat, aes(x = group, y = y.fake, fill=n)) +
  geom_boxplot() +
  scale_fill_viridis_c(option = "magma") +
  theme_bw() +
  geom_hline(yintercept=g.0, linetype="dashed",
             color = "red", linewidth=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  labs(title = "", y = "Response variable (Y)", x = "Group #", fill = 'Sample\n size')

# Stan fit
m1a <- stan_glm(formula = y.fake ~  -1 + as.factor(group),
                 family = gaussian,
                 data = dat,
                 seed = 349, iter = 3000, chains = 3)
# summary(m1)
# Grab posterior samples
posteriora <- as.data.frame(m1a)

# Grab random intercepts
rand.intsa <- posteriora[, grep("as.", colnames(posteriora)) ]


rand.ints2a  <- rand.intsa %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name)) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            upper = quantile(value, 0.975),
            lower = quantile(value, 0.025))

rand.ints2a <- data.frame(rand.ints2a, samp.size)


ggplot(rand.ints2a, aes(site, mean)) +
  geom_pointrange(aes(ymax=upper, ymin=lower, color=n),size=0.5) +
  scale_colour_viridis_c(option = "magma") +
  geom_hline(yintercept=g.0, linetype="dashed",
             color = "red", linewidth=1) +
  labs(title = "", y = "Estimate", x = "Group #", color = 'Sample\n size') +
  ylim(0.5, 2)



# Stan fit - hierarchical model

m1 <- stan_glmer(formula = y.fake ~  1 +
                   (1 | group),
                 family = gaussian,
                 data = dat,
                 seed = 349, iter = 2000, chains = 3)
# summary(m1)
# Grab posterior samples
posterior <- as.data.frame(m1)

# Grab random intercepts
rand.ints <- posterior[, grep("b", colnames(posterior)) ]
# Grab overall intercept
int <- posterior[, "(Intercept)"]
# Calculate group means
rand.ints <- rand.ints + int

rand.ints2  <- rand.ints %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name)) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
         upper = quantile(value, 0.975),
         lower = quantile(value, 0.025))

rand.ints2 <- data.frame(rand.ints2, samp.size)


ggplot(rand.ints2, aes(site, mean)) +
  geom_pointrange(aes(ymax=upper, ymin=lower, color=n),size=0.5) +
  scale_colour_viridis_c(option = "magma") +
  geom_hline(yintercept=g.0, linetype="dashed",
             color = "red", linewidth=1) +
  labs(title = "", y = "Estimate", x = "Group #", color = 'Sample\n size') +
  ylim(0.5,2)

