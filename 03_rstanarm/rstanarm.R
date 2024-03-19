library(dplyr) # data management
library(tidyverse)
library(ggplot2) # plot
library(lubridate) # dates
library(bslib) # a modern UI toolkit for Shiny and R Markdown
library(stringr) # manipulate character stings
library(sf) # map creation
library(spData) # spatial data
library(car)  # logit function
library(ggmap) # Use Google Map services

# Stan
library(rstanarm)
library(bayesplot)
library(kableExtra)
library(loo)
# Let's simulate some data
# Number of observations
N <- 50
# Intercept
B0 <- 1
# Slope
B1 <- 0.50
# Predictor variable
X <- rnorm(N, 0, 1)
# Expected value
mu <- B0 + B1 * X
# Residual standard deviation
sigma <- 0.50
# Response with error
Y <- rnorm(n = N, mean = mu, sd = sigma)

data <- data.frame(Y, X)

# Fit regression model
summary(lm(Y~X, data=data))

# Fit Bayesian regression model
stan_m1 <- stan_glm(Y ~ X, data = data,
  family = gaussian)

summary(stan_m1, digits = 3)

# Check trace plots
# Convert stan object to array for plotting
# using bayesplot functions
posterior <- as.array(stan_m1)

mcmc_trace(posterior, pars = c("(Intercept)", "X"))


# Extract posterior samples for summarizing, plotting, etc.
samples <- stan_m1 %>%
  as.data.frame
head(samples)

# Plot posterior distributions
# Go from wide to long for facet plotting
post_plot <- samples %>%
  gather(parameter, value, `(Intercept)`:sigma, factor_key=TRUE)
head(post_plot)

# Plot
ggplot(data = post_plot, aes(x = value)) +
  facet_wrap(~parameter)+
  geom_histogram(bins = 50, color = 'blue', alpha = 0.5) +
  theme_bw() +
  theme(strip.text.y = element_text(size=11, face="bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))+
  ylab("Count") + xlab("Value")

# Look at credible interval
posterior_interval(stan_m1, prob = 0.95)

# Can do this by hand too
post_plot %>%
  group_by(parameter) %>%
  summarise(mean = mean(value),
            lower = quantile(value,probs = 0.025),
            upper = quantile(value,probs = 0.975)) %>%
  kbl(caption = "Posterior estimtes", digits=3,
    col.names = c("Parameter", "Mean", "Lower CI", "Upper CI")) %>%
  kable_classic(full_width = F, html_font = "Cambria")

# Plot estimates
# Shows 80% and 90% credible intervals
plot(stan_m1)

#--- Plot model fit
# Create a fake predictor variable corresponding to our X
fake.x <- seq(min(X), max(X), length=50)
# Container to hold predicted values
preds <- array(NA, c(dim(samples)[1], length(fake.x)))

for(i in 1:length(fake.x)){
  preds[,i] <- samples$`(Intercept)` + samples$X * fake.x[i]
}

head(preds)
mean.fit <- apply(preds, 2, mean)
lower <- apply(preds, 2, quantile, 0.025)
upper <- apply(preds, 2, quantile, 0.975)

plot_fit_m1 <- data.frame(mean.fit, lower, upper, fake.x)

ggplot() +
  geom_line(data=plot_fit_m1, aes(x=fake.x, y=mean.fit), linewidth=0.8) +
  geom_ribbon(data=plot_fit_m1, aes(x=fake.x, ymax=upper, ymin=lower),
              fill="green", alpha=0.3) +
  geom_point(data=data, aes(x=X, y=Y)) +
  xlab('X') +
  ylab('Y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12),
        axis.text = element_text(size=12))

# Stan arguments (options)

# Fit Bayesian regression model
stan_m2 <- stan_glm(Y ~ X, data = data,
                    family = gaussian,
                    iter = 2000, # of iterations per chain
                    warmup = 1000, # of iter to use in warmup
                    chains = 1, # of chains to run
                    cores = 1, # of parallel cores to run chains on
                    prior = normal(),
                    control =
                      list(max_treedepth = 15, # controls how far the model will look for a new proposal before giving up <!-- - higher values allow model to explore a flatter posterior -->
                           adapt_delta = 0.8) # he target rate that new proposals are accepted <!-- - higher means the model takes smaller steps -->
                            )

summary(stan_m2, digits = 3)

rstanarm::prior_summary(stan_m2)


## Informative prior for X

stan_m3 <- stan_glm(Y ~ X, data = data,
                    family = gaussian,
                    prior = normal(2,0.025), # prior on the model coefficients
                    prior_intercept = normal(autoscale = TRUE), # default prior for intercept
                    prior_aux = exponential(autoscale = TRUE) ) # default prior for sigma


summary(stan_m3, digits = 3)

samples3 <- stan_m3 %>%
  as.data.frame

post_plot3 <- samples3 %>%
  gather(parameter, value, `(Intercept)`:sigma, factor_key=TRUE)
head(post_plot3)

post_plot3 %>%
  group_by(parameter) %>%
  summarise(mean = mean(value),
            lower = quantile(value,probs = 0.025),
            upper = quantile(value,probs = 0.975)) %>%
  kbl(caption = "Posterior estimtes", digits=3,
      col.names = c("Parameter", "Mean", "Lower CI", "Upper CI")) %>%
  kable_classic(full_width = F, html_font = "Cambria")


rstan::Rhat(posterior[,,1])

rstan::check_hmc_diagnostics(stan_m1$stanfit)

## Compare samples from the posterior and prior
posterior_vs_prior(stan_m1)
posterior_vs_prior(stan_m3)


# Extract MCMC samples

# samples2 <- rstan::extract(stan_m1$stanfit, permute = TRUE)

# rstan::extract(stan_m1$stanfit, permute = TRUE) %>%
#   listviewer::jsonedit()



# Statistical inferences
mean(samples$X > 0)


# Compare two or more models

# Let's simulate some data
# Number of observations
N <- 100
# Intercept
B0 <- 1
# Slope 1
B1 <- -0.30
# Slope 2
B2 <- 0.15
# Predictor variables
X1 <- rnorm(N, 0, 1)
X2 <- rnorm(N, 0, 1)
# Expected value
mu <- B0 + B1 * X1 + B2 * X2
# Residual standard deviation
sigma <- 0.50
# Response with error
Y <- rnorm(n = N, mean = mu, sd = sigma)

data <- data.frame(Y, X1, X2)

# Fit regression model
summary(lm(Y~X1 + X2, data=data))

stan_mX1 <- stan_glm(Y ~ X1, data = data,
                    family = gaussian)

stan_mX1X2 <- stan_glm(Y ~ X1 + X2, data = data,
                       family = gaussian)

# Calculate LOO-IC, smaller is better
looX1 <- loo(stan_mX1)
looX1$estimates[3,]
looX1X2 <- loo(stan_mX1X2)
looX1X2$estimates[3,]

