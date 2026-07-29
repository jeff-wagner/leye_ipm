####Clean Annotated code
### Week 3 only of surveys from 2014-2025. Week 3 is most representative of a closed population of Lesser Yellowlegs

#generalized linear mixed model, looking for trend in abundance, not actual abundance
#install relevant packages
install.packages("rjags")
install.packages("runjags")
install.packages("ggplot")
library(runjags)
library(rjags)
library(MCMCvis)
library(ggplot2)
load.module("glm")

d <- read.csv("data_2025.csv", header = TRUE) #input data 2014-2025

# label years from 1 to 12, 2014-2025
d$year <- d$year - min(d$year) + 1

# mcmc settings
niters <- 10^4 #number of iterations to monitor
nburn <- niters / 10
nchains <- 3

# fit model to week 3
d3 <- d[d$week == 3, ]

#covariate key below
#et is time surveyed (min)
#ed is distance traveled by volunteer during survey (miles)
#wc is distance to coast of each wetland
#wt is wetland type, 3 options.1=freshwater ponds / lakes / riverine, 2=Estuarine and marine deepwater / estuarine and marine wetland, 3=freshwater emergent wetland / freshwater forested / shrub wetland
#wu is urban or rural. 1=urban, 2=rural. Labelled as urban if a 500m radius circle around wetland had more than 3 residential buildings in it
#wh is highway presence. 1= yes, 2= no. Yes if one of 3 major highways in Southcentral Alaska were within a 500 m radius of wetland

## Week 3 data only
n <- dim(d3)[1] #281 observations for week 3

#check data
count <- d3$count #LEYE count
w <- d3$wetland #wetland id for each observation in week 3, 32 wetlands
y <- d3$year #12 years for each week for each wetland
nw <- length(unique(w)) #32, number of wetlands, extracts unique elements (removes duplicates)
ny <- length(unique(y)) #12 years of data collection

#standardize data
et <- as.numeric(scale(d3$et)) #center and scale numeric values of time to improve model fitting. Used function scale so et recognized as numeric
ed <- as.numeric(scale(d3$ed)) #center and scale numeric values of distance to improve model fitting
wc <- as.numeric(scale(aggregate(wc ~ wetland, unique, data = d3)[, 2])) #allow distance to coast to be comparable
wa <- as.numeric(scale(aggregate(wa ~ wetland, unique, data = d3)[, 2])) #allow wetland shape to be comparable
wh <- aggregate(wh ~ wetland, unique, data = d3)[, 2] #highway presence. Aggregate splits data into subsets, computes summary stats, returns in convenient form, Simplifies highway data per wetland
wt <- aggregate(wt ~ wetland, unique, data = d3)[, 2] #wetland type. The [,2] command says we only want second column input as data into Rjags
wu <- aggregate(wu ~ wetland, unique, data = d3)[, 2] #urban/rural
nwt <- length(unique(wt)) #There are 3 options for wetland type
nwh <- length(unique(wh)) #There are 2 options for highway presence
nwu <- length(unique(wu)) #There are 2 options for urban/rural

# set up and run model
data <- list(
  n = n,
  count = count,
  w = w,
  y = y,
  nw = nw,
  ny = ny,
  et = et,
  ed = ed,
  wc = wc,
  wa = wa,
  wh = wh,
  wt = wt,
  wu = wu, #defining elements in model
  nwt = nwt,
  nwh = nwh,
  nwu = nwu
)
monitor <- c(
  "b0",
  "b1",
  "bet",
  "bed", #variable names we care about #monitor regression coefficients, error
  "bwc",
  "bwa",
  "bwt",
  "bwh",
  "bwu",
  "sd.e",
  "sd.ew",
  "sd.ey",
  "countsim"
) #add countsim for zero inflation check
mod <- jags.model(
  "model david_added overdispersion.R",
  data = data,
  n.chains = nchains
) #run model
samp <- coda.samples(mod, monitor, n.iter = niters) #Generate posterior samples in mcmc.list format

##plot the posterior distributions for the model parameters
samp <- samp[(nburn + 1):niters, ]
res1 <- as.matrix(samp[[1]]) #res is a list of MCMC samples, each res is the first, second, and third chain
res2 <- as.matrix(samp[[2]])
res3 <- as.matrix(samp[[3]])
res <- rbind(res1, res2, res3)

#probability that the following elements have a positive effect on bird count, these will be manually added to MCMCplot
#Repeat for all covariates, use colnames(res) to find which column matches each covariate
b1 <- res[, 2] #b1 is the trend, which is 10% chance of increase, 90% chance of decrease over 12 years
quantile(res[, 2], probs = c(0.025, 0.975))
mean(b1 > 0) #0.102, 10% chance of increase
bed <- res[, 3]
mean(bed > 0)
bet <- res[, 4]
mean(bet > 0)
bwa <- res[, 5]
mean(b1 > 0)
bwc <- res[, 6]
mean(bwc > 0)
bwh <- res[, 8] #bwh[1] is just going to be zero effect, need to check the other categories
mean(bwh > 0)
bwt <- res[, 10] #This is bwt[2], [1] is zero effect
mean(bwt > 0)
bwt <- res[, 11] #This is bwt[3], [1] is zero effect
mean(bwt > 0)
bwu <- res[, 13] #This is bwu[2], [1] is zero effect
mean(bwu > 0)

#make smaller plots showing probability of + effect of each covariate on bird count
mcmc.mod <- as.mcmc(res)
summary(mcmc.mod) #probability of positive effect on covariate. This should match the above line 78 for each element
par(mfrow = c(1, 1))
MCMCplot(
  mcmc.mod,
  params = c(
    "b1",
    "bed",
    "bet",
    "bwa",
    "bwc",
    "bwh[2]",
    "bwt[2]",
    "bwt[3]",
    "bwu[2]"
  ),
  ISB = FALSE,
  labels = c(
    "Year",
    "Distance traveled",
    "Time surveyed",
    "Wetland area",
    "Distance to coast",
    "No highway nearby",
    "Estuarine and marine deepwater/wetland",
    "Freshwater shrub wetland",
    "Rural"
  )
)
prob.pos <- MCMCvis::MCMCsummary(mcmc.mod, pg0 = TRUE)
#manually add mean p>0 to the right side of each line.
mtext("0.90", side = 4, line = -2, at = c(9), las = 2, font = 2) #in summary, There is a 0.1034074 chance that trend is increasing over time, 10% chance, which means 90% chance of decrease
mtext("0.62", side = 4, line = -2, at = c(8), las = 2, font = 2) #distance traveled, 8 from bottom of graph
mtext("0.99", side = 4, line = -2, at = c(7), las = 2, font = 2) #Time surveyed
mtext("0.62", side = 4, line = -2, at = c(6), las = 2, font = 2) #wetland area
mtext("0.13", side = 4, line = -2, at = c(5), las = 2, font = 2) #distance to coast
mtext("0.54", side = 4, line = -2, at = c(4), las = 2, font = 2) #no highway nearby
mtext("0.90", side = 4, line = -2, at = c(3), las = 2, font = 2) #estuarine/marine deepwater/wetland
mtext("0.76", side = 4, line = -2, at = c(2), las = 2, font = 2) #freshwater shrub wetland
mtext("0.85", side = 4, line = -2, at = c(1), las = 2, font = 2) #rural environment

year_beta_chains <- MCMCchains(mcmc.mod, params = c("b1"))
mean(year_beta_chains)
sd(year_beta_chains)
sd(mu_beta_chains)
mean(mu_beta_chains)
mu_beta_chains <- MCMCchains(mcmc.mod, params = c("b0"))
n_samps <- length(year_beta_chains)
yr_pred <- seq(from = 1, to = 12, length.out = 100)
n_x_vals <- 100
y_pred_matrix <- matrix(NA, nrow = n_samps, ncol = n_x_vals)
for (i in 1:n_samps) {
  y_pred_matrix[i, ] <- exp(mu_beta_chains[i] + year_beta_chains[i] * yr_pred)
}
y_pred_mean <- apply(y_pred_matrix, 2, mean)
y_pred_ci <- apply(y_pred_matrix, 2, quantile, probs = c(0.025, 0.975))


# 4. Plot the results (using ggplot2)
plot_df <- data.frame(
  x = yr_pred,
  mean = y_pred_mean,
  lower_ci = y_pred_ci[1, ],
  upper_ci = y_pred_ci[2, ]
)
plot_df$untrans_x <- seq(2014, 2025, length.out = 100)
ggplot(plot_df, aes(x = x)) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    fill = "blue",
    alpha = 0.3
  ) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  labs(x = "Year", y = "Lesser Yellowlegs Count", ) +
  scale_x_continuous(
    breaks = seq(1, 12, length.out = 12),
    labels = seq(2014, 2025, length.out = 12)
  ) +
  theme_minimal(base_size = 18)

###questions: why is CI so large on graph but so small on MCMCplot? Shouldn't they be the same?, if it is because
#the GGplot uses mu_beta_chains[i], what would be the best way to reflect what the MCMC plot shows on the graph?
