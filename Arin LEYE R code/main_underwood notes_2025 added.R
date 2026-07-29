#generalized linear mixed model, more robust, as it is not trying to estimate actual abundance, only a trend in relative abundance
#Oct 2025, add 2025 data to this code from original file: main_underwood notes_2024 added

install.packages("rjags")
install.packages("runjags")
library(runjags)
library(rjags)
load.module("glm")

# NB wetlands with all zero counts removed from csv file

d <- read.csv("data_2025.csv", header = TRUE) #input data

# label years from 1 to 12, 2014-2025

d$year <- d$year - min(d$year) + 1

# mcmc settings

niters <- 10^4 #number of iterations to monitor
nburn <- niters / 10
nchains <- 3

# fit model separately for each week

d1 <- d[d$week == 1, ] #week 1
d2 <- d[d$week == 2, ] #week 2
d3 <- d[d$week == 3, ] #week 3

# only showing code for week 1 below, then week 2, then week 3. We end up only using week 3 for results as that is most representative of a closed population.
#LEYE GPS return dates not known when project started in 2013.

n <- dim(d1)[1] #283 observations for week 1

count <- d1$count #LEYE count for week 1
#z<-count>0        #check for zero-inflation? don't think this is right, do that in another code and there isn't any zero inflation or overdispersion
w <- d1$wetland #wetland id for each observation in week 1, 32 wetlands, 3 wetlands removed for having no sightings
y <- d1$year #12 years for each week for each wetland

nw <- length(unique(w)) #32, number of wetlands, extracts unique elements (removes duplicates)
ny <- length(unique(y)) #12 years of data collection

#scaling employed when using vectors or columns in a data frame. Allows you to compare data that isn't measured in the same way
#The normalizing of a dataset using the mean value and standard deviation is known as scaling.
#QUESTION: diff between this process and standardization? Same purpose of comparing data? How to know which to use?
#answer: same thing as standardization, subtracting mean and diving by SD. Only applies to numerical predictors.

et <- as.numeric(scale(d1$et)) #center and scale numeric values of time to improve model fitting. Used function scale so et recognized as numeric
ed <- as.numeric(scale(d1$ed)) #center and scale numeric values of distance to improve model fitting

wc <- as.numeric(scale(aggregate(wc ~ wetland, unique, data = d1)[, 2])) #allow distance to coast to be comparable
wa <- as.numeric(scale(aggregate(wa ~ wetland, unique, data = d1)[, 2])) #allow wetland shape to be comparable

#categorical covariates. Not numerical, dont need numerical standardization. aggregate command finds out unique values of wh for each wetland.Wetland 1 is 2 for highway, etc
wh <- aggregate(wh ~ wetland, unique, data = d1)[, 2] #highway presence. Aggregate splits data into subsets, computes summary stats, returns in convenient form, Simplifies highway data per wetland
wt <- aggregate(wt ~ wetland, unique, data = d1)[, 2] #wetland type. The [,2] command says we only want second column input as data into Rjags
wu <- aggregate(wu ~ wetland, unique, data = d1)[, 2] #urban/rural

#research aggregate command for data summarization in future.

nwt <- length(unique(wt)) #There are 3 options for wetland type
nwh <- length(unique(wh)) #There are 2 options for highway presence
nwu <- length(unique(wu)) #There are 2 options for urban/rural

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
  "sd.ey"
)
mod <- jags.model(
  "model david_added overdispersion.R",
  data = data,
  n.chains = nchains
) #run model
samp <- coda.samples(mod, monitor, n.iter = niters) #Generate posterior samples in mcmc.list format


colnames(samp[[1]]) #look at names of variables so we can define indicator
#bed=distance traveled, bet=time surveyed, bwa=wetland shape (area), bwc= distance to coast, bwh= highway presence (2 options), bwt=wetland type (3 options), bwu=urban/rural (2 options)
ind <- 1:16 #we are focusing on 16 regression coefficients, b0, b1, bed, bet, bwc...look at the convergence on these
ind <- ind[-c(7, 9, 12)] #1-16 minus 7,9,12 (categorical covariate options)

#Gelman and Rubin's convergence diagnostic. Values substantially above 1 indicate lack of convergence
gelman.diag(samp[, ind], multivariate = T) #take a look at convergence, looks good if close to 1. Diagnostics . 1.02 we good

##plot the posterior distributions for the model parameters
samp <- samp[(nburn + 1):niters, ]

res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])
res3 <- as.matrix(samp[[3]])

res <- rbind(res1, res2, res3)

# plot densities
par(mfrow = c(4, 4)) #change 4,4 to 1,1 to zoom in on each plot
for (i in ind) {
  plot(density(res[, i]), ylab = "", xlab = "", main = colnames(res)[i])
}
#what are we looking for? graph one of b1 is trend on log scale. 0 is a stable population
#type in colnames(res) to get column names
#if density includes 0 then not necessarily significant,
#there's just a chance it isn't But with Bayes, ##b1<-res[,2]  #quantile(res[,2],probs=c(0.025, 0.975)
#mean(b1>0), 0.0439, 4% chance trend  is increasing
#mean(b1<0) 0.956, 96% chance pop is going down.

#bed <-res[,3] type into console, then mean (bed>0), what are chances its over 0? 99? chance that distance is significant
#median(bed), 0.42

#bwa, 80% chance that wetland area has a positive impact on bird count. larger the wetland is the larger the count. some evidence that larger areas mean larger counts, but not totally conclusive, 20% chancethat larger areas have lower counts.
#bwh[2]: where you've got a highway, there is a certain chance count i slower. type head(res) into console to see values. then type bwh2 <-res[,8]. mean (bwh2<0), 72% chance that if by highway, lower count then if theres no highway.bwh measures the difference with and without highway
#bwt[1] is 0. any of them could be 0. then estimate other 2. bswt[2] tells diff between 2 and 1. bwt[3] tells diff between 3 and 1. to learn diff between 2 and 3, type colnames(res) into console, column 10 and 11. bwt<-res[,10], bwt<-res[,11]. bwt[3] is kind o fin middle. probability that. typ eino console, bwt3vs2<-bwt3-bwt2. calculate diff between 2 columns. then mean(bwt3vs2>0). 0.267. chance that diff is positive is 27%.. if looking at lower than 0. 73% chance that wetlany 3 has fewer birds than wetand 2.

# check for correlations between parameters
#we are looking for very strong correlation. mild is expected
#looking for very narrow plots, that means bad correlation. cloud is good, slightly narrow is okay. it is a little subjective though, how narrow is too narrow?

np <- 10^3
pairs(res[1:np, ind])
corr <- cor(res[, ind])
max(abs(corr - diag(diag(corr)))) #max correlation, dont care if neg or pos, if it's close to 1? if it's below 0.9 then its probs fine, above 0.9 then model might be funny
#if this number is less than 0.9 and graphs aren't very narrow, everything is probably fine

save.image("results.RData")


############Repeat for week 2
n <- dim(d2)[1] #291 observations for week 2

count <- d2$count #LEYE count
w <- d2$wetland #wetland id for each observation in week 2, 29 wetlands, 3 wetlands removed for having no sightings
y <- d2$year #10 years for each week for each wetland

nw <- length(unique(w)) #29, number of wetlands, extracts unique elements (removes duplicates)
ny <- length(unique(y)) #10 years of data collection

et <- as.numeric(scale(d2$et)) #center and scale numeric values of time to improve model fitting. Used function scale so et recognized as numeric
ed <- as.numeric(scale(d2$ed)) #center and scale numeric values of distance to improve model fitting

wc <- as.numeric(scale(aggregate(wc ~ wetland, unique, data = d2)[, 2])) #allow distance to coast to be comparable
wa <- as.numeric(scale(aggregate(wa ~ wetland, unique, data = d2)[, 2])) #allow wetland shape to be comparable

wh <- aggregate(wh ~ wetland, unique, data = d2)[, 2] #highway presence. Aggregate splits data into subsets, computes summary stats, returns in convenient form, Simplifies highway data per wetland
wt <- aggregate(wt ~ wetland, unique, data = d2)[, 2] #wetland type. The [,2] command says we only want second column input as data into Rjags
wu <- aggregate(wu ~ wetland, unique, data = d2)[, 2] #urban/rural

nwt <- length(unique(wt)) #There are 3 options for wetland type
nwh <- length(unique(wh)) #There are 2 options for highway presence
nwu <- length(unique(wu)) #There are 2 options for urban/rural

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
  "sd.ey"
)
mod <- jags.model(
  "model david_added overdispersion.R",
  data = data,
  n.chains = nchains
) #run model
samp <- coda.samples(mod, monitor, n.iter = niters) #Generate posterior samples in mcmc.list format
summary(mod)

colnames(samp[[1]]) #look at names of variables so we can define indicator
ind <- 1:16 #we are focusing on 16 regression coefficients, b0, b1, bed, bet, bwc...look at the convergence on these
ind <- ind[-c(7, 9, 12)] #1-16 minus 7,9,12 (categorical covariate options)

#Gelman and Rubin's convergence diagnostic. Values substantially above 1 indicate lack of convergence
gelman.diag(samp[, ind], multivariate = T) #take a look at convergence, looks good if close to 1. Diagnostics

##plot the posterior distributions for the model parameters
samp <- samp[(nburn + 1):niters, ]

res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])
res3 <- as.matrix(samp[[3]])

res <- rbind(res1, res2, res3)

# plot densities
par(mfrow = c(4, 4)) #chnage 4,4 to 1,1 to zoom in on each plot
for (i in ind) {
  plot(density(res[, i]), ylab = "", xlab = "", main = colnames(res)[i])
}

np <- 10^3
pairs(res[1:np, ind])
corr <- cor(res[, ind])
max(abs(corr - diag(diag(corr)))) #max correlation, dont care if neg or pos, if it's close to 1? if it's below 0.9 then its probs fine, above 0.9 then model might be funny
#if this number is less than 0.9 and graphs arent very narrow, everything is probs fine

save.image("results.RData")

###########
#######################Repeat for week 3    ################## Week 3 is most representative of population, use this for conclusions

n <- dim(d3)[1] #281 observations for week 3

count <- d3$count #LEYE count
w <- d3$wetland #wetland id for each observation in week 2, 29 wetlands, 3 wetlands removed for having no sightings
y <- d3$year #10 years for each week for each wetland
nw <- length(unique(w)) #32, number of wetlands, extracts unique elements (removes duplicates)
ny <- length(unique(y)) #12 years of data collection

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
summary(samp)
colnames(samp[[1]]) #look at names of variables so we can define indicator
ind <- 1:16 #we are focusing on 16 regression coefficients, b0, b1, bed, bet, bwc...look at the convergence on these
ind <- ind[-c(7, 9, 12)] #1-16 minus 7,9,12 (categorical covariate options)

#Gelman and Rubin's convergence diagnostic. Values substantially above 1 indicate lack of convergence
gelman.diag(samp[, ind], multivariate = T) #take a look at convergence, looks good if close to 1. Diagnostics

##plot the posterior distributions for the model parameters
samp <- samp[(nburn + 1):niters, ]

res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])
res3 <- as.matrix(samp[[3]])

res <- rbind(res1, res2, res3)

mu <- res[, 1:n] #is this needed? 1 to 255 observations for week 3
countsim <- res[, (n + 1):(2 * n)] #zero inflation check, error, subscript out of bounds. But it's supposed to be isn't it? We're simulating extra data?
colnames(countsim)

# plot densities
par(mfrow = c(4, 4)) #change 4,4 to 1,1 to zoom in on each plot
for (i in ind) {
  plot(density(res[, i]), ylab = "", xlab = "", main = colnames(res)[i])
}

np <- 10^3
pairs(res[1:np, ind])
corr <- cor(res[, ind])
max(abs(corr - diag(diag(corr)))) #max correlation, dont care if neg or pos, if it's close to 1? if it's below 0.9 then its probs fine, above 0.9 then model might be funny
#if this number is less than 0.9 and graphs arent very narrow, everything is probs fine

save.image("results.RData")

#what are we looking for? graph one of b1 is trend on log scale. 0 is a stable population
#type in colnames(res) to get column names
#if density includes 0 then not necessarily significant,
#there's just a chance it isn't But with Bayes, ##b1<-res[,2]  #quantile(res[,2],probs=c(0.025, 0.975)
#mean(b1>0), 0.099, 2% chance trend  is increasing
#mean(b1<0) 0.9005, 90% chance pop is going down.

#bed <-res[,3] type into console, then mean (bed>0), what are chances its over 0? 99? chance that distance is significant
#median(bed), 0.6296

#bwa <-res[,5] and mean (bwa>0)
#bwa, 62% chance that wetland area has a positive impact on bird count. larger the wetland is the larger the count. some evidence that larger areas mean larger counts, but not totally conclusive, 20% chancethat larger areas have lower counts.
#bwh <-res[,7] and mean (bwh>0)
#bwh[2]: where you've got a highway, there is a certain chance count i slower. type head(res) into console to see values. then type bwh2 <-res[,8]. mean (bwh2<0), 72% chance that if by highway, lower count then if theres no highway.bwh measures the difference with and without highway
#bwt <-res[,4] and mean (bwt>0)
#bwt[1] is 0. any of them could be 0. then estimate other 2. bswt[2] tells diff between 2 and 1. bwt[3] tells diff between 3 and 1. to learn diff between 2 and 3, type colnames(res) into console, column 10 and 11. bwt<-res[,10], bwt<-res[,11]. bwt[3] is kind o fin middle. probability that. typ eino console, bwt3vs2<-bwt3-bwt2. calculate diff between 2 columns. then mean(bwt3vs2>0). 0.267. chance that diff is positive is 27%.. if looking at lower than 0. 73% chance that wetlany 3 has fewer birds than wetand 2.

##95% posterior credible interval for the trend by working out the appropriate quantiles of the posterior (2.5th and 97.5th)
#using the quantile function
b1 <- res[, 2]
round(t(apply(b1, 2, quantile, probs = c(0.025, 0.975))), 2) #error so far,
#the transpose function t() simply makes the output match the format of output from the summarise function

b1 <- res[, 2] #quantile(res[,2],probs=c(0.025, 0.975)

#############
##############checking for zero-inflation

# observed proportion of observations that are zero
p0 <- mean(count > 0) #this will only check d3 since that is the last time it got defined, if this works, input into other weeks

# proportion of zeros we get in data simulated from model
p0sim <- apply(countsim > 0, 1, mean) #this doesn't work since countsim not successfully defined

# compare these (can discuss details if needed)
pvalue <- mean(p0sim > p0)
p0max <- max(hist(p0sim, plot = F)$counts)
hist(
  p0sim,
  main = paste0("Proportion zero (p = ", round(pvalue, 2), ")"),
  ylab = "",
  xlab = "",
  col = "white"
)
segments(p0, 0, p0, p0max, col = "red", lwd = 2)

# check for overdispersion relative to the model
#won't work until countsim is defined  #david will do another better version

nsim <- 10^3
ind <- sample(1:niters, nsim)
Dobs <- Dsim <- vector()
for (i in 1:nsim) {
  Dobs[i] <- mean((count - mu[ind[i], ])^2 / mu[ind[i], ])
  Dsim[i] <- mean((countsim[ind[i], ] - mu[ind[i], ])^2 / mu[ind[i], ])
}
par(mfrow = c(1, 2))
xylims <- range(Dsim, Dobs)
plot(Dsim, Dobs, xlim = xylims, ylim = xylims)
abline(0, 1)
hist(Dobs / Dsim, main = "")
summarise(Dobs / Dsim)
mean(Dobs / Dsim > 1)


##Find abundance and plot time vs abundance with CI
###########################################################

# equation 2 from Nichols: N=C/p  N is the estimate of true abundance, C is the count, p is the detection probability
# p = 0.98 from above calculations. C is the count data from historical data set, we want N, abundance
#resulting population estimates from this equation can be used to draw inferences about changes in abundance over time/space
#focus only on data from week 3, most representative of closed population

d <- read.csv("data_2024_chrono order.csv", header = TRUE) #input data
library(dplyr) #be able to organize data more
library(ggplot2)

week_3_data <- d %>%
  filter(week == 3) #filters out only rows of data from week 3, the most representative of the results

year3 <- week_3_data$year #pull out years from week 3
LEYEcount3 <- week_3_data$count #pull out count data from week 3 only
N3 <- (LEYEcount3) / 0.98 #0.98=prob of det. Abundance equation 2 from Nichols
Time <- year3
Abundance <- N3

#attempt 4
data <- data.frame(
  Time = year3,
  Abundance = N3
)

jags_data <- list(
  N = nrow(data),
  Time = data$Time,
  Abundance = data$Abundance
)

model_string <- "
model {
  for (i in 1:N) {
    Abundance[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * Time[i]
  }
  
  alpha ~ dnorm(0, 0.001)
  beta ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
}
"
# Write the model to a file
writeLines(model_string, con = "model.txt")

# Fit the model using JAGS
jags_model <- jags.model(
  "model.txt",
  data = jags_data,
  n.chains = 3,
  n.adapt = 1000
)

# Generate posterior samples
samples <- coda.samples(
  jags_model,
  variable.names = c("alpha", "beta", "tau"),
  n.iter = 2000
)


#falcon method
n <- length(year3)
mean.mu <- apply(mu, 2, mean)
low.mu <- apply(mu, 2, quantile, probs = 0.025) #then take quantiles of mu. 10,000 roles and n columns #type dim(mu) into console to confirm
upp.mu <- apply(mu, 2, quantile, probs = 0.975)
ylims <- range(N3, low.mu, upp.mu)
par(mfrow = c(1, 1))
plot(year3, N3, ylab = "LEYE count", xlab = "Year", ylim = ylims)
lines(year3, mean.mu) #posterior mean for mu                    #get crazy lines
lines(year3, low.mu, lty = 2) #quantiles
lines(year3, upp.mu, lty = 2)

####falcon example
# posterior mean and 95% credible bands for population mean
mean.mu <- apply(mu, 2, mean) #take mean of mu first
low.mu <- apply(mu, 2, quantile, probs = 0.025) #then take quantiles of mu. 10,000 roles and n columns #type dim(mu) into console to confirm
upp.mu <- apply(mu, 2, quantile, probs = 0.975) #these are 95% credible intervals for mu for each year
ylims <- range(fledglings, low.mu, upp.mu)
plot(year, fledglings, ylab = "fledglings", ylim = ylims)
lines(year, mean.mu) #posterior mean for mu
lines(year, low.mu, lty = 2) #quantiles
lines(year, upp.mu, lty = 2)


##Find abundance and plot time vs abundance with CI  ################ Nov 2025 attempt
###########################################################
###Email from David Fletcher 5/8/25:
## 1.	In the "trend datafile" the values for time surveyed (et) and distance travelled (ed) vary, which makes allowing for the detection probability
#within the trend code non-trivial. If you had more money available we could look into this, but what I suggest below would be quicker, less complicated,
#and potentially all you need.

## 2.For the "general public" plot you want, I would simply use a relative abundance axis, e.g. just start at 100 (or whatever total count you had in the
#first year, as a rough guide) and show a graph in which this relative abundance decreases at the estimated rate. This exponential trend line can be
#calculated as exp(mb0+mb1*year), where mb0 and mb1 are the posterior means or medians of b0 and b1, and year starts at 1 (as in the trend code).

## 3.For a scientific report or paper, you would report the estimate and 95% credible interval for exp(b1), plus evidence that the detection probability
#does not change from year to year (and is very close to 1 anyway). exp(b1) is the ratio of the relative abundance (after allowing for the effort and
#wetland covariates) from one year to the next, e.g. if exp(b1) = 0.98, we have a 2% decline per year, if exp(b1) = 1.05, we have a 5% increase per year.

trend_line <- exp(mb0 + mb1 * year) #mb0 and mb1 are the posterior means or medians of b0 and b1, and year starts at 1 (as in the trend code)
