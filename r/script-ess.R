#### Required packages ####

# # IMPORTANT:
# # Next line shoud be run once to install the latest development version from GitHub, enabling the use of HMC with dcar_leroux
# remotes::install_github("nimble-dev/nimble", subdir = "packages/nimble")

# install.packages("pacman", dep = TRUE)

pacman::p_load(sf, spdep, ggplot2, RColorBrewer, patchwork, 
               nimble, nimbleHMC, MCMCvis, scales, install = FALSE)

rm(list = ls())

# # Loading functions for calling Nimble
# source("https://raw.githubusercontent.com/MigueBeneito/pNimble/refs/heads/main/RutinasNimble.0.2.R")
# load.leroux()

#### Data loading ####

# ------------------------------------------------------------- #
# Data source: European Social Survey (ESS)
#
# To replicate the analysis, please download the dataset from:
# https://ess.sikt.no/en/datafile/242aaa39-3bbb-40f5-98bf-bfb1ce53d8ef
#
# Access to the ESS data requires a free user registration.
# Once downloaded, extract the .zip file and place the extracted
# CSV file in the "data" folder before running this script.
# ------------------------------------------------------------- #

# ESS data
ESSData <- read.csv(file.path("data", "ESS11e04_1.csv"))
# Germany data
ESSData <- ESSData[ESSData$cntry == "DE", ]

# Definition of the categorical covariates
# Age groups: 1 = (14, 34]; 2 = (34, 54]; 3 = (54, 100]
ESSData$agecut <- cut(ESSData$agea, c(14, 34, 54, 100))
levels(ESSData$agecut) <- c("15-34", "35-54", "55...")
# Education groups: 1 = (-1, 225]; 2 = (225, 500]; 3 = (500, 1000]
ESSData$educut <- cut(ESSData$edulvlb,c(-1, 225, 500, 1000))
levels(ESSData$educut) <- c("edulow", "edumid", "eduhigh")

# ESS cleaned
ESS <- ESSData[, c("health", "region", "gndr", "agecut", "educut", 
                   "dweight", "pspwght", "anweight", "cntry", 
                   "domain", "psu", "stratum", "prob")]
# Removing NA rows for age and education groups for simplicity
ESS <- ESS[!is.na(ESS$agecut), ]
ESS <- ESS[!is.na(ESS$educut), ]
str(ESS)

# Transforming gender and region into factors
ESS$gndr <- factor(ESS$gndr, labels = c("Man", "Woman"))
ESS$region <- factor(ESS$region, levels = c("DE1", "DE2", "DE3", "DE4", "DE5", 
                                            "DE6", "DE7", "DE8", "DE9", "DEA", 
                                            "DEB", "DEC", "DED", "DEE", "DEF", "DEG"))

# # Design weight
# weight <- ESS$dweight/mean(ESS$dweight)
# # Post-stratification weight including design weight
# weight <- ESS$pspwght/mean(ESS$pspwght)
# Analysis weight
weight <- ESS$anweight/mean(ESS$anweight)

# Gender of each respondent
gnd <- as.numeric(ESS$gndr)
# Age group of each respondent
age <- as.numeric(ESS$agecut)
# Education group of each respondent
edu <- as.numeric(ESS$educut)
# NUTS of each respondent
nuts <- as.numeric(ESS$region)

# How is your health in general? Would you say it is...
# 1 = Very good, 2 = Good, 3 = Fair, 4 = Bad, 5 = Very Bad
# 7 = Refusal, 8 = Don't know, 9 = No answer
levels(ESS$health)
table(ESS$health)
y <- as.numeric(ESS$health)
y[y > 5] <- NA

# survey: data frame containing the variables collected for each respondent (ESS)
# - nuts: small-area (NUTS1) from 1 to number of areas
# - gnd: 1 = Man, 2 = Woman
# - age: 1 = (14, 34]; 2 = (34, 54]; 3 = (54, 100]
# - edu: 1 = (-1, 225]; 2 = (225, 500]; 3 = (500, 1000]
# - weight: design/post-stratification/analysis weight
survey <- data.frame("y" = y, "nuts" = nuts, "gnd" = gnd, 
                     "age" = age, "edu" = edu, "weight" = weight)

#### Preparation of maps ####

# Loading Germany states GeoJSON
cartography <- st_read(file.path("data", "1_sehr_hoch.geo.json"))

# Neighborhood structure by contiguity
Neigh <- poly2nb(cartography)

# Adjacency matrix
W <- nb2mat(Neigh, style = "B")
# D - W matrix 
Q <- diag(apply(W, 1, sum)) - W
# Number of areas
NNUTS <- nrow(W)
# Number of neighbors of each area
nadj <- apply(W, 1, sum)
# Neighbors of each area
map <- unlist(apply(W, 1, function(x) which(x != 0)))
# Sum of all the neighbor numbers of all areas
nadj.tot <- length(map)
# Cumulative sums of the number of neighbors of each area
index <- c(0, cumsum(nadj))
# All the neighborhoods j ~ i where i < j
from.to <- cbind(rep(1:NNUTS, times = nadj), map); colnames(from.to) <- c("from", "to")
from.to <- from.to[which(from.to[, 1] < from.to[, 2]), ]
NDist <- nrow(from.to)
# Eigenvalues of D - W
Lambda <- eigen(Q)$values

#### Descriptive analysis ####

# Gender selection
Gender <- 1
Gender_Cat <- c("Man", "Woman")[Gender]

ordinal_survey <- survey[survey$gnd == Gender, ]
y <- as.numeric(ordinal_survey$y)
y[y > 5] <- NA
nuts_all <- factor(ordinal_survey$nuts, levels = 1:NNUTS)
table(nuts_all)
cartography$y_mean_ordinal <- as.numeric(by(y, nuts_all, mean, na.rm = TRUE))

bernoulli1_survey <- survey[survey$gnd == Gender, ]
y <- as.numeric(bernoulli1_survey$y)
y[y > 5] <- NA
y[y == 1 | y == 2] <- 0
y[y == 3 | y == 4 | y == 5] <- 1
nuts_all <- factor(bernoulli1_survey$nuts, levels = 1:NNUTS)
table(nuts_all)
cartography$y_mean_bernoulli1 <- as.numeric(by(y, nuts_all, mean, na.rm = TRUE))

bernoulli2_survey <- survey[survey$gnd == Gender, ]
y <- as.numeric(bernoulli2_survey$y)
y[y > 5] <- NA
y[y == 1 | y == 2 | y == 3] <- 0
y[y == 4 | y == 5] <- 1
nuts_all <- factor(bernoulli2_survey$nuts, levels = 1:NNUTS)
table(nuts_all)
cartography$y_mean_bernoulli2 <- as.numeric(by(y, nuts_all, mean, na.rm = TRUE))

# Ordinal
limit <- c(min(cartography$y_mean_ordinal, na.rm = TRUE),
           max(cartography$y_mean_ordinal, na.rm = TRUE))
p_y_mean_ordinal <- ggplot(cartography) + 
  geom_sf(aes(fill = y_mean_ordinal), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), 
                       limits = c(limit[1], limit[2]),
                       name = "Ordinal") +
  theme_void()

# Bernoulli (O1)
limit <- c(min(cartography$y_mean_bernoulli1, na.rm = TRUE),
           max(cartography$y_mean_bernoulli1, na.rm = TRUE))
p_y_mean_bernoulli1 <- ggplot(cartography) + 
  geom_sf(aes(fill = y_mean_bernoulli1), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), 
                       limits = c(limit[1], limit[2]),
                       name = "Bernoulli (D+)") +
  theme_void()

# Bernoulli (O2)
limit <- c(min(cartography$y_mean_bernoulli2, na.rm = TRUE),
           max(cartography$y_mean_bernoulli2, na.rm = TRUE))
p_y_mean_bernoulli2 <- ggplot(cartography) + 
  geom_sf(aes(fill = y_mean_bernoulli2), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), 
                       limits = c(limit[1], limit[2]),
                       name = "Bernoulli (D-)") +
  theme_void()

(p_y_mean_ordinal + p_y_mean_bernoulli1 + p_y_mean_bernoulli2)

ggsave(file.path("images", paste0("descriptive_", Gender_Cat, ".png")), 
       device = "png", width = 8, height = 6, dpi = 600)

#### Data preparation ####

# Gender selection
Gender <- 1
survey <- survey[survey$gnd == Gender, ]

# Number of respondents
NResp <- nrow(survey)
# Vector of zeros for the zero-trick
zero <- rep(0, NResp)

# NUTS of each respondent
nuts <- survey$nuts
# Age group of each respondent
age <- survey$age
# Education group of each respondent
edu <- survey$edu
# Sampling weight
weight <- survey$weight

# Response variable
y <- survey$y
# Number of levels
NCats <- length(table(y))
# Vector of ones for Dirichlet prior
ones <- rep(1, NCats)

# Number of levels of each (categorical) covariate
NAge <- length(table(survey$age))
NEdu <- length(table(survey$edu))

#### Population loading ####

# population_eu: four-dimensional array containing population counts by 
# nuts, gender, age and education group.
# The population data were obtained from Eurostat:
# https://ec.europa.eu/eurostat/databrowser/explore/all/popul?sort=category&lang=en&subtheme=cens.cens_21.cens_21dc&display=list
population_eu <- readRDS(file = file.path("data", "population-germany-eu.rds"))
population_eu <- population_eu[, Gender, , ]

# # population_wp: four-dimensional array containing population counts by
# # nuts, gender, age and education group. 
# # The population data were obtained from WorldPop:
# # https://hub.worldpop.org/geodata/summary?id=96779
# population_wp <- readRDS(file = file.path("data", "population-germany-wp.rds"))
# population_wp <- population_wp[, Gender, , ]

#### Ordinal model ####

### Model code ###

modelCode <- nimbleCode(
  {
    # Likelihood
    for (Resp in 1:NResp) {
      y[Resp] ~ dcat(prlevels[Resp, 1:NCats])
      
      # Definition of the probabilities of each category as a function of the
      # cumulative probabilities
      prlevels[Resp, 1] <- p.gamma[Resp, 1]
      for (Cat in 2:(NCats-1)) {
        prlevels[Resp, Cat] <- p.gamma[Resp, Cat] - p.gamma[Resp, Cat-1]
      }
      prlevels[Resp, NCats] <- 1 - p.gamma[Resp, NCats-1]
      
      # Linear predictor
      for (Cat in 1:(NCats-1)) {
        logit(p.gamma[Resp, Cat]) <- kappa[Cat] - 
          beta_age[age[Resp]] - beta_edu[edu[Resp]] - 
          sd.theta * theta[nuts[Resp]]
      }
    }
    
    # Prior distributions
    
    # kappa[1:(NCats-1)] cut points
    # Monotonic transformation
    for (Cat in 1:(NCats-1)) {
      kappa[Cat] <- logit(sum(delta[1:Cat]))
    }
    # delta[1:NCats] Dirichlet prior
    delta[1:NCats] ~ ddirch(ones[1:NCats])
    
    # beta_age[1:NAge] age group fixed effect (corner constraint)
    beta_age[1] <- 0
    for (AgeGroup in 2:NAge) {
      beta_age[AgeGroup] ~ dflat()
    }
    
    # beta_edu[1:NEdu] education group fixed effect (corner constraint)
    beta_edu[1] <- 0
    for (EduGroup in 2:NEdu) {
      beta_edu[EduGroup] ~ dflat()
    }
    
    # theta[1:NNUTS] spatial random effect
    # LCAR distribution
    theta[1:NNUTS] ~ dcar_leroux(rho = rho,
                                 sd.theta = 1,
                                 Lambda = Lambda[1:NNUTS],
                                 from.to = from.to[1:NDist, 1:2])
    
    # Hyperparameters of the spatial random effect
    rho ~ dunif(0, 1)
    sd.theta ~ dhalfflat()
    
    # Stochastic restrictions in order to avoid confounding problems
    # Required vectors
    for (Resp in 1:NResp) {
      theta.Resp[Resp] <- theta[nuts[Resp]]
    }
    
    # Weighted constraint
    # Zero-mean constraint for theta.Resp[1:NResp]
    zero.theta.resp ~ dnorm(mean.thetas.resp, 10000)
    mean.thetas.resp <- mean(theta.Resp[1:NResp])
    
  }
)

### Data to be loaded ###

modelData <- list(y = y,
                  zero.theta.resp = 0)

modelConstants <- list(age = age, edu = edu, nuts = nuts,
                       NResp = NResp, NCats = NCats, NAge = NAge,
                       NEdu = NEdu, ones = ones, NNUTS = NNUTS,
                       NDist = NDist, Lambda = Lambda, from.to = from.to)

### Parameters to be saved ###

modelParameters <- c("kappa", "beta_age", "beta_edu",
                     "theta", "sd.theta", "rho")

### Loading functions for calling Nimble ###

source("https://raw.githubusercontent.com/MigueBeneito/pNimble/refs/heads/main/RutinasNimble.0.2.R")
load.leroux()

# Inital values
modelInits <- local({
  constants <- modelConstants
  function() {
    library(extraDistr)
    NAge <- constants$NAge
    NEdu <- constants$NEdu; NNUTS <- constants$NNUTS
    NCats <- constants$NCats; ones <- constants$ones
    list(delta = as.numeric(rdirichlet(1, ones)),
         beta_age = c(NA, rnorm(NAge - 1)),
         beta_edu = c(NA, rnorm(NEdu - 1)),
         rho = runif(1), sd.theta = runif(1),
         theta = rnorm(NNUTS, sd = 0.1))
  }
})

# Number of chains to run in parallel
nchains <- 5
# pNimble call
salnimble <- pNimble(code = modelCode, data = modelData, constants = modelConstants, 
                     inits = modelInits, nchains = nchains, seeds = 1:nchains, 
                     niter = 2000, nburnin = 1000, thin = 5, 
                     summary = TRUE, WAIC = TRUE, monitors = modelParameters, 
                     # ntfyAccount = "MigueBeneito", 
                     HMC = TRUE, parallel = TRUE)

saveRDS(salnimble, file = file.path("results", "ordinal-results-srh-leroux-hmc-waic-man.rds"))

#### Bernoulli model ####

# Option 1 (D+): 0 = Very good and Good; 1 = Fair, Bad, and Very bad.
y[y > 5] <- NA
y[y == 1 | y == 2] <- 0
y[y == 3 | y == 4 | y == 5] <- 1
# # Option 2 (D-): 0 = Very good, Good, and Fair; 1 = Bad and Very bad.
# y[y > 5] <- NA
# y[y == 1 | y == 2 | y == 3] <- 0
# y[y == 4 | y == 5] <- 1

### Model code ###

modelCode <- nimbleCode(
  {
    # Likelihood
    for (Resp in 1:NResp) {
      y[Resp] ~ dbern(prsuccess[Resp])
      
      # Linear predictor  
      logit(prsuccess[Resp]) <- beta_0 + 
        beta_age[age[Resp]] + beta_edu[edu[Resp]] + 
        sd.theta * theta[nuts[Resp]]
    }
    
    # Prior distributions
    
    # beta_0 intercept
    beta_0 ~ dflat()
    
    # beta_age[1:NAge] age group fixed effect (corner constraint)
    beta_age[1] <- 0
    for (AgeGroup in 2:NAge) {
      beta_age[AgeGroup] ~ dflat()
    }
    
    # beta_edu[1:NEdu] education group fixed effect (corner constraint)
    beta_edu[1] <- 0
    for (EduGroup in 2:NEdu) {
      beta_edu[EduGroup] ~ dflat()
    }
    
    # theta[1:NNUTS] spatial random effect
    # LCAR distribution
    theta[1:NNUTS] ~ dcar_leroux(rho = rho,
                                 sd.theta = 1,
                                 Lambda = Lambda[1:NNUTS],
                                 from.to = from.to[1:NDist, 1:2])
    
    # Hyperparameters of the spatial random effect
    rho ~ dunif(0, 1)
    sd.theta ~ dhalfflat()
    
    # Stochastic restrictions in order to avoid confounding problems
    # Required vectors
    for (Resp in 1:NResp) {
      theta.Resp[Resp] <- theta[nuts[Resp]]
    }
    
    # Weighted constraint
    # Zero-mean constraint for theta.Resp[1:NResp]
    zero.theta.resp ~ dnorm(mean.thetas.resp, 10000)
    mean.thetas.resp <- mean(theta.Resp[1:NResp])
    
  }
)

### Data to be loaded ###

modelData <- list(y = y,
                  zero.theta.resp = 0)

modelConstants <- list(age = age, edu = edu, nuts = nuts,
                       NResp = NResp, NAge = NAge,
                       NEdu = NEdu, NNUTS = NNUTS,
                       NDist = NDist, Lambda = Lambda, from.to = from.to)

### Parameters to be saved ###

modelParameters <- c("beta_0", "beta_age", "beta_edu",
                     "theta", "sd.theta", "rho")

### Loading functions for calling Nimble ###

source("https://raw.githubusercontent.com/MigueBeneito/pNimble/refs/heads/main/RutinasNimble.0.2.R")
load.leroux()

# Inital values
modelInits <- local({
  constants <- modelConstants
  function() {
    NAge <- constants$NAge
    NEdu <- constants$NEdu; NNUTS <- constants$NNUTS
    list(beta_0 = rnorm(1),
         beta_age = c(NA, rnorm(NAge - 1)),
         beta_edu = c(NA, rnorm(NEdu - 1)),
         rho = runif(1), sd.theta = runif(1),
         theta = rnorm(NNUTS, sd = 0.1))
  }
})

# Number of chains to run in parallel
nchains <- 5
# pNimble call
salnimble <- pNimble(code = modelCode, data = modelData, constants = modelConstants, 
                     inits = modelInits, nchains = nchains, seeds = 1:nchains, 
                     niter = 2000, nburnin = 1000, thin = 5, 
                     summary = TRUE, WAIC = TRUE, monitors = modelParameters, 
                     # ntfyAccount = "MigueBeneito", 
                     HMC = TRUE, parallel = TRUE)

saveRDS(salnimble, file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-man-o1.rds"))

#### Model results ####

### Man ###

# Ordinal regression
ordinal_results <- readRDS(file = file.path("results", "ordinal-results-srh-leroux-hmc-waic-man.rds"))
# Logistic regression (Option 1): 0 = Very good and Good; 1 = Fair, Bad, and Very bad.
bernoulli_results1 <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-man-o1.rds"))
# Logistic regression (Option 2): 0 = Very good, Good, and Fair; 1 = Bad and Very bad.
bernoulli_results2 <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-man-o2.rds"))

### Woman ###

# Ordinal regression
ordinal_results <- readRDS(file = file.path("results", "ordinal-results-srh-leroux-hmc-waic-woman.rds"))
# Logistic regression (Option 1): 0 = Very good and Good; 1 = Fair, Bad, and Very bad.
bernoulli_results1 <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-woman-o1.rds"))
# Logistic regression (Option 2): 0 = Very good, Good, and Fair; 1 = Bad and Very bad.
bernoulli_results2 <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-woman-o2.rds"))

#### Convergence assessment ####

### Ordinal model ###

ordinal_results$summary

# c("kappa", "theta", "sd.theta", "rho")

MCMCsummary(object = ordinal_results$samples, params = "rho",
            # exact = TRUE,
            # ISB = FALSE,
            round = 4)

MCMCtrace(object = ordinal_results$samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE,
          n.eff = TRUE,
          params = "kappa")

test <- ordinal_results$samples

which((MCMCsummary(object = test, params = "kappa", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "kappa", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "beta_age", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "beta_age", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "beta_edu", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "beta_edu", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "theta", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "theta", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "sd.theta", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "sd.theta", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "rho", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "rho", round = 4)[, 7] < 400))

### Bernoulli model ###

bernoulli_results1$summary

# c("kappa", "theta", "sd.theta", "rho")

MCMCsummary(object = bernoulli_results1$samples, params = "rho",
            # exact = TRUE,
            # ISB = FALSE,
            round = 4)

MCMCtrace(object = bernoulli_results2$samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE,
          n.eff = TRUE,
          params = "beta_0")

test <- bernoulli_results1$samples

which((MCMCsummary(object = test, params = "beta_0", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "beta_0", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "beta_age", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "beta_age", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "beta_edu", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "beta_edu", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "theta", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "theta", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "sd.theta", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "sd.theta", round = 4)[, 7] < 400))
which((MCMCsummary(object = test, params = "rho", round = 4)[, 6]) > 1.02 | (MCMCsummary(object = test, params = "rho", round = 4)[, 7] < 400))

#### NimToWin function ####

nchains <- 5
nsims <- nchains * nrow(ordinal_results$samples[[1]])

# NimToWin: 
# - transforms NIMBLE output to WinBUGS sims.list output format

NimToWin <- function(salnimble) {
  
  kappa <- matrix(nrow = nsims, ncol = NCats - 1)
  beta_age <- matrix(nrow = nsims, ncol = NAge)
  beta_edu <- matrix(nrow = nsims, ncol = NEdu)
  theta <- matrix(nrow = nsims, ncol = NNUTS)
  sd.theta <- numeric(length = nsims)
  rho <- numeric(length = nsims)
  
  for (Cat in 1:(NCats - 1)) {
    kappa[, Cat] <- c(salnimble[[1]][,  paste0("kappa[", Cat, "]")],
                      salnimble[[2]][,  paste0("kappa[", Cat, "]")], 
                      salnimble[[3]][,  paste0("kappa[", Cat, "]")],
                      salnimble[[4]][,  paste0("kappa[", Cat, "]")], 
                      salnimble[[5]][,  paste0("kappa[", Cat, "]")])
  }
  
  for (AgeGroup in 1:NAge) {
    beta_age[, AgeGroup] <- c(salnimble[[1]][, paste0("beta_age[", AgeGroup, "]")], 
                              salnimble[[2]][, paste0("beta_age[", AgeGroup, "]")], 
                              salnimble[[3]][, paste0("beta_age[", AgeGroup, "]")],
                              salnimble[[4]][, paste0("beta_age[", AgeGroup, "]")],
                              salnimble[[5]][, paste0("beta_age[", AgeGroup, "]")])
  }
  
  for (EduGroup in 1:NEdu) {
    beta_edu[, EduGroup] <- c(salnimble[[1]][, paste0("beta_edu[", EduGroup, "]")], 
                              salnimble[[2]][, paste0("beta_edu[", EduGroup, "]")], 
                              salnimble[[3]][, paste0("beta_edu[", EduGroup, "]")],
                              salnimble[[4]][, paste0("beta_edu[", EduGroup, "]")],
                              salnimble[[5]][, paste0("beta_edu[", EduGroup, "]")])
  }
  
  for (NUTS in 1:NNUTS) {
    theta[, NUTS] <- c(salnimble[[1]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[2]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[3]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[4]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[5]][, paste0("theta[", NUTS, "]")])
  }
  
  sd.theta <- c(salnimble[[1]][, "sd.theta"], salnimble[[2]][, "sd.theta"], 
                salnimble[[3]][, "sd.theta"], salnimble[[4]][, "sd.theta"], 
                salnimble[[5]][, "sd.theta"])
  
  rho <- c(salnimble[[1]][, "rho"], salnimble[[2]][, "rho"], 
           salnimble[[3]][, "rho"], salnimble[[4]][, "rho"], 
           salnimble[[5]][, "rho"])
  
  summary <- MCMCsummary(object = salnimble, round = 4)
  sims.list <- list("kappa" = kappa, "beta_age" = beta_age, "beta_edu" = beta_edu,
                    "theta" = theta, "sd.theta" = sd.theta, "rho" = rho)
  
  salwinbugs <- list("summary" = summary, "sims.list" = sims.list,
                     "nchains" = nchains, "nsims" = nsims)
  
  return(salwinbugs)
}

ordinal_salwinbugs <- NimToWin(salnimble = ordinal_results$samples)

NimToWin <- function(salnimble) {
  
  beta_0 <- numeric(length = nsims)
  beta_age <- matrix(nrow = nsims, ncol = NAge)
  beta_edu <- matrix(nrow = nsims, ncol = NEdu)
  theta <- matrix(nrow = nsims, ncol = NNUTS)
  sd.theta <- numeric(length = nsims)
  rho <- numeric(length = nsims)
  
  beta_0 <- c(salnimble[[1]][, "beta_0"], salnimble[[2]][, "beta_0"], 
              salnimble[[3]][, "beta_0"], salnimble[[4]][, "beta_0"], 
              salnimble[[5]][, "beta_0"])
  
  for (AgeGroup in 1:NAge) {
    beta_age[, AgeGroup] <- c(salnimble[[1]][, paste0("beta_age[", AgeGroup, "]")], 
                              salnimble[[2]][, paste0("beta_age[", AgeGroup, "]")], 
                              salnimble[[3]][, paste0("beta_age[", AgeGroup, "]")],
                              salnimble[[4]][, paste0("beta_age[", AgeGroup, "]")],
                              salnimble[[5]][, paste0("beta_age[", AgeGroup, "]")])
  }
  
  for (EduGroup in 1:NEdu) {
    beta_edu[, EduGroup] <- c(salnimble[[1]][, paste0("beta_edu[", EduGroup, "]")], 
                              salnimble[[2]][, paste0("beta_edu[", EduGroup, "]")], 
                              salnimble[[3]][, paste0("beta_edu[", EduGroup, "]")],
                              salnimble[[4]][, paste0("beta_edu[", EduGroup, "]")],
                              salnimble[[5]][, paste0("beta_edu[", EduGroup, "]")])
  }
  
  for (NUTS in 1:NNUTS) {
    theta[, NUTS] <- c(salnimble[[1]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[2]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[3]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[4]][, paste0("theta[", NUTS, "]")], 
                       salnimble[[5]][, paste0("theta[", NUTS, "]")])
  }
  
  sd.theta <- c(salnimble[[1]][, "sd.theta"], salnimble[[2]][, "sd.theta"], 
                salnimble[[3]][, "sd.theta"], salnimble[[4]][, "sd.theta"], 
                salnimble[[5]][, "sd.theta"])
  
  rho <- c(salnimble[[1]][, "rho"], salnimble[[2]][, "rho"], 
           salnimble[[3]][, "rho"], salnimble[[4]][, "rho"], 
           salnimble[[5]][, "rho"])
  
  summary <- MCMCsummary(object = salnimble, round = 4)
  sims.list <- list("beta_0" = beta_0, "beta_age" = beta_age, "beta_edu" = beta_edu,
                    "theta" = theta, "sd.theta" = sd.theta, "rho" = rho)
  
  salwinbugs <- list("summary" = summary, "sims.list" = sims.list,
                     "nchains" = nchains, "nsims" = nsims)
  
  return(salwinbugs)
}

bernoulli1_salwinbugs <- NimToWin(salnimble = bernoulli_results1$samples)
bernoulli2_salwinbugs <- NimToWin(salnimble = bernoulli_results2$samples)

#### Fixed effects ####

df_plot_ord <- data.frame("param" = rownames(ordinal_results$summary)[1:6], 
                          "mean" = ordinal_results$summary[1:6, 1],
                          "lower" = ordinal_results$summary[1:6, 3],
                          "upper" = ordinal_results$summary[1:6, 5],
                          "model" = "Ordinal")
df_plot_bern1 <- data.frame("param" = rownames(bernoulli_results1$summary)[2:7], 
                            "mean" = bernoulli_results1$summary[2:7, 1],
                            "lower" = bernoulli_results1$summary[2:7, 3],
                            "upper" = bernoulli_results1$summary[2:7, 5],
                            "model" = "Bernoulli (D+)")
df_plot_bern2 <- data.frame("param" = rownames(bernoulli_results2$summary)[2:7], 
                            "mean" = bernoulli_results2$summary[2:7, 1],
                            "lower" = bernoulli_results2$summary[2:7, 3],
                            "upper" = bernoulli_results2$summary[2:7, 5],
                            "model" = "Bernoulli (D-)")
df_plot <- rbind(df_plot_ord, df_plot_bern1, df_plot_bern2)
df_plot$param <- factor(df_plot$param, levels = rev(df_plot$param[1:6]))
df_plot$model <- factor(df_plot$model, levels = c("Ordinal", "Bernoulli (D+)", "Bernoulli (D-)"))

ggplot(df_plot, aes(y = param, x = mean, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 position = position_dodge(width = 0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Ordinal" = "tomato", 
                                "Bernoulli (D+)" = "steelblue", 
                                "Bernoulli (D-)" = "seagreen")) +
  theme_minimal() +
  labs(x = "Posterior mean and 95% CI", y = "", color = "")

df_plot <- rbind(df_plot_ord, df_plot_bern1, df_plot_bern2)
df_plot$model <- factor(df_plot$model, levels = c("Ordinal", "Bernoulli (D+)", "Bernoulli (D-)"))

ggplot(df_plot, aes(y = param, x = mean, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 position = position_dodge(width = 0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Ordinal" = "tomato", 
                                "Bernoulli (D+)" = "steelblue", 
                                "Bernoulli (D-)" = "seagreen")) +
  theme_minimal() +
  labs(x = "Posterior mean and 95% CI", y = "", color = "") + coord_flip()

ggsave(file.path("images", paste0("fixed_", Gender_Cat, ".png")), 
       device = "png", width = 8, height = 6, dpi = 600)

#### Spatial effect ####

# Posterior mean of theta's
cartography$thetamean_ordinal <- ordinal_results$summary[startsWith(rownames(ordinal_results$summary), "theta"), 1]
cartography$thetamean_bernoulli1 <- bernoulli_results1$summary[startsWith(rownames(bernoulli_results1$summary), "theta"), 1]
cartography$thetamean_bernoulli2 <- bernoulli_results2$summary[startsWith(rownames(bernoulli_results2$summary), "theta"), 1]

# Posterior sd of theta's
cartography$thetasd_ordinal <- ordinal_results$summary[startsWith(rownames(ordinal_results$summary), "theta"), 2]
cartography$thetasd_bernoulli1 <- bernoulli_results1$summary[startsWith(rownames(bernoulli_results1$summary), "theta"), 2]
cartography$thetasd_bernoulli2 <- bernoulli_results2$summary[startsWith(rownames(bernoulli_results2$summary), "theta"), 2]

# Mean
limit <- max(abs(c(cartography$thetamean_ordinal, 
                   cartography$thetamean_bernoulli1,
                   cartography$thetamean_bernoulli2)), na.rm = TRUE)
p_thetamean_ordinal <- ggplot(cartography) + 
  geom_sf(aes(fill = thetamean_ordinal), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "BrBG")[9:1], 
                       limits = c(-limit, limit),     
                       values = rescale(c(-limit, 0, limit)),
                       name = NULL) + theme_void()
p_thetamean_bernoulli1 <- ggplot(cartography) + 
  geom_sf(aes(fill = thetamean_bernoulli1), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "BrBG")[9:1], 
                       limits = c(-limit, limit),     
                       values = rescale(c(-limit, 0, limit)),
                       name = NULL) + theme_void()
p_thetamean_bernoulli2 <- ggplot(cartography) + 
  geom_sf(aes(fill = thetamean_bernoulli2), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "BrBG")[9:1], 
                       limits = c(-limit, limit),     
                       values = rescale(c(-limit, 0, limit)),
                       name = NULL) + theme_void()

# Sd
limit <- c(min(c(cartography$thetasd_ordinal, 
                 cartography$thetasd_bernoulli1,
                 cartography$thetasd_bernoulli2), na.rm = TRUE),
           max(c(cartography$thetasd_ordinal, 
                 cartography$thetasd_bernoulli1,
                 cartography$thetasd_bernoulli2), na.rm = TRUE))
p_thetasd_ordinal <- ggplot(cartography) + 
  geom_sf(aes(fill = thetasd_ordinal), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), 
                       limits = c(limit[1], limit[2]), 
                       name = NULL) + theme_void()
p_thetasd_bernoulli1 <- ggplot(cartography) + 
  geom_sf(aes(fill = thetasd_bernoulli1), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), 
                       limits = c(limit[1], limit[2]), 
                       name = NULL) + theme_void()
p_thetasd_bernoulli2 <- ggplot(cartography) + 
  geom_sf(aes(fill = thetasd_bernoulli2), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), 
                       limits = c(limit[1], limit[2]), 
                       name = NULL) + theme_void()

p_thetamean_ordinal <- p_thetamean_ordinal + ggtitle("Ordinal")
p_thetamean_bernoulli1 <- p_thetamean_bernoulli1 + ggtitle("Bernoulli (D+)")
p_thetamean_bernoulli2 <- p_thetamean_bernoulli2 + ggtitle("Bernoulli (D-)")

p_thetamean_ordinal <- p_thetamean_ordinal + labs(tag = "Mean")
p_thetasd_ordinal <- p_thetasd_ordinal + labs(tag = "Sd")

tema_mapas <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                    plot.tag = element_text(face = "bold", size = 13, angle = 90),
                    plot.tag.position = c(-0.08, 0.5),
                    plot.margin = margin(5.5, 5.5, 5.5, 20))

final_plot <- wrap_plots(p_thetamean_ordinal + tema_mapas,
                         p_thetamean_bernoulli1 + tema_mapas,
                         p_thetamean_bernoulli2 + tema_mapas,
                         p_thetasd_ordinal + tema_mapas,
                         p_thetasd_bernoulli1 + tema_mapas,
                         p_thetasd_bernoulli2 + tema_mapas, ncol = 3)
final_plot

ggsave(file.path("images", paste0("spatialEffect_", Gender_Cat, ".png")), 
       device = "png", width = 8, height = 6, dpi = 600)

#### Post-stratification ####

# poststratify: 
# - (1) computes the nsims simulated probabilities for each age, education and nuts
# - (2) post-stratifies these probabilities for each nuts and category

poststratify <- function(salwinbugs) {
  
  p.gamma <- array(dim = c(nsims, NNUTS, NAge, NEdu, NCats - 1))
  prlevels <- array(dim = c(nsims, NNUTS, NAge, NEdu, NCats))
  prlevels_post <- array(dim = c(nsims, NNUTS, NCats))
  population <- population_eu
  
  # Probabilities for each age, education and nuts
  for (sim in 1:nsims) {
    for (AgeGroup in 1:NAge) {
      for (EduGroup in 1:NEdu) {
        for (NUTS in 1:NNUTS) {
          for (Cat in 1:(NCats - 1)) {
            p.gamma[sim, NUTS, AgeGroup, EduGroup, Cat] <- 
              ilogit(salwinbugs$sims.list$kappa[sim, Cat] - 
                       salwinbugs$sims.list$beta_age[sim, AgeGroup] - 
                       salwinbugs$sims.list$beta_edu[sim, EduGroup] - 
                       salwinbugs$sims.list$sd.theta[sim] * salwinbugs$sims.list$theta[sim, NUTS])
          }
          
          prlevels[sim, NUTS, AgeGroup, EduGroup, 1] <- p.gamma[sim, NUTS, AgeGroup, EduGroup, 1]
          prlevels[sim, NUTS, AgeGroup, EduGroup, NCats] <- 1 - p.gamma[sim, NUTS, AgeGroup, EduGroup, NCats - 1]
          
          for (Cat in 2:(NCats - 1)) {
            prlevels[sim, NUTS, AgeGroup, EduGroup, Cat] <- 
              p.gamma[sim, NUTS, AgeGroup, EduGroup, Cat] - p.gamma[sim, NUTS, AgeGroup, EduGroup, Cat-1]
          }
        }
      }
    }
    
    # Post-stratification
    for (NUTS in 1:NNUTS) {
      for (Cat in 1:NCats) {
        prlevels_post[sim, NUTS, Cat] <- sum(population[NUTS, , ])^(-1) * sum(prlevels[sim, NUTS, , , Cat] * population[NUTS, , ])
      }
    }
    if (sim %in% c(1, seq(nsims/nchains, nsims, nsims/nchains))) {
      cat(sim, "of", nsims, "simulations", "\n")
    } else {}
  }
  
  return(prlevels_post)
}

ordinal_post <- poststratify(salwinbugs = ordinal_salwinbugs)

poststratify <- function(salwinbugs) {
  
  prsuccess <- array(dim = c(nsims, NNUTS, NAge, NEdu))
  prsuccess_post <- matrix(nrow = nsims, ncol = NNUTS)
  population <- population_eu
  
  # Probabilities for each age, education and nuts
  for (sim in 1:nsims) {
    for (AgeGroup in 1:NAge) {
      for (EduGroup in 1:NEdu) {
        for (NUTS in 1:NNUTS) {
          prsuccess[sim, NUTS, AgeGroup, EduGroup] <- 
            ilogit(salwinbugs$sims.list$beta_0[sim] + 
                     salwinbugs$sims.list$beta_age[sim, AgeGroup] + 
                     salwinbugs$sims.list$beta_edu[sim, EduGroup] + 
                     salwinbugs$sims.list$sd.theta[sim] * salwinbugs$sims.list$theta[sim, NUTS])
        }
      }
    }
    
    # Post-stratification
    for (NUTS in 1:NNUTS) {
      prsuccess_post[sim, NUTS] <- sum(population[NUTS, , ])^(-1) * sum(prsuccess[sim, NUTS, , ] * population[NUTS, , ])
    }
    if (sim %in% c(1, seq(nsims/nchains, nsims, nsims/nchains))) {
      cat(sim, "of", nsims, "simulations", "\n")
    } else {}
  }
  
  return(prsuccess_post)
}

bernoulli1_post <- poststratify(salwinbugs = bernoulli1_salwinbugs)
bernoulli2_post <- poststratify(salwinbugs = bernoulli2_salwinbugs)

# Post-stratified posterior means and sd's of the percentages for each nuts and category

### Ordinal model ###

# Mean
cartography$percentage_mean1 <- apply(ordinal_post, 2:3, mean)[, 1] * 100
p_percentage_mean1 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean1), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = NULL) +
  theme_void()
cartography$percentage_mean2 <- apply(ordinal_post, 2:3, mean)[, 2] * 100
p_percentage_mean2 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean2), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = NULL) +
  theme_void()
cartography$percentage_mean3 <- apply(ordinal_post, 2:3, mean)[, 3] * 100
p_percentage_mean3 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean3), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = NULL) +
  theme_void()
cartography$percentage_mean4 <- apply(ordinal_post, 2:3, mean)[, 4] * 100
p_percentage_mean4 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean4), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = NULL) +
  theme_void()
cartography$percentage_mean5 <- apply(ordinal_post, 2:3, mean)[, 5] * 100
p_percentage_mean5 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean5), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = NULL) +
  theme_void()

# Sd
cartography$percentage_sd1 <- apply(ordinal_post, 2:3, sd)[, 1] * 100
p_percentage_sd1 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd1), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = NULL) +
  theme_void()
cartography$percentage_sd2 <- apply(ordinal_post, 2:3, sd)[, 2] * 100
p_percentage_sd2 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd2), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = NULL) +
  theme_void()
cartography$percentage_sd3 <- apply(ordinal_post, 2:3, sd)[, 3] * 100
p_percentage_sd3 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd3), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = NULL) +
  theme_void()
cartography$percentage_sd4 <- apply(ordinal_post, 2:3, sd)[, 4] * 100
p_percentage_sd4 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd4), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = NULL) +
  theme_void()
cartography$percentage_sd5 <- apply(ordinal_post, 2:3, sd)[, 5] * 100
p_percentage_sd5 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd5), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = NULL) +
  theme_void()

# CV
cartography$percentage_CV1 <- 100 * cartography$percentage_sd1/cartography$percentage_mean1
p_percentage_CV1 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV1), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = NULL) + 
  theme_void()
cartography$percentage_CV2 <- 100 * cartography$percentage_sd2/cartography$percentage_mean2
p_percentage_CV2 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV2), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = NULL) + 
  theme_void()
cartography$percentage_CV3 <- 100 * cartography$percentage_sd3/cartography$percentage_mean3
p_percentage_CV3 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV3), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = NULL) + 
  theme_void()
cartography$percentage_CV4 <- 100 * cartography$percentage_sd4/cartography$percentage_mean4
p_percentage_CV4 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV4), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = NULL) + 
  theme_void()
cartography$percentage_CV5 <- 100 * cartography$percentage_sd5/cartography$percentage_mean5
p_percentage_CV5 <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV5), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = NULL) + 
  theme_void()

p_percentage_mean1 <- p_percentage_mean1 + ggtitle("Very good")
p_percentage_mean2 <- p_percentage_mean2 + ggtitle("Good")
p_percentage_mean3 <- p_percentage_mean3 + ggtitle("Fair")
p_percentage_mean4 <- p_percentage_mean4 + ggtitle("Bad")
p_percentage_mean5 <- p_percentage_mean5 + ggtitle("Very bad")

p_percentage_mean1 <- p_percentage_mean1 + labs(tag = "Mean")
p_percentage_sd1   <- p_percentage_sd1 + labs(tag = "Sd")
p_percentage_CV1   <- p_percentage_CV1 + labs(tag = "CV (%)")

tema_mapas <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                    plot.tag = element_text(face = "bold", size = 12, angle = 90),
                    plot.tag.position = c(-0.08, 0.5),
                    plot.margin = margin(5.5, 5.5, 5.5, 20))

final_plot <- wrap_plots(p_percentage_mean1 + tema_mapas, 
                         p_percentage_mean2 + tema_mapas,
                         p_percentage_mean3 + tema_mapas, 
                         p_percentage_mean4 + tema_mapas,
                         p_percentage_mean5 + tema_mapas,
                         p_percentage_sd1 + tema_mapas, 
                         p_percentage_sd2 + tema_mapas,
                         p_percentage_sd3 + tema_mapas, 
                         p_percentage_sd4 + tema_mapas,
                         p_percentage_sd5 + tema_mapas,
                         p_percentage_CV1 + tema_mapas, 
                         p_percentage_CV2 + tema_mapas,
                         p_percentage_CV3 + tema_mapas, 
                         p_percentage_CV4 + tema_mapas, 
                         p_percentage_CV5 + tema_mapas, ncol = 5)
final_plot

ggsave(file.path("images", paste0("prevalenceOrdinal_", Gender_Cat, ".png")), 
       device = "png", width = 15, height = 10, dpi = 600)

### Bernoulli (O1) ###

# Mean
cartography$percentage_mean <- apply(bernoulli1_post, 2, mean) * 100
p_percentage_mean <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = "Mean") + 
  theme_void()

# Sd
cartography$percentage_sd <- apply(bernoulli1_post, 2,  sd) * 100
p_percentage_sd <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = "Sd") + 
  theme_void()

# CV
cartography$percentage_CV <- 100 * cartography$percentage_sd/cartography$percentage_mean
p_percentage_CV <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = "CV (%)") + 
  theme_void()

wrap_plots(p_percentage_mean, p_percentage_sd, p_percentage_CV, ncol = 3)

ggsave(file.path("images", paste0("prevalenceDplus_", Gender_Cat, ".png")), 
       device = "png", width = 8, height = 6, dpi = 600)

### Bernoulli (O2) ###

# Mean
cartography$percentage_mean <- apply(bernoulli2_post, 2, mean) * 100
p_percentage_mean <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_mean), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrBr"), name = "Mean") + 
  theme_void()

# Sd
cartography$percentage_sd <- apply(bernoulli2_post, 2,  sd) * 100
p_percentage_sd <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_sd), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Blues"), name = "Sd") + 
  theme_void()

# CV
cartography$percentage_CV <- 100 * cartography$percentage_sd/cartography$percentage_mean
p_percentage_CV <- ggplot(cartography) + 
  geom_sf(aes(fill = percentage_CV), color = "grey30", linewidth = 0.2) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Greens"), name = "CV (%)") + 
  theme_void()

wrap_plots(p_percentage_mean, p_percentage_sd, p_percentage_CV, ncol = 3)

ggsave(file.path("images", paste0("prevalenceDminus_", Gender_Cat, ".png")), 
       device = "png", width = 8, height = 6, dpi = 600)

#### Additional work: models with interaction between age and edu ####

### Ordinal model ###

### Model code ###

modelCode <- nimbleCode(
  {
    # Likelihood
    for (Resp in 1:NResp) {
      y[Resp] ~ dcat(prlevels[Resp, 1:NCats])
      
      # Definition of the probabilities of each category as a function of the
      # cumulative probabilities
      prlevels[Resp, 1] <- p.gamma[Resp, 1]
      for (Cat in 2:(NCats-1)) {
        prlevels[Resp, Cat] <- p.gamma[Resp, Cat] - p.gamma[Resp, Cat-1]
      }
      prlevels[Resp, NCats] <- 1 - p.gamma[Resp, NCats-1]
      
      # Linear predictor
      for (Cat in 1:(NCats-1)) {
        logit(p.gamma[Resp, Cat]) <- kappa[Cat] - 
          beta_age[age[Resp]] - beta_edu[edu[Resp]] - 
          beta_inter[age[Resp], edu[Resp]] -
          sd.theta * theta[nuts[Resp]]
      }
    }
    
    # Prior distributions
    
    # kappa[1:(NCats-1)] cut points
    # Monotonic transformation
    for (Cat in 1:(NCats-1)) {
      kappa[Cat] <- logit(sum(delta[1:Cat]))
    }
    # delta[1:NCats] Dirichlet prior
    delta[1:NCats] ~ ddirch(ones[1:NCats])
    
    # beta_age[1:NAge] age group fixed effect (corner constraint)
    beta_age[1] <- 0
    for (AgeGroup in 2:NAge) {
      beta_age[AgeGroup] ~ dflat()
    }
    
    # beta_edu[1:NEdu] education group fixed effect (corner constraint)
    beta_edu[1] <- 0
    for (EduGroup in 2:NEdu) {
      beta_edu[EduGroup] ~ dflat()
    }
    
    # beta_inter[1:NAge, 1:NEdu] interaction group fixed effect (corner constraint)
    beta_inter[1, 1] <- 0
    for (EduGroup in 2:NEdu) {
      beta_inter[1, EduGroup] <- 0
    }
    for (AgeGroup in 2:NAge) {
      beta_inter[AgeGroup, 1] <- 0
      for (EduGroup in 2:NEdu) {
        beta_inter[AgeGroup, EduGroup] ~ dflat()
      }
    }
    
    # theta[1:NNUTS] spatial random effect
    # LCAR distribution
    theta[1:NNUTS] ~ dcar_leroux(rho = rho,
                                 sd.theta = 1,
                                 Lambda = Lambda[1:NNUTS],
                                 from.to = from.to[1:NDist, 1:2])
    
    # Hyperparameters of the spatial random effect
    rho ~ dunif(0, 1)
    sd.theta ~ dhalfflat()
    
    # Stochastic restrictions in order to avoid confounding problems
    # Required vectors
    for (Resp in 1:NResp) {
      theta.Resp[Resp] <- theta[nuts[Resp]]
    }
    
    # Weighted constraint
    # Zero-mean constraint for theta.Resp[1:NResp]
    zero.theta.resp ~ dnorm(mean.thetas.resp, 10000)
    mean.thetas.resp <- mean(theta.Resp[1:NResp])
    
  }
)

### Data to be loaded ###

modelData <- list(y = y,
                  zero.theta.resp = 0)

modelConstants <- list(age = age, edu = edu, nuts = nuts,
                       NResp = NResp, NCats = NCats, NAge = NAge,
                       NEdu = NEdu, ones = ones, NNUTS = NNUTS,
                       NDist = NDist, Lambda = Lambda, from.to = from.to)

### Parameters to be saved ###

modelParameters <- c("kappa", "beta_age", "beta_edu", "beta_inter",
                     "theta", "sd.theta", "rho")

### Loading functions for calling Nimble ###

source("https://raw.githubusercontent.com/MigueBeneito/pNimble/refs/heads/main/RutinasNimble.0.2.R")
load.leroux()

# Inital values
modelInits <- local({
  constants <- modelConstants
  function() {
    library(extraDistr)
    NAge <- constants$NAge
    NEdu <- constants$NEdu; NNUTS <- constants$NNUTS
    NCats <- constants$NCats; ones <- constants$ones
    list(delta = as.numeric(rdirichlet(1, ones)),
         beta_age = c(NA, rnorm(NAge - 1)),
         beta_edu = c(NA, rnorm(NEdu - 1)),
         beta_inter = {
           mat <- matrix(NA, nrow = NAge, ncol = NEdu)
           mat[2:NAge, 2:NEdu] <- rnorm((NAge - 1) * (NEdu - 1), 0, 0.1)
           mat
         }, 
         rho = runif(1), sd.theta = runif(1),
         theta = rnorm(NNUTS, sd = 0.1))
  }
})

# Number of chains
nchains <- 5
salnimble <- pNimble(code = modelCode, data = modelData, constants = modelConstants, 
                     inits = modelInits, nchains = nchains, seeds = 1:nchains, 
                     niter = 2000, nburnin = 1000, thin = 5, 
                     summary = TRUE, WAIC = TRUE, monitors = modelParameters, 
                     # ntfyAccount = "MigueBeneito", 
                     HMC = TRUE, parallel = TRUE)

saveRDS(salnimble, file = file.path("results", "ordinal-results-srh-leroux-hmc-waic-man-inter.rds"))

### Bernoulli model ###

# Option 1: 0 = Very good and Good; 1 = Fair, Bad, and Very bad.
y[y > 5] <- NA
y[y == 1 | y == 2] <- 0
y[y == 3 | y == 4 | y == 5] <- 1
# # Option 2: 0 = Very good, Good, and Fair; 1 = Bad and Very bad.
# y[y > 5] <- NA
# y[y == 1 | y == 2 | y == 3] <- 0
# y[y == 4 | y == 5] <- 1

### Model code ###

modelCode <- nimbleCode(
  {
    # Likelihood
    for (Resp in 1:NResp) {
      y[Resp] ~ dbern(prsuccess[Resp])
      
      # Linear predictor  
      logit(prsuccess[Resp]) <- beta_0 + 
        beta_age[age[Resp]] + beta_edu[edu[Resp]] + 
        beta_inter[age[Resp], edu[Resp]] +
        sd.theta * theta[nuts[Resp]]
    }
    
    # Prior distributions
    
    # beta_0 intercept
    beta_0 ~ dflat()
    
    # beta_age[1:NAge] age group fixed effect (corner constraint)
    beta_age[1] <- 0
    for (AgeGroup in 2:NAge) {
      beta_age[AgeGroup] ~ dflat()
    }
    
    # beta_edu[1:NEdu] education group fixed effect (corner constraint)
    beta_edu[1] <- 0
    for (EduGroup in 2:NEdu) {
      beta_edu[EduGroup] ~ dflat()
    }
    
    # beta_inter[1:NAge, 1:NEdu] interaction group fixed effect (corner constraint)
    beta_inter[1, 1] <- 0
    for (EduGroup in 2:NEdu) {
      beta_inter[1, EduGroup] <- 0
    }
    for (AgeGroup in 2:NAge) {
      beta_inter[AgeGroup, 1] <- 0
      for (EduGroup in 2:NEdu) {
        beta_inter[AgeGroup, EduGroup] ~ dflat()
      }
    }
    
    # theta[1:NNUTS] spatial random effect
    # LCAR distribution
    theta[1:NNUTS] ~ dcar_leroux(rho = rho,
                                 sd.theta = 1,
                                 Lambda = Lambda[1:NNUTS],
                                 from.to = from.to[1:NDist, 1:2])
    
    # Hyperparameters of the spatial random effect
    rho ~ dunif(0, 1)
    sd.theta ~ dhalfflat()
    
    # Stochastic restrictions in order to avoid confounding problems
    # Required vectors
    for (Resp in 1:NResp) {
      theta.Resp[Resp] <- theta[nuts[Resp]]
    }
    
    # Weighted constraint
    # Zero-mean constraint for theta.Resp[1:NResp]
    zero.theta.resp ~ dnorm(mean.thetas.resp, 10000)
    mean.thetas.resp <- mean(theta.Resp[1:NResp])
    
  }
)

### Data to be loaded ###

modelData <- list(y = y,
                  zero.theta.resp = 0)

modelConstants <- list(age = age, edu = edu, nuts = nuts,
                       NResp = NResp, NAge = NAge,
                       NEdu = NEdu, NNUTS = NNUTS,
                       NDist = NDist, Lambda = Lambda, from.to = from.to)

### Parameters to be saved ###

modelParameters <- c("beta_0", "beta_age", "beta_edu", "beta_inter",
                     "theta", "sd.theta", "rho")

### Loading functions for calling Nimble ###

source("https://raw.githubusercontent.com/MigueBeneito/pNimble/refs/heads/main/RutinasNimble.0.2.R")
load.leroux()

# Inital values
modelInits <- local({
  constants <- modelConstants
  function() {
    NAge <- constants$NAge
    NEdu <- constants$NEdu; NNUTS <- constants$NNUTS
    list(beta_0 = rnorm(1),
         beta_age = c(NA, rnorm(NAge - 1)),
         beta_edu = c(NA, rnorm(NEdu - 1)),
         beta_inter = {
           mat <- matrix(NA, nrow = NAge, ncol = NEdu)
           mat[2:NAge, 2:NEdu] <- rnorm((NAge - 1) * (NEdu - 1), 0, 0.1)
           mat
         }, 
         rho = runif(1), sd.theta = runif(1),
         theta = rnorm(NNUTS, sd = 0.1))
  }
})

# Number of chains
nchains <- 5
salnimble <- pNimble(code = modelCode, data = modelData, constants = modelConstants, 
                     inits = modelInits, nchains = nchains, seeds = 1:nchains, 
                     niter = 2000, nburnin = 1000, thin = 5, 
                     summary = TRUE, WAIC = TRUE, monitors = modelParameters, 
                     # ntfyAccount = "MigueBeneito", 
                     HMC = TRUE, parallel = TRUE)

saveRDS(salnimble, file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-man-o1-inter.rds"))

#### Model results (interaction) ####

### Man ###

# Ordinal regression
ordinal_results_inter <- readRDS(file = file.path("results", "ordinal-results-srh-leroux-hmc-waic-man-inter.rds"))
# Logistic regression (Option 1): 0 = Very good and Good; 1 = Fair, Bad, and Very bad.
bernoulli_results1_inter <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-man-o1-inter.rds"))
# Logistic regression (Option 2): 0 = Very good, Good, and Fair; 1 = Bad and Very bad.
bernoulli_results2_inter <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-man-o2-inter.rds"))

### Woman ###

# Ordinal regression
ordinal_results_inter <- readRDS(file = file.path("results", "ordinal-results-srh-leroux-hmc-waic-woman-inter.rds"))
# Logistic regression (Option 1): 0 = Very good and Good; 1 = Fair, Bad, and Very bad.
bernoulli_results1_inter <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-woman-o1-inter.rds"))
# Logistic regression (Option 2): 0 = Very good, Good, and Fair; 1 = Bad and Very bad.
bernoulli_results2_inter <- readRDS(file = file.path("results", "bernoulli-results-srh-leroux-hmc-waic-woman-o2-inter.rds"))

#### Covariate effects (interaction) ####

df_plot_ord <- data.frame("param" = rownames(ordinal_results_inter$summary)[1:15], 
                          "mean" = ordinal_results_inter$summary[1:15, 1],
                          "lower" = ordinal_results_inter$summary[1:15, 3],
                          "upper" = ordinal_results_inter$summary[1:15, 5],
                          "model" = "Ordinal")
df_plot_bern1 <- data.frame("param" = rownames(bernoulli_results1_inter$summary)[2:16], 
                            "mean" = bernoulli_results1_inter$summary[2:16, 1],
                            "lower" = bernoulli_results1_inter$summary[2:16, 3],
                            "upper" = bernoulli_results1_inter$summary[2:16, 5],
                            "model" = "Bernoulli (D+)")
df_plot_bern2 <- data.frame("param" = rownames(bernoulli_results2_inter$summary)[2:16], 
                            "mean" = bernoulli_results2_inter$summary[2:16, 1],
                            "lower" = bernoulli_results2_inter$summary[2:16, 3],
                            "upper" = bernoulli_results2_inter$summary[2:16, 5],
                            "model" = "Bernoulli (D-)")
df_plot <- rbind(df_plot_ord, df_plot_bern1, df_plot_bern2)
df_plot$param <- factor(df_plot$param, levels = rev(df_plot$param[1:15]))
df_plot$model <- factor(df_plot$model, levels = c("Ordinal", "Bernoulli (D+)", "Bernoulli (D-)"))

ggplot(df_plot, aes(y = param, x = mean, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 position = position_dodge(width = 0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Ordinal" = "tomato", 
                                "Bernoulli (D+)" = "steelblue", 
                                "Bernoulli (D-)" = "seagreen")) +
  theme_minimal() +
  labs(x = "Posterior mean and 95% CI", y = "", color = "")

#### WAIC ####

c(ordinal_results$WAIC$WAIC, ordinal_results_inter$WAIC$WAIC)
c(bernoulli_results1$WAIC$WAIC, bernoulli_results1_inter$WAIC$WAIC)
c(bernoulli_results2$WAIC$WAIC, bernoulli_results2_inter$WAIC$WAIC)
