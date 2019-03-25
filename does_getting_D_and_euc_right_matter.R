library(tidyverse)
#library(devtools)
#install_github("rachelphillip/SCR-Book/scrmlebook")
library(scrmlebook)
library(pals)

# be sure to check you're using the right kind of LC distance calculations!!
source("noneuc-utils.R")

dat = readRDS("data/TNN.Rds")
Tost = dat$Tost

aoiR <- secr::raster(Tost$mask, "stdGC")
aoiR <- raster::aggregate(aoiR, fact = 3, fun = mean)
aoi <- as.data.frame(aoiR) %>% rename(stdGC = layer)
aoi_df <- cbind(as.data.frame(coordinates(aoiR)), aoi, D = 1) %>% filter(!is.nan(stdGC),!is.na(stdGC))
aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = stdGC))
aoi_mesh <- read.mask(data = aoi_df)

# turn the matrix of trap co-ordinates into a trap object
detectors_df <- as.data.frame(traps(Tost$capthist))
colnames(detectors_df) <- c("x", "y")
rownames(detectors_df) <- rownames(traps(Tost$capthist))
detectors <- read.traps(data = detectors_df, detector = "count", binary.usage = FALSE)
#detectors <- read.traps(data = traps(Tost$capthist), detector = "count", binary.usage = FALSE)

#######################################
#### Simulate activity centres
#######################################

# desired number of points to generate
n_pts <- 500

# uniform AC density 
D0_for_sim = n_pts / attr(aoi_mesh, "area") * (1 / nrow(aoi_df))
simulated_points_D0 <- sim.popn(D = D0_for_sim, 
                                core = aoi_mesh, 
                                model2D = "IHP",
                                Ndist = "fixed",
                                seed = 123)

# non-uniform AC density (fiddle n_pts a bit to get exactly 200, not sure why)
# ac density is potentially a function of distance to stream, so add this as a mesh covariate.
# introduce an alpha3 parameter controlling strength of relationship
alpha3 <- 1.5 
covariates(aoi_mesh)$Dac <- exp(alpha3 * covariates(aoi_mesh)$stdGC)
Dcov_for_sim = n_pts / attr(aoi_mesh, "area") * (covariates(aoi_mesh)$Dac / sum(covariates(aoi_mesh)$Dac)) 
simulated_points_Dcov <- sim.popn(D = Dcov_for_sim, 
                                  core = aoi_mesh, 
                                  model2D = "IHP",
                                  Ndist = "fixed",
                                  seed = 123)

# plot the ACs generated from uniform density
aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = stdGC)) + 
  geom_point(data = simulated_points_D0, color = "red") +
  geom_point(data = detectors_df, color = "green", size = 2, shape = 3) +
  ggtitle("Flat density")

# plot the ACs generated from D ~ dist_to_stream
aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = stdGC)) + 
  geom_point(data = simulated_points_Dcov, color = "red") +
  geom_point(data = detectors_df, color = "green", size = 2, shape = 3) +
  ggtitle("ACs closer to streams")

#######################################
#### Generate capture histories
#######################################

# there are 4 kinds of simulated capture histories: D~0 / D~cov x euclidean / noneuclidean

# parameter values from Sutherland (2014)
alpha0 <- -1 # implies lambda0 = invlogit(-1) = 0.2689414
sigma <- 1700
alpha1 <- 1 / (2 * sigma^2)
alpha2 <- 0.3 # true non-euc parameter -- cost/friction ~ exp(alpha2 * cov) (see below)
K <- 10 # sampling over 10 occasions, collapsed to 1 occasion

# create pixel-specific cost/friction ~ distance to stream, and assign to the simulated popn objects 
covariates(aoi_mesh)$noneuc <- exp(alpha2 * covariates(aoi_mesh)$stdGC)
attr(simulated_points_D0, "mask") <- aoi_mesh
attr(simulated_points_Dcov, "mask") <- aoi_mesh

# choose the function used to compute the path weight (non-Euc distance) between two points from 
# the covariate values i.e. path_weight([x1,y1],[x2,y2]) = mymean(cov([x1,y1]), cov[x2,y2])
# this is where you specify arithmetic or geometric mean, and 1-upon if inverse reln between covariate
# and conductance (higher distance from river = LESS conductance / MORE resistance)
mymean <- function(x){mean(x)}
#mymean <- function(x){(prod(x)^(1/length(x)))}

# uniform AC density, Euclidean distance
ch_flatD_euc <- sim.capthist(detectors, pop = simulated_points_D0, 
                             noccasions = K,
                             detectpar = list(lambda0 = invlogit(alpha0), sigma = sigma), 
                             detectfn = "HHN")

# non-uniform AC density, Euclidean distance
ch_Dcov_euc <- sim.capthist(detectors, pop = simulated_points_Dcov, 
                            noccasions = K,
                            detectpar = list(lambda0 = invlogit(alpha0), sigma = sigma), 
                            detectfn = "HHN")

# uniform AC density, non-Euclidean distance
ch_flatD_noneuc <- sim.capthist(detectors, pop = simulated_points_D0, userdist = myLCdist, 
                                noccasions = K,
                                detectpar = list(lambda0 = invlogit(alpha0), sigma = sigma), 
                                detectfn = "HHN")

# non-uniform AC density, non-Euclidean distance
ch_Dcov_noneuc <- sim.capthist(detectors, pop = simulated_points_Dcov, userdist = myLCdist, 
                               noccasions = K,
                               detectpar = list(lambda0 = invlogit(alpha0), sigma = sigma), 
                               detectfn = "HHN")

#######################################
#### Fit models
#######################################

# there are four kinds of models: D~0 / D~cov x euclidean / noneuclidean
# each of these models can be fitted to any of the four kinds of capture histories
# Sutherland et al. (2014) simulated one kind of CH (D~0, non-euc) and two kinds of 
# model (D~0, euc & D~0, non-euc).

# set the capture history you want to use here
ch <- ch_Dcov_noneuc  

# set starting values 
startvals = list(D = exp(-5.8952), lambda0 = exp(-1.3335), sigma = exp(7.9654), noneuc = 1)

# model 1: uniform AC density, Euclidean distance, use to get startvals
mod_1 <-secr.fit(ch, detectfn = "HHN", mask = aoi_mesh,
                 model = list(D ~ 1))
coefficients(mod_1)

# model 2: non-uniform AC density, Euclidean distance
mod_2 <-secr.fit(ch, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ stdGC),
                 start = startvals)
coefficients(mod_2)

# model 3: uniform AC density, non-Euclidean distance
mod_3 <-secr.fit(ch, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ 1, noneuc ~ stdGC -1),
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))
coefficients(mod_3)

# model 4: non-uniform AC density, non-Euclidean distance
mod_4 <-secr.fit(ch, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ stdGC, noneuc ~ stdGC -1),
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))
coefficients(mod_4)

# model 5: as for model 4 but exponentiate inside LCdist and use identity link
mymean = function(x) mean(exp(x))
mod_5l <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                        model= list(D ~ stdGC, noneuc ~ stdGC -1),
                        details = list(userdist = myLCdist),
                        start = startvals,
                        link = list(noneuc = "identity"))
# reset to what it was
mymean <- function(x){mean(x)}
#mymean <- function(x){(prod(x)^(1/length(x)))}

#######################################
#### With actual capture history
#######################################

ch_actual <- Tost$capthist

# model 1: uniform AC density, Euclidean distance, use to get startvals
mod_1_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                 model = list(D ~ 1))

startvals = list(D = exp(-9.6371), lambda0 = exp(-4.3387), sigma = exp(8.8508), noneuc = 1)

# model 2: non-uniform AC density, Euclidean distance
mod_2_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ stdGC),
                 start = startvals)
coefficients(mod_2)

# model 3: uniform AC density, non-Euclidean distance
mod_3_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ 1, noneuc ~ stdGC -1),
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))
coefficients(mod_3)

# model 4: non-uniform AC density, non-Euclidean distance
mod_4_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ stdGC, noneuc ~ stdGC -1),
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))

# model 5: as for model 4 but exponentiate inside LCdist and use identity link
mymean = function(x) mean(exp(x))
mod_5_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                        model= list(D ~ stdGC, noneuc ~ stdGC -1),
                        details = list(userdist = myLCdist),
                        start = startvals,
                        link = list(noneuc = "identity"))

# review results
coefficients(mod_1_actual)
coefficients(mod_2_actual)
coefficients(mod_3_actual)
coefficients(mod_4_actual)
coefficients(mod_5_actual)

region.N(mod_1_actual)
region.N(mod_2_actual)
region.N(mod_3_actual)
region.N(mod_4_actual)
region.N(mod_5_actual)

plotcovariate(predictDsurface(mod_2_actual),"D.0",col=parula(40))
plotcovariate(predictDsurface(mod_4_actual),"D.0",col=parula(40))
plotcovariate(predictDsurface(mod_5_actual),"D.0",col=parula(40))

