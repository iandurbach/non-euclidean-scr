library(tidyverse)
#library(devtools)
#install_github("rachelphillip/SCR-Book/scrmlebook")
library(scrmlebook)
library(pals)

# general LCdist function that allows you to call your own "mymean" fn (this is the only diff
# between the arithLCdist and geomLCdist in noneuc-utils.R)

myLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = mymean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

dat = readRDS("./Analysis4paper/TNN.Rds")
Tost = dat$Tost

# reduce by some factor for faster processing (leave at 1 for as is)
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

# recreate the flipping problem
ch_actual <- Tost$capthist
mymean = function(x){mean(exp(x))}
startvals = list(D = exp(-9.6371), lambda0 = exp(-4.3387), sigma = exp(8.8508), noneuc = 1)

mod_5_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                        model= list(D ~ stdGC, noneuc ~ stdGC -1),
                        details = list(userdist = myLCdist),
                        start = startvals,
                        link = list(noneuc = "identity"))

plotcovariate(predictDsurface(mod_5_actual),"D.0",col=parula(40))

# but CIs include zero
coefficients(mod_5_actual)

# rather use a log link and don't exponentiate inside LCdist
mymean = function(x){mean(x)}
mod_4_actual <-secr.fit(ch_actual, detectfn = "HHN", mask = aoi_mesh,
                        model= list(D ~ stdGC, noneuc ~ stdGC -1),
                        details = list(userdist = myLCdist),
                        start = startvals,
                        link = list(noneuc = "log"))

plotcovariate(predictDsurface(mod_4_actual),"D.0",col=parula(40))
coefficients(mod_4_actual) # still includes zero

# above is for arithmetic means, same for geometric means

###########################
# check if can reproduce for a simulated dataset with similar properties
###########################

# desired number of points to generate
n_pts <- 15

# non-uniform AC density
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

# plot the ACs generated from non-uniform D
aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = stdGC)) + 
  geom_point(data = simulated_points_Dcov, color = "red") +
  geom_point(data = detectors_df, color = "green", size = 2, shape = 3) +
  ggtitle("ACs with non-uniform density")

# some parameter values
alpha0 <- -1 # implies lambda0 = invlogit(-1) = 0.2689414
sigma <- 4000
alpha1 <- 1 / (2 * sigma^2)
alpha2 <- 0.3 
K <- 10 # sampling over 10 occasion, with 1 occasion variability is HIGH

# create pixel-specific cost/friction and assign to the simulated popn objects 
covariates(aoi_mesh)$noneuc <- exp(alpha2 * covariates(aoi_mesh)$stdGC)
attr(simulated_points_Dcov, "mask") <- aoi_mesh

# choose the function used to compute non-Euc distance (transitionFunction)
mymean <- function(x){mean(x)}
#mymean <- function(x){(prod(x)^(1/length(x)))}

# simulate capture history for non-uniform AC density, non-Euclidean distance
ch_Dcov_noneuc <- sim.capthist(detectors, pop = simulated_points_Dcov, userdist = myLCdist, 
                               noccasions = K,
                               detectpar = list(lambda0 = invlogit(alpha0), sigma = sigma), 
                               detectfn = "HHN")

# fit a constant density Euc model to get starting values
mod_0 <-secr.fit(ch_Dcov_noneuc, detectfn = "HHN", mask = aoi_mesh, model = D ~ 1)
# set starting values 
startvals = list(D = exp(-9.536054), lambda0 = exp(-1.339195), sigma = exp(8.704232), noneuc = 1)

# model 4: non-uniform AC density, non-Euclidean distance
mod_4 <-secr.fit(ch_Dcov_noneuc, detectfn = "HHN", mask = aoi_mesh,
                 model= list(D ~ stdGC, noneuc ~ stdGC -1),
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))
coefficients(mod_4)

# model 5: as for model 4 but exponentiate inside LCdist and use identity link
mymean = function(x) mean(exp(x))
mod_5 <-secr.fit(ch_Dcov_noneuc, detectfn = "HHN", mask = aoi_mesh,
                  model= list(D ~ stdGC, noneuc ~ stdGC -1),
                  details = list(userdist = myLCdist),
                  start = startvals,
                  link = list(noneuc = "identity"))
# reset to what it was
mymean <- function(x){mean(x)}
#mymean <- function(x){(prod(x)^(1/length(x)))}

# seems like things are consistent and correct simulated values are returned
plotcovariate(predictDsurface(mod_4),"D.0",col=parula(40))
plotcovariate(predictDsurface(mod_5),"D.0",col=parula(40))

coefficients(mod_4) 
coefficients(mod_5) 
