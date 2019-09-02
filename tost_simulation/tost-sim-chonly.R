library(tidyverse)
#library(devtools)
#install_github("rachelphillip/SCR-Book/scrmlebook")
library(scrmlebook)
library(secrdesign)
library(pals)
library(furrr)

# mean function for non-euc distance calcs
mymean = function(x) exp(mean(log(x))) 

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

# load Tost data
dat = readRDS("data/TNN.Rds")
Tost = dat$Tost

# load model to be used to estimate lambda0 needed to give desired n recaps
load("output/Tost_Enrm_calcs.RData")

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

## simulate a capture history with desired properties

sim_noneuc_ch_only <- function(n_pts, b_ac, b_con, mean_r, sigma, n_occasions, mesh, detectors, mod0){
  
  # get value of lambda0 that produces r recaps in expectation
  # comes from a previous model, see Enrm_for_noneuc_Tost.R
  lambda0 <- predict(mod0, newdata = data.frame(b_ac = b_ac, b_con = b_con, mean_r = mean_r))
  
  # desired number of points to generate
  n_pts <- n_pts
  
  ## generate non-uniform AC density
  # ac density ~ ruggedness, so add this as a mesh covariate.
  # introduce an b_ac parameter controlling strength of relationship
  b_ac <- b_ac 
  covariates(mesh)$Dac <- exp(b_ac * covariates(mesh)$stdGC)
  Dcov_for_sim = n_pts / attr(mesh, "area") * (covariates(mesh)$Dac / sum(covariates(mesh)$Dac)) 
  simulated_points_Dcov <- sim.popn(D = Dcov_for_sim, 
                                    core = mesh, 
                                    model2D = "IHP",
                                    Ndist = "fixed")
  
  ## generate non-Euclidean conductance
  # conductance ~ ruggedness, so add this as a mesh covariate.
  # create pixel-specific cost/friction and assign to the simulated popn objects 
  # introduce an b_con parameter controlling strength of relationship
  lambda0 <- lambda0
  sigma <- sigma
  b_con <- b_con 
  
  covariates(mesh)$noneuc <- exp(b_con * covariates(mesh)$stdGC)
  attr(simulated_points_Dcov, "mask") <- mesh
  
  ## simulate a capture history 
  ch <- sim.capthist(traps = detectors, pop = simulated_points_Dcov, 
                     userdist = myLCdist, 
                     noccasions = n_occasions,
                     detectpar = list(lambda0 = lambda0, sigma = sigma), 
                     detectfn = "HHN")
  
  return(ch)
}

chh <- sim_noneuc_ch_only(n_pts = 20, b_ac = 0.5, b_con = 2, mean_r = 1000, sigma = 4000, 
            n_occasions = 1, mesh = aoi_mesh, detectors = detectors, mod0 = mod0c)
