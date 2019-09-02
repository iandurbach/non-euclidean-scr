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

dat = readRDS("tost_simulation/data/TNN.Rds")
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

# # plot the ACs generated from non-uniform D
# aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = stdGC)) + 
#   geom_point(data = simulated_points_Dcov, color = "red") +
#   geom_point(data = detectors_df, color = "green", size = 2, shape = 3) +
#   ggtitle("ACs with non-uniform density")

### find N recaps as function of lambda

# function that simulates many capture histories and takes

calc_recaps <- function(n_pts, b_ac, b_con, lambda0, sigma, n_occasions, mesh, detectors){
  
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

  ## number of unique animals detected
  n <- summary(ch)$counts$Total[4]
  
  ## number of recaptures
  r <- summary(ch)$counts$Total[6]
  
  res <- data.frame(n_pts, b_ac, b_con, lambda0, sigma, n_occasions, n, r)
  
  return(res)
}

# calc_recaps(n_pts = 15, b_ac = 1.5, b_con = 0.3, lambda0 = .26, sigma = 4000, 
#             n_occasions = 1, mesh = aoi_mesh, detectors = detectors)

pars <- expand.grid(b_ac = c(0.3,0.5,1,2), b_con = c(0.3,0.5), lambda0 = seq(0.5, 2, by = 0.5))

future::plan(multiprocess)

f_allres <- future_map_dfr(1:50, ~ future_pmap_dfr(pars, calc_recaps, n_pts = 20, sigma = 4000, n_occasions = 1, mesh = aoi_mesh, detectors = detectors),
                           .progress = TRUE)

f_allres %>% 
  ggplot(aes(x = lambda0, y = r)) + 
  facet_grid(b_ac ~ b_con) +
  stat_summary(fun.y = "mean", geom = "point") + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange") + 
  ylab("Number of recaps") + xlab("lambda0")

newdata <- f_allres %>% group_by(b_ac, b_con, lambda0) %>%
  summarize(mean_r = mean(r), mean_n = mean(n)) %>% ungroup() %>% data.frame()

mod0 <- lm(lambda0 ~ -1 + factor(b_ac)*factor(b_con)*mean_r, data = newdata)

# best model uses b_ac and b_con as factor variables (fits line to each combo) but 
# this means you can only predict with same values of b_ac and b_con (only mean_r changes)
summary(mod0)
plot(newdata$lambda0,predict(mod0))
predict(mod0, newdata = data.frame(b_ac = 0.5, b_con = 2, mean_r = 100))

# model with continuous b_ac and b_con, useful if want to predict with different values
# of these
mod0c <- lm(lambda0 ~ -1 + b_ac*b_con*mean_r, data = newdata)
plot(newdata$lambda0,predict(mod0c))
predict(mod0c, newdata = data.frame(b_ac = 0.3, b_con = 1, mean_r = 100))

save(mod0, mod0c, f_allres, file = "tost_simulation/output/Tost_Enrm_calcs.RData")
