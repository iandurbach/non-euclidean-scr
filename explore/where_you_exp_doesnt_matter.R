library(tidyverse)
#library(devtools)
#install_github("rachelphillip/SCR-Book/scrmlebook")
library(scrmlebook)

# be sure to check you're using the right kind of LC distance calculations!!
source("noneuc-utils.R")

load("data/mee312316-sup-0002-appendix_s2.rdata") # from supplement zip file!!

aoiR <- raster::aggregate(aoiR, fact = 5, fun = mean)
aoi <- as.data.frame(aoiR) %>% rename(dist_stream = layer)
aoi_df <- cbind(as.data.frame(coordinates(aoiR)), aoi, D = 1)
aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = dist_stream))
aoi_mesh <- read.mask(data = aoi_df)

# turn the matrix of trap co-ordinates into a trap object
detectors_df <- as.data.frame(gridTrapsXY)
colnames(detectors_df) <- c("x", "y")
rownames(detectors_df) <- 1:nrow(detectors_df)
detectors <- read.traps(data = detectors_df, detector = "count", binary.usage = FALSE)

#### Simulate activity centers 

# these can be from a flat density or D ~ dist_to_stream

# desired number of points to generate
n_pts <- 200  

# using secr, flat density
D0_for_sim = n_pts / attr(aoi_mesh, "area") * (1 / nrow(aoi_df))
simulated_points_D0 <- sim.popn(D = D0_for_sim, 
                                core = aoi_mesh, 
                                model2D = "IHP",
                                Ndist = "fixed",
                                seed = 123)

# plot the ACs generated from D ~ dist_to_stream
aoi_df %>% ggplot(aes(x, y)) + geom_raster(aes(fill = dist_stream)) + 
  geom_point(data = simulated_points_D0, color = "red") +
  geom_point(data = detectors_df, color = "green", size = 2, shape = 3) +
  ggtitle("ACs closer to streams")

#### Generate capture histories

# parameter values from Sutherland (2014)
alpha0 <- -1 # implies lambda0 = invlogit(-1) = 0.2689414
sigma <- 1.4
alpha1 <- 1 / (2 * sigma^2)
alpha2 <- 5 # just one scenario from the range 0..10
K <- 2 # sampling over 10 occasions, collapsed to 1 occasion

# create pixel-specific cost/friction ~ distance to stream, and assign to the simulated popn objects 
covariates(aoi_mesh)$noneuc <- exp(alpha2 * covariates(aoi_mesh)$dist_stream)
attr(simulated_points_D0, "mask") <- aoi_mesh

## simulate capture histories from flat AC density with non-Euclidean distance function
ch_flatD_noneuc <- sim.capthist(detectors, pop = simulated_points_D0, userdist = arithLCdist, 
                                noccasions = K,
                                detectpar = list(lambda0 = invlogit(alpha0), sigma = sigma), 
                                detectfn = "HHN")

# fit model (arithmetic mean, not exponentiating in LCdist, log link)
mymean = function(x){1/mean(x)}
startvals = list(D = exp(8.0156), lambda0 = exp(-1.6518), sigma = exp(-0.2381), noneuc = 1)
mod_1 <-secr.fit(ch_flatD_noneuc, detectfn = "HHN", mask = aoi_mesh,
                 model= noneuc ~ dist_stream -1,
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))
coefficients(mod_1)

# fit model (arithmetic mean, exponentiating in myLCdist, identity link)
mymean = function(x){1/mean(exp(x))}
startvals = list(D = exp(8.0156), lambda0 = exp(-1.6518), sigma = exp(-0.2381), noneuc = exp(1))
mod_2 <-secr.fit(ch_flatD_noneuc, detectfn = "HHN", mask = aoi_mesh,
                 model= noneuc ~ dist_stream -1,
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "identity"))
coefficients(mod_2)

# fit model (geometric mean, not exponentiating in LCdist, log link)
mymean = function(x){1 / (prod(x)^(1/length(x)))}
startvals = list(D = exp(8.0156), lambda0 = exp(-1.6518), sigma = exp(-0.2381), noneuc = 1)
mod_3 <-secr.fit(ch_flatD_noneuc, detectfn = "HHN", mask = aoi_mesh,
                 model= noneuc ~ dist_stream -1,
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "log"))
coefficients(mod_3)

# fit model (geometric mean, exponentiating in myLCdist, identity link)
mymean = function(x){exp(-mean(x))}
startvals = list(D = exp(8.0156), lambda0 = exp(-1.6518), sigma = exp(-0.2381), noneuc = exp(1))
mod_4 <-secr.fit(ch_flatD_noneuc, detectfn = "HHN", mask = aoi_mesh,
                 model= noneuc ~ dist_stream -1,
                 details = list(userdist = myLCdist),
                 start = startvals,
                 link = list(noneuc = "identity"))
coefficients(mod_4)


