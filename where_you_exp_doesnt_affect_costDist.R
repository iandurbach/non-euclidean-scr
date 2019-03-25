library(secr)
library(pals)
library(raster)
library(spatstat)
library(maptools)
library(gdistance) 
library(dplyr)

load("data/mee312316-sup-0002-appendix_s2.rdata") # from supplement zip file!!

aoiR <- raster::aggregate(aoiR, fact = 8, fun = mean)
aoi <- as.data.frame(aoiR) %>% rename(dist_stream = layer)
aoi_df <- cbind(as.data.frame(coordinates(aoiR)), aoi, D = 1)
aoi_mesh <- read.mask(data = aoi_df)
covariates(aoi_mesh)$exp_dist_stream <- exp(covariates(aoi_mesh)$dist_stream)

# simulate
# desired number of points to generate
n_pts <- 200  

# simulate some activity centers
D0_for_sim = n_pts / attr(aoi_mesh, "area") * (1 / nrow(aoi_df))
acs <- sim.popn(D = D0_for_sim, core = aoi_mesh, model2D = "IHP", Ndist = "fixed", seed = 123)

Sraster1 <- raster(aoi_mesh, "dist_stream") 
trans1 <- transition(Sraster1, transitionFunction = function(x) exp(mean(x)), directions = 16)
trans1 <- geoCorrection(trans1)
cd1 <- costDistance(trans1, as.matrix(gridTrapsXY), as.matrix(acs))

Sraster2 <- raster(aoi_mesh, "exp_dist_stream") 
trans2 <- transition(Sraster2, transitionFunction = function(x) prod(x)^(1/length(x)), directions = 16)
trans2 <- geoCorrection(trans2)
cd2 <- costDistance(trans2, as.matrix(gridTrapsXY), as.matrix(acs))

trans3 <- transition(Sraster1, transitionFunction = function(x) prod(exp(x))^(1/length(x)), directions = 16)
trans3 <- geoCorrection(trans3)
cd3 <- costDistance(trans3, as.matrix(gridTrapsXY), as.matrix(acs))

cd_all <- data.frame(cd1 = as.vector(cd1), cd2 = as.vector(cd2), cd3 = as.vector(cd3))
cd_all <- cd_all %>% mutate(diff_cd12 = abs(cd1 - cd2),
                            diff_cd13 = abs(cd1 - cd3),
                            diff_cd23 = abs(cd2 - cd3))
apply(cd_all, 2, max)