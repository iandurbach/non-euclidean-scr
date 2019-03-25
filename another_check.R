library(secr)
library(pals)
library(raster)
library(spatstat)
library(maptools)
library(gdistance) 
library(dplyr)

# Get some user functions for using noneuc in secr:
# ================================================
source("Design_simulator/noneuc-utils.R")

# Get the Tost, Noyon and Nemegt data
# ===================================
dat = readRDS("./Analysis4paper/TNN.Rds")
Tost = dat$Tost


# Consider Tost data
# ==================
ch = dat$Tost$capthist
mask = dat$Tost$mask
boundary = dat$Tost$boundary
cams = traps(ch)

## Fit a LC distance model with flat D:
## ===================================
## (startvals from previous fit)
startvals = list(D=exp(-9.4515904),lambda0=exp(-4.2505931),sigma=exp(8.6914951),noneuc=0.3314881)
sl.ne1 <-secr.fit(ch, detectfn="HHN", mask=mask,
                model=list(D~1, lambda0~1, sigma~1, noneuc~stdGC-1),
                details = list(userdist = geomLCdist),
                start=startvals,
                link = list(noneuc="identity"))

# this uses exp(stdGC) as covariate but then takes the geom mean directly
covariates(mask)$expstdGC <- exp(covariates(mask)$stdGC)
geommean = function(x) prod(x)^(1/length(x))

startvals = list(D=exp(-9.5227750),lambda0=exp(-4.4012107),sigma=exp(8.8538228))
sl.nuD2 <-secr.fit(ch, detectfn="HHN", mask=mask, start=startvals,
                   model=list(D~expstdGC, lambda0~1, sigma~1))
nuD2 = predictDsurface(sl.nuD2)
plotcovariate(nuD2,"D.0",col=parula(40))
region.N(sl.nuD2)

# this uses stdGC as covariate but uses an unsimplified expression for the geom mean
geommean = function(x) prod(exp(x))^(1/length(x))

startvals = list(D=exp(-9.5227750),lambda0=exp(-4.4012107),sigma=exp(8.8538228))
sl.nuD3 <-secr.fit(ch, detectfn="HHN", mask=mask, start=startvals,
                   model=list(D~stdGC, lambda0~1, sigma~1))
nuD3 = predictDsurface(sl.nuD3)
plotcovariate(nuD3,"D.0",col=parula(40))
region.N(sl.nuD3)

# D1 and D3 give the same Dsurface, D2 is flipped

## Fit varying density model with LC distance:
sl.ne.fit.D <-secr.fit(ch, detectfn="HHN", mask=mask,
                      model=list(D~expstdGC, lambda0~1, sigma~1, noneuc~stdGC-1),
                      details = list(userdist = geomLCdist),
                      link = list(noneuc="log"))
sl.noneuc.D = predictDsurface(sl.ne.fit.D)
plotcovariate(sl.noneuc.D,"D.0",col=parula(40))

geommean = function(x) mean(x)
sl.ne.fit.D2 <-secr.fit(ch, detectfn="HHN", mask=mask,
                       model=list(D~stdGC, lambda0~1, sigma~1, noneuc~stdGC-1),
                       details = list(userdist = geomLCdist),
                       link = list(noneuc="log"))
sl.noneuc.D2 = predictDsurface(sl.ne.fit.D2)
plotcovariate(sl.noneuc.D2,"D.0",col=parula(40))

sl.ne.fit.D3 <-secr.fit(ch, detectfn="HHN", mask=mask,
                        model=list(D~stdGC, lambda0~1, sigma~1, noneuc~expstdGC-1),
                        details = list(userdist = geomLCdist),
                        link = list(noneuc="identity"))
sl.noneuc.D3 = predictDsurface(sl.ne.fit.D3)
plotcovariate(sl.noneuc.D3,"D.0",col=parula(40))


## Fit a Euclidian distance model with varying D
## ==================================================
## (startvals from previous fit)
startvals = list(D=exp(-9.5227750),lambda0=exp(-4.4012107),sigma=exp(8.8538228))
sl.nuD <-secr.fit(ch, detectfn="HHN", mask=mask, start=startvals,
                  model=list(D~stdGC, lambda0~1, sigma~1),
                  link=list(D="log"))
nuD = predictDsurface(sl.nuD)
plotcovariate(nuD,"D.0",col=parula(40))
region.N(sl.nuD)

# this uses exp(stdGC) as covariate but then takes the geom mean directly
covariates(mask)$expstdGC <- exp(covariates(mask)$stdGC)
geommean = function(x) prod(x)^(1/length(x))

startvals = list(D=0.0002,lambda0=exp(-4.4012107),sigma=exp(8.8538228))
sl.nuD2 <-secr.fit(ch, detectfn="HHN", mask=mask, 
                  model=list(D~expstdGC, lambda0~1, sigma~1),
                  link=list(D="identity"))
nuD2 = predictDsurface(sl.nuD2)
plotcovariate(nuD2,"D.0",col=parula(40))
region.N(sl.nuD2)

# this uses stdGC as covariate but uses an unsimplified expression for the geom mean
geommean = function(x) prod(exp(x))^(1/length(x))

startvals = list(D=exp(-9.5227750),lambda0=exp(-4.4012107),sigma=exp(8.8538228))
sl.nuD3 <-secr.fit(ch, detectfn="HHN", mask=mask, start=startvals,
                   model=list(D~stdGC, lambda0~1, sigma~1))
nuD3 = predictDsurface(sl.nuD3)
plotcovariate(nuD3,"D.0",col=parula(40))
region.N(sl.nuD3)

# D1 and D3 give the same Dsurface, D2 is flipped







# check the outputs of costDistance and transition to see if that's where the difference is
# answer: its not

# simulate some activity centers
acs <- sim.popn(D = 1/4000, core = mask, model2D = "IHP", Ndist = "fixed", seed = 123)

Sraster1 <- raster(mask, "stdGC") 
trans1 <- transition(Sraster1, transitionFunction = function(x) exp(mean(x)), directions = 16)
trans1 <- geoCorrection(trans1)
cd1 <- costDistance(trans1, as.matrix(cams), as.matrix(acs))

Sraster2 <- raster(mask, "expstdGC") 
trans2 <- transition(Sraster2, transitionFunction = function(x) prod(x)^(1/length(x)), directions = 16)
trans2 <- geoCorrection(trans2)
cd2 <- costDistance(trans2, as.matrix(cams), as.matrix(acs))

trans3 <- transition(Sraster1, transitionFunction = function(x) prod(exp(x))^(1/length(x)), directions = 16)
trans3 <- geoCorrection(trans3)
cd3 <- costDistance(trans3, as.matrix(cams), as.matrix(acs))

cd_all <- data.frame(cd1 = as.vector(cd1), cd2 = as.vector(cd2), cd3 = as.vector(cd3))
cd_all <- cd_all %>% mutate(diff_cd12 = abs(cd1 - cd2),
                            diff_cd13 = abs(cd1 - cd3),
                            diff_cd23 = abs(cd2 - cd3))
apply(cd_all, 2, max)
