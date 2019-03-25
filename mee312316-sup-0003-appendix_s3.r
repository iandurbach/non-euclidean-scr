
################################################################################
# load the R workspace and libraries:
  wd <- getwd()
  load("connectivitySpace.Rdata") # from supplement zip file!!
  library("gdistance")
  library("raster")
  greyRamp <- colorRampPalette(c("gray50", "white"))( 100 )

# Which contains:
  ls()
  # aoiR: A raster with 'distance to stream' values in kms
  # e2dist: A function to calculate Euclidean diastance between 2 sets of points
  # gridTrapsXY: Trap locations
  # intlik4: Function to fit SCR model using Euclidean distance
  # intlik4ed: Function to fit SCR model estimating ecological distance


################################################################################
# function to simulate a single data set
  simCH.fn <- function(x){
    out <- NULL
    N  <- 200
    a0 <- -1
    sigma <- 1.4
    a1 <- 1/(2*sigma^2)
    K  <- 10
    a2vec <- seq(0,10,1)
    scen <- length(a2vec)
    r <- aoiR
    aoiXY <- as.matrix(coordinates(r))
    for(s in 1:scen){
# generate individuals
      acXY <- as.matrix(aoiXY[sample(1:nrow(aoiXY),N,replace=T),c("x","y")])
# generate cost surface
      a2 <- a2vec[s]
      trLayer <- NULL
      cost <- exp(a2 * r)
      trFn <- function(x)(1/(mean(x)))
      tr <- transition(cost,transitionFunction=trFn, direction = 8)
      trLayer <- geoCorrection(tr, scl = F)
# generate capture histories -- GRID
      nt <- nrow(gridTrapsXY)
      d <- costDistance(trLayer,acXY,gridTrapsXY)
      probcap <- plogis(a0) * exp(-a1 * d^2)
      Y <- matrix(NA,nrow=N,ncol=nt)
      for(i in 1:nrow(Y)){
        Y[i,] <- rbinom(nt, K, probcap[i,])
      }
      Y <- Y[apply(Y,1,sum)>0,]
# fit the models
      mink1 <- nlm(intlik4,p=c(a0,log(a1),log(N-nrow(Y))), hessian=F, 
                   y=Y, K=K,X=gridTrapsXY,G=aoiXY)
      mink2 <- nlm(intlik4ed,p=c(a0,log(a1),log(N-nrow(Y)),log(2)), hessian=F, 
                   y=Y, K=K,X=gridTrapsXY,G=aoiXY,cov=r,directions=8)
      out <- rbind(out,
                  c(mink1$estimate[1], mink1$estimate[2], mink1$estimate[3], 
                    NA, nrow(Y), exp(mink1$estimate[3])+nrow(Y),x,s,1),
                  c(mink2$estimate[1], mink2$estimate[2], mink2$estimate[3], 
                    mink2$estimate[4], nrow(Y), exp(mink2$estimate[3])+nrow(Y),x,s,2))
    }
    return(out)
  }
  plot(aoiR,col=greyRamp)
  points(gridTrapsXY,pch=16,col=2)
  # we used snowfall to parellize our processing so this looping version will:
  #  a) take ages
  #  b) differ very slightly due to different random numbers (monte carlo error)
  sim.out <- list()
  nsim <- 1
  set.seed(04091984)
  for(i in 1:nsim){
    sim.out[[i]] <- simCH.fn()
  }
  