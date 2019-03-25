################################################################################
## Integrated likelihood function (intlik4 above) to calculate cost distance

## Required objects for running the model:
#
# 'start':  starting parameter values
# 'y':      individual x trap (row x col) counts of detections across K occasions
# 'K':      number of sampling occassions
# 'delta':




SCRed <- function(start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs,
                      cov=cov,G = NULL, ssbuffer = 2,directions=16, dist="ecol"){

    if (is.null(G)) {
        Xl <- min(X[, 1]) - ssbuffer
        Xu <- max(X[, 1]) + ssbuffer
        Yu <- max(X[, 2]) + ssbuffer
        Yl <- min(X[, 2]) - ssbuffer
        SSarea <- (Xu - Xl) * (Yu - Yl)
        if (is.null(K)) 
            return("need sample size")
        xg <- seq(Xl + delta/2, Xu - delta/2, delta)
        yg <- seq(Yl + delta/2, Yu - delta/2, delta)
        npix.x <- length(xg)
        npix.y <- length(yg)
        area <- (Xu - Xl) * (Yu - Yl)/((npix.x) * (npix.y))
        G <- cbind(rep(xg, npix.y), sort(rep(yg, npix.x)))
    }
    else {
        G <- G
        SSarea <- nrow(G)
    }
    nG <- nrow(G)
                                                                          
    ## cost distance -- This (#~#) differs from Euclidean distance model  #~#
    alpha2 <- exp(start[4])                                               #~#
    cost <- exp(alpha2 * cov)                                             #~#
    tr <- transition(cost, transitionFunction=function(x) (1/(mean(x))),  #~#
                     direction = directions)                              #~#
    trLayer <- geoCorrection(tr, scl = F)                                 #~#
    D <- costDistance(trLayer,as.matrix(X),as.matrix(G))                  #~#
    
    alpha0 <- start[1]
    alpha1 <- exp(start[2])
    n0 <- exp(start[3])
    probcap <- plogis(alpha0) * exp(-alpha1 * D * D)
    Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
    ymat <- y
    ymat <- rbind(y, rep(0, ncol(y)))
    lik.marg <- rep(NA, nrow(ymat))
    for (i in 1:nrow(ymat)) {
        Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, 
            nG), probcap[1:length(Pm)], log = TRUE))
        lik.cond <- exp(colSums(Pm))
        lik.marg[i] <- sum(lik.cond * (1/nG))
    }
    nv <- c(rep(1, length(lik.marg) - 1), n0)
    part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
    part2 <- sum(nv * log(lik.marg))
    out <- -1 * (part1 + part2)
    attr(out, "SSarea") <- SSarea
    out
}

