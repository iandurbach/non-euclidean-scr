# general LCdist function that allows you to call your own "mymean" fn (this is the only diff
# between the arithLCdist and geomLCdist below)

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



# geommean = function(x){exp(mean(x))}   ## if conductance increases with covariate
geommean = function(x){exp(-mean(x))}    ## if conductance decreases with covariate
#geommean = function(x){mean(x)}    ## if you already exponentiate outside the fn and con PT cov

geomLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = geommean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

#arithmean = function(x){mean(exp(x))}  ## if conductance increases with covariate
#arithmean = function(x){1/mean(exp(x))}  ## if conductance decreases with covariate
arithmean = function(x){1/mean(x)}      ## if you already exponentiate outside the fn and con PT cov

arithLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = arithmean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

expmax = function(x) exp(max(x))

expmaxLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = expmax, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

lcusageplot = function(fit,n=512,mask=NULL,base="noneuc.0",lcdfun="geomLCdist",maskcol=parula(40),...) {
  require(pals)
  
  if(is.null(mask)) mask = fit$mask
  
  lcd.fun = match.fun(lcdfun)
  
  if(!is.element("noneuc.0",names(covariates(mask))))
    stop("Must have 'noneuc.0' as one of the mask covariates. It is not there.")
  if(!is.element(base,names(covariates(mask)))) {
    warning(paste("mask does not have a covariate called ",base,"; noneuc.0 being used instead."))
  }
  covariates(mask)$base = covariates(mask)$noneuc.0
  
  for(i in 1:n){
    plotcovariate(mask,"base",col=maskcol,what="image")
    fromind = nearesttrap(unlist(locator(1)), mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    from=matrix(from,ncol=2)
    dists = lcd.fun(from,mask,mask)
    dfn <- secr:::getdfn(fit$detectfn)
    p = dfn(dists[1,],unlist(detectpar(fit)))
    covariates(mask)$p = p/sum(p)
    plotcovariate(mask,"p",what="image",...)
    points(from,col="white",pch=19)
    points(from,cex=0.8,pch=19)
    waitclick = unlist(locator(1))
  }
  
  invisible(mask)
}

lcpathplot = function(mask,transitionFunction,type="noneuc",n=512,background="d2.river",linecol="white", 
                      directions=16, symm=TRUE,directed=FALSE,lwd=1,...) {
  
  if(!is.element("noneuc.0",names(covariates(mask))))
    stop("Must have 'noneuc.0' as one of the mask covariates. It is not there.")
  if(!is.element(type,c("noneuc","sigma")))
    stop(paste("Invalid type: '",type,"'"," It must be `noneuc` or `sigma`.",sep=""))
  rastermask = raster(mask,"noneuc.0") # make raster with covariates(mask)$noneuc.0 as values of pixels
  
  transfun=match.fun(transitionFunction)
  
  coords = coordinates(rastermask) # lookup table for vertex coordinates
  # secr models conductance (sigma) with a log link. Line below assumes that $noneuc.0 is on linear predictor scale 
  # tr1<-transition(rastermask,transitionFunction=function(x) exp(lp(x)),directions=directions,symm=symm)
  tr1<-transition(rastermask,transitionFunction=transfun,directions=directions,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  
  if(type=="noneuc") plotcovariate(mask,"noneuc.0",...)
  else {
    covariates(mask)$sigma = exp(covariates(mask)$noneuc.0)
    plotcovariate(mask,"sigma",...)
  }
  
  dists = rep(NA,n) # to keep least-cost path distances in
  
  for(i in 1:n) {
    fromind = nearesttrap(unlist(locator(1)), mask)
    toind = nearesttrap(unlist(locator(1)), mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    to = c(x=mask$x[toind],y=mask$y[toind])
    from=matrix(from,ncol=2)
    to=matrix(to,ncol=2)
    #    npts=dim(from)[1]
    #    nptsto=dim(to)[1]
    #    if(nptsto != npts) stop("Must have same number of points in 'from' as in 'to'.")
    #    if(npts>1) pts = closest_coords(from[i,],to[i,],rastermask)
    #    else pts = closest_coords(from,to,rastermask)
    pts = closest_coords(from,to,rastermask)
    vpts = get_vertices(pts,rastermask)
    
    trmat=summary(transitionMatrix(tr1CorrC))
    #cbind(trmat,1/trmat$x)
    #    rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
    rel=data.frame(from=trmat$i,to=trmat$j,weight=trmat$x)
    #rel
    g = graph_from_data_frame(rel,directed=directed,vertices=NULL)
    #    attributes(g)$noneuc.0=1/trmat$x
    #    E(g)$weight=1/trmat$x
    attributes(g)$noneuc.0=trmat$x
    E(g)$weight=trmat$x
    #vertices = as_data_frame(g, what="vertices")
    #edges = as_data_frame(g, what="edges")
    svert=which(names(V(g))==vpts[1])
    evert=which(names(V(g))==vpts[2])
    # NB: Need to invert E(g) so that higher values lead to shorter distances:
    spath=as.numeric(names(shortest_paths(g,from=svert,to=evert,weights=1/E(g)$weight)$vpath[[1]]))
    dists[i]=igraph:::distances(g,v=svert,to=evert,weights=attributes(g)$noneuc.0)
    
    nppts=length(spath)
    segments(coords[spath[-nppts],1],coords[spath[-nppts],2],coords[spath[-1],1],coords[spath[-1],2],col=linecol,lwd=lwd)
    points(coords[spath[c(1,nppts)],],pch=19,col="white",cex=1.5)
    points(coords[spath[c(1,nppts)],],pch=19,col=c("green","red"),cex=0.75)
  }
  
  invisible(dists)
}

closest_coords=function(from,to,rastermask){
  ends=SpatialPoints(rbind(from,to))
  grid=as(rastermask, 'SpatialGrid') 
  xy=over(ends,grid)
  return(coordinates(grid)[xy,])
}




plotcovariate = function(mask,covariate, ...) {
  require(sp)
  
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  pts = cbind(mask$x,mask$y)
  spdfdat = data.frame(covariate=covariates(mask)[[cnum]])
  spdf = SpatialPixelsDataFrame(pts, spdfdat)
  plot(spdf, ...)
  
  invisible(spdf)
}


#' @title Finds closest coordinates on raster to two points
#'
#' @description
#'  Uses function over() from package sp to overlay points on raster and return closest raster coordinates
#'  
#' @param from pair of coordinates (x,y) from which to start
#' @param to pair of coordinates (x,y) to get to
#' @param rastermask Raster object (typically created from mask by something like 
#' rastermask = raster(mask,"noneuc"))
#' 
#' @return Returns the coordinates of the closest point on the raster, as a matrix with two columns (x,y), 
#' named s1 and s2, with first row corresponding to 'from' coordinates, and second row corresponding to 'to' 
#' coordinates.
#' @export closest_coords
#' 
closest_coords=function(from,to,rastermask){
  ends=SpatialPoints(rbind(from,to))
  grid=as(rastermask, 'SpatialGrid') 
  xy=over(ends,grid)
  return(coordinates(grid)[xy,])
}


#' @title Finds vertex index on graph made from raster
#'
#' @description
#'  Finds vertex index on graph made from raster, of points at coordinates pts. Vertex index is just the row of 
#'  the point in the raster object.
#'  
#' @param pts Matrix whose rows are (x,y) coordinates of points on raster
#' @param raster Raster object.
#' 
#' @return Returns the row numbers of raster that correpond to pts. Note that pts must match exactly some 
#' coordinates of raster (use \code{closest_coords} to find closest coordinates if necessary).
#' 
#' @export get_vertices
#' 
get_vertices = function(pts,rastermask){
  target = nearest.raster.point(pts[,1],pts[,2],as.im(rastermask),indices=FALSE)
  coords = coordinates(rastermask) # lookup table from index produced by transition() to coordinates
  npts = dim(pts)[1]
  vert = rep(NA,npts)
  for(i in 1:npts){
    #    vert[i] = which(coords[,1]==target$x[i] & coords[,2]==target$y[i])
    dst = sqrt((coords[,1]-target$x[i])^2 + (coords[,2]-target$y[i])^2)
    vert[i] = which(dst == min(dst))[1]
  }
  return(vert)
}

#' @title Creates the igraph of a mask object
#'
#' @description
#'  Creates an igraph object with a vertex for each mask point and edges to neighbours, weighted according 
#'  to the cost function \code{costfun}, using the mask covariate \code{costname}.
#'  
#'  Requires packages raster, gdistance, igraph
#'  
#' @param mask \code{secr} mask object. Must have covariate called 'noneuc' containing cost
#' @param costname Name of variable to use in cost calculation
#' @param costfun Cost function name
#' @param directed If TRUE, use directed graph for transition between cells, else undirected 
#' @param symm If TRUE, cost is same in both directions 
#' 
#' @examples 
#' ny=4; nx=4 # dimensions of mask
#' # set costs (NA is "hole" - nothing there & can't go there):
#' costs=c(100,100,100,100,NA,100,100,NA,1,NA,100,1,1,1,1,1) 
#' rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),noneuc=costs) # make data frame with coords and costs
#' 
#' rmask=read.mask(data=rmesh,columns="noneuc")  # make mask with covariate 'noneuc' containing cost
#' ig=make_igraph(rmask,"noneuc")
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' cfun=function(x) exp(diff(x)) # asymmetric cost function
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun,directed=TRUE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun,directed=TRUE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun,directed=TRUE,symm=FALSE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' @export make_igraph
#' 
make_igraph = function(mask,costname,costfun="mean",directed=FALSE,symm=TRUE) {
  
  
  if(!is.element(costname,names(covariates(mask))))
    stop(paste("'",costname,"'"," is not the name of one of the mask covariates.",sep=""))
  rastermask = raster(mask,costname) # make raster with covariates(mask)$costname as values of pixels
  
  f=match.fun(costfun)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/mean(x),directions=8)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/exp(diff(x)),directions=8,symm=FALSE)
  tr1<-transition(rastermask,transitionFunction=function(x) 1/f(x),directions=8,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  #costs1<-costDistance(tr1CorrC,pts)
  
  pts = closest_coords(from,to,rastermask)
  vpts = get_vertices(pts,rastermask)
  
  trmat=summary(tr1CorrC)
  rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
  if(directed) g = graph_from_data_frame(rel,directed=TRUE,vertices=NULL)
  else g = graph_from_data_frame(rel,directed=FALSE,vertices=NULL)
  attributes(g)$noneuc=1/trmat$x
  E(g)$weight=1/trmat$x
  
  return(g)  
}




plot.eland.detprob = function(fit,mask=NULL,occ="all",...) {
  cams = traps(fit$capthist)
  ch = fit$capthist
  if(is.null(mask)) mask = fit$mask
  
  betas = coefficients(fit) # beta paameters
  
  # Make linear predictor and calculate g0
  g0s = which(substr(row.names(betas),1,2)=="g0")
  beta.g0 = betas$beta[g0s]
  X.g0 = model.matrix(fit$model$g0,data=elandmask)
  masklp.g0 = X.g0%*%beta.g0
  maskg0.hat=invlogit(masklp.g0)[,1]
  
  # Make linear predictor and calculate sigma
  sigmas = which(substr(row.names(betas),1,5)=="sigma")
  beta.sigma = betas$beta[sigmas]
  X.sigma = model.matrix(fit$model$sigma,data=elandmask)
  masklp.sigma = X.sigma%*%beta.sigma
  masksigma.hat=exp(masklp.sigma)[,1]
  
  #function to calculate Halfnormal dectection function:
  pdet = function(d,g0,sigma,nocc) {
    npar = length(g0)
    ncam = dim(d)[1]
    p = matrix(rep(NA,npar*ncam),nrow=ncam)
    for(i in 1:ncam) p[i,] = 1 - (1-g0*exp(-d[i,]^2/(2*sigma^2)))^nocc
    return(p)
  }
  # Function to combine across independen detect probs
  combdet = function(ps) 1 - prod(1-ps)
  
  ncam = dim(cams)[1]
  nmask = dim(elandmask)[1]
  dists = matrix(rep(NA,nmask*ncam),nrow=ncam)
  for(i in 1:ncam) {
    dists[i,] = distancetotrap(mask,cams[i,])
  }
  nocc = dim(ch)[2]
  #  pest = rbind(pest,pestot)
  
  if(occ=="all") {
    pest = pdet(dists,maskg0.hat,masksigma.hat,nocc=nocc)
    pestot = apply(pest,2,combdet)
  } else if(occ=="single") {
    pest = pdet(dists,maskg0.hat,masksigma.hat,nocc=1)
    pestot = apply(pest,2,combdet)
  } else stop("Invalid occ passed")
  
  covariates(mask)$p = pestot
  plotcovariate(mask,"p",...)
  
}


plotcovariate = function(mask,covariate, ...) {
  require(sp)
  
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  pts = cbind(mask$x,mask$y)
  spdfdat = data.frame(covariate=covariates(mask)[[cnum]])
  spdf = SpatialPixelsDataFrame(pts, spdfdat)
  plot(spdf, ...)
  
  invisible(spdf)
}


# Plot encounter rate in 3D
plotER3d = function(capthist, mask=NULL, occasions=1, 
                    add.points=TRUE,add.text=FALSE, add.maskedge=FALSE,add.mask=TRUE, add.traps=TRUE) {
  require(rgl)
  require(pals)
  
  if(length(session)>1) stop("Only one session at a time: capthist must be a matrix, not a list.")
  if(max(occasions)>dim(capthist)[2]) stop("occasion bigger than number of occasions in cpature history.")
  if(min(occasions)<1) stop("occasion must be greater than zero.")
  
  occasions = as.integer(occasions)
  
  dets = traps(capthist)
  asp= c(1,diff(range(dets$y))/diff(range(dets$x)),1)
  effort = usage(dets)[,occasions]
  ndets = apply(capthist[,occasions,],2,"sum")
  er = ndets/effort
  zlim=range(0,er)
  if(!is.null(mask)){
    xlim = range(mask$x)
    ylim = range(mask$y)
    plot3d(dets$x,dets$y,er, size=10,type="h",lwd=1,xlim=xlim,ylim=ylim,zlim=zlim,
           xlab="Easting",ylab="Northing",zlab="Encounter Rate",aspect=asp)
  } else {
    plot3d(dets$x,dets$y,er, size=10,type="h",lwd=1,zlim=zlim,
           xlab="Easting",ylab="Northing",zlab="Encounter Rate",aspect=asp)
  }
  rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = 'white')
  if(add.points) {
    pointcolr = parula(max(ndets[er>0])+1)[ndets[er>0]+1]
    points3d(dets$x[er>0],dets$y[er>0],er[er>0],col=pointcolr,size=10,add=TRUE)
  }
  if(add.text) {
    textcolr = parula(max(ndets)+1)[ndets+1]
    text3d(dets$x,dets$y,1.05*er, texts=signif(er,2),col=textcolr,add=TRUE, cex=0.75)
  }
  
  if(!is.null(mask)){
    if(!add.maskedge & !add.mask) add.mask=TRUE # if pass mask, plot it!
    # add mask edge if asked to
    if(add.maskedge) {
      edge = plotMaskEdge(mask,plt=FALSE,add=TRUE)
      xedge = edge[c(1,3),]
      yedge = edge[c(2,4),]
      segments3d(xedge,yedge,z=-0.01*diff(zlim))
    }
    
    # add mask if asked to
    if(add.mask) {
      points3d(mask$x,mask$y,z=-0.01*diff(zlim),col="gray")
    }
  }
  # add traps
  if(add.traps) points3d(dets$x,dets$y,z=0,pch=1,col="red")
  
}

