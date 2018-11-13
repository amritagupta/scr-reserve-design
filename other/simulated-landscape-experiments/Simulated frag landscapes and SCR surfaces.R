#20160829
#D. MORIN

#simulating continuous to binary landscape surfaces and getting fragstat metrics
#to select high and low frag surfaces for Amrita's optimization experiments

library(raster)
library(scrbook)
library(rgeos)
library(rgdal)
library(gdistance)
library(igraph)
library(SDMTools)
library(RColorBrewer)

#source modified patchy covariate surface (allows for different dimensions and seeds)
patch.fun<-function(cells,seed){ ####function to change raster size (cells=sqrt of total)
  library("gdistance")
  library("raster")
  par(mfrow = c(1, 1))
  set.seed(seed)     #67 is the seed for DWC paper
  r <- raster(nrows = cells, ncols = cells)
  projection(r) <- "+proj=utm +zone=12 +datum=WGS84"
  extent(r) <- c(0.5, 4.5, 0.5, 4.5)
  delta <- (4.5 - 0.5)/cells
  gx <- seq(0.5 + delta/2, 4.5 - delta/2, , cells)
  gy <- rev(gx)
  gx <- sort(rep(gx, cells))
  gy <- rep(gy, cells)
  grid <- cbind(gx, gy)
  Dmat <- as.matrix(dist(grid))
  V <- exp(-Dmat/0.5)
  z <- t(chol(V)) %*% rnorm(cells*cells)
  z <- (z - mean(z))/sqrt(var(as.vector(z))) #
  values(r) <- matrix(z, cells, cells, byrow = FALSE)
  plot(r)
  covariate.patchy <- r
  class(covariate.patchy)  
  return(covariate.patchy)
  list(covariate.patchy)
}

#function to simulate continuous rasters, convert to binary, and produce fragstat metrics
#cells is the sqrt of total number of cells (number of rows = columns)
sim.fun<-function(cells,seed){
  out<-list()
covariate<-patch.fun(cells=cells,seed=seed) #cells X cells raster (covariate.patchy)
#convert to a binary surface
z.cont<-values(covariate)
z.tmp<-NULL
for (i in 1:length(z.cont)){
  if(z.cont[i]>median(z.cont)){
    z.tmp[i]<-0 #high resistance = bad habitat
  }
  else {
    z.tmp[i]<-1
  }
}
rast.tmp<-rasterFromXYZ(cbind(coordinates(covariate),z.tmp))
class.tmp<-ClassStat(mat=rast.tmp, cellsize = 0.1, bkgd = 0, latlon = FALSE) #bkgd means 0 values are the background (1 is the focal)
patch.tmp<-PatchStat(mat=rast.tmp, cellsize = 0.1, latlon = FALSE)
stats<-c(sum(z.tmp),patch.tmp[2,]$perim.area.ratio,
         patch.tmp[2,]$shape.index,patch.tmp[2,]$frac.dim.index,
         patch.tmp[2,]$core.area.index,
         class.tmp$aggregation.index,class.tmp$n.patches,class.tmp$patch.density,
         class.tmp$edge.density,class.tmp$mean.patch.area,
         class.tmp$perimeter.area.frac.dim)
out<-c(seed,stats,z.cont)
return(out)
}

nsim<-10000
seed.v<-c(1:nsim)

out.mat<-matrix(NA,nrow=nsim,ncol=1612) #seed + 11 metrics + 1600 pixel values
for (i in 1:nsim){
  out.mat[i,]<-sim.fun(seed=seed.v[i],cells=40)
}

out.df<-as.data.frame(out.mat[,1:12])#separate from z values for the corresponding raster
colnames(out.df)<-c("seed","ncell","p.perim.area.ratio","p.shape.index",
                    "p.frac.dim.index","p.core.area.index","c.aggregation.index",
                    "c.n.patches","c.patch.density","c.edge.density",
                    "c.mean.patch.area","c.perimeter.area.frac.dim"
                    )
head(out.df)
which(out.df$c.aggregation.index<(min(out.df$c.aggregation.index)+1))
which(out.df$c.aggregation.index>(max(out.df$c.aggregation.index)-2))
which(out.df$c.n.patches<(min(out.df$c.n.patches)+2))
which(out.df$c.n.patches>max(out.df$c.n.patches)-4)
#which(out.df$c.n.patches==min(out.df$p.frac.dim.index))
#which(out.df$c.n.patches==max(out.df$p.frac.dim.index))

#create vector of seeds want maps for
z.seed<-c(which(out.df$c.aggregation.index<(min(out.df$c.aggregation.index)+1)),
          which(out.df$c.aggregation.index>(max(out.df$c.aggregation.index)-2)),
          which(out.df$c.n.patches<(min(out.df$c.n.patches)+2)),
          which(out.df$c.n.patches>max(out.df$c.n.patches)-4))#,
          #which(out.df$c.n.patches==min(out.df$p.frac.dim.index)),
          #which(out.df$c.n.patches==max(out.df$p.frac.dim.index)))
length(z.seed)

#get coordinates so can create rasters with z-values
coor.fun<-function(cells){ #coordniates for 40 by 40 raster with same extent
  r <- raster(nrows = cells, ncols = cells)
  projection(r) <- "+proj=utm +zone=12 +datum=WGS84"
  extent(r) <- c(0.5, 4.5, 0.5, 4.5)
  return(coordinates(r))
}

xy<-coor.fun(cells=40)

#running again with the median as cut-off for good and bad hab, so need to
#rename previous output so don't lose it in WS.
#out.mat.0<-out.mat
#out.df.0<-out.df

####function to produce frag figures
frag.fig.fun<-function(z.seed){
  z.cont<-out.mat[z.seed,13:1612]
  z.tmp<-NULL
  for (i in 1:length(z.cont)){
    if(z.cont[i]>median(z.cont)){
      z.tmp[i]<-0 #high resistance = bad habitat
    }
    else {
      z.tmp[i]<-1
    }
  }
  z.bin<-z.tmp
  library(RColorBrewer)
  col<-brewer.pal(10,"Spectral")
  col.bin<-c("grey","forestgreen")
  par(mfrow=c(1,2))
  plot(rasterFromXYZ(cbind(xy,z.cont)),col=rev(col),
       main=paste("seed =",z.seed,", agg.index=",round(out.df$c.aggregation.index[z.seed],2)))
  plot(rasterFromXYZ(cbind(xy,z.bin)),col=col.bin, 
       main=paste(round(out.df$ncell[z.seed]/1600,2),"% good habitat,", out.df$c.n.patches[z.seed], "patches"))
  
} #end figure function

frag.fig.fun(z.seed=9799)

for (i in 1:length(z.seed)){
  frag.fig.fun(z.seed=z.seed[i])
}

#create rasters for each seed
rast.list<-list()
for (i in 1:length(z.seed)){
  rast.list[[i]]<-rasterFromXYZ(cbind(xy,out.mat[z.seed[i],13:1612]))
}


#Table of selected seeds for candidate surfaces

cand.seed<-c(75, 6208, 6270, 8901,9799,245,3809,6530,8606)
cand.df<-rbind(out.df[75,],out.df[6208,],out.df[6270,],out.df[8901,],out.df[9799,],out.df[245,],
               out.df[3809,],out.df[6530,],out.df[8606,])
cand.low<-rbind(out.df[2175,], out.df[5125,],out.df[8435,],out.df[602,],out.df[1467,],
                out.df[3115,], out.df[8336,])

low.frag<-rasterFromXYZ(cbind(xy,out.mat[5125,13:1612]))
high.frag<-rasterFromXYZ(cbind(xy,out.mat[9799,13:1612]))
save(low.frag,file="low.frag.RData")
save(high.frag,file="high.frag.RData")



##################################################
#simulating SCR data

#First, need to find appropriate theta value (effective sigma) so mean home range area is consistent 
#across a2 values 
mean.HR.fun<-function(sigma,a2,covariate){
  a1 <- 1/(2*sigma*sigma) 
  xy<-coordinates(covariate)
  cost.out <- exp(a2 * covariate) #estimate cost surface based on covariate
  trFn <- function(x)(1/(mean(x))) #transition function. By taking 1/mean, transforms cost to conductance
  tr1.out <- transition(cost.out, transitionFunction=trFn, directions=8) #apply transition function
  tr1.out <- geoCorrection(tr1.out,type="c",multpl=F,scl=F) #correct for curvature and diagonals
  dmat <- costDistance(tr1.out, xy,xy) #create a matrix with a row for each raster cell and a column for every other cell and the conductance to that cell.
  potmat <- 1*exp(-a1*dmat*dmat) #encouter model with a0=1 and cost-weighted distance
  area.p<-NULL
  for (i in 1:ncell(covariate)){
    area.p[i]<-length(which(potmat[i,]>(0.95*(sum(potmat[i,]))/ncell(covariate))))
  }
  area.r<-rasterFromXYZ(cbind(xy,area.p))
  x.box<-c(2,3) 
  y.box<-c(2,3)
  box<-raster(nrow=1,ncol=1)
  extent(box)<-c(x.box,y.box)
  HR<-crop(area.r,box) #this creastes a raster with the home range area for an individual with an activity center in a pixel, for 100 pixels in the middle of the covariate surface.
  data.row<-c(a2,mean(HR@data@values),sigma)
  out<-list(data.row,HR)
  return(data.row)
}

sig.lcp<-c(seq(from=0.1,to=1,by=0.1))
a2.lcp<-c(seq(from=0,to=2.25,by=0.25))
scen<-expand.grid(sig.lcp,a2.lcp)
data.lcp.high<-matrix(NA,nrow=nrow(scen),ncol=3)
for (i in 1:nrow(scen)){
  rbind(data.lcp.high[i,]<-mean.HR.fun(sigma=scen[i,1],a2=scen[i,2],covariate=high))
}

data.lcp.low<-matrix(NA,nrow=nrow(scen),ncol=3)
for (i in 1:nrow(scen)){
  rbind(data.lcp.low[i,]<-mean.HR.fun(sigma=scen[i,1],a2=scen[i,2],covariate=low))
}
#

colnames(data.lcp.low)<-c("a2","area","sigma")
colnames(data.lcp.high)<-c("a2","area","sigma")
data.lcp.low<-as.data.frame(data.lcp.low)
data.lcp.high<-as.data.frame(data.lcp.high)

#fit linear model and predict x (sigma) with y held constant for each a2 value

#Y = intecept + coefficient*(X)
# to solve for X, X = (Y - intercept)/coefficient

#we solve for px when euc.sig = 0.15 to allow for constant trap spacing

#make this a function
theta.fun<-function(data,px){
  lm<-lm(data$area~data$sigma) 
  coef(lm)
  intercept<-coef(lm)[1] #intercept
  slope.coef<-coef(lm)[2] #slope coefficient
  sigx<-(px-(intercept))/slope.coef
  d.row<-c(data$a2[1],intercept,slope.coef,sigx)
  return(d.row)
}

#make a list of the data frames then use a for loop
data.low<-list()
data.high<-list()

for (i in 1:length(a2.lcp)){
  data.low[[i]]<-subset(data.lcp.low,a2==a2.lcp[i])
}

for (i in 1:length(a2.lcp)){
  data.high[[i]]<-subset(data.lcp.high,a2==a2.lcp[i])
}
#solve for px when sigma = 0.15

px.fun<-function(data,euc.sig){
  lm.tmp<-lm(data[[1]]$area~data[[1]]$sigma)
  intercept<-coef(lm.tmp)[1] #intercept
  slope.coef<-coef(lm.tmp)[2] #slope coefficient
  px<-intercept+slope.coef*euc.sig
  return(px)
}
px.low<-px.fun(data=data.low,euc.sig=0.15)
px.high<-px.fun(data=data.high,euc.sig=0.15) # should be the same, high or low 
#because it's the area when theta = euc.sigma
# > px.low
# (Intercept) 
# 78.28533 
# > px.high
# (Intercept) 
# 78.28533 

data.low.df<-matrix(NA,nrow=length(data.low),ncol=4)
data.high.df<-matrix(NA,nrow=length(data.high),ncol=4)

for (i in 1:length(data.low)){
  data.low.df[i,]<-theta.fun(data=data.low[[i]],px=px.low)
}

for (i in 1:length(data.high)){
  data.high.df[i,]<-theta.fun(data=data.high[[i]],px=px.high)
}
# > data.low.df
# [,1]      [,2]     [,3]      [,4]
# [1,] 0.00 -31.93867 734.8267 0.1500000
# [2,] 0.25 -41.76400 678.3727 0.1769666
# [3,] 0.50 -48.66867 622.3739 0.2039835
# [4,] 0.75 -50.24933 562.5697 0.2284778
# [5,] 1.00 -50.27133 504.9879 0.2545738
# [6,] 1.25 -48.98600 450.7618 0.2823472
# [7,] 1.50 -46.32800 400.1400 0.3114243
# [8,] 1.75 -43.39467 353.9921 0.3437365
# [9,] 2.00 -40.05533 312.2897 0.3789452
# [10,] 2.25 -36.56667 274.2230 0.4188270

# > data.high.df
# [,1]      [,2]     [,3]      [,4]
# [1,] 0.00 -31.93867 734.8267 0.1500000
# [2,] 0.25 -47.34133 674.9297 0.1861330
# [3,] 0.50 -55.01200 606.8073 0.2196700
# [4,] 0.75 -58.79333 536.1533 0.2556706
# [5,] 1.00 -57.64533 463.6442 0.2931788
# [6,] 1.25 -53.42333 393.5970 0.3346283
# [7,] 1.50 -47.88400 329.1255 0.3833472
# [8,] 1.75 -41.45733 270.9842 0.4418806
# [9,] 2.00 -35.33533 220.6879 0.5148478
# [10,] 2.25 -29.55933 177.8461 0.6063933

########################################
#so for the low frag surface
theta.075.l<-data.low.df[4,4]
theta.100.l<-data.low.df[5,4]
theta.150.l<-data.low.df[7,4]
theta.200.l<-data.low.df[9,4]
theta.225.l<-data.low.df[10,4]

theta.075.h<-data.high.df[4,4]
theta.100.h<-data.high.df[5,4]
theta.150.h<-data.high.df[7,4]
theta.200.h<-data.high.df[9,4]
theta.225.h<-data.high.df[10,4]

############################################
#simulate IPP for N= 100 and N=250, a2 = 0.75 and 2.25, and high and low frag

IPP.fun<-function(exp.N=NULL,cov=NULL){
  xy<-coordinates(cov)
  #exp.N<-N #specify in function
  area<-nrow(xy) #number of pixels
  exp.D<-exp.N/area #density per pixel
  #simulate the s matrix first using the covariate
  set.seed(31025) 
  #library(R2Cuba)
  beta0<-log(exp.D)
  beta1<- -3 #this determines the intensity of the ipp
  cov.fun <- function(s,cov){
    s <- matrix(s,ncol=2)
    return(extract(cov,s))
  }
  ######################
  #intesity function
  mu<-function(s,beta0,beta1) exp(beta0+beta1*cov.fun(s=s, cov=cov))
  
  s<-matrix(NA,exp.N,2)
  cov.min<-min(values(cov))
  cov.max<-max(values(cov))
  Q<-max(c(exp(beta0+beta1*cov.min),exp(beta0+beta1*cov.max)))
  counter<-1
  while(counter<=exp.N){
    x.c<-runif(1,xmin(cov),xmax(cov));y.c<-runif(1,ymin(cov),ymax(cov))
    s.cand<-c(x.c,y.c)
    pr<-mu(s.cand,beta0,beta1) #/ EN
    if(runif(1) < pr/Q) {
      s[counter,]<-s.cand
      counter<-counter+1
    }
  }
 return(s) 
  
}

s.high.100<-IPP.fun(exp.N=100,cov=high)
s.high.250<-IPP.fun(exp.N=250,cov=high)
s.low.100<-IPP.fun(exp.N=100,cov=low)
s.low.250<-IPP.fun(exp.N=250,cov=low)

#hist.plot<-hist(extract(covariate,s))
#summ<-summary(extract(covariate,s))
#library(RColorBrewer)
IPP.plot.fun<-function(cov=NULL,s=NULL){
  plot(cov, col=rev(brewer.pal(9,"YlOrBr")))
  points(s,pch=20,col="black")
}
par(mfrow=c(2,2))
IPP.plot.fun(cov=high,s=s.high.100)
IPP.plot.fun(cov=high,s=s.high.250)
IPP.plot.fun(cov=low,s=s.low.100)
IPP.plot.fun(cov=low,s=s.low.250)
par(mfrow=c(1,1))

#simulate truth
tru.fn<-function(cov,alpha2,theta,directions,s){
  cell.trap<-coordinates(cov)#to get a point for each cell, just use the raster coordinates
  ntraps.1<-nrow(cell.trap) #this should be the number of pixels
  K.1<-1 #detection is perfect so only need one occasion 
  cost<-exp(alpha2*cov)
  alpha1<-1/(2*theta*theta)
  S.tmp<-SpatialPoints(s)
  Nx.tmp<-rasterize(S.tmp,cov,fun='count')
  Nx.tmp[is.na(Nx.tmp[])] <- 0
  D.z<-values(Nx.tmp)
    ## compute least-cost path distance
    tr1<-transition(cost, transitionFunction=function(x) 1/mean(x),directions=directions)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  D.1<-costDistance(tr1CorrC,cell.trap,cell.trap)
  potmat<-1*exp(-alpha1*D.1*D.1) 
  #pc.z<-rowSums(potmat) #vector of pc values for raster
  #pc.r<-rasterFromXYZ(cbind(cell.trap,pc.z)) #pc raster
  realmat <- t(D.z *t(exp(-(alpha1)*D.1^2)))
  out<-list(D.z=D.z,potmat=potmat,realmat=realmat,s=s,cov=cov) #if want if want the z values as a vector, switch to return(pc.z), the matrix, switch to return(potmat)
  return(out)
}

#theta estimated above
th.low.075<-0.229
th.low.225<-0.419
th.high.075<-0.256
th.high.225<-0.606

tru.low.N100.a2.075<-tru.fn(cov=low,alpha2=0.75,theta=th.low.075,directions=8,s=s.low.100)
tru.high.N100.a2.075<-tru.fn(cov=high,alpha2=0.75,theta=th.high.075,directions=8,s=s.high.100)
tru.low.N250.a2.075<-tru.fn(cov=low,alpha2=0.75,theta=th.low.075,directions=8,s=s.low.250)
tru.high.N250.a2.075<-tru.fn(cov=high,alpha2=0.75,theta=th.high.075,directions=8,s=s.high.250)
tru.low.N100.a2.225<-tru.fn(cov=low,alpha2=2.25,theta=th.low.225,directions=8,s=s.low.100)
tru.high.N100.a2.225<-tru.fn(cov=high,alpha2=2.25,theta=th.high.225,directions=8,s=s.high.100)
tru.low.N250.a2.225<-tru.fn(cov=low,alpha2=2.25,theta=th.low.225,directions=8,s=s.low.250)
tru.high.N250.a2.225<-tru.fn(cov=high,alpha2=2.25,theta=th.high.225,directions=8,s=s.high.250)


TRU<-list(tru.low.N100.a2.075,
          tru.high.N100.a2.075,
          tru.low.N250.a2.075,
          tru.high.N250.a2.075,
          tru.low.N100.a2.225,
          tru.high.N100.a2.225,
          tru.low.N250.a2.225,
          tru.high.N250.a2.225)
save(TRU,file="TRU.list.RData") #this is a workspace with only the requested true data

#for estimated surfaces need to simulate trap locations
trap.fun<-function(){
  tr<-1.5 #trap spacing relative to euc.sig
  cen<-2.5 #center of the covariate surface
  euc.sig<-0.15
  xg.10<-c(seq(cen-9.5*(tr*euc.sig),cen+9.5*(tr*euc.sig),(tr*euc.sig))); yg.10<-c(seq(cen+9.5*(tr*euc.sig),cen-9.5*(tr*euc.sig),-(tr*euc.sig)))    
  xg.sub<-xg.10[xg.10>1]      #creates buffer around trap locs 
  xg.sub<-xg.sub[xg.sub<4]    #so as trap spcing increases
  yg.sub<-yg.10[yg.10>1]      #overall number of traps decreases
  yg.sub<-yg.sub[yg.sub<4]    #area sampled remains the same
  X<-cbind(sort(rep(xg.sub,length(xg.sub))),rep(yg.sub,length(yg.sub)))
  #ntraps<-nrow(X)
  return(X)
}
X<-trap.fun()
ntraps<-nrow(X)

#simulate detecting inds as they move around and then are detected in traps
#simulating with IPP, but fitting using homogeneous point process

Yij.fun<-function(theta = NULL,S=NULL,alpha2=NULL,cov=NULL,X=NULL, directions=NULL){
  N<-nrow(S)
  sigma<-theta    
  alpha0<- -2 #keep constant
  alpha1<-1/(2*sigma*sigma) 
  K<- 10
  ntraps<-nrow(X)
  ## compute least-cost path distance
  cost<-exp(alpha2*cov)
  tr1<-transition(cost, transitionFunction=function(x) 1/mean(x),directions=directions)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  # generate cost surface
  D<-costDistance(tr1CorrC,S,X)
  probcap<-plogis(alpha0)*exp(-alpha1*D*D) #for estimation
  # generate capture histories
  Y<-matrix(NA,nrow=N,ncol=ntraps)
  for(i in 1:nrow(Y)) {
  Y[i,]<-rbinom(ntraps,K,probcap[i,])
  }
  Y<-Y[apply(Y,1,sum)>0,] #remove those not detected
  return(Y)
}

Y.low.N100.a2.075<-Yij.fun(theta = th.low.075,S=s.low.100,alpha2=0.75,cov=low,X=X, 
                           directions=8)
Y.low.N100.a2.225<-Yij.fun(theta = th.low.225,S=s.low.100,alpha2=2.25,cov=low,X=X, 
                           directions=8)
Y.low.N250.a2.075<-Yij.fun(theta = th.low.075,S=s.low.250,alpha2=0.75,cov=low,X=X, 
                           directions=8)
Y.low.N250.a2.225<-Yij.fun(theta = th.low.225,S=s.low.250,alpha2=2.25,cov=low,X=X, 
                           directions=8)
Y.high.N100.a2.075<-Yij.fun(theta = th.high.075,S=s.high.100,alpha2=0.75,cov=high,X=X, 
                           directions=8)
Y.high.N100.a2.225<-Yij.fun(theta = th.high.225,S=s.high.100,alpha2=2.25,cov=high,X=X, 
                           directions=8)
Y.high.N250.a2.075<-Yij.fun(theta = th.high.075,S=s.high.250,alpha2=0.75,cov=high,X=X, 
                           directions=8)
Y.high.N250.a2.225<-Yij.fun(theta = th.high.225,S=s.high.250,alpha2=2.25,cov=high,X=X, 
                           directions=8)

#likelihood function

############
SCRed.DM<-function(start = NULL, y = NULL, K = NULL, X = NULL,
                   cov=NULL, directions=NULL, model="B",dist="ecol",predict=FALSE){
  G = coordinates(cov)
  SSarea <- nrow(G)
  nG <- nrow(G)
  
  ## cost distance -- This (#~#) differs from Euclidean distance model  #~#
  alpha2 <- exp(start[4])                                               #~# start value provided 4th
  cost <- exp(alpha2 * cov)                                             #~# have to provide a covariate layer
  tr <- transition(cost, transitionFunction=function(x) (1/(mean(x))),  #~#
                   direction = directions)                              #~#
  trLayer <- geoCorrection(tr, scl = F)                                 #~#
  D <- costDistance(trLayer,as.matrix(X),as.matrix(G))                  #~# calculate ecol dist between trap locs (X) and pixels
  
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
  if(predict==FALSE){                             #this gives model estimates
    nv <- c(rep(1, length(lik.marg) - 1), n0)
    part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
    part2 <- sum(nv * log(lik.marg))
    out <- -1 * (part1 + part2)
    attr(out, "SSarea") <- SSarea
    return(out)
  }
  if(predict==TRUE){                              #this uses estimates to predict realized density surface
    posterior<-matrix(NA,nrow=nG,ncol=nrow(ymat))
    for(i in 1:nrow(ymat)){
      if(model=="B")
        Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap[1:length(Pm)], log = TRUE))
      if(model=="P")
        Pm[1:length(Pm)] <- (dpois(rep(ymat[i, ], nG), rep(K, nG)*probcap[1:length(Pm)], log = TRUE))
      lik.cond <- exp(colSums(Pm))
      posterior[,i]<- (lik.cond*(1/nG))/lik.marg[i]      # pr(s) = pr(y|s)*pr(s) / pr(y)
    }
    return(cbind(G,posterior))
  }  
}

#fit and predict
fit.fn<-function(theta=NULL,alpha2=NULL,N=NULL,Y=NULL,
                 S=NULL,K=NULL,cov=NULL,directions=NULL,X=NULL){
  xy<-coordinates(cov)
  alpha0<- -2 #keep constant
  alpha1<-1/(2*theta*theta)
  #estimate parameters
  outEcol<-nlm(SCRed.DM,c(-alpha0,log(alpha1),log(N-nrow(Y)),log(alpha2)),hessian=TRUE,y=Y,K=K,
               X=X,cov=cov,
               directions=directions, model = "B", dist="ecol",predict=FALSE)
  #predict = realized density
  pred<-SCRed.DM(start=c(outEcol$estimate[1],outEcol$estimate[2],outEcol$estimate[3],
                  outEcol$estimate[4]), y=Y, K=K,X=X,
                  cov=cov,directions=directions, model = "B", dist="ecol",predict=TRUE)
  #estimate surfaces
  n0 <- exp(outEcol$estimate[3]) # *******Set this to your estimate of n0
  dees <- pred[,-c(1,2)] #just the values without the columns with coordinates
  dees[,ncol(dees)] <- dees[,ncol(dees)] * n0
  D.z<-rowSums(dees) #realized density surface
  
  #for potential connectivity,
  #evaluate the SCR encounter model using the mle
  #estimates for alpha1 and alpha2
  a1<-exp(outEcol$estimate[2])
  a2<-exp(outEcol$estimate[4])
  cost.out <- exp(a2 * cov) #estimated cost, as opposed to TRUE cost above
  tr2.out <- transition(cost.out, transitionFunction=function(x)1/mean(x), directions=directions)
  tr2.out <- geoCorrection(tr2.out,type="c",multpl=F,scl=F)
  
  ###
  dmat <- costDistance(tr2.out, xy, coordinates(cov)) #lcp matrix
  potmat <- exp(-(a1)*dmat^2)
  #pc<-rowSums(potmat) #potential connectivity surface
  realmat <- t(D.z *t(exp(-(a1)*dmat^2)))
  #dwc<-rowSums(realmat) #estimated density weighted connectivity surface
  #compile output
  pred.out<-list(D.z=D.z, potmat=potmat, realmat=realmat,S=S,cov=cov)

  #compile output
  out.2<-list()
  out.1 <- as.vector(c(outEcol$estimate[1], outEcol$estimate[2], outEcol$estimate[3], 
                       outEcol$estimate[4], nrow(Y), exp(outEcol$estimate[3])+nrow(Y), 
                       alpha2, theta))
  
  out.2<-list(pred.out=pred.out,pred=pred,out.1=out.1)
  return(out.2)
}

out.low.N100.a2.075<-fit.fn(theta=th.low.075,alpha2=0.75,N=100,S=s.low.100,
                            Y=Y.low.N100.a2.075,K=10,cov=low,
                            directions=8,X=X)
out.low.N100.a2.225<-fit.fn(theta=th.low.225,alpha2=2.25,N=100,S=s.low.100,
                            Y=Y.low.N100.a2.225,K=10,cov=low,
                            directions=8,X=X)
out.high.N100.a2.075<-fit.fn(theta=th.high.075,alpha2=0.75,N=100,S=s.high.100,
                            Y=Y.high.N100.a2.075,K=10,cov=high,
                            directions=8,X=X)
out.high.N100.a2.225<-fit.fn(theta=th.high.225,alpha2=2.25,N=100,S=s.high.100,
                            Y=Y.high.N100.a2.225,K=10,cov=high,
                            directions=8,X=X)
out.low.N250.a2.075<-fit.fn(theta=th.low.075,alpha2=0.75,N=250,S=s.low.250,
                            Y=Y.low.N250.a2.075,K=10,cov=low,
                            directions=8,X=X)
out.low.N250.a2.225<-fit.fn(theta=th.low.225,alpha2=2.25,N=250,S=s.low.250,
                            Y=Y.low.N250.a2.225,K=10,cov=low,
                            directions=8,X=X)
out.high.N250.a2.075<-fit.fn(theta=th.high.075,alpha2=0.75,N=250,S=s.high.250,
                             Y=Y.high.N250.a2.075,K=10,cov=high,
                             directions=8,X=X)
out.high.N250.a2.225<-fit.fn(theta=th.high.225,alpha2=2.25,N=250,S=s.high.250,
                             Y=Y.high.N250.a2.225,K=10,cov=high,
                             directions=8,X=X)

EST<-list(est.low.N100.a2.075=out.low.N100.a2.075$pred.out,
          est.low.N100.a2.225=out.low.N100.a2.225$pred.out,
          est.high.N100.a2.075=out.high.N100.a2.075$pred.out,
          est.high.N100.a2.225=out.high.N100.a2.225$pred.out,
          est.low.N250.a2.075=out.low.N250.a2.075$pred.out,
          est.low.N250.a2.225=out.low.N250.a2.225$pred.out,
          est.high.N250.a2.075=out.high.N250.a2.075$pred.out,
          est.high.N250.a2.225=out.high.N250.a2.225$pred.out)
save(EST,file="EST.list.RData")
