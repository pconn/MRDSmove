### function to produce horvitz-thompson like estimators of abundance that
# integrate over true distance bin using Bayes rule to get the mixing distribution

expit<-function(x)(1/(1+exp(-x)))


#' Calculate Horvitz-Thompson like estimator from detection, movement, and measurement error parameters by integrating
#' over latent distance using Bayes rule.  Currently conditioning on animal being available to be detected (i.e. within Obs.bins) for the 1st observer
#' @param Par Parameter vector, including detection parameters, movement error SDs (left and right tail), measurement error SD
#' @param Data A design.matrix with the following column names: "match" indicates which records match with which (there should be two records
#'        for each detection, one for each observer), "observer","species" (provides species or other grouping variable: abundance estimates will be provided separately for each), 
#'       "obs.dist" (observed distance; NA if missing), "g_size" (group/cluster size), "moving" (binary indicator for moving/not moving), 
#'       "detected" (binary detection/nondetection). Finally, the optional column "count" gives the count of histories observed a particular type (to make likelihood calculations faster).
#'        Additional covariates may also be provided and used in formula (e.g. "other"; see below)
#' @param G A vector of group sizes, one for each detected animal cluster
#' @param mod.formula Formula object giving formula for detection probability for each observer.  Variables can be linked to column names in Data, 
#'        but need to use "distance" to represent distance and "distance2" to represent distance^2.  Adding in a variable "other" into Data (whether the other observer detected it or not)
#'        can be used to implement symmetric point/limiting independence as described in discussion of MacKenzie and Clement 2016.  Note that the
#'        interaction between "other" and "distance" can be used to implement point indendence.
#' @param Bin.widths Vector of distance bin widths
#' @param Obs.bins A vector giving which bins are observed (e.g. 1:3 if bins 1-3 are observed)
#' @param Move.fix A indicator vector giving which movement/measurement error parameters to fix to 0 (omit if all are estimated)
#' @param gaussian.move If TRUE, uses two half-normals for movement; if FALSE (default), uses exponential / half exponential (Laplace dist)
#' @param gaussian.meas If TRUE, uses normal distribution for measurement error; if FALSE (default), uses double exponential (Laplace dist)
#' @return a community mrds dataset
#' @export
#' @keywords simulation, mrds
#' @author Paul B. Conn
ht_mrds <- function(Par,Data,G,mod.formula,Bin.widths,Obs.bins,Move.fix=Move.fix,gaussian.move=FALSE,gaussian.meas=FALSE){
  Cur.par = Par
  n.par = length(Par)
  if(is.null(Move.fix)==FALSE | sum(Move.fix)>0){
    n.est = 3-sum(Move.fix)
    Est.ind = c(rep(1,n.par-n.est),1-Move.fix)
    Which.est = which(Est.ind==1)
    Cur.par=rep(-3,n.par+sum(Move.fix))
    if(Cur.par[length(Cur.par)]==-3 & gaussian.meas==FALSE)Cur.par[length(Cur.par)]=3  #large value of this parameter equates to no measurement error
    if(Cur.par[length(Cur.par)-1]==-3 & gaussian.move==FALSE)Cur.par[length(Cur.par)-1]=3
    if(Cur.par[length(Cur.par)-2]==-3 & gaussian.move==FALSE)Cur.par[length(Cur.par)-2]=3
  }  
  Temp=Bin.midpoints=0*Bin.widths
  n.obs.bins=length(Obs.bins)
  n.bins = length(Bin.widths)
  n.hists= nrow(Data)/2
  n.par=length(Cur.par)
  if(n.bins<2)cat('ERROR: function not set up for <2 distance bins')
  if(is.null(Data$count))Count=rep(1,n.hists)
  else Count = matrix(Data$count,2,n.hists)[1,]

  #calculate bin midpoints given bin widths (could be passed into function in to make it faster)
  Temp[1]=Bin.widths[1]
  Bin.midpoints=Bin.widths/2
  Obs.dists = matrix(Data$obs.dist,2,n.hists)
  for(ibin in 2:n.bins){
    Temp[ibin]=Temp[ibin-1]+Bin.widths[ibin]
    Bin.midpoints[ibin] = Temp[ibin-1]+0.5*(Temp[ibin]-Temp[ibin-1])
  }
  Dists = Bin.midpoints+0.5
  Dist.probs = Bin.widths/sum(Bin.widths)  #probability of groups being in each distance bin given a uniform distribution
  
  #categorize histories into 11, 10, 01
  Detected = matrix(Data$detected,2,n.hists)
  Move = Data$moving[2*(1:n.hists)-1]
  Hist.type = Detected[1,]-Detected[2,]  #so 01 is a -1, 11 is a 0, 10 is a 1
  Which.move = which(Move==1)
  Which.11 = which(Hist.type==0)
  Which.01 = which(Hist.type==-1)
  Which.10 = which(Hist.type==1)

  
  #parameterize Psi transition matrix, Measurement error matrix
  Psi = Meas = matrix(0,n.bins,n.bins)
  if(gaussian.move==TRUE){
    for(ibin1 in 1:n.bins)Psi[ibin1,ibin1:n.bins]=dnorm(Bin.midpoints[ibin1:n.bins],Bin.midpoints[ibin1],exp(Cur.par[n.par-1]))+0.00001
    for(ibin1 in 2:n.bins)Psi[ibin1,1:(ibin1-1)]=dnorm(Bin.midpoints[1:(ibin1-1)],Bin.midpoints[ibin1],exp(Cur.par[n.par-2]))
  }
  if(gaussian.move==FALSE){  #double exponential
    for(ibin1 in 1:n.bins)Psi[ibin1,ibin1:n.bins]=dexp(Bin.midpoints[ibin1:n.bins]-Bin.midpoints[ibin1],exp(Cur.par[n.par-1]))+0.00001
    for(ibin1 in 2:n.bins)Psi[ibin1,1:(ibin1-1)]=dexp(Bin.midpoints[ibin1]-Bin.midpoints[1:(ibin1-1)],exp(Cur.par[n.par-2]))
  }
  if(gaussian.meas==TRUE){
    for(ibin1 in 1:n.bins)Meas[ibin1,] = dnorm(Bin.midpoints,Bin.midpoints[ibin1],exp(Cur.par[n.par]))+0.00001
  }
  if(gaussian.meas==FALSE){
    for(ibin1 in 1:n.bins)Meas[ibin1,] = dexp(abs(Bin.midpoints-Bin.midpoints[ibin1]),exp(Cur.par[n.par]))+0.00001
  }
  Psi = Psi/rowSums(Psi)
  Psi.nomove = diag(n.bins)
  Meas = Meas/rowSums(Meas)
  Pin = rowSums(Meas[,Obs.bins])  #probability of being observed in the strip width given location
  
  P.hist = rep(0,n.hists)
  P = array(0,dim=c(n.bins,n.bins,2,n.hists))  #to hold (conditional) detection probabilities for each observer
  Pdot = rep(0,n.hists)
  
  # #1) calculate probability of being in d1, d2 given detection history
  P.dist = Pdot = array(1,dim=c(n.bins,n.bins,n.hists))  
  for(ibin1 in 1:n.bins){  #observer 1 latent distance
    for(ibin2 in 1:n.bins){  #observer 2 latent distance
      Data$distance = rep(c(Dists[ibin1],Dists[ibin2]),n.hists)
      #Data$distance = rep(c(Bin.midpoints[ibin1],Bin.midpoints[ibin2]),n.hists)
      Data$distance2 = Data$distance^2
      X = model.matrix(mod.formula,Data)
      P[ibin1,ibin2,,]=expit(X%*%Par[1:ncol(X)])
      #1) probability of d1, d2 given detected
      Pdot[ibin1,ibin2,] = (1-(1-P[ibin1,ibin2,1,]*Pin[ibin1])*(1-P[ibin1,ibin2,2,]*Pin[ibin2]))
      P.dist[ibin1,ibin2,Which.move] = Dist.probs[ibin1]*Psi[ibin1,ibin2]
      P.dist[ibin1,ibin2,-Which.move] = Dist.probs[ibin1]*Psi.nomove[ibin1,ibin2]
      P.dist[ibin1,ibin2,Which.11] = P.dist[ibin1,ibin2,Which.11]*P[ibin1,ibin2,1,Which.11]*P[ibin1,ibin2,2,Which.11]*Meas[ibin1,Obs.dists[1,Which.11]]*Meas[ibin2,Obs.dists[2,Which.11]]
      #components of Pr(10 hists)
      P.dist[ibin1,ibin2,Which.10] = P.dist[ibin1,ibin2,Which.10]*P[ibin1,ibin2,1,Which.10]*(1-P[ibin1,ibin2,2,Which.10]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10]]
      #components of Pr(01 hists)
      P.dist[ibin1,ibin2,Which.01] = P.dist[ibin1,ibin2,Which.01]*P[ibin1,ibin2,2,Which.01]*(1-P[ibin1,ibin2,1,Which.01]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01]]
    }
  }
  N.hat = rep(0,n.hists)
  for(ihist in 1:n.hists){
    P.dist[,,ihist]=P.dist[,,ihist]/sum(P.dist[,,ihist])
    N.hat[ihist] = sum(G[ihist] * Count[ihist] * P.dist[c(1:n.obs.bins),,ihist] / Pdot[c(1:n.obs.bins),,ihist])
  }
  n.hat = sum(N.hat)
  n.hat
}
