#  integrated likelihood approach for movement

#simulate data
n.species=1
n.bins=5
n.obs.bins=5
source('./MRDSmove/R/simulate_mrds.R')
Sim_data <- simulate_mrds(n_species=n.species,n_bins=n.bins,n_obs_bins=n.obs.bins,measure_par=0.6,move_par=c(0.01,0.01),seed=12345,gaussian=TRUE)
Obs_data=Sim_data$Obs_data

expit<-function(x)(1/(1+exp(-x)))

#' Estimate abundance from MRDS double observer data subject to movement and measurement error
#' @param Parameter vector, including detection parameters, movement error SDs (left and right tail), measurement error SD
#' @param Data A design.matrix with the following column names: "match" indicates which records match with which (there should be two records
#'        for each detection, one for each observer), "observer","species" (provides species or other grouping variable: abundance estimates will be provided separately for each), 
#'       "obs.dist" (observed distance; NA if missing), "g_size" (group/cluster size), "moving" (binary indicator for moving/not moving), 
#'       "detected" (binary detection/nondetection). Additional covariates may also be provided and used in formula (e.g. "other"; see below)
#' @param mod.formula Formula object giving formula for detection probability for each observer.  Variables can be linked to column names in Data, 
#'        but need to use "distance" to represent distance and "distance2" to represent distance^2.  Adding in a variable "other" into Data (whether the other observer detected it or not)
#'        can be used to implement symmetric point/limiting independence as described in discussion of MacKenzie and Clement 2016.  Note that the
#'        interaction between "other" and "distance" can be used to implement point indendence.
#' @param Bin.widths Vector of distance bin widths
#' @param Obs.bins A vector giving which bins are observed (e.g. 1:3 if bins 1-3 are observed)
#' @param gaussian If TRUE, uses a Gaussian distribution for measurement error and a half-normal for movement; if FALSE (default), uses exponential / half exponential (Laplace dist)
#' @return a community mrds dataset
#' @export
#' @keywords simulation, mrds
#' @author Paul B. Conn
MRDSmove_IntLik <- function(Par,Data,mod.formula,Bin.widths,Obs.bins,gaussian=FALSE){
  #cat(Par)
  Temp=Bin.midpoints=0*Bin.widths
  n.bins = length(Bin.widths)
  n.hists= nrow(Data)/2
  n.par=length(Par)
  if(n.bins<2)cat('ERROR: function not set up for <2 distance bins')
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

  #Dist = Obs.dists[1,]
  #Dist[is.na(Dist)]=Obs.dists[2,is.na(Dist)]
  
  #categorize histories into 11, 10, 01
  Detected = matrix(Data$detected,2,n.hists)
  Move = Data$moving[2*(1:n.hists)-1]
  Hist.type = Detected[1,]-Detected[2,]  #so 01 is a -1, 11 is a 0, 10 is a 1
  Which.move = which(Move==1)
  Which.11 = which(Hist.type==0)
  Which.01 = which(Hist.type==-1)
  Which.10 = which(Hist.type==1)
  # Which.11.move = which(Hist.type==0 & Move==1)
  # Which.01.move = which(Hist.type==-1 & Move==1)
  # Which.10.move = which(Hist.type==1 & Move==1)
  # Which.11.nomove = which(Hist.type==0 & Move==0)
  # Which.01.nomove = which(Hist.type==-1 & Move==0)
  # Which.10.nomove = which(Hist.type==1 & Move==0)
  
  #parameterize Psi transition matrix, Measurement error matrix
  Psi = Meas = matrix(0,n.bins,n.bins)
  if(gaussian==TRUE){
    for(ibin1 in 1:n.bins){
      Psi[ibin1,ibin1:n.bins]=dnorm(Bin.midpoints[ibin1:n.bins],Bin.midpoints[ibin1],exp(Par[n.par-2]))
      Meas[ibin1,] = dnorm(Bin.midpoints,Bin.midpoints[ibin1],exp(Par[n.par]))
    }
    for(ibin1 in 2:n.bins)Psi[ibin1,1:(ibin1-1)]=dnorm(Bin.midpoints[1:(ibin1-1)],Bin.midpoints[ibin1],exp(Par[n.par-1]))
  }
  if(gaussian==FALSE){  #double exponential
    for(ibin1 in 1:n.bins){
      Psi[ibin1,ibin1:n.bins]=dexp(Bin.midpoints[ibin1:n.bins]-Bin.midpoints[ibin1],exp(Par[n.par-2]))
      Meas[ibin1,] = dexp(abs(Bin.midpoints-Bin.midpoints[ibin1]),exp(Par[n.par]))
    }
    for(ibin1 in 2:n.bins)Psi[ibin1,1:(ibin1-1)]=dexp(Bin.midpoints[ibin1]-Bin.midpoints[1:(ibin1-1)],exp(Par[n.par-1]))
  }
  Psi = Psi/rowSums(Psi)
  Psi.nomove = diag(n.bins)
  Meas = Meas/rowSums(Meas)
  Pin = rowSums(Meas[,Obs.bins])  #probability of being observed in the strip width given location
  
  P.hist = rep(0,n.hists)
  P = array(0,dim=c(n.bins,n.bins,2,n.hists))  #to hold (conditional) detection probabilities for each observer
  Pdot = rep(0,n.hists)

  # #1) calculate probability of being in d1, d2 given detected
   P.dist = array(1,dim=c(n.bins,n.bins,n.hists))  
   for(ibin1 in 1:n.bins){  #observer 1 latent distance
     for(ibin2 in 1:n.bins){  #observer 2 latent distance
       Data$distance = rep(c(Dists[ibin1],Dists[ibin2]),n.hists)
       #Data$distance = rep(c(Bin.midpoints[ibin1],Bin.midpoints[ibin2]),n.hists)
       Data$distance2 = Data$distance^2
       X = model.matrix(mod.formula,Data)
       P[ibin1,ibin2,,]=expit(X%*%Par[1:ncol(X)])
       #1) probability of d1, d2 given detected
       #P.dist[ibin1,ibin2,Which.move] = (1-(1-P[ibin1,ibin2,1,Which.move]*Pin[ibin1])*(1-P[ibin1,ibin2,2,Which.move]*Pin[ibin2]))
       #P.dist[ibin1,ibin2,Which.move] = P.dist[ibin1,ibin2,Which.move]*Psi[ibin1,ibin2]
       #if(ibin1==ibin2)P.dist[ibin1,ibin2,-Which.move] = (1-(1-P[ibin1,ibin2,1,-Which.move]*Pin[ibin1])*(1-P[ibin1,ibin2,2,-Which.move]*Pin[ibin2]))
       P.dist[ibin1,ibin2,Which.move] = Psi[ibin1,ibin2]
       P.dist[ibin1,ibin2,-Which.move] = Psi.nomove[ibin1,ibin2]
       P.dist[ibin1,ibin2,Which.11] = P.dist[ibin1,ibin2,Which.11]*P[ibin1,ibin2,1,Which.11]*P[ibin1,ibin2,2,Which.11]*Meas[ibin1,Obs.dists[1,Which.11]]*Meas[ibin2,Obs.dists[2,Which.11]]
       #components of Pr(10 hists)
       P.dist[ibin1,ibin2,Which.10] = P.dist[ibin1,ibin2,Which.10]*P[ibin1,ibin2,1,Which.10]*(1-P[ibin1,ibin2,2,Which.10]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10]]
       #components of Pr(01 hists)
       P.dist[ibin1,ibin2,Which.01] = P.dist[ibin1,ibin2,Which.01]*P[ibin1,ibin2,2,Which.01]*(1-P[ibin1,ibin2,1,Which.01]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01]]
     }
   }
   for(ihist in 1:n.hists){
     P.dist[,,ihist]=P.dist[,,ihist]/sum(P.dist[,,ihist])
   }
   
   for(ibin1 in 1:n.bins){  #observer 1 latent distance
     for(ibin2 in 1:n.bins){  #observer 2 latent distance
       #2) probability of being detected at least once
       #Pdot=Pdot+ P.dist[ibin1,ibin2,]*(1-(1-Pin[ibin1]*P[ibin1,ibin2,1,])*(1-Pin[ibin2]*P[ibin1,ibin2,2,]))
       P.dot = (1-(1-Pin[ibin1]*P[ibin1,ibin2,1,])*(1-Pin[ibin2]*P[ibin1,ibin2,2,]))
       #probability of Pr(11 hists)
       P.hist[Which.11] = P.hist[Which.11]+P.dist[ibin1,ibin2,Which.11]*P[ibin1,ibin2,1,Which.11]*P[ibin1,ibin2,2,Which.11]*Meas[ibin1,Obs.dists[1,Which.11]]*Meas[ibin2,Obs.dists[2,Which.11]]/P.dot[Which.11]
       #components of Pr(10 hists)
       P.hist[Which.10] = P.hist[Which.10]+P.dist[ibin1,ibin2,Which.10]*P[ibin1,ibin2,1,Which.10]*(1-P[ibin1,ibin2,2,Which.10]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10]]/P.dot[Which.10]
       #components of Pr(01 hists)
       P.hist[Which.01] = P.hist[Which.01]+P.dist[ibin1,ibin2,Which.01]*P[ibin1,ibin2,2,Which.01]*(1-P[ibin1,ibin2,1,Which.01]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01]]/P.dot[Which.01]
     }
   }
  
   #straight huggins-alho - works!!
   # Data$distance = rep(Dist,each=2)
   # Data$distance2 = Data$distance^2
   # X = model.matrix(mod.formula,Data)
   # P = matrix(0,2,n.hists)  #to hold (conditional) detection probabilities for each observer
   # P[]=expit(X%*%Par[1:ncol(X)])
   # P.dot = (1-(1-P[1,])*(1-P[2,]))
   # #probability of Pr(11 hists)
   # P.hist[Which.11] = P[1,Which.11]*P[2,Which.11]/P.dot[Which.11]
   # # #components of Pr(10 hists)
   # P.hist[Which.10] = P[1,Which.10]*(1-P[2,Which.10])/P.dot[Which.10]
   # # #components of Pr(01 hists)
   # P.hist[Which.01] = P[2,Which.01]*(1-P[1,Which.01])/P.dot[Which.01]
   

  #Dist.probs=rep(1,5)    
  # for(ibin1 in 1:n.bins){  #observer 1 latent distance
  #   for(ibin2 in 1:n.bins){  #observer 2 latent distance
  #     Data$distance = rep(c(Bin.midpoints[ibin1],Bin.midpoints[ibin2]),n.hists)
  #     Data$distance2 = Data$distance^2
  #     X = model.matrix(mod.formula,Data)
  #     P[] = expit(X%*%Par[1:ncol(X)])
  #     
  #     # #compute components of log(probability) of being seen at least once (need to sum over s1, s2)
  #     # Pdot[Which.move]=Pdot[Which.move]+Dist.probs[ibin1]*Psi[ibin1,ibin2]*(1-(1-Pin[ibin1]*P[1,Which.move])*(1-Pin[ibin2]*P[2,Which.move]))
  #     # if(ibin1==ibin2)Pdot[-Which.move]=Pdot[-Which.move]+Dist.probs[ibin1]*(1-(1-Pin[ibin1]*P[1,-Which.move])*(1-Pin[ibin2]*P[2,-Which.move]))
  #     # #probability of Pr(11 hists)
  #     # P.hist[Which.11.move] = P.hist[Which.11.move]+Dist.probs[ibin1]*Psi[ibin1,ibin2]*P[1,Which.11.move]*P[2,Which.11.move]*Meas[ibin1,Obs.dists[1,Which.11.move]]*Meas[ibin2,Obs.dists[2,Which.11.move]]
  #     # if(ibin1==ibin2)P.hist[Which.11.nomove] = P.hist[Which.11.nomove]+Dist.probs[ibin1]*P[1,Which.11.nomove]*P[2,Which.11.nomove]*Meas[ibin1,Obs.dists[1,Which.11.nomove]]*Meas[ibin2,Obs.dists[2,Which.11.nomove]]
  #     # #components of Pr(10 hists)
  #     # P.hist[Which.10.move] = P.hist[Which.10.move]+Dist.probs[ibin1]*Psi[ibin1,ibin2]*P[1,Which.10.move]*(1-P[2,Which.10.move]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10.move]]
  #     # if(ibin1==ibin2)P.hist[Which.10.nomove] = P.hist[Which.10.nomove]+Dist.probs[ibin1]*P[1,Which.10.nomove]*(1-P[2,Which.10.nomove]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10.nomove]]
  #     # #components of Pr(01 hists)
  #     # P.hist[Which.01.move] = P.hist[Which.01.move]+Dist.probs[ibin1]*Psi[ibin1,ibin2]*P[2,Which.01.move]*(1-P[1,Which.01.move]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01.move]]
  #     # if(ibin1==ibin2)P.hist[Which.01.nomove] = P.hist[Which.01.nomove]+Dist.probs[ibin1]*P[2,Which.01.nomove]*(1-P[1,Which.01.nomove]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01.nomove]]
  # 
  #     # #compute components of log(probability) of being seen at least once (need to sum over s1, s2)
  #     # Pdot[Which.move]=Pdot[Which.move]+ P.init[ibin1,1,Which.move]*Psi[ibin1,ibin2]*(1-(1-Pin[ibin1]*P[1,Which.move])*(1-Pin[ibin2]*P[2,Which.move]))
  #     # if(ibin1==ibin2)Pdot[-Which.move]=Pdot[-Which.move]+P.init[ibin1,2,-Which.move]*(1-(1-Pin[ibin1]*P[1,-Which.move])*(1-Pin[ibin2]*P[2,-Which.move]))
  #     # #probability of Pr(11 hists)
  #     # P.hist[Which.11.move] = P.hist[Which.11.move]+P.init[ibin1,1,Which.11.move]*Psi[ibin1,ibin2]*P[1,Which.11.move]*P[2,Which.11.move]*Meas[ibin1,Obs.dists[1,Which.11.move]]*Meas[ibin2,Obs.dists[2,Which.11.move]]
  #     # if(ibin1==ibin2)P.hist[Which.11.nomove] = P.hist[Which.11.nomove]+P.init[ibin1,2,Which.11.nomove]*P[1,Which.11.nomove]*P[2,Which.11.nomove]*Meas[ibin1,Obs.dists[1,Which.11.nomove]]*Meas[ibin2,Obs.dists[2,Which.11.nomove]]
  #     # #components of Pr(10 hists)
  #     # P.hist[Which.10.move] = P.hist[Which.10.move]+P.init[ibin1,1,Which.10.move]*Psi[ibin1,ibin2]*P[1,Which.10.move]*(1-P[2,Which.10.move]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10.move]]
  #     # if(ibin1==ibin2)P.hist[Which.10.nomove] = P.hist[Which.10.nomove]+P.init[ibin1,2,Which.10.nomove]*P[1,Which.10.nomove]*(1-P[2,Which.10.nomove]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10.nomove]]
  #     # #components of Pr(01 hists)
  #     # P.hist[Which.01.move] = P.hist[Which.01.move]+P.init[ibin1,1,Which.01.move]*Psi[ibin1,ibin2]*P[2,Which.01.move]*(1-P[1,Which.01.move]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01.move]]
  #     # if(ibin1==ibin2)P.hist[Which.01.nomove] = P.hist[Which.01.nomove]+P.init[ibin1,2,Which.01.nomove]*P[2,Which.01.nomove]*(1-P[1,Which.01.nomove]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01.nomove]]
  # 
  #     #compute components of log(probability) of being seen at least once (need to sum over s1, s2)
  #     
  #     Pdot[Which.move]=Pdot[Which.move]+ P.dist[ibin1,ibin2,Which.move]*(1-(1-Pin[ibin1]*P[1,Which.move])*(1-Pin[ibin2]*P[2,Which.move]))
  #     if(ibin1==ibin2)Pdot[-Which.move]=Pdot[-Which.move]+P.dist[ibin1,2,-Which.move]*(1-(1-Pin[ibin1]*P[1,-Which.move])*(1-Pin[ibin2]*P[2,-Which.move]))
  #     #probability of Pr(11 hists)
  #     P.hist[Which.11.move] = P.hist[Which.11.move]+P.init[ibin1,1,Which.11.move]*Psi[ibin1,ibin2]*P[1,Which.11.move]*P[2,Which.11.move]*Meas[ibin1,Obs.dists[1,Which.11.move]]*Meas[ibin2,Obs.dists[2,Which.11.move]]
  #     if(ibin1==ibin2)P.hist[Which.11.nomove] = P.hist[Which.11.nomove]+P.init[ibin1,2,Which.11.nomove]*P[1,Which.11.nomove]*P[2,Which.11.nomove]*Meas[ibin1,Obs.dists[1,Which.11.nomove]]*Meas[ibin2,Obs.dists[2,Which.11.nomove]]
  #     #components of Pr(10 hists)
  #     P.hist[Which.10.move] = P.hist[Which.10.move]+P.init[ibin1,1,Which.10.move]*Psi[ibin1,ibin2]*P[1,Which.10.move]*(1-P[2,Which.10.move]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10.move]]
  #     if(ibin1==ibin2)P.hist[Which.10.nomove] = P.hist[Which.10.nomove]+P.init[ibin1,2,Which.10.nomove]*P[1,Which.10.nomove]*(1-P[2,Which.10.nomove]*sum(Meas[ibin2,Obs.bins]))*Meas[ibin1,Obs.dists[1,Which.10.nomove]]
  #     #components of Pr(01 hists)
  #     P.hist[Which.01.move] = P.hist[Which.01.move]+P.init[ibin1,1,Which.01.move]*Psi[ibin1,ibin2]*P[2,Which.01.move]*(1-P[1,Which.01.move]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01.move]]
  #     if(ibin1==ibin2)P.hist[Which.01.nomove] = P.hist[Which.01.nomove]+P.init[ibin1,2,Which.01.nomove]*P[2,Which.01.nomove]*(1-P[1,Which.01.nomove]*sum(Meas[ibin1,Obs.bins]))*Meas[ibin2,Obs.dists[2,Which.01.nomove]]
  #     
  #   }
  # }
  #if(sum(P.hist<=0 | P.hist>1 | Pdot<=0 | Pdot>1 | Pdot<P.hist))cat("\n BUG detected!!, probability of history is not in (0,1) \n")  
  
  log.lik = sum(log(P.hist)) #-log(Pdot))
  -log.lik
}

####   run through Huggins-Alho in MARK
Dat_MARK = data.frame(Obs_data)
Dat_MARK$ch = paste0(as.character(Dat_MARK$det1),as.character(Dat_MARK$det2))
Dat_MARK$dist = Dat_MARK$d1_obs
Dat_MARK$dist[is.na(Dat_MARK$d1_obs)]=Dat_MARK$d2_obs[is.na(Dat_MARK$d1_obs)]
Dat_MARK$dist2 = Dat_MARK$dist^2

library('RMark')
HA_proc = process.data(Dat_MARK,model="Huggins")
HA_ddl = make.design.data(HA_proc)
p_model = list(formula=~dist+dist2+time+g_size+fly,share=TRUE)
HA_mod= mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))
#####

#convert simulated data into input format for MRDSmove
n.hists=nrow(Obs_data)
Obs_data=data.frame(Obs_data)
Data = data.frame(match=rep(c(1:n.hists),each=2),observer=rep(c(0,1),n.hists),species=rep(Obs_data$species,each=2),
                     obs.dist=as.vector(rbind(Obs_data$d1_obs,Obs_data$d2_obs)),g_size=rep(Obs_data$g_size,each=2),
                     moving = rep(Obs_data$fly,each=2),detected=as.vector(rbind(Obs_data$det1,Obs_data$det2)))

my.formula=~distance+distance2+observer+g_size+moving
my.par = c(0,.03,-.04,-.5,.02,.5,-4,-4,-.51)
my.par = c(0,.03,-.04,-.5,.02,.5,-4,-4,-.51)

MRDSmove_IntLik(Par=my.par,Data,mod.formula=my.formula,Bin.widths=rep(1,n.bins),Obs.bins=c(1:n.obs.bins),gaussian=TRUE)

glm_out = optim(par=my.par,MRDSmove_IntLik,hessian=TRUE,method="BFGS",Data=Data,mod.formula=my.formula,Bin.widths=rep(1,n.bins),Obs.bins=c(1:n.obs.bins),gaussian=TRUE)
solve(glm_out$hessian)
glm_out$par
sqrt(diag(solve(glm_out$hessian[c(1:6,9),c(1:6,9)])))

DM.pred = matrix(1,5,6)
DM.pred[,2]=c(1:5)
DM.pred[,3]=DM.pred[,2]^2
DM.pred[,4]=0

cat('\n real scale true values:')
expit(DM.pred%*%my.par[1:6])
cat('\n real estimates, my code:')
expit(DM.pred%*%glm_out$par[1:6])
cat('\n real estimates, Rmark:')
expit(DM.pred%*%HA_mod$results$beta$estimate)