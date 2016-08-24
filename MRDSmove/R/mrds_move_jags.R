mrds.move.jags <- function(){
  for (isp in 1:n.species) {
    for(iind in 1:(M[isp])){
      #detection model
      Det1[isp,iind] ~ dbern(P.obs1[isp,iind]*Z[isp,iind]*I.obs1[isp,iind])
      Det2[isp,iind] ~ dbern(P.obs2[isp,iind]*Z[isp,iind]*I.obs2[isp,iind])
      #crap[isp,iind] <- P.obs1[isp,iind]*Z[isp,iind]*I.obs1[isp,iind]
      P.obs1[isp,iind] <- (.0000001+exp(-(Dist1.true[isp,iind]-1)*(Dist1.true[isp,iind]-1)/(2*Sigma1[isp,iind]*Sigma1[isp,iind])))*0.99999  #include 0.99999 because numerical values = 1.0 can be problematic for dbern
      P.obs2[isp,iind] <- (.0000001+exp(-(Dist2.obs[isp,iind]-1)*(Dist2.obs[isp,iind]-1)/(2*Sigma2[isp,iind]*Sigma2[isp,iind])))*0.99999
      Sigma1[isp,iind] <- exp(beta.det.0+beta.det.sp[isp]+beta.det.fly*Fly[isp,iind]+beta.det.group*Group[isp,iind])
      Sigma2[isp,iind] <- exp(beta.det.0+beta.det.obs+beta.det.sp[isp]+beta.det.fly*Fly[isp,iind]+beta.det.group*Group[isp,iind])
      I.obs2[isp,iind] <- ObsBins2.ext[Dist2.obs[isp,iind]] #perceived distance within strip width?
      I.obs1[isp,iind] <- ObsBins1.ext[Dist1[isp,iind]]  #perceived distance within strip width?
      
      #P.obs1[isp,iind] <- P0.obs1[isp,iind]*dnorm1[isp,iind]/dnorm(0,0,Tau.obs1[isp,iind]) 
      #dnorm1[isp,iind] <- dnorm(Dist1.true[isp,iind]-1,0.0,Tau.obs1[isp,iind])
      #P.obs2[isp,iind] <- P0.obs2[isp,iind]*dnorm(Dist2.true[isp,iind]-1,0,Tau.obs2[isp,iind])/dnorm(0,0,Tau.obs2[isp,iind]) 
      #Tau.obs1[isp,iind]=1/exp(beta.det.0+beta.det.sp[isp]+beta.det.fly*Fly[isp,iind]+beta.det.group*Group[isp,iind])
      #Tau.obs2[isp,iind]=1/exp(beta.det.0+beta.det.obs+beta.det.sp[isp]+beta.det.fly*Fly[isp,iind]+beta.det.group*Group[isp,iind])
      #P0.obs1[isp,iind]=1/(1+exp(-(beta.p0.0+beta.p0.sp[isp]+beta.p0.fly*Fly[isp,iind]+beta.p0.group*Group[isp,iind])))  
      #P0.obs2[isp,iind]=1/(1+exp(-(beta.p0.0+beta.p0.obs+beta.p0.sp[isp]+beta.p0.fly*Fly[isp,iind]+beta.p0.group*Group[isp,iind])))  
      #distance model
      Dist2.obs[isp,iind] ~ dcat(BinProb.obs2.adj[isp,iind,1:(n.obs.bins+1)])
      BinProb.obs2.adj[isp,iind,1:(n.obs.bins+1)]<-BinProb.obs2[isp,iind,1:(n.obs.bins+1)]+Small #Small needed to prevent probabilities identically zero (dcat throws an error otherwise)
      BinProb.obs2[isp,iind,n.obs.bins+1] <- 1-sum(BinProb.obs2[isp,iind,1:n.obs.bins])  #n.bins + 1 is equivalent to 'observed distance not in strip'
      BinProb.obs2[isp,iind,1:n.obs.bins] <- up[isp,iind,1:n.obs.bins]-low[isp,iind,1:n.obs.bins]
      up[isp,iind,1:n.obs.bins]=pnorm(BinVal+0.5*BinWidth[1:n.obs.bins],Dist2.true[isp,iind],tau.measure)
      low[isp,iind,1]=0  #larger probability of being seen in first bin because assuming not possible to have measurement error that crosses the center line
      low[isp,iind,2:n.obs.bins]=pnorm(BinVal[2:n.obs.bins]-0.5*BinWidth[2:n.obs.bins],Dist2.true[isp,iind],tau.measure)
      
      Dist2.true[isp,iind] ~ dcat(BinProb.true[isp,iind,1:n.bins])
      BinProb.true[isp,iind,1:n.bins] <- dnorm(Dist1.true[isp,iind],BinVal2[Dist1[isp,iind],],Tau.move[isp,iind])  #need to update to use Rate
      Dist1.true[isp,iind]<-Dist1[isp,iind]-first.bin+1  #first distance bin here is '1'
      Dist1[isp,iind] ~ dcat(BinWidth/sum(BinWidth))  
      Tau.move[isp,iind]=1000000*(1-Fly[isp,iind])+tau.move*Fly[isp,iind]
       
      #proportion flying model
      Fly[isp,iind] ~ dbern(Fly.sp.expit[isp])
      
      #abundance
      ZG[isp,iind]<-Z[isp,iind]*Group[isp,iind]*I.obs1[isp,iind]  #only counting individuals within strip width
      Z.in[isp,iind]<-Z[isp,iind]*I.obs1[isp,iind]
      Z[isp,iind] ~ dbern(Psi[isp])
      
      #cluster/group model
      Group[isp,iind] <- Group.min1[isp,iind] + 1
      Group.min1[isp,iind] ~ dpois(Mu.grp.ind[isp,iind])
      Mu.grp.ind[isp,iind] <- exp(Mu.grp[isp]+RE.grp[isp,iind])
    }
    #abundance model
    #Psi[isp] ~ dbeta(0.01,1)  #Link (2013) scale prior approx
    Psi[isp] ~ dbeta(1.0,1.0)
    G[isp] <- sum(Z.in[isp,1:(M[isp])])
    N[isp] <- sum(ZG[isp,1:(M[isp])])
    
    #cluster/group size model
    Mu.grp[isp] ~ dnorm(mu.grp,tau.grp)
    tau.grp.ind.exp[isp] <- exp(tau.grp.ind[isp])
    tau.grp.ind[isp] ~ dnorm(tau.grp.mu,10000)   #dnorm(tau.grp.mu,0.1)
    #proportion flying model
    Fly.sp.expit[isp] <- 1/(1+exp(-Fly.sp[isp]))
    Fly.sp[isp] ~ dnorm(fly.mean,tau.fly)
    #movement
    Rate.exp[isp]= exp(Rate[isp])
    Rate[isp] ~ dnorm(rate.mu,tau.move)
    beta.det.sp[isp] ~ dnorm(0.0,tau.sp.det)  
    beta.p0.sp[isp] ~ dnorm(0.0,tau.sp.p0)
  }
 
  for(isp2 in 1:n.species){
    for(iind2 in 1:M.max){
      RE.grp[isp2,iind2] ~ dnorm(0,tau.grp.ind.exp[isp2])
    }
  }
  
  #More priors
  mu.grp ~ dnorm(0,0.01)
  tau.grp ~ dgamma(1.0,0.01)
  tau.grp.mu ~ dgamma(1.0,0.01)
  fly.mean ~ dnorm(0,0.01)
  tau.fly ~ dgamma(1.0,0.01)
  tau.sp.det ~ dgamma(1.0,0.01)
  tau.sp.p0 ~ dgamma(1.0,0.01)
  rate.mu ~ dnorm(0,0.01)
  tau.move ~ dgamma(1.0,0.01) 
  tau.measure ~ dgamma(1.0,0.01)
  beta.det.0 ~ dnorm(0.0,10)
  beta.det.obs ~ dnorm(0.0,10)
  beta.det.fly~ dnorm(0.0,10)
  beta.det.group~ dnorm(0.0,10)
  #beta.det.0 ~ dnorm(0.0,0.01)
  #beta.det.obs ~ dnorm(0.0,0.01)
  #beta.det.fly~ dnorm(0.0,0.01)
  #beta.det.group~ dnorm(0.0,0.01)
  beta.p0.0 ~ dnorm(0.0,0.01)
  beta.p0.obs ~ dnorm(0.0,0.01)
  beta.p0.fly~ dnorm(0.0,0.01)
  beta.p0.group~ dnorm(0.0,0.01)
}

n.species=10
n.bins=9
n.obs.bins=5
source('./MRDSmove/R/simulate_mrds.R')
Data <- simulate_mrds(n_species=n.species,n_bins=n.bins,n_obs_bins=n.obs.bins,seed=12345,p0=FALSE)
M=M.max=500
Det1 = Det2 =  matrix(0,n.species,M)
Group = Fly = Z = Dist1 = Dist2.obs = Dist2.true = matrix(NA,n.species,M)
Nobs = rep(0,n.species)
for(isp in 1:n.species){
  Nobs[isp]=length(which(Data$Obs_data[,"species"]==isp))
  Cur_data=Data$Obs_data[which(Data$Obs_data[,"species"]==isp),]
  Z[isp,1:nrow(Cur_data)]=1
  Group[isp,1:nrow(Cur_data)]=Cur_data[,"g_size"]
  #Group[isp,(nrow(Cur_data)+1):M]=sample(Cur_data[,"g_size"],M-nrow(Cur_data),replace=TRUE)
  Fly[isp,1:nrow(Cur_data)]=Cur_data[,"fly"]
  #Fly[isp,(nrow(Cur_data)+1):M]=sample(Cur_data[,"fly"],M-nrow(Cur_data),replace=TRUE)
  Det1[isp,1:nrow(Cur_data)]=Cur_data[,"det1"]
  Det2[isp,1:nrow(Cur_data)]=Cur_data[,"det2"]
  Dist1[isp,1:nrow(Cur_data)]=Cur_data[,"d1_obs"]
  Dist2.obs[isp,1:nrow(Cur_data)]=Cur_data[,"d2_obs"]
  #Dist1_true[isp,1:nrow(Cur_data)]=Cur_data[,"d1_obs"]
  #Missing = which(is.na(Cur_data[,"d1_obs"]))
  #if(length(Missing)>0)Dist1_true[isp,Missing]=sample(c(1:n_bins),length(Missing),replace=TRUE)
  #Missing = which(is.na(Cur_data[,"d2_obs"]))
  #if(length(Missing)>0)Dist2_true[isp,Missing]=sample(c(1:n_bins),length(Missing),replace=TRUE)
}

#change detections for observed distance > 5 or < 1 to zero
Det1[which(Dist1>n.obs.bins | Dist1<1)]=0
Det2[which(Dist2.obs>n.obs.bins | Dist2.obs<1)]=0
Dist1[which(Dist1>n.obs.bins  | Dist1<1)]=NA
Dist2.obs[which(Dist2.obs>n.obs.bins | Dist2.obs<1)]=NA
Which.0=which(Dist1==0)
if(length(Which.0)>0){
  Dist1[Which.0]=NA
  Det1[Which.0]=0
}
Which.0=which(Dist2.obs==0)
if(length(Which.0)>0){
  Dist2.obs[Which.0]=NA
  Det2[Which.0]=0
}

#remove sequences of detections that were both 0 
for(isp in n.species){
  Cur.sum=Det1[isp,1:Nobs[isp]]+Det2[isp,1:Nobs[isp]]
  Which.0 = which(Cur.sum==0)
  if(length(Which.0)>0){
    Tmp=c(1:Nobs[isp])
    Det1[isp,1:(Nobs[isp]-length(Which.0))]=Det1[isp,Tmp[-Which.0]] 
    Dist1[isp,1:(Nobs[isp]-length(Which.0))]=Dist1[isp,Tmp[-Which.0]]
    Det2[isp,1:(Nobs[isp]-length(Which.0))]=Det2[isp,Tmp[-Which.0]] 
    Dist2.obs[isp,1:(Nobs[isp]-length(Which.0))]=Dist2.obs[isp,Tmp[-Which.0]]
    Nobs[isp]=Nobs[isp]-length(Which.0)
  }
}



library(rjags)
library(R2jags)
set.factory("bugs::Conjugate",FALSE, type="sampler")  #FALSE recommended by R. Sollman 
set.factory("bugs::Slice",TRUE, type="sampler") 
BinWidth=rep(1,n.bins)
BinVal=c(1:5)
BinVal.ext = c(0:8)
ObsBins1.ext = c(1,1,1,1,1,0,0,0,0)
ObsBins2.ext = c(1,1,1,1,1,0)  #probabilities consolidated to one "unobservable" bin at end

BinVal2=matrix(-100,n.bins,n.bins)
BinVal2[1,]=c(1:n.bins)
for(ibin in 2:n.bins)BinVal2[ibin,c(ibin:n.bins)]=c(ibin:n.bins)
Dist2.true=matrix(NA,n.species,M)


M=Nobs*5
n_iter=1000
n_burnin=500
n_thin=1
n_chains=1
first.bin=1
Small=rep(0.00001,n.obs.bins+1)
Group.min1 = Group-1
jags_data = list("Dist1","Dist2.true","Z","n.species","n.bins","n.obs.bins","Group.min1","Fly","Det1","Det2","Dist2.obs","BinWidth","BinVal","BinVal.ext","ObsBins1.ext","ObsBins2.ext","BinVal2","M","Small","first.bin","M.max")
jags_params = c("Psi","beta.p0.0","beta.p0.obs","beta.p0.fly","beta.p0.group","beta.p0.sp",
                "beta.det.0","beta.det.obs","beta.det.fly","beta.det.group","beta.det.sp",
                "mu.grp","tau.grp","tau.grp.mu","fly.mean","Fly.sp","tau.fly","tau.sp.det","tau.sp.p0","rate.mu","Rate","tau.move","tau.measure",
                "Mu.grp","RE.grp")
jags_save =c("G","N","tau.move","tau.measure","Psi","beta.det.0","beta.det.obs","beta.det.fly","beta.det.group","Z[1,68]","Dist1[1,68]","Dist2.obs[1,68]","P.obs1[1,68]","P.obs2[1,68]","Group[1,68]")

jags.inits = function(){
  list("Psi"=runif(n.species,0.2,0.4),"beta.p0.0"=rnorm(1,1,.5),"beta.p0.obs"=rnorm(1,0,0.25),"beta.p0.fly"=rnorm(1,0,0.25),"beta.p0.group"=rnorm(1,0,0.25),beta.p0.sp=rnorm(n.species,0,0.25),
       "beta.det.0"=rnorm(1,log(2),.1),"beta.det.obs"=rnorm(1,0,.1),"beta.det.fly"=rnorm(1,0,.1),"beta.det.group"=rnorm(1,0,.1),beta.det.sp=rnorm(n.species,0,.1),
       "mu.grp" = rpois(1,3), "tau.grp"=runif(1,0.5,5.0), "tau.grp.mu"=runif(1,0.5,5.0),"fly.mean"=runif(1,0.2,0.8),"Fly.sp"=runif(n.species,0.2,0.8),"tau.fly"=runif(1,0.5,5.0),"tau.sp.det"=runif(1,0.5,5.0),"tau.sp.p0"=runif(1,0.5,5.0),
       "rate.mu"=runif(1,0.1,2.0),"Rate"=runif(n.species,0.1,2.0),"tau.move"=runif(1,0.5,5.0),"tau.measure"=runif(1,10.0,11.0),"Mu.grp"=rpois(n.species,3),"RE.grp"=matrix(rnorm(M.max*n.species,0,0.1),n.species,M.max))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,
                n.iter=n_iter,
                model.file=mrds.move.jags,
                DIC=FALSE,
                parameters.to.save=jags_save,
                n.chains=n_chains,
                n.burnin=n_burnin,
                n.thin=n_thin,
                working.directory=getwd())

jags_fit$BUGSoutput$sims.matrix[1,]