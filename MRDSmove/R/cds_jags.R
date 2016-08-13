mrds.cds.jags <- function(){
  for (isp in 1:n.species) {
    for(iind in 1:(M[isp])){
      #detection model
      Det1[isp,iind] ~ dbern(P.obs1[isp,iind]*Z[isp,iind])
      P.obs1[isp,iind] <- (.0000001+exp(-(Dist1[isp,iind]-1)*(Dist1[isp,iind]-1)/(2*Sigma1[isp,iind]*Sigma1[isp,iind])))*0.99999  #include 0.99999 because numerical values = 1.0 can be problematic for dbern
      Sigma1[isp,iind] <- exp(beta.det.0+beta.det.sp[isp]+beta.det.fly*Fly[isp,iind]+beta.det.group*Group[isp,iind])
      #distance model
      Dist1[isp,iind] ~ dcat(BinWidth)  

      #proportion flying model
      Fly[isp,iind] ~ dbern(Fly.sp.expit[isp])
      
      #abundance
      ZG[isp,iind]<-Z[isp,iind]*Group[isp,iind]      
      Z[isp,iind] ~ dbern(Psi[isp])
      
      #cluster/group model
      Group[isp,iind] <- Group.min1[isp,iind] + 1
      Group.min1[isp,iind] ~ dpois(exp(Mu.grp.ind[isp,iind]))
      Mu.grp.ind[isp,iind] ~ dnorm(Mu.grp[isp],tau.grp.ind.exp[isp])
      
    }
    #abundance model
    #Psi[isp] ~ dbeta(0.01,1)  #Link (2013) scale prior approx
    Psi[isp] ~ dbeta(1.0,1.0)
    G[isp] <- sum(Z[isp,1:(M[isp])])
    N[isp] <- sum(ZG[isp,1:(M[isp])])
    
    #cluster/group size model
    Mu.grp[isp] ~ dnorm(mu.grp,tau.grp)
    tau.grp.ind.exp[isp] <- exp(tau.grp.ind[isp])
    tau.grp.ind[isp] ~ dnorm(tau.grp.mu,10000)   #dnorm(tau.grp.mu,0.1)
    #proportion flying model
    Fly.sp.expit[isp] <- 1/(1+exp(-Fly.sp[isp]))
    Fly.sp[isp] ~ dnorm(fly.mean,tau.fly)
    
    beta.det.sp[isp]~ dnorm(0.0,10)
    
  }
 
  
  #More priors
  mu.grp ~ dnorm(0,0.01)
  tau.grp ~ dgamma(1.0,0.01)
  tau.grp.mu ~ dgamma(1.0,0.01)
  fly.mean ~ dnorm(0,0.01)
  tau.fly ~ dgamma(1.0,0.01)
  tau.sp.det ~ dgamma(1.0,0.01)
  beta.det.0 ~ dnorm(0.0,10)
  beta.det.obs ~ dnorm(0.0,10)
  beta.det.fly~ dnorm(0.0,10)
  beta.det.group~ dnorm(0.0,10)
}

n.species=10
n.bins=9
n.obs.bins=5
source('./MRDSmove/R/simulate_mrds.R')
Data <- simulate_mrds(n_species=n.species,n_bins=n.bins,n_obs_bins=n.obs.bins,seed=12345,p0=FALSE)
Obs1_data=Data$Obs_data[which(Data$Obs_data[,"det1"]==1),]

M=500
Det1 = matrix(0,n.species,M)
Group = Fly = Z = Dist1 =  matrix(NA,n.species,M)
Nobs = rep(0,n.species)
for(isp in 1:n.species){
  Nobs[isp]=length(which(Obs1_data[,"species"]==isp))
  Cur_data=Obs1_data[which(Obs1_data[,"species"]==isp),]
  Z[isp,1:nrow(Cur_data)]=1
  Group[isp,1:nrow(Cur_data)]=Cur_data[,"g_size"]
  Fly[isp,1:nrow(Cur_data)]=Cur_data[,"fly"]
  Det1[isp,1:nrow(Cur_data)]=Cur_data[,"det1"]
  Dist1[isp,1:nrow(Cur_data)]=Cur_data[,"d1_obs"]
}





library(rjags)
library(R2jags)
set.factory("bugs::Conjugate",FALSE, type="sampler")  #FALSE recommended by R. Sollman 
set.factory("bugs::Slice",TRUE, type="sampler") 
BinWidth=rep(1,n.bins)
BinVal=c(1:5)
ObsBins1.ext = c(1,1,1,1,1,0,0,0,0)
ObsBins2.ext = c(1,1,1,1,1,0)  #probabilities consolidated to one "unobservable" bin at end


M=Nobs*5
n_iter=1000
n_burnin=500
n_thin=1
n_chains=1
first.bin=1
Small=rep(0.00001,n.obs.bins+1)
Group.min1 = Group-1
jags_data = list("Dist1","Z","n.species","n.bins","n.obs.bins","Group.min1","Fly","Det1","BinWidth","BinVal","ObsBins1.ext","M","Small","first.bin")
jags_params = c("Psi","beta.det.0","beta.det.fly","beta.det.group","beta.det.sp",
                "mu.grp","tau.grp","tau.grp.mu","fly.mean","Fly.sp","tau.fly","tau.sp.det",
                "Mu.grp")
jags_save =c("G","N","Psi","beta.det.0","beta.det.obs","beta.det.fly","beta.det.group","Z[1,68]","Dist1[1,68]","Dist2.obs[1,68]","P.obs1[1,68]","P.obs2[1,68]","Group[1,68]")

jags.inits = function(){
  list("Psi"=runif(n.species,0.2,0.4),
       "beta.det.0"=rnorm(1,log(2),.1),"beta.det.fly"=rnorm(1,0,.1),"beta.det.group"=rnorm(1,0,.1),"beta.det.sp"=rnorm(n.species,0,.1),
       "mu.grp" = rpois(1,3), "tau.grp"=runif(1,0.5,5.0), "tau.grp.mu"=runif(1,0.5,5.0),"fly.mean"=runif(1,0.2,0.8),"Fly.sp"=runif(n.species,0.2,0.8),"tau.fly"=runif(1,0.5,5.0),"tau.sp.det"=runif(1,0.5,5.0),
       "Mu.grp"=rpois(n.species,3))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,
                n.iter=n_iter,
                model.file=mrds.cds.jags,
                DIC=FALSE,
                parameters.to.save=jags_save,
                n.chains=n_chains,
                n.burnin=n_burnin,
                n.thin=n_thin,
                working.directory=getwd())

jags_fit$BUGSoutput$sims.matrix[1,]