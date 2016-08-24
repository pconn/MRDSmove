
n.indiv=100
n.sp = n.species = 10
Group=matrix(0,n.sp,n.indiv)

log_grp_mean = log(5)
log_grp_se = log(2)

Log_grp_means = rnorm(n.sp,log_grp_mean,log_grp_se)
Grp_means = exp(Log_grp_means)

for(isp in 1:n.sp)Group[isp,]=rpois(n.indiv,Grp_means[isp])+1

Group.min1 = Group-1

group.jags <- function(){
  for (isp in 1:n.species) {
    for(iind in 1:n.indiv){
      #cluster/group model
      Group[isp,iind] <- Group.min1[isp,iind]+1
      Group.min1[isp,iind] ~ dpois(max(0,Mu.grp.ind[isp,iind]))
      #Group.min1[isp,iind] ~ dpois(max(0,Mu.grp[isp]))
      #Group.min1[isp,iind] ~ dpois(exp(Mu.grp[isp]))
      Mu.grp.ind[isp,iind] <-  exp(Mu.grp[isp] + Grp.ind.RE[isp,iind])
      Grp.ind.RE[isp,iind] ~ dnorm(0.0,tau.grp.ind)
       
    }
    #abundance model
    Mu.grp[isp] ~ dnorm(mu.grp,tau.grp)
  }

  #More priors
  mu.grp ~ dnorm(0.0,0.01)
  tau.grp ~ dgamma(1.0,0.01)
  tau.grp.ind ~ dnorm(1000,1)
  #tau.grp.ind ~ dgamma(1.0,0.01)   #dnorm(tau.grp.mu,0.1)
}

library(rjags)
library(R2jags)
set.factory("bugs::Conjugate",TRUE, type="sampler")  #FALSE recommended by R. Sollman 
#set.factory("bugs::Slice",TRUE, type="sampler") 
n_iter=1000
n_burnin=500
n_thin=1
n_chains=1

jags_data = list("n.species","Group.min1","n.indiv")
#jags_params = c("mu.grp","tau.grp","Mu.grp","Mu.grp.ind","tau.grp.ind")
jags_params = c("mu.grp","tau.grp","Mu.grp","Grp.ind.RE","tau.grp.ind")

jags_save =c("Mu.grp","tau.grp","tau.grp.ind","mu.grp","Mu.grp.ind[1,1]")

jags.inits = function(){
  list("Psi"=runif(n.species,0.2,0.4),
       "beta.det.0"=rnorm(1,log(2),.1),"beta.det.fly"=rnorm(1,0,.1),"beta.det.group"=rnorm(1,0,.1),"beta.det.sp"=rnorm(n.species,0,.1),
       "mu.grp" = rpois(1,5), "tau.grp"=runif(1,0.5,5.0),"fly.mean"=runif(1,0.2,0.8),"Fly.sp"=runif(n.species,0.2,0.8),"tau.fly"=runif(1,0.5,5.0),"tau.sp.det"=runif(1,0.5,5.0),
       "Mu.grp"=runif(n.species,5,10),"Grp.ind.RE"=matrix(rnorm(n.species*n.indiv,0,0.1),n.species,n.indiv)) 
       #"Mu.grp.ind"=matrix(rnorm(n.species*n.indiv,5,1),n.species,n.indiv))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,
                n.iter=n_iter,
                model.file=group.jags,
                DIC=FALSE,
                parameters.to.save=jags_save,
                n.chains=n_chains,
                n.burnin=n_burnin,
                n.thin=n_thin,
                working.directory=getwd())

head(jags_fit$BUGSoutput$sims.matrix)
summary(sqrt(1/jags_fit$BUGSoutput$sims.matrix[,"tau.grp"]))

Grp_means
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[1]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[2]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[3]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[4]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[5]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[6]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[7]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[8]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[9]"]))
summary(exp(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[10]"]))

Grp_means
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[1]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[2]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[3]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[4]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[5]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[6]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[7]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[8]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[9]"])
summary(jags_fit$BUGSoutput$sims.matrix[,"Mu.grp[10]"])


