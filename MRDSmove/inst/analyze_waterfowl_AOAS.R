
#analyze waterfowl data using movement & measurement error models, HA models
library(MRDSmove)
library(mvtnorm)
library(RMark)
#source('./MRDSmove/R/MRDSmove_IntLik.R')
#source('./MRDSmove/R/MRDSmove_IntLik_ObsDep.R')
#source('./MRDSmove/R/sim_mrds_GOF.R')

#Data=read.csv('QMG 2014 Detections Distance Double-observer.csv')
data(waterfowl)
Data=waterfowl

#source('./MRDSmove/R/ht_mrds_ObsDep.R')
#source('./MRDSmove/R/ht_mrds.R')

SMALL = 0.0000001

Data$Species=as.character(Data$Species)
Sp.list = c("CAGO","KIEI","NOPI","ROPT","WFGO","LTDU","SACR") #,"TUSW","UK LOON")
Data=Data[Data$Species%in%Sp.list,]
n.species=length(Sp.list)

Data$DISTFRONT[Data$DISTFRONT==6]=NA
#Data$DISTBACK[Data$DISTBACK==6]=NA
Data$Dist = Data$DISTFRONT
Data$Dist[is.na(Data$Dist)]=Data$DISTBACK[is.na(Data$Dist)]
Data$Dist[Data$Dist==6]=NA  # So, only include observations for which distance of first observer was in bins 1-5 OR first observer missed it but second observer saw it in bins 1-5; i.e. observer 2 distance bin 6 is still included as a detection of first observer saw it in 1-5

Data$DISTBACK[Data$DISTBACK==6]=NA
Data$Distmean=(Data$DISTFRONT+Data$DISTBACK)/2
Data$Distmean[is.na(Data$Distmean)]=Data$DISTFRONT[is.na(Data$Distmean)]
Data$Distmean[is.na(Data$Distmean)]=Data$DISTBACK[is.na(Data$Distmean)]
DataHA1 = Data[-which(is.na(Data$Dist)==1),]
DataHA2 = Data[-which(is.na(Data$Distmean)==1),]
DataHA2$Dist=DataHA2$Distmean

Data = Data[-which(is.na(Data$DISTBACK) & is.na(Data$DISTFRONT)),]
n.hists=nrow(Data)
DataMML = data.frame(match=rep(c(1:n.hists),each=2),observer=rep(c(0,1),n.hists),
                     obs.dist=as.vector(rbind(Data$DISTFRONT,Data$DISTBACK)),g_size=rep(Data$GROUPSZ,each=2),species=rep(Data$Species,each=2),
                     moving = rep(Data$groupFLY,each=2),detected=as.vector(rbind(1-is.na(Data$DISTFRONT),1-is.na(Data$DISTBACK))))
DataMML$det_other=as.vector(rbind(DataMML$detected[c(1:n.hists)*2],DataMML$detected[c(1:n.hists)*2-1]))


G.obs = Data$GROUPSZ
p.formula=~distance+distance2+moving+g_size+observer+observer:distance+observer:distance2
pi.formula=~0+det_other:distance
#li.formula=~0+det_other+det_other:distance
fi.formula=~0
Move.fix=c(0,0,0)

#n.bootstraps=1000
Obs.bins=c(1:5)

Results = data.frame(Model=rep("a",14),k=rep(0,14),LogL=rep(0,14),AIC=rep(0,14),Nhat=rep(0,14),SE.N=rep(0,14),Ghat=rep(0,14),SE.G=rep(0,14))
Results$Model = rep("a",14) #for some reason this defaults to factor in the above

imod=1
# gaussian, gaussian, fully indep
my.par = c(1,.07,-.09,.5,0.2,-0.5,.07,-.07,log(.5),log(1.5),log(.5))  
MML_gg_fi = glm_out = optim(par=my.par,MRDSmove_IntLik,hessian=TRUE,method="BFGS",Data=DataMML,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
HT.N = ht_mrds(Par=glm_out$par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_gg_fi"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1

# gaussian, gaussian, obs point indep
my.par = c(glm_out$par[1:8],0.2,log(.5),log(1.5),log(.5))  
MML_gg_pi = glm_out = optim(par=my.par,MRDSmove_IntLik_ObsDep,hessian=TRUE,method="BFGS",Data=DataMML,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
HT.N = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_gg_pi"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1


# laplace,laplace, fully indep
my.par = c(glm_out$par[1:8],log(5),log(1),log(5))  
MML_ll_fi = glm_out = optim(par=my.par,MRDSmove_IntLik,hessian=TRUE,method="BFGS",Data=DataMML,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=FALSE,gaussian.meas=FALSE)
Sigma = solve(glm_out$hessian)
HT.G = ht_mrds(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
HT.N = ht_mrds(Par=glm_out$par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_ll_fi"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1

# laplace, laplace, obs point indep
my.par = c(glm_out$par[1:8],0.2,glm_out$par[9:11])  
MML_ll_pi = glm_out = optim(par=my.par,MRDSmove_IntLik_ObsDep,hessian=TRUE,method="BFGS",Data=DataMML,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=FALSE,gaussian.meas=FALSE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
HT.N = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_ll_pi"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1




# gaussian, gaussian, fully indep, distance*flying
p.formula=~distance+distance2+moving+g_size+observer+observer:distance+observer:distance2+moving:distance+moving:distance2

my.par = c(1,.07,-.09,.5,0.2,-0.5,0,0,0,0,log(.5),log(1.5),log(.5))  
MML_gg_fi_fly = glm_out = optim(par=my.par,MRDSmove_IntLik,hessian=TRUE,method="BFGS",Data=DataMML,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
HT.N = ht_mrds(Par=glm_out$par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_gg_fi_fly"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1


# gaussian, gaussian, obs point indep
my.par = c(glm_out$par[1:10],0.2,log(.5),log(1.5),log(.5))  
MML_gg_pi_fly = glm_out = optim(par=my.par,MRDSmove_IntLik_ObsDep,hessian=TRUE,method="BFGS",Data=DataMML,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
HT.N = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_gg_pi_fly"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1



# laplace,laplace, fully indep
my.par = c(MML_gg_fi_fly$par[1:10],log(5),log(1),log(5))  
MML_ll_fi_fly = glm_out = optim(par=my.par,MRDSmove_IntLik,hessian=TRUE,method="BFGS",Data=DataMML,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=FALSE,gaussian.meas=FALSE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
HT.N = ht_mrds(Par=glm_out$par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds(Par=Temp.par,Data=DataMML,G=G.obs,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_ll_fi_fly"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1

# laplace, laplace, obs point indep
my.par = c(glm_out$par[1:10],0.2,glm_out$par[11:13])  
MML_ll_pi_fly = glm_out = optim(par=my.par,MRDSmove_IntLik_ObsDep,hessian=TRUE,method="BFGS",Data=DataMML,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=FALSE,gaussian.meas=FALSE)
Sigma = solve(glm_out$hessian)
#ht estimate and SE
HT.G = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
HT.N = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)
#additional variance 
Grad.N = Grad.G = rep(0,length(my.par))
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.G[ipar]=(New-g_est)/SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var=NULL,gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.N[ipar]=(New-n_est)/SMALL
}
E.Var.G = crossprod(Grad.G,Sigma%*%Grad.G)
E.Var.N = crossprod(Grad.N,Sigma%*%Grad.N)
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N
Results$Model[imod]="MML_ll_pi_fly"
Results$k[imod]=length(my.par)
Results$SE.G[imod] = sqrt(Var.G)
Results$SE.N[imod] = sqrt(Var.N)
Results$LogL[imod] = -glm_out$value
Results$AIC[imod] = 2*Results$k[imod]-2*Results$LogL[imod]
Results$Nhat[imod] = n_est
Results$Ghat[imod] = g_est
imod=imod+1





n.gof=1000
Data.sim=vector("list",n.gof)
for(i in 1:n.gof){
  Data.sim[[i]]=sim_mrds_GOF_ObsDep(Par=MML_ll_pi_fly$par,Data=DataMML,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=c(1:5),Move.fix=Move.fix,gaussian.move=FALSE,gaussian.meas=FALSE)
}

##make some plots
#1) plot detected distances as a function of observer and moving/not moving
Dist.total=array(0,c(2,2,5,n.gof))
Obs.total=array(0,c(2,2,5))
Hist.freq = array(0,c(2,3,n.gof))
for(imove in 1:2){
  for(iobs in 1:2){
    for(ibin in 1:5){
      Obs.total[imove,iobs,ibin]=length(which(DataMML$observer==(iobs-1) & DataMML$moving==(imove-1) & DataMML$obs.dist==ibin))
      for(igof in 1:n.gof){
        Dist.total[imove,iobs,ibin,igof]=length(which(Data.sim[[igof]]$observer==(iobs-1) & Data.sim[[igof]]$moving==(imove-1) & Data.sim[[igof]]$obs.dist==ibin))
      }
    }
  }
}

for(imove in 1:2){
  for(igof in 1:n.gof){
    Hist.freq[imove,1,igof]=length(which(Data.sim[[igof]][Data.sim[[igof]]$observer==0 & Data.sim[[igof]]$moving==(imove-1),"detected"]==1 & Data.sim[[igof]][Data.sim[[igof]]$observer==1 & Data.sim[[igof]]$moving==(imove-1),"detected"]==1))
    Hist.freq[imove,2,igof]=length(which(Data.sim[[igof]][Data.sim[[igof]]$observer==0 & Data.sim[[igof]]$moving==(imove-1),"detected"]==1 & Data.sim[[igof]][Data.sim[[igof]]$observer==1 & Data.sim[[igof]]$moving==(imove-1),"detected"]==0))
    Hist.freq[imove,3,igof]=length(which(Data.sim[[igof]][Data.sim[[igof]]$observer==0 & Data.sim[[igof]]$moving==(imove-1),"detected"]==0 & Data.sim[[igof]][Data.sim[[igof]]$observer==1 & Data.sim[[igof]]$moving==(imove-1),"detected"]==1))
  }
}


Quantile = array(0,dim=c(2,2,5,2))
Mean = array(0,dim=c(2,2,5))
for(imove in 1:2){
  for(iobs in 1:2){
    for(ibin in 1:5){
      Quantile[imove,iobs,ibin,]=quantile(Dist.total[imove,iobs,ibin,],c(0.025,0.975))
      Mean[imove,iobs,ibin] = mean(Dist.total[imove,iobs,ibin,])
    }
  }
}
#format for ggplot
Plot_df = data.frame(Count=as.vector(Quantile),Move=rep(c("Not moving","Moving"),2*5*2),
                     Observer=rep(c("Observer 1","Observer 1","Observer 2","Observer 2"),5*2),Distance=rep(rep(c(1:5),each=4),2),
                     Quantity=rep(c("SimQuantileLower","SimQuantileUpper"),each=20))
Plot_df_mean = data.frame(Count=as.vector(Mean),Move=rep(c("Not moving","Moving"),2*5),
                     Observer=rep(c("Observer 1","Observer 1","Observer 2","Observer 2"),5),Distance=rep(rep(c(1:5),each=4)),
                     Quantity="Mean")


library(ggplot2)
GOF_plot1 = ggplot(Plot_df)+geom_line(aes(x=Distance,y=Count,group=Quantity),linetype=2)+facet_grid(Move~Observer)

GOF_plot1 = GOF_plot1 + geom_line(data=Plot_df_mean,linetype=2,aes(x=Distance,y=Count),size=1.2)+theme_grey(base_size=16)


Plot_df2=data.frame(Count=as.vector(Obs.total),Move=rep(c("Not moving","Moving"),2*5),
                    Observer=rep(c("Observer 1","Observer 1","Observer 2","Observer 2"),5),Distance=rep(c(1:5),each=4),
                    Quantity="Observed")
GOF_plot1 = GOF_plot1 + geom_line(data=Plot_df2,aes(x=Distance,y=Count),size=1.2)+theme_grey(base_size=16)
                     
pdf("bootstrapGOF.pdf")
  GOF_plot1
dev.off()

## calculate frequency of 11, 10, 01 hists
# move = 0, 11
quantile(Hist.freq[1,1,],c(0.025,0.975)) #11
quantile(Hist.freq[1,2,],c(0.025,0.975)) #10
quantile(Hist.freq[1,3,],c(0.025,0.975)) #01
# move = 1
quantile(Hist.freq[2,1,],c(0.025,0.975)) #11
quantile(Hist.freq[2,2,],c(0.025,0.975)) #10
quantile(Hist.freq[2,3,],c(0.025,0.975)) #01

#actual numbers
length(which(DataMML$detected[DataMML$observer==0 & DataMML$moving==0]==1 & DataMML$detected[DataMML$observer==1 & DataMML$moving==0]==1))
length(which(DataMML$detected[DataMML$observer==0 & DataMML$moving==0]==1 & DataMML$detected[DataMML$observer==1 & DataMML$moving==0]==0))
length(which(DataMML$detected[DataMML$observer==0 & DataMML$moving==0]==0 & DataMML$detected[DataMML$observer==1 & DataMML$moving==0]==1))
length(which(DataMML$detected[DataMML$observer==0 & DataMML$moving==1]==1 & DataMML$detected[DataMML$observer==1 & DataMML$moving==1]==1))
length(which(DataMML$detected[DataMML$observer==0 & DataMML$moving==1]==1 & DataMML$detected[DataMML$observer==1 & DataMML$moving==1]==0))
length(which(DataMML$detected[DataMML$observer==0 & DataMML$moving==1]==0 & DataMML$detected[DataMML$observer==1 & DataMML$moving==1]==1))


##plot kernel estimates

MML_ll_pi_fly
Movement = Measure = Distance = c(-4:4)
Movement[1:4]=dexp(abs(Distance[1:4]),exp(2.67))
Movement[5:9]=dexp(abs(Distance[5:9]),exp(.42))
Measure=dexp(abs(Distance),exp(1.24))
Movement=Movement/sum(Movement)
Measure=Measure/sum(Measure)
Plot_df = data.frame(Distance=rep(c(-4:4),2),Kernel=rep(c("Movement","Measure"),each=9),Probability=c(Movement[1:9],Measure[1:9]))

library(ggplot2)
kernel_plot = ggplot(Plot_df)+geom_line(aes(x=Distance,y=Probability,lty=Kernel),size=1.2)+theme_gray(base_size=16)+xlab(expression(paste("Distance bin discrepancy (",delta,")")))

pdf("Kernel_plot.pdf")
 kernel_plot
dev.off()


####### estimate species-specific density for "best" model
#ht estimate and SE
glm_out=MML_ll_pi_fly
my.par=glm_out$par
HT.G = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var="species",gaussian.move=FALSE,gaussian.meas=FALSE)
HT.N = ht_mrds_ObsDep(Par=glm_out$par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var="species",gaussian.move=FALSE,gaussian.meas=FALSE)
#additional variance 
Grad.N = Grad.G = matrix(0,length(my.par),n.species)
g_est=HT.G$n.hat
n_est=HT.N$n.hat
for(ipar in 1:length(my.par)){
  Temp.par=glm_out$par
  Temp.par[ipar]=Temp.par[ipar]+SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=rep(1,n.hists),p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var="species",gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.G[ipar,]=(New-g_est)/SMALL
  New = ht_mrds_ObsDep(Par=Temp.par,Data=DataMML,G=G.obs,p.formula=p.formula,dep.formula=pi.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,Agg.var="species",gaussian.move=FALSE,gaussian.meas=FALSE)$n.hat
  Grad.N[ipar,]=(New-n_est)/SMALL
}
E.Var.G = E.Var.N = rep(0,n.species)
for(isp in 1:n.species){
  E.Var.G[isp] = crossprod(Grad.G[,isp],Sigma%*%Grad.G[,isp])
  E.Var.N[isp] = crossprod(Grad.N[,isp],Sigma%*%Grad.N[,isp])
}
Var.G = HT.G$var.E + E.Var.G
Var.N = HT.N$var.E + E.Var.N

SE.N = sqrt(Var.N)
C = exp(1.96*sqrt(log(1+Var.N/n_est^2)))
CI.low = n_est/C  #log-based CI for N (see e.g. Buckland et al. Distance sampling book, pg 88)
CI.hi = n_est*C

#corresponding density estimates
D = n_est/189.6
D.low = CI.low/189.6
D.hi = CI.hi/189.6


CI_DF = data.frame(Density=D,D.low=D.low,D.hi=D.hi,Species=names(D),x=c(1:length(D)))
CI_plot = ggplot(CI_DF,aes(x=x))+geom_errorbar(size=1,aes(ymin=D.low,ymax=D.hi),width=0.1)+geom_point(size=2,aes(y=Density))
CI_plot = CI_plot + scale_x_discrete(limits=c(1:7),labels=as.character(CI_DF$Species)) + xlab('Species') + ylab(expression(paste("Density (individuals/",km^2,")")))+theme(text=element_text(size=14))
CI_plot

pdf("SpeciesN.pdf")
  CI_plot
dev.off()

