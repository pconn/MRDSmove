#  run simulation study examining performance of integrated likelihood approach for modeling
# movement and measurement error in MRDS surveys with binned distances
library('MRDSmove')
library('RMark')
library('mvtnorm')
#source('c:/users/paul.conn/git/Alisauskas/MRDSmove/R/ht_mrds.R')
#source('c:/users/paul.conn/git/Alisauskas/MRDSmove/R/ht_mrds_ObsDep.R')
#source('c:/users/paul.conn/git/Alisauskas/MRDSmove/R/MRDSmove_IntLik.R')
#source('c:/users/paul.conn/git/Alisauskas/MRDSmove/R/MRDSmove_IntLik_ObsDep.R')

n.obs.bins=5


expit<-function(x)1/(1+exp(-x))
n.sims = 500
set.seed(123456)

start.time = Sys.time()

M.pars = matrix(0.0001,2,3)
#M.pars[1,] = c(.7,.7,.3)
M.pars[2,] = c(0.0001,1.5,.5)
Obs.bins=c(1:n.obs.bins)

N.vals = c(400,2000)
N.est= SE = Cov = array(0,dim=c(n.sims,2,2,5))  # 2 N options * 2 M options * 4 estimation models
N.true = array(0,dim=c(n.sims,2,2))
p.formula=~distance+distance2+moving
dep.formula.li=~0+det_other+det_other:distance
dep.formula=~0+det_other:distance

SMALL = 0.0000001 #for gradient calculation

for(isim in 1:n.sims){
  cat(paste("\n simulation ",isim,"\n"))
  for(iN in 1:2){   #1:2
    for(iM in 1:2){  #1:2
      for(idep in 1:1){
        Sim_data <- simulate_mrds_hn(n_species=1,n_bins=10,n_obs_bins=5,N=N.vals[iN],measure_par=M.pars[iM,3],move_par=M.pars[iM,1:2],gaussian=TRUE,point_indep=idep)
        Obs_data=data.frame(Sim_data$Obs_data)
        Char_hist = factor(apply(cbind(Obs_data$d1_obs,Obs_data$d2_obs,Obs_data$det1,Obs_data$det2,Obs_data$fly),1,'paste',collapse=''))
        Counts = tabulate(Char_hist)
        Duplicated = duplicated(Char_hist)
        Which.unique = which(Duplicated==0)
        Obs_data2 = Obs_data[Which.unique,]
        n.hists=nrow(Obs_data2)
        Order = rep(0,n.hists)
        for(ihist in 1:n.hists)Order[ihist]=which(levels(Char_hist)==Char_hist[Which.unique[ihist]])
        Data = data.frame(match=rep(c(1:n.hists),each=2),observer=rep(c(0,1),n.hists),species=rep(Obs_data2$species,each=2),
                          obs.dist=as.vector(rbind(Obs_data2$d1_obs,Obs_data2$d2_obs)),g_size=rep(Obs_data2$g_size,each=2),
                          moving = rep(Obs_data2$fly,each=2),detected=as.vector(rbind(Obs_data2$det1,Obs_data2$det2)),det_other=as.vector(rbind(Obs_data2$det2,Obs_data2$det1)))
        Data$count = rep(Counts[Order],each=2)
        
        N.true[isim,iN,iM] = sum(Sim_data$Complete_data[,"d1_true"]%in%Obs.bins)
        Which.fix = which(M.pars[iM,] < 0.01)
        Move.fix = rep(0,3)
        my.par = c(1,.07,-.09,.5,0)
        if(length(Which.fix)==0)my.par = c(my.par,log(M.pars[iM,]))
        if(length(Which.fix)<3 & length(Which.fix)>0)my.par = c(my.par,log(M.pars[iM,-Which.fix]))
        if(length(Which.fix)>0)Move.fix[Which.fix]=1
        Grad=rep(0,length(my.par))
        
        #### MODELS WITH PI PARAM 
        # (1) 8 bin integrated likelihood
        glm_out = optim(par=my.par,MRDSmove_IntLik_ObsDep,hessian=TRUE,method="BFGS",Data=Data,p.formula=p.formula,dep.formula=dep.formula,Bin.widths=rep(1,8),Obs.bins=c(1:n.obs.bins),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
        Sigma = try(solve(glm_out$hessian))
        if(class(Sigma)!="matrix")N.est[isim,iN,iM,1] = SE[isim,iN,iM,1] = Cov[isim,iN,iM,1] = NA
        else{
          if(sum(diag(Sigma)<0)>0 | max(diag(Sigma))>10)N.est[isim,iN,iM,1] = SE[isim,iN,iM,1] = Cov[isim,iN,iM,1] = NA
          else{
            #ht estimate and SE
            HT = ht_mrds_ObsDep(Par=glm_out$par,Data=Data,G=rep(1,n.hists),p.formula=p.formula,dep.formula=dep.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
            g_est = HT$n.hat
            for(ipar in 1:length(my.par)){
              Temp.par=glm_out$par
              Temp.par[ipar]=Temp.par[ipar]+SMALL
              new_g = ht_mrds_ObsDep(Par=Temp.par,Data=Data,G=rep(1,n.hists),p.formula=p.formula,dep.formula=dep.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
              Grad[ipar]=(new_g-g_est)/SMALL
            }
            E.Var = crossprod(Grad,solve(glm_out$hessian,Grad))
            N.est[isim,iN,iM,1]=g_est
            Var.E = HT$var.E
            Var = E.Var+Var.E
            SE[isim,iN,iM,1] = sqrt(Var)
            C = exp(1.96*sqrt(log(1+Var/N.est[isim,iN,iM,1]^2)))
            CI.low = N.est[isim,iN,iM,1]/C  #log-based CI for N (see e.g. Buckland et al. Distance sampling book, pg 88)
            CI.hi = N.est[isim,iN,iM,1]*C
            Cov[isim,iN,iM,1] = ((N.true[isim,iN,iM] > CI.low) & (N.true[isim,iN,iM]<= CI.hi))
          }
        }
        
        ##### LI Param
        
        my.par = c(1,.07,-.09,.5,0,0)
        if(length(Which.fix)==0)my.par = c(my.par,log(M.pars[iM,]))
        if(length(Which.fix)<3 & length(Which.fix)>0)my.par = c(my.par,log(M.pars[iM,-Which.fix]))
        if(length(Which.fix)>0)Move.fix[Which.fix]=1
        Grad=rep(0,length(my.par))
        
        # (1) 8 bin integrated likelihood
        glm_out = optim(par=my.par,MRDSmove_IntLik_ObsDep,hessian=TRUE,method="BFGS",Data=Data,p.formula=p.formula,dep.formula=dep.formula.li,Bin.widths=rep(1,8),Obs.bins=c(1:n.obs.bins),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
        #Sigma = solve(hessian(MRDSmove_IntLik_ObsDep,glm_out$par,Data=Data,p.formula=p.formula,dep.formula=dep.formula,Bin.widths=rep(1,8),Obs.bins=c(1:n.obs.bins),Move.fix=Move.fix,gaussian=TRUE))
        Sigma = try(solve(glm_out$hessian))
        if(class(Sigma)!="matrix")N.est[isim,iN,iM,5] = SE[isim,iN,iM,5] = Cov[isim,iN,iM,5] = NA
        else{
          if(sum(diag(Sigma)<0)>0 | max(diag(Sigma))>10)N.est[isim,iN,iM,5] = SE[isim,iN,iM,5] = Cov[isim,iN,iM,5] = NA
          else{
            #ht estimate and SE
            HT = ht_mrds_ObsDep(Par=glm_out$par,Data=Data,G=rep(1,n.hists),p.formula=p.formula,dep.formula=dep.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
            g_est = HT$n.hat
            for(ipar in 1:length(my.par)){
              Temp.par=glm_out$par
              Temp.par[ipar]=Temp.par[ipar]+SMALL
              new_g = ht_mrds_ObsDep(Par=Temp.par,Data=Data,G=rep(1,n.hists),p.formula=p.formula,dep.formula=dep.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
              Grad[ipar]=(new_g-g_est)/SMALL
            }
            E.Var = crossprod(Grad,solve(glm_out$hessian,Grad))
            N.est[isim,iN,iM,5]=g_est
            ####
            Var.E = HT$var.E
            Var = E.Var+Var.E
            SE[isim,iN,iM,5] = sqrt(Var)
            C = exp(1.96*sqrt(log(1+Var/N.est[isim,iN,iM,5]^2)))
            CI.low = N.est[isim,iN,iM,5]/C  #log-based CI for N (see e.g. Buckland et al. Distance sampling book, pg 88)
            CI.hi = N.est[isim,iN,iM,5]*C
            Cov[isim,iN,iM,5] = ((N.true[isim,iN,iM] > CI.low) & (N.true[isim,iN,iM]<= CI.hi))
          }
        }       
        
        my.formula=~distance+distance2+moving
        #### FI models
        # (2) 8 bin integrated likelihood
        my.par = c(1,.07,-.09,.5)
        if(length(Which.fix)==0)my.par = c(my.par,log(M.pars[iM,]))
        if(length(Which.fix)<3 & length(Which.fix)>0)my.par = c(my.par,log(M.pars[iM,-Which.fix]))
        Grad=rep(0,length(my.par))
        glm_out = optim(par=my.par,MRDSmove_IntLik,hessian=TRUE,method="BFGS",Data=Data,mod.formula=p.formula,Bin.widths=rep(1,8),Obs.bins=c(1:n.obs.bins),Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
        Sigma = try(solve(glm_out$hessian))
        if(class(Sigma)!="matrix")N.est[isim,iN,iM,2] = SE[isim,iN,iM,2] = Cov[isim,iN,iM,2] = NA
        else{
          if(sum(diag(Sigma)<0)>0 | max(diag(Sigma))>10)N.est[isim,iN,iM,2] = SE[isim,iN,iM,2] = Cov[isim,iN,iM,2] = NA
          else{
            #ht estimate and SE
            HT = ht_mrds(Par=glm_out$par,Data=Data,G=rep(1,n.hists),mod.formula=my.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)
            g_est = HT$n.hat           
            for(ipar in 1:length(my.par)){
              Temp.par=glm_out$par
              Temp.par[ipar]=Temp.par[ipar]+SMALL
              new_g = ht_mrds(Par=Temp.par,Data=Data,G=rep(1,n.hists),mod.formula=my.formula,Bin.widths=rep(1,8),Obs.bins=Obs.bins,Move.fix=Move.fix,gaussian.move=TRUE,gaussian.meas=TRUE)$n.hat
              Grad[ipar]=(new_g-g_est)/SMALL
            }
            E.Var = crossprod(Grad,solve(glm_out$hessian,Grad))
            N.est[isim,iN,iM,2]=g_est
            Var.E = HT$var.E
            Var = E.Var+Var.E
            SE[isim,iN,iM,2] = sqrt(Var)
            C = exp(1.96*sqrt(log(1+Var/N.est[isim,iN,iM,2]^2)))
            CI.low = N.est[isim,iN,iM,2]/C  #log-based CI for N (see e.g. Buckland et al. Distance sampling book, pg 88)
            CI.hi = N.est[isim,iN,iM,2]*C
            Cov[isim,iN,iM,2] = ((N.true[isim,iN,iM] > CI.low) & (N.true[isim,iN,iM]<= CI.hi))
          }
        }
        # (3) FI model with Huggins-Alho in RMark - priority to first observer distance
        Dat_MARK = data.frame(Obs_data)
        Dat_MARK$ch = paste0(as.character(Dat_MARK$det1),as.character(Dat_MARK$det2))
        Dat_MARK$dist = Dat_MARK$d1_obs
        Dat_MARK$dist[is.na(Dat_MARK$d1_obs)]=Dat_MARK$d2_obs[is.na(Dat_MARK$d1_obs)]
        Dat_MARK$dist2 = Dat_MARK$dist^2
        HA_proc = process.data(Dat_MARK,model="Huggins")
        HA_ddl = make.design.data(HA_proc)
        p_model = list(formula=~dist+dist2+fly,share=TRUE)
        HA_mod= try(mark(HA_proc,HA_ddl,model.parameters=list(p=p_model)),TRUE)
        if(class(HA_mod)[1]=="mark"){
          N.est[isim,iN,iM,3]=unlist(HA_mod$results$derived)[1]
          SE[isim,iN,iM,3] = unlist(HA_mod$results$derived)[2]
          Cov[isim,iN,iM,3] = ((N.true[isim,iN,iM] > unlist(HA_mod$results$derived)[3]) & (N.true[isim,iN,iM]<= unlist(HA_mod$results$derived)[4]))
        }
               
        # (4) FI model with Huggins-Alho in RMark - take average distance if both observers see an animal cluster
        Dat_MARK = data.frame(Obs_data)
        Dat_MARK$ch = paste0(as.character(Dat_MARK$det1),as.character(Dat_MARK$det2))
        Dat_MARK$dist = (Dat_MARK$d2_obs+Dat_MARK$d1_obs)/2
        Dat_MARK$dist[is.na(Dat_MARK$dist)]=Dat_MARK$d1_obs[is.na(Dat_MARK$dist)]
        Dat_MARK$dist[is.na(Dat_MARK$dist)]=Dat_MARK$d2_obs[is.na(Dat_MARK$dist)]
        Dat_MARK$dist2 = Dat_MARK$dist^2
        HA_proc = process.data(Dat_MARK,model="Huggins")
        HA_ddl = make.design.data(HA_proc)
        p_model = list(formula=~dist+dist2+fly,share=TRUE)
        HA_mod= try(mark(HA_proc,HA_ddl,model.parameters=list(p=p_model)),TRUE)
        if(class(HA_mod)[1]=="mark"){
          N.est[isim,iN,iM,4]=unlist(HA_mod$results$derived)[1]
          SE[isim,iN,iM,4] = unlist(HA_mod$results$derived)[2]
          Cov[isim,iN,iM,4] = ((N.true[isim,iN,iM] > unlist(HA_mod$results$derived)[3]) & (N.true[isim,iN,iM]<= unlist(HA_mod$results$derived)[4]))
        }
               
        #clean up directory 
        shell('del *.vcv')
        shell('del *.res')
        shell('del *.inp')
        shell('del *.out')
        shell('del *.tmp')
      }
    } 
    if(isim%%50==0){
      Out = list(N.true=N.true,N.est=N.est,SE=SE,Cov=Cov)
      save(Out,file="SimResultsPILI_revision.Rdata")
    }
  }
}

Sys.time()-start.time
