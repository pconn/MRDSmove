#' simulate waterfowl community MRDS double observer data subject to movement and measurement error.  This version uses a half-normal detection function.
#' @param n_species Number of species to simulate data for
#' @param n_bins Number of bins for complete data (note that measurement error can result in some of these being detected)
#' @param n_obs_bins Number of bins observed
#' @param N true abundance within the complete data area
#' @param measure_par Gives measurement error sd (if gaussian = TRUE) or measurement error rate (if gaussian=FALSE)
#' @param move_par 2 parameter vector movement sd (if gaussian = TRUE) or rate (if gaussian=FALSE) for responseive movement 
#' @param seed Optional; set the random # seed
#' @param gaussian If TRUE, uses a Gaussian distribution for measurement error and a half-normal for movement; if FALSE (default), uses exponential / half exponential (Laplace dist)
#' @param point_indep if TRUE (default is FALSE), implements point independence via a random effect
#' @return a community mrds dataset
#' @export
#' @keywords simulation, mrds
#' @author Paul B. Conn
simulate_mrds_hn <- function(n_species,n_bins,n_obs_bins,N,measure_par=0.6,move_par=c(0.1,1.2),seed=NULL,gaussian=FALSE,point_indep=FALSE){
  if(is.null(seed)==FALSE)set.seed(seed)
  expit<-function(x)1/(1+exp(-x))
  logit<-function(x)log(x/(1-x))
 
  ###  model for detection probability
  #1) sigma
  beta_sigma_0 = 1
  beta_sigma_obs2 = 0 #-0.5
  beta_sigma_sp = 0
  if(n_species>1)beta_sigma_sp=c(0,rnorm(n_species-1,0,0.1))
  beta_sigma_fly = 0.2
  RE=0
  if(point_indep)RE=runif(N,-0.7,0.7)
  beta_sigma_group = 0 #0.1

  #2) g0
  beta_g0_0 = 1
  beta_g0_obs2 = 0 #-0.5
  beta_g0_sp = 0
  if(n_species>1)beta_g0_sp=c(0,rnorm(n_species-1,0,0.25))
  beta_g0_fly = 0.5
  beta_g0_group = 0 #0.1
  
  Bin.midpoints = c(0:(n_bins-1))+0.5
  
  ### group model - poisson-normal mixture
  Mu_group=rnorm(n_species,log(1),0.1)
  sigma_grp=0.000001
  #sigma_grp = 1.2
  
  ### measurement error model
  rdexp <- function(n,my_par) ifelse(runif(n) > 0.5, 1, -1) * rexp(n,my_par) #double exponential random variates
  #measure_par = 1.5  #par parameter for exponential distribution

  ### proportion flying
  Prop_fly = 0.75
  
  Complete_data = matrix(0,sum(N),10)
  colnames(Complete_data) = c("species","d1_true","d2_true","d1_obs","d2_obs","g_size","fly","det1","det2","pi_re")
  counter=1
  
  Psi = Meas = matrix(0,n_bins,n_bins)
  if(gaussian==TRUE){
    for(ibin1 in 1:n_bins){
      Psi[ibin1,ibin1:n_bins]=dnorm(Bin.midpoints[ibin1:n_bins],Bin.midpoints[ibin1],move_par[2])
      Meas[ibin1,] = dnorm(Bin.midpoints,Bin.midpoints[ibin1],measure_par)
    }
    for(ibin1 in 2:n_bins)Psi[ibin1,1:(ibin1-1)]=dnorm(Bin.midpoints[1:(ibin1-1)],Bin.midpoints[ibin1],move_par[1])
  }
  if(gaussian==FALSE){  #double exponential
    for(ibin1 in 1:n_bins){
      Psi[ibin1,ibin1:n_bins]=dexp(Bin.midpoints[ibin1:n_bins]-Bin.midpoints[ibin1],move_par[2])
      Meas[ibin1,] = dexp(abs(Bin.midpoints-Bin.midpoints[ibin1]),measure_par)
    }
    for(ibin1 in 2:n_bins)Psi[ibin1,1:(ibin1-1)]=dexp(Bin.midpoints[ibin1]-Bin.midpoints[1:(ibin1-1)],move_par[1])
  }
  Psi = Psi/rowSums(Psi)
  Meas = Meas/rowSums(Meas)
  
  
  for(isp in 1:n_species){
    Cur_rows=counter:(counter+N[isp]-1)
    Complete_data[Cur_rows,"species"]=isp
    Complete_data[Cur_rows,"fly"] = rbinom(N[isp],1,Prop_fly[isp])
    Complete_data[Cur_rows,"g_size"] = rpois(N[isp],exp(rnorm(N[isp],Mu_group[isp],sigma_grp)))+1  
    Complete_data[Cur_rows,"d1_true"] = round(runif(N[isp],0.5,n_bins+0.5))
    Complete_data[Cur_rows,"d2_true"] = Complete_data[Cur_rows,"d1_true"]
    
    for(irow in 1:(length(Cur_rows))){
      if(Complete_data[Cur_rows[irow],"fly"]==1)Complete_data[Cur_rows[irow],"d2_true"]=sample(c(1:n_bins),1,prob=Psi[Complete_data[Cur_rows[irow],"d1_true"],],replace=TRUE)
      Complete_data[Cur_rows[irow],"d1_obs"]=sample(c(1:n_bins),1,prob=Meas[Complete_data[Cur_rows[irow],"d1_true"],])
      Complete_data[Cur_rows[irow],"d2_obs"]=sample(c(1:n_bins),1,prob=Meas[Complete_data[Cur_rows[irow],"d2_true"],])
    }

    #detection | true distances, species, fly, group

    ### polynomial detection model for g0, sigma
    G0_det_obs1 = expit(beta_g0_0 + beta_g0_sp[Complete_data[Cur_rows,"species"]] + beta_g0_fly*Complete_data[Cur_rows,"fly"] + beta_g0_group*Complete_data[Cur_rows,"g_size"])
    G0_det_obs2 = expit(beta_g0_0 + beta_g0_obs2 + beta_g0_sp[Complete_data[Cur_rows,"species"]] + beta_g0_fly*Complete_data[Cur_rows,"fly"] + beta_g0_group*Complete_data[Cur_rows,"g_size"]) 
    Sigma_det_obs1 = exp(beta_sigma_0 + beta_sigma_sp[Complete_data[Cur_rows,"species"]] + beta_sigma_fly*Complete_data[Cur_rows,"fly"] + beta_sigma_group*Complete_data[Cur_rows,"g_size"]+RE)
    Sigma_det_obs2 = exp(beta_sigma_0 + beta_sigma_obs2 + beta_sigma_sp[Complete_data[Cur_rows,"species"]] + beta_sigma_fly*Complete_data[Cur_rows,"fly"] + beta_sigma_group*Complete_data[Cur_rows,"g_size"]+RE) 
    P_det_obs1 = G0_det_obs1*dnorm(Complete_data[Cur_rows,"d1_true"],1,Sigma_det_obs1)/dnorm(1,1,Sigma_det_obs1)
    P_det_obs2 = G0_det_obs2*dnorm(Complete_data[Cur_rows,"d2_true"],1,Sigma_det_obs2)/dnorm(1,1,Sigma_det_obs2)

    
    
    ### determine detections
    Complete_data[Cur_rows,"det1"]=rbinom(N[isp],rep(1,N[isp]),P_det_obs1)
    Complete_data[Cur_rows,"det2"]=rbinom(N[isp],rep(1,N[isp]),P_det_obs2)
    
    counter=counter+N[isp]
    
  }
  Complete_data=Complete_data[1:(counter-1),]
  Obs_data=Complete_data[,-c(2,3,10)]
  Obs_data[which(Obs_data[,"d1_obs"]<1 | Obs_data[,"d1_obs"]>n.obs.bins),"det1"]=0
  Obs_data[which(Obs_data[,"d2_obs"]<1 | Obs_data[,"d2_obs"]>n.obs.bins),"det2"]=0
  Obs_data[which(Obs_data[,"det1"]==0),"d1_obs"]=NA
  Obs_data[which(Obs_data[,"det2"]==0),"d2_obs"]=NA
  Obs_data=Obs_data[-which(Obs_data[,"det1"]==0 & Obs_data[,"det2"]==0),]
  
  N_true=G_true=0*N
  for(isp in 1:n_species){
    Which_sp_obs=which(Complete_data[,"species"]==isp & Complete_data[,"d1_true"] %in% c(1:n.obs.bins))
    G_true[isp]=length(Which_sp_obs)
    N_true[isp]=sum(Complete_data[Which_sp_obs,"g_size"])
  }
  return(list(Obs_data=Obs_data,Complete_data=Complete_data,G_true=G_true,N_true=N_true,Mu_group=Mu_group,Beta=list(beta_sigma_0,beta_sigma_sp,beta_sigma_obs2,beta_sigma_fly,beta_sigma_group)))
}