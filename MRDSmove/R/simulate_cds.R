#' simulate waterfowl community MRDS double observer data subject to movement and measurement error
#' @param n_species Number of species to simulate data for
#' @param n_bins Number of bins for complete data (note that measurement error can result in some of these being detected)
#' @param n_obs_bins Number of detection bins
#' @param seed Optional; set the random # seed
#' @param p0 If TRUE (default), detection probability is less than 1.0 in the first bin.  If FALSE, detection is 1.0 in first bin
#' @param gaussian If TRUE, uses a Gaussian distribution for measurement error and a half-normal for movement; if FALSE (default), uses exponential / half exponential (Laplace dist)
#' @return a community mrds dataset
#' @export
#' @keywords simulation, mrds
#' @author Paul B. Conn
simulate_mrds <- function(n_species,n_bins,n_obs_bins,seed,p0=TRUE,gaussian=FALSE){
  set.seed(seed)
  expit<-function(x)1/(1+exp(-x))
  logit<-function(x)log(x/(1-x))
  N<-round(runif(n_species,9.5,100.5))
  
  ###  model for maximum detection probability (logit scale)
  beta_pmax_0 = logit(0.8)
  beta_pmax_obs2 = -0.5
  beta_pmax_sp = rnorm(n_species,0,0.25)
  beta_pmax_fly = 0.5
  beta_pmax_group = 0.02
  
  ### model for sigma (half normal detection fall off; log link)
  beta_sigma_0 = log(1)
  beta_sigma_obs2 = -0.2
  beta_sigma_sp = rnorm(n_species,0,0.25)
  beta_sigma_fly = 0.5
  beta_sigma_group=0.02
  
  ### group model - poisson-normal mixture
  Mu_group=rnorm(n_species,log(1),0.1)
  sigma_grp = 1.2
  
  ### measurement error model
  rdexp <- function(n,my_rate) ifelse(runif(n) > 0.5, 1, -1) * rexp(n,my_rate) #double exponential random variates
  #measure_rate = 1.5  #rate parameter for exponential distribution
  measure_rate=0.6
  
  ### movement model
  move_rate = 1.2
  
  ### proportion flying
  Prop_fly = runif(n_species,0.4,0.8)
  
  Complete_data = matrix(0,2000,9)
  colnames(Complete_data) = c("species","d1_true","d2_true","d1_obs","d2_obs","g_size","fly","det1","det2")
  counter=1
  
  if(gaussian==TRUE){
    Diff_cell = c(0,1,2,3,4)  #number of cells possible for movement
    Prob_move = dnorm(Diff_cell,0,1/move_rate)
    Prob_move=Prob_move/sum(Prob_move)
    Prob_measure = matrix(0,n_bins,n_bins)
    for(ibin in 1:n_bins){
      Low=c(0,pnorm(c(2:n_bins)-0.5,ibin,1/measure_rate))
      Up = pnorm(c(1:n_bins)+0.5,ibin,1/measure_rate)
      Prob_measure[ibin,]=Up-Low
    }           
  }
  for(isp in 1:n_species){
    Cur_rows=counter:(counter+N[isp]-1)
    Complete_data[Cur_rows,"species"]=isp
    Complete_data[Cur_rows,"fly"] = rbinom(N[isp],1,Prop_fly[isp])
    Complete_data[Cur_rows,"g_size"] = rpois(N[isp],exp(rnorm(N[isp],Mu_group[isp],sigma_grp)))+1  
    Complete_data[Cur_rows,"d1_true"] = round(runif(N[isp],0.5,n_bins+0.5))
    Complete_data[Cur_rows,"d1_obs"] =  Complete_data[Cur_rows,"d1_true"]
    if(gaussian){
      Complete_data[Cur_rows,"d2_true"]=Complete_data[Cur_rows,"d1_true"] + Complete_data[Cur_rows,"fly"]*rmultinom(length(Cur_rows),1,Prob_move)
      for(irow in 1:(length(Cur_rows)))Complete_data[Cur_rows[irow],"d2_obs"]=rmultinom(1,Prob_measure[Complete_data[Cur_rows[irow],"d1_true"],])
    }
    else{
      Complete_data[Cur_rows,"d2_true"] = floor(Complete_data[Cur_rows,"d1_true"] +  Complete_data[Cur_rows,"fly"]*rexp(N[isp],move_rate))
      Complete_data[Cur_rows,"d2_obs"] = Complete_data[Cur_rows,"d2_true"]+round(rdexp(N[isp],measure_rate))
    }
    #detection | true distances, species, fly, group
    if(p0==TRUE){
      Temp=beta_pmax_0 + beta_pmax_sp[Complete_data[Cur_rows,"species"]] + beta_pmax_fly*Complete_data[Cur_rows,"fly"] + beta_pmax_group*Complete_data[Cur_rows,"g_size"]
      P_max_obs1 = expit(Temp)
      P_max_obs2 = expit(Temp + beta_pmax_obs2)
    }
    else{
      P_max_obs1 = 1.0
      P_max.obs2 = 1.0
    }
    
    ### model for sigma (half normal detection fall off; log link)
    Temp=beta_sigma_0 + beta_sigma_sp[Complete_data[Cur_rows,"species"]] + beta_sigma_fly*Complete_data[Cur_rows,"fly"] + beta_sigma_group*Complete_data[Cur_rows,"g_size"]
    Sigma_obs1 = exp(Temp)
    Sigma_obs2 = exp(Temp+beta_sigma_obs2)
    P_det_obs1 = dnorm(Complete_data[Cur_rows,"d1_true"]-1,0,Sigma_obs1)/dnorm(0,0,Sigma_obs1)
    P_det_obs2 = dnorm(Complete_data[Cur_rows,"d2_true"]-1,0,Sigma_obs2)/dnorm(0,0,Sigma_obs2)
    
    ### determine detections
    Complete_data[Cur_rows,"det1"]=rbinom(N[isp],rep(1,N[isp]),P_max_obs1*P_det_obs1)
    Complete_data[Cur_rows,"det2"]=rbinom(N[isp],rep(1,N[isp]),P_max_obs1*P_det_obs2)
    
    
    counter=counter+N[isp]
    
  }
  Complete_data=Complete_data[1:(counter-1),]
  Obs_data=Complete_data[,-c(2,3)]
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
  return(list(Obs_data=Obs_data,Complete_data=Complete_data,G_true=G_true,N_true=N_true))
}