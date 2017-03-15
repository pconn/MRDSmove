
### analyze Canada waterfowl data for the Alisauskas and Conn Euring paper - note this does not use the movement
# and measurement error model

data(waterfowl)
Data=waterfowl

Delta_measure=-(Data$DISTFRONT[which(Data$groupFLY==0)]-Data$DISTBACK[which(Data$groupFLY==0)])
Delta_move = -(Data$DISTFRONT[which(Data$groupFLY==1)]-Data$DISTBACK[which(Data$groupFLY==1)])
library(ggplot2)
Plot_df = data.frame(Type = c(rep("Flying",length(Delta_move)),rep("Not flying",length(Delta_measure))), Error = c(Delta_move,Delta_measure))
Dist_error_hist = ggplot(Plot_df)+geom_histogram(aes(Error))+facet_grid(.~Type)+ylab('Number of waterfowl groups')+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16),strip.text.x = element_text(size = 16))
png("Dist_error_hists.png")
Dist_error_hist
dev.off()

Data$DISTFRONT[Data$DISTFRONT==6]=NA
#Data$DISTBACK[Data$DISTBACK==6]=NA
Data$Dist = Data$DISTFRONT
Data$Dist[is.na(Data$Dist)]=Data$DISTBACK[is.na(Data$Dist)]
Data$Dist[Data$Dist==6]=NA  # So, only include observations for which distance of first observer was in bins 1-5 OR first observer missed it but second observer saw it in bins 1-5; i.e. observer 2 distance bin 6 is still included as a detection of first observer saw it in 1-5

#Data$Dist[Data$Dist==1]=NA  # take out first bin because it looks like 2nd observer doesn't detect these very well
Data=Data[-which(is.na(Data$Dist)),]
Data$Species=as.character(Data$Species)
Data$Species[Data$Species=="UK LOON"]="UKLOON"
Data$Species = factor(Data$Species)
Which.sp = unique(Data$Species)[which(tabulate(Data$Species)>=20)]  #which species seen 20 or more times
Data=Data[Data$Species %in% Which.sp,]  
Data$Species =  factor(as.character(Data$Species))
Species.list = unique(Data$Species)
n.species=length(Species.list)
Species=Data$Species
Data$Species=as.character(Data$Species)
for(isp in 1:n.species)Data$Species[Data$Species==Species.list[isp]]=isp
Data$Species=as.numeric(Data$Species)
#Data$Species = factor(as.character(Data$Species))  #reset species levels for factor var


n.transects = 1
n.obs = nrow(Data)
ch_split <- function(x){
  as.numeric(c(substr(x,1,1),substr(x,2,2)))
}
Data$CH = as.character(Data$CH)
Data$CH[Data$CH=="1"]="01"
Obs = as.vector(sapply(Data$CH,"ch_split"))
Dat = data.frame(Transect = rep(1,n.obs*2),Match = rep(c(1:n.obs),each=2),Observer=rep(c(1,2),n.obs),Obs=Obs,Species=rep(Data$Species,each=2),Distance = rep(Data$Dist,each=2),Group=rep(Data$GROUPSZ,each=2),Fly=rep(Data$groupFLY,each=2))

Dat_LP = Dat
Dat_LP[,"Distance"]=1

#conduct double observer analysis using L-P and stratified L-P w Huggins in RMARK- note
#estimates are only for # of groups; estimates of # of individuals would need to be made using a 
#bootstrap  with group/cluster size in the numerator of the H-T procedure
Dat_MARK = Data
Dat_MARK$ch = Dat_MARK$CH
Dat_MARK$Cluster = Dat_MARK$GROUPSZ
Dat_MARK$Species = Species
Dat_MARK$Obs11 = 0
Dat_MARK$Obs12 = as.numeric(substr(Dat_MARK$CH,1,1))
Dat_MARK$Dist21 = (Dat_MARK$Dist==2)*1
Dat_MARK$Dist31 = (Dat_MARK$Dist==3)*1
Dat_MARK$Dist41 = (Dat_MARK$Dist==4)*1
Dat_MARK$Dist51 = (Dat_MARK$Dist==5)*1
Dat_MARK$Dist22 = (Dat_MARK$Dist==2)*1
Dat_MARK$Dist32 = (Dat_MARK$Dist==3)*1
Dat_MARK$Dist42 = (Dat_MARK$Dist==4)*1
Dat_MARK$Dist52 = (Dat_MARK$Dist==5)*1
Dat_MARK$FlyDist21 = (Dat_MARK$Dist==2)*Dat_MARK$groupFLY
Dat_MARK$FlyDist31 = (Dat_MARK$Dist==3)*Dat_MARK$groupFLY
Dat_MARK$FlyDist41 = (Dat_MARK$Dist==4)*Dat_MARK$groupFLY
Dat_MARK$FlyDist51 = (Dat_MARK$Dist==5)*Dat_MARK$groupFLY
Dat_MARK$FlyDist22 = (Dat_MARK$Dist==2)*Dat_MARK$groupFLY
Dat_MARK$FlyDist32 = (Dat_MARK$Dist==3)*Dat_MARK$groupFLY
Dat_MARK$FlyDist42 = (Dat_MARK$Dist==4)*Dat_MARK$groupFLY
Dat_MARK$FlyDist52 = (Dat_MARK$Dist==5)*Dat_MARK$groupFLY
Dat_MARK$Ddist = abs(Dat_MARK$DISTBACK-2)
Dat_MARK$Ddist[is.na(Dat_MARK$Ddist)]=0

library('RMark')
HA_proc = process.data(Dat_MARK,model="Huggins",groups="Species")
HA_ddl = make.design.data(HA_proc)

p_model = list(formula=~1,share=TRUE)
HA_mod_dot = mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))

p_model = list(formula=~Species+Cluster+time,share=TRUE)
HA_mod_SOG = mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))

p_model = list(formula=~Species+Cluster+time+groupFLY,share=TRUE)
HA_mod_SOGF = mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))

#p_model = list(formula=~Species+Cluster+time+groupFLY+Dist2:time+Dist3:time+Dist4:time+Dist5:time,share=TRUE)
#HA_mod_SOGFD = mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))

p_model = list(formula=~Species+Cluster+time+groupFLY+Dist2:time+Dist3:time+Dist4:time+Dist5:time,share=TRUE)
HA_mod_SOGFD = mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))

p_model = list(formula=~Species+Cluster+time+groupFLY+Dist2:time+Dist3:time+Dist4:time+Dist5:time+FlyDist2:time+FlyDist3:time+FlyDist4:time+FlyDist5:time,share=TRUE)
HA_mod_SOGFDint = mark(HA_proc,HA_ddl,model.parameters=list(p=p_model))

#X = matrix(0,9,17)  #transformation matrix for computing estimate of group size per species combining flying and not flying.  Note: no loons were ever flying
#diag(X[1:9,1:9]) = 1
#diag(X[1:7,10:16]) = 1
#X[9,17] = 1
#Est_comb = X%*%LP_mod_fly$results$derived[["N Population Size"]][,"estimate"]
#SE_comb = sqrt(diag(X %*% tcrossprod(LP_mod_fly$results$derived.vcv[["N Population Size"]],X)))

######   mrds analysis
library(Distance)

Dat_MRDS = Dat[Dat$Fly==0,]
Dat_MRDS = Dat_MRDS[Dat_MRDS$Species %in% c(1,2,3,4,5,7,8,9),]  #take out
Dat_MRDS$Species[Dat_MRDS$Species>5]=Dat_MRDS$Species[Dat_MRDS$Species>5]-1
Dat_MRDS = Dat_MRDS[Dat_MRDS$Distance>1,]
Dat_MRDS$Distance = Dat_MRDS$Distance-1  #limit to distance bins 2:5 so we can use the same distance key
n.ind = nrow(Dat_MRDS)/2
Dat_MRDS$Match = rep(c(1:n.ind),each=2) #reorder Match to go from 1:n.ind
Dat_MRDS[,"detected"]=Dat_MRDS$Obs
Dat_MRDS[,"distance"]=Dat_MRDS$Distance/4-0.125
Dat_MRDS[,"distbegin"]=Dat_MRDS$distance-0.125
Dat_MRDS[,"distend"]=Dat_MRDS$distance+0.125
Dat_MRDS[,"object"]=Dat_MRDS$Match
Dat_MRDS[,"observer"]=Dat_MRDS$Observer
Dat_MRDS$Dist2 = 1.0*(Dat_MRDS[,"Distance"]==2)
Dat_MRDS$Dist3 = 1.0*(Dat_MRDS[,"Distance"]==3)
Dat_MRDS$Dist4 = 1.0*(Dat_MRDS[,"Distance"]==4)
Dat_MRDS$Species = factor(Dat_MRDS$Species)
#Dat_MRDS$Dist5 = 1.0*(Dat_MRDS[,"distance"]==5)

#Dat_MRDS$Sample.Label=Dat_MRDS$Species
#Dat_MRDS$Region.Label=factor(Dat_MRDS$Species)
fi.mr.dist <- ddf(method='io.fi',mrmodel=~glm(link='logit',formula=~observer*distance+Group+Species), # +observer+Fly+Species+Group), #+observer+Fly+Species+Group),
                  data=Dat_MRDS,meta.data=list(width=1,binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))
fi.mr.dist <- ddf(method='io.fi',mrmodel=~glm(link='logit',formula=~observer+Dist2+Dist3+Dist4+Dist2:observer+Dist3:observer+Dist4:observer+Group+Species), # +observer+Fly+Species+Group), #+observer+Fly+Species+Group),
                  data=Dat_MRDS,meta.data=list(width=1,binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))

#fi.mr.dist <- ddf(method='io.fi',mrmodel=~glm(link='logit',formula=~observer+Dist2+Dist3+Dist4), # +observer+Fly+Species+Group), #+observer+Fly+Species+Group),
#                 data=Dat_MRDS,meta.data=list(width=1,binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))

Region_table = data.frame(Region.Label=factor(c(1:8)),Area=rep(1,8))
Sample_table = data.frame(Region.Label=factor(c(1:8)),Sample.Label=factor(c(1:8)),Effort=rep(1,8))
Cur_dat=Dat_MRDS[c(1:(nrow(Dat_MRDS)/2))*2,] 
Obs_table = data.frame(object=Cur_dat$Match,Region.Label=as.integer(Cur_dat$Species),Sample.Label=as.integer(Cur_dat$Species))
Strata_ests = dht(fi.mr.dist,region.table=Region_table,sample.table=Sample_table,obs.table=Obs_table) #function for computing strata-specific ht estimates
#NEED TO MULTIPLY STRATA ESTS AND SEs by 2 since mrds assumes both sides of line are searched

#mr.dist <- ddf(method='io',dsmodel=~mcds(key="hn",formula=~Group+Species),mrmodel=~glm(link='logit',formula=~observer+Dist2+Dist3+Dist4+Dist2:observer+Dist3:observer+Dist4:observer+Group+Species), #+observer+Fly+Species+Group),
#               data=Dat_MRDS,meta.data=list(binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))
mr.dist.hr <- ddf(method='io',dsmodel=~mcds(key="hr",formula=~Group+Species),mrmodel=~glm(link='logit',formula=~observer+Dist2+Dist3+Dist4+Dist2:observer+Dist3:observer+Dist4:observer+Group+Species), #+observer+Fly+Species+Group),
               data=Dat_MRDS,meta.data=list(binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))
mr.dist.hr.sp <- ddf(method='io',dsmodel=~mcds(key="hr",formula=~Species),mrmodel=~glm(link='logit',formula=~observer+Dist2+Dist3+Dist4+Dist2:observer+Dist3:observer+Dist4:observer+Group+Species), #+observer+Fly+Species+Group),
                  data=Dat_MRDS,meta.data=list(binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))
mr.dist.hr.gr <- ddf(method='io',dsmodel=~mcds(key="hr",formula=~Group),mrmodel=~glm(link='logit',formula=~observer+Dist2+Dist3+Dist4+Dist2:observer+Dist3:observer+Dist4:observer+Group+Species), #+observer+Fly+Species+Group),
                     data=Dat_MRDS,meta.data=list(binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))

Strata_ests_mrds = dht(mr.dist.hr.gr,region.table=Region_table,sample.table=Sample_table,obs.table=Obs_table) #function for computing strata-specific ht estimates

#plot conditional and absolute detection functions - point independence vs. full independence, histogram of all distances
#empirical conditional detection prob
Which.1 = which(Dat_MRDS$detected==1 & Dat_MRDS$observer==1)
Cur_dat=Dat_MRDS[Which.1+1,]
P2giv1=rep(0,4)
for(i in 1:4)P2giv1[i]=sum(Cur_dat$detected==1 & Cur_dat$Distance==i)/sum(Cur_dat$Distance==i)
Which.2 = which(Dat_MRDS$detected==1 & Dat_MRDS$observer==2)
Cur_dat=Dat_MRDS[Which.1-1,]
P1giv2=rep(0,4)
for(i in 1:4)P1giv2[i]=sum(Cur_dat$detected==1 & Cur_dat$Distance==i)/sum(Cur_dat$Distance==i)
Pdet = 1-(1-P2giv1)*(1-P1giv2)

expit<-function(x){1/(1+exp(-x))}
fi.dot = ddf(method='io.fi',mrmodel=~glm(link='logit',formula=~observer*Distance), # +observer+Fly+Species+Group), #+observer+Fly+Species+Group),
                    data=Dat_MRDS,meta.data=list(width=1,binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))
pi.dot= ddf(method='io',dsmodel=~mcds(key="hr",formula=~1),mrmodel=~glm(link='logit',formula=~observer*Distance), #+observer+Fly+Species+Group),
            data=Dat_MRDS,meta.data=list(binned=TRUE,breaks=c(0,.25,0.5,0.75,1)))
DM1 = matrix(0,4,8)
DM1[,1]=1
DM1[2,3]=1
DM1[3,4]=1
DM1[4,5]=1
DM2 = DM1
DM2[,2]=1
DM2[2,6]=1
DM2[3,7]=1
DM2[4,8]=1
P1.fi = expit(DM1 %*% fi.dot$par)
P2.fi = expit(DM2 %*% fi.dot$par)

#realized delta
Hist_scaled = tabulate(Dat_MRDS$Distance[1:(nrow(Dat_MRDS)/1)])
Hist_scaled = Hist_scaled/Hist_scaled[1] * Pdet[1]
Delta = Pdet[1]/Hist_scaled 

Plot_df = data.frame(Distance=c(2,3,4,5),Freq=tabulate(Dat_MRDS$Distance[1:(nrow(Dat_MRDS)/1)]))
Line_df = data.frame(Distance=rep(c(2:5),4),Detection=100*c(P1giv2,P2giv1,Pdet,Delta),Type=c(rep("P(1|2)",4),rep("P(2|1)",4),rep("P(Either)",4),rep("delta",4)))
NotFlying_hists = ggplot() + geom_bar(data=Plot_df,stat="identity",fill="#999999",aes(x=Distance,y=Freq)) + 
  geom_line(data=Line_df,aes(x=Distance,y=Detection,linetype=Type),size=1.0)+
  geom_point(data=Line_df,aes(x=Distance,y=Detection,shape=Type),size=2)+theme(legend.key.width=unit(1.8,"cm"))+
  scale_linetype_manual(values=c("solid","dotted","dashed","dotdash"),labels=expression(delta,paste("P(1|2)"),paste("P(2|1)"),paste("P(1",union(),"2)")))+theme(text = element_text(size=16))+
  scale_shape_manual(values=c(15:18),labels=expression(delta,paste("P(1|2)"),paste("P(2|1)"),paste("P(1",union(),"2)")))
  #scale_colour_manual(values=c('black','black','red','red')) #,labels=expression(delta,paste("P(1|2)"),paste("P(2|1)"),union(1,2)))
  #scale_colour_manual(values=c('P(1|2)'='red','P(2|1)'='black','P(Either)'='black','delta'='black'),name='',labels=c('1','2','3','4'))
pdf('NotFlying_hists.pdf')
  NotFlying_hists
dev.off()
png('NotFlying_hists.png')
NotFlying_hists
dev.off()
NotFlying_hists

tabulate(Dat_MRDS$Distance[1:(nrow(Dat_MRDS)/1)])


#plot various unconditional and conditional detection functions; partition by flying/not flying
Dat_plot=Dat
Dat_plot$Observer=Dat_plot$Observer-1
Dist1 = tabulate(Dat_plot[which(Dat_plot$Observer==0 & Dat_plot$Obs==1),"Distance"])
Dist1Fly1 = tabulate(Dat_plot[which(Dat_plot$Observer==0 & Dat_plot$Obs==1 & Dat_plot$Fly==1),"Distance"])
Dist1Fly0 = tabulate(Dat_plot[which(Dat_plot$Observer==0 & Dat_plot$Obs==1 & Dat_plot$Fly==0),"Distance"])
#Dist1Fly1 = Dist1Fly1/max(Dist1Fly1)
#Dist1Fly0 = Dist1Fly0/max(Dist1Fly0)
plot(Dist1Fly1/max(Dist1Fly1),type="b",col="blue",ylim=c(0,1))
lines(Dist1Fly0/max(Dist1Fly0),type="b",col="green")

Dist2=tabulate(Dat_plot[which(Dat_plot$Observer==1 & Dat_plot$Obs==1),"Distance"])  #all observer 2 distances
Dist2Fly1=tabulate(Dat_plot[which(Dat_plot$Observer==1 & Dat_plot$Obs==1 & Dat_plot$Fly==1),"Distance"])  #all observer 2 distances
Dist2Fly0=tabulate(Dat_plot[which(Dat_plot$Observer==1 & Dat_plot$Obs==1 & Dat_plot$Fly==0),"Distance"])  #all observer 2 distances
#Dist2Fly1 = Dist2Fly1/max(Dist2Fly1)
#Dist2Fly0 = Dist2Fly0/max(Dist2Fly0)
plot(Dist2Fly1/max(Dist2Fly1),type="b",col="blue",ylim=c(0,1))
lines(Dist2Fly0/max(Dist2Fly0),type="b",col="green")

#plot conditional detections 
Temp = Dat_plot[WhichObs+1,]  #condition on 1st observer detection
Dist2giv1Fly0=tabulate(Temp[which(Temp$Obs==1 & Temp$Fly==0),"Distance"])  #observer 2 distances for all also detected by Obs 1
#Dist2giv1Fly0=Dist2giv1Fly0/max(Dist2giv1Fly0)
Dist2giv1Fly1=tabulate(Temp[which(Temp$Obs==1 & Temp$Fly==1),"Distance"])  #observer 2 distances for all also detected by Obs 1
#Dist2giv1Fly1=Dist2giv1Fly1/max(Dist2giv1Fly1)

#
WhichObs2 = which(Dat_plot$Obs==1 & Dat_plot$Observer==1)
Temp = Dat_plot[WhichObs-1,]  #condition on 1st observer detection
Dist1giv2Fly0=tabulate(Temp[which(Temp$Obs==1 & Temp$Fly==0),"Distance"])  #observer 2 distances for all also detected by Obs 1
#Dist1giv2Fly0=Dist1giv2Fly0/max(Dist1giv2Fly0)
Dist1giv2Fly1=tabulate(Temp[which(Temp$Obs==1 & Temp$Fly==1),"Distance"])  #observer 2 distances for all also detected by Obs 1
#Dist1giv2Fly1=Dist1giv2Fly1/max(Dist1giv2Fly1)

Plot_df = data.frame(
  Freq=(c(Dist1Fly1,Dist1Fly0,Dist2Fly1,Dist2Fly0,Dist1giv2Fly1,Dist1giv2Fly0,Dist2giv1Fly1,Dist2giv1Fly0)),
  Observer = rep(c("Observer 1","Observer 1","Observer 2","Observer 2","Observer 1","Observer 1","Observer 2","Observer 2"),each=5),
  Fly = factor(rep(c("Flying","Not flying","Flying","Not flying","Flying","Not flying","Flying","Not flying"),each=5)),
  Distance = rep(c(1:5),8),
  Conditional = factor(rep(c("Unconditional","Unconditional","Unconditional","Unconditional","Conditional","Conditional","Conditional","Conditional"),each=5))
)
Plot_df=Plot_df[-which(Plot_df[,"Conditional"]=="Conditional"),]

#ggplot(Plot_df,aes(x=Distance,y=Freq,linetype=Fly,colour=Conditional))+geom_line(size=1.2)+facet_wrap(~Observer)

#calculate (empirical) conditional detection probability
library(doBy)
WhichObs = which(Dat_plot$Obs==1 & Dat_plot$Observer==0 & Dat_plot$Fly==1)
Temp = data.frame(Obs=Dat_plot[WhichObs+1,"Obs"],Dist=Dat_plot[WhichObs,"Distance"])
Det_obs2_fly1 = summaryBy(Obs~Dist,data=Temp,FUN=mean)$Obs.mean
WhichObs = which(Dat_plot$Obs==1 & Dat_plot$Observer==0 & Dat_plot$Fly==0)
Temp = data.frame(Obs=Dat_plot[WhichObs+1,"Obs"],Dist=Dat_plot[WhichObs,"Distance"])
Det_obs2_fly0 = summaryBy(Obs~Dist,data=Temp,FUN=mean)$Obs.mean
WhichObs = which(Dat_plot$Obs==1 & Dat_plot$Observer==1 & Dat_plot$Fly==1)
Temp = data.frame(Obs=Dat_plot[WhichObs-1,"Obs"],Dist=Dat_plot[WhichObs,"Distance"])
Det_obs1_fly1 = summaryBy(Obs~Dist,data=Temp,FUN=mean)$Obs.mean
WhichObs = which(Dat_plot$Obs==1 & Dat_plot$Observer==1 & Dat_plot$Fly==0)
Temp = data.frame(Obs=Dat_plot[WhichObs+1,"Obs"],Dist=Dat_plot[WhichObs,"Distance"])
Det_obs1_fly0 = summaryBy(Obs~Dist,data=Temp,FUN=mean)$Obs.mean

Line_df = data.frame(Detection=c(Det_obs1_fly0*100,Det_obs1_fly1*100,Det_obs2_fly0*100,Det_obs2_fly1*100),
                     Fly=rep(c("Not flying","Flying","Not flying","Flying"),each=5),
                     Observer=rep(c("Observer 1","Observer 2"),each=10),
                     Distance=rep(c(1:5),4)
)

myPalette =  c("#999999", "#E69F00")
DF_expanded <- Plot_df[rep(row.names(Plot_df), Plot_df$Freq), ]
Dist_hists = ggplot() + geom_bar(data=DF_expanded,fill="#999999",aes(x=Distance)) + 
  geom_line(data=Line_df,aes(x=Distance,y=Detection),size=1)+
  geom_point(data=Line_df,aes(x=Distance,y=Detection),size=1.5)+
  facet_grid(Fly~Observer) #+ labs(fill="") + scale_fill_manual(values=myPalette)
png("Distance_hists.png")
  Dist_hists
dev.off()



