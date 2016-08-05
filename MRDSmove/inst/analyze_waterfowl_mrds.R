Data=read.csv('Alisauskas_MRDS_2014.csv')

hist(Data$DISTFRONT[which(Data$groupFLY==0)]-Data$DISTBACK[which(Data$groupFLY==0)])