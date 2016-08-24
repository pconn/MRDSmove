Data=read.csv('Alisauskas_MRDS_2014.csv')

hist(Data$DISTFRONT[which(Data$groupFLY==0)]-Data$DISTBACK[which(Data$groupFLY==0)])

hist(Data$DISTFRONT[which(Data$groupFLY==0 & Data$DISTFRONT==2)]-Data$DISTBACK[which(Data$groupFLY==0 & Data$DISTFRONT==2)])

hist(Data$DISTFRONT)
hist(Data$DISTBACK)