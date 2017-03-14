### tabulate simualtion results

#Simulation study 1
load('c:/users/paul.conn/git/alisauskas/output/SimResultsPI.Rdata')

#produce results by 6 scenarios, 3 estimation models

RelBias = RelBias.median = CV = Cover = CV.median = RMSE = rep(0,18)

irow=1 
for(iscen in 1:6){
  for(imod in 1:3){
    Which.conv = which(is.na(Out$N.est[,iscen,imod])==FALSE & Out$N.est[,iscen,imod]>0)
    RelBias[irow]=mean((Out$N.est[Which.conv,iscen,imod]-Out$N.true[Which.conv,iscen])/Out$N.true[Which.conv,iscen])
    RelBias.median[irow]=median((Out$N.est[Which.conv,iscen,imod]-Out$N.true[Which.conv,iscen])/Out$N.true[Which.conv,iscen])
    CV[irow] = mean(Out$SE[Which.conv,iscen,imod]/Out$N.est[Which.conv,iscen,imod])
    CV.median[irow] = median(Out$SE[Which.conv,iscen,imod]/Out$N.est[Which.conv,iscen,imod])
    Cover[irow] = mean(Out$Cov[Which.conv,iscen,imod])
    RMSE[irow] = mean((Out$N.est[Which.conv,iscen,imod]-Out$N.true[Which.conv,iscen])^2)
    irow=irow+1
  }
}

cbind(RelBias,CV,Cover,RMSE)

#Simulation study 2
load('c:/users/paul.conn/git/alisauskas/output/SimResultsPI.Rdata')
RelBias = RelBias.median = CV = Cover = CV.median = RMSE = rep(0,16)

irow=1 
for(iN in 1:2){
  for(iscen in 1:2){
    for(imod in 1:4){
      Which.conv = which(is.na(Out$N.est[,iN,iscen,imod])==FALSE & Out$N.est[,iN,iscen,imod]>0)
      RelBias[irow]=mean((Out$N.est[Which.conv,iN,iscen,imod]-Out$N.true[Which.conv,iN,iscen])/Out$N.true[Which.conv,iN,iscen])
      RelBias.median[irow]=median((Out$N.est[Which.conv,iN,iscen,imod]-Out$N.true[Which.conv,iN,iscen])/Out$N.true[Which.conv,iN,iscen])
      CV[irow] = mean(Out$SE[Which.conv,iN,iscen,imod]/Out$N.est[Which.conv,iN,iscen,imod])
      CV.median[irow] = median(Out$SE[Which.conv,iN,iscen,imod]/Out$N.est[Which.conv,iN,iscen,imod])
      Cover[irow] = mean(Out$Cov[Which.conv,iN,iscen,imod])
      RMSE[irow] = mean((Out$N.est[Which.conv,iN,iscen,imod]-Out$N.true[Which.conv,iN,iscen])^2)
      irow=irow+1
    }
  }
}

cbind(RelBias,CV,Cover,RMSE)



