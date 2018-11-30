rm(list = ls())
library(gtools)
library(GREMLIN)


#
# nFG <- 4;
 namesfg <- c('famers','plants')
# vK <- c(3,2)
# vNQ <- c(30,37)
# E = matrix(1,2,2)
# E[1,] = c(1,1)
# E[2,] = c(1,2)
# nnet <- 2
# typeinter = c('diradj','inc')
#

#---------------------- parameters
# lpi = list(c(0.31,0.42,0.27),c(0.65,0.35))
# ltheta   <- lapply(1:nnet,function(i){matrix(0,vK[E[i,1]],vK[E[i,2]])})
# ltheta[[1]][1,] = c(0.025,0.123,0.053)
# ltheta[[1]][2,]  = c(0.159, 0.300,0.070)
# ltheta[[1]][3,] =  c(0.374, 0.585,  0.357)
#
# ltheta[[2]][1,] = c(0.186 ,0.653)
# ltheta[[2]][2,] = c(0.559,0.905 )
# ltheta[[2]][3,] = c( 0.390,0.696)
#
# save(ltheta,lpi,file = '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/paramSim1.Rdata')
#------------------------SIMULATION -----------------------------
nSimu = 100 ;
vdistrib = rep('bernoulli',2)
dirSaveSimuData <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/data'
#
# for (i in 1:nSimu) {
#   datasim <- rMBM(vNQ ,E , typeinter, vdistrib, lpi, ltheta, seed = NULL, namesfg = namesfg,keepclassif = TRUE)
#   namedata  <- paste('datasim',i,'.Rdata',sep = "")
#   save(datasim,file = paste(dirSaveSimuData,namedata,sep='/' ))
# }
#

#---------------- ESTIMATION ------------------------------------
dirSaveSimuRes <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/res_MBM'
for (i in 1:nSimu)
{
  print(i)
  #load data
  namedata  <- paste('datasim',i,'.Rdata',sep = "")
  load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
  listNet <- datasim$list_net
  # estim
  res_estim <- MultipartiteBM(listNet, namesfg = namesfg, vdistrib = vdistrib , vKmin = 1 , vKmax = 10 , vKinit = c(1,1), init.BM = TRUE, save = FALSE , verbose = FALSE,nb_cores = 10)
  nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
  save(res_estim, file = paste(dirSaveSimuRes,nameres,sep = '/' ))
}
#

#-----------------------------------------------------
#------------------- EXPLOITATION -------------------
#------------------------------------------------------
#
# truevK <- vK_estim <- resCompar <- matrix(0,nSimu,2)
# for (i in 1:nSimu)
# {
#   print(i)
#   #load data
#   namedata  <- paste('datasim',i,'.Rdata',sep = "")
#   load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
#   #listNet <- datasim$list_net
#   trueZ <- datasim$classif
#   truevK[i,] <- vapply(1:2,function(q){length(unique(trueZ[[q]]))},1)
#   # estim
#   #res_estim <- MultipartiteBM(listNet, namesfg = namesfg, vdistrib = vdistrib , vKmin = 1 , vKmax = 10 , vKinit = c(1,1,1,1), init.BM = TRUE, save = FALSE , verbose = FALSE,nb_cores = 10)
#   nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
#   load(file = paste(dirSaveSimuRes,nameres,sep = '/' ))
#   estimZ <- res_estim$fitted.model[[1]]$param_estim$Z
#   resCompar[i,] <-  comparClassif(trueZ,estimZ)
#   vK_estim[i,] <- res_estim$fitted.model[[1]]$param_estim$vK
# }
#
# vK_estim - truevK
# summary(resCompar)
# par(mfrow=c(2,2))
# for (q in 1:2){plot(density(resCompar[,q])); abline(v=1,col='red')}

