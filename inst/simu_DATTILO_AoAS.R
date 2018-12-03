rm(list = ls())
library(gtools)
library(GREMLIN)



nFG <- 4;
namesfg <- c('plants','pollinators','birds','ants')
vK <- c(7,2,2,1)
vNQ <- c(141,173,46,30)
E = matrix(1,nFG - 1,2)
E[,2] = c(2,3,4)
nnet <- 3
typeinter = c('inc','inc','inc')

#---------------------- parameters
lpi = list(c(0.4675,0.1606, 0.1351,0.0784,0.1061, 0.0142,0.0381),c(0.06,0.94),c(0.1,0.9),c(1))
lpi[[2]] = c(0.1,0.9)
lpi[[1]] = c(0.46,0.16, 0.15,0.184,0.1061, 0.1,0.1)
lpi[[1]] = lpi[[1]]/sum(lpi[[1]])
ltheta   <- lapply(1:nnet,function(i){matrix(0,vK[E[i,1]],vK[E[i,2]])})
ltheta[[1]][1,] = c(0.0957,0.0075)
ltheta[[1]][2,1]  = 0.01
ltheta[[1]][3,2] = 0.0003
ltheta[[1]][4,] = c(0.1652, 0.0343)
ltheta[[1]][5,] = c(0.2018,0.138)

ltheta[[2]][1,2] = 0.0006
ltheta[[2]][2,] = c(0.5431, 0.0589)
ltheta[[2]][4,] = c(0.6620,0.1542)
ltheta[[2]][7,] = c(0.8492,0.3565)

ltheta[[3]][1,] = c(0.0013)
ltheta[[3]][3,] = 0.0753
ltheta[[3]][5,] = 0.0163
ltheta[[3]][6,] = 0.5108

save(ltheta,lpi,file ='/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param2/paramSim2.Rdata')
#------------------------SIMULATION -----------------------------
nSimu = 100 ;
vdistrib = rep('bernoulli',3)
dirSaveSimuData <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param2/data'

# for (i in 1:nSimu){
#   datasim <- rMBM(vNQ ,E , typeinter, vdistrib, lpi, ltheta, seed = NULL, namesfg = namesfg,keepclassif = TRUE)
#   namedata  <- paste('datasim',i,'.Rdata',sep = "")
#   save(datasim,file=paste(dirSaveSimuData,namedata,sep='/' ))
# }


#---------------- ESTIMATION ------------------------------------
dirSaveSimuRes <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param2/res_MBM'
for (i in 1:nSimu)
  {
  print(i)
  #load data
  namedata  <- paste('datasim',i,'.Rdata',sep = "")
  load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
  list_Net <- datasim$list_Net
  # estim
  res_estim <- MultipartiteBM(list_Net, namesfg = namesfg, vdistrib = vdistrib , vKmin = 1 , vKmax = 10 , vKinit = c(1,1,1,1), init.BM = TRUE, save = FALSE , verbose = FALSE,nb_cores = 10)
  nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
  save(res_estim, file = paste(dirSaveSimuRes,nameres,sep = '/' ))
}
#

#-----------------------------------------------------
#------------------- EXPLOITATION -------------------
#------------------------------------------------------

truevK <- vK_estim <- resCompar <- matrix(0,nSimu,4)
for (i in 1:nSimu)
{
  print(i)
  #load data
  namedata  <- paste('datasim',i,'.Rdata',sep = "")
  load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
  #listNet <- datasim$list_net
  trueZ <- datasim$classif
  truevK[i,] <- vapply(1:4,function(q){length(unique(trueZ[[q]]))},1)
  # estim
  #res_estim <- MultipartiteBM(listNet, namesfg = namesfg, vdistrib = vdistrib , vKmin = 1 , vKmax = 10 , vKinit = c(1,1,1,1), init.BM = TRUE, save = FALSE , verbose = FALSE,nb_cores = 10)
  nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
  load(file = paste(dirSaveSimuRes,nameres,sep = '/' ))
  estimZ <- res_estim$fitted.model[[1]]$param_estim$Z
  resCompar[i,] <-  comparClassif(trueZ,estimZ)
  vK_estim[i,] <- res_estim$fitted.model[[1]]$param_estim$vK
}

vK_estim - truevK
summary(resCompar)
par(mfrow=c(2,2))
for (q in 1:4){plot(density(resCompar[,q])); abline(v=1,col='red')}

